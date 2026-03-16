#+tests(in-package #:clasp-tests)

;;; LK solvation regression test
;;;
;;; Tests that the C++ LK solvation energy computation matches an independent
;;; implementation of the mathematical model derived from:
;;;
;;;   Gaussian desolvation function:
;;;     Lazaridis & Karplus (1999) "Effective energy function for proteins in
;;;     solution." Proteins 35(2):133-152.  Appendix, Equation A3.
;;;
;;;   5-region piecewise form with cubic spline transitions:
;;;     Alford et al. (2017) "The Rosetta All-Atom Energy Function for
;;;     Macromolecular Modeling and Design." J. Chem. Theory Comput. 13:3031-3048.
;;;     Equations 11-13.
;;;
;;;   Cubic spline coefficients:
;;;     Standard cubic Hermite interpolation (see derivation in comments below).
;;;
;;; The Lisp model is derived entirely from these published sources and basic
;;; calculus, NOT translated from the C++ implementation.

;;; ============================================================================
;;; Rosetta smoothing constants (Alford 2017, text below eq 12)
;;; ============================================================================

(defparameter *lk-c0* 0.3d0)       ; contact offset (Angstroms)
(defparameter *lk-c1* 0.2d0)       ; Gaussian onset offset (Angstroms)
(defparameter *lk-r-solv-low* 4.5d0)  ; high-spline onset (Angstroms)
(defparameter *lk-r-solv-high* 6.0d0) ; cutoff distance (Angstroms)

;;; ============================================================================
;;; Per-atom test parameters (arbitrary but physically reasonable values)
;;; ============================================================================

(defparameter *lk-test-type* :CH2)
(defparameter *lk-test-lj-radius* 2.0d0)    ; Angstroms
(defparameter *lk-test-lk-dgfree* -1.7d0)   ; kcal/mol
(defparameter *lk-test-lk-lambda* 3.5d0)    ; Angstroms
(defparameter *lk-test-lk-volume* 14.7d0)   ; Angstroms^3

;;; ============================================================================
;;; Force field and aggregate setup (same pattern as rosetta-nonbond.lisp)
;;; ============================================================================

(defun lksolvation--ensure-force-field ()
  "Create a force field with known LK solvation parameters and register it as :rosetta."
  (let* ((force-field (chem:force-field/make))
         (lkdb (chem:get-lksolvation-db force-field))
         (fflk (chem:make-fflksolvation *lk-test-type*
                                        :lj-radius *lk-test-lj-radius*
                                        :lk-dgfree *lk-test-lk-dgfree*
                                        :lk-lambda *lk-test-lk-lambda*
                                        :lk-volume *lk-test-lk-volume*))
         ;; Also need nonbond entry so the energy function can be built
         (nonbond-db (chem:get-nonbond-db force-field))
         (ffnonbond (chem:make-ffnonbond *lk-test-type*)))
    (chem:fflksolvation-db-add lkdb fflk)
    (chem:ffnonbond/set-radius-angstroms ffnonbond *lk-test-lj-radius*)
    (chem:ffnonbond/set-epsilon-kcal ffnonbond 0.2d0)
    (chem:ffnonbond-db-add nonbond-db ffnonbond)
    (when (chem:find-force-field :rosetta nil)
      (leap:clear-force-field :rosetta))
    (leap.core:add-force-field-or-modification
     force-field
     :force-field-name :rosetta
     :combined-force-field-class-name 'chem:combined-force-field)
    fflk))

(defun lksolvation--make-two-atom-aggregate (atom-type)
  "Create a two-atom aggregate in separate molecules with the given LK solvation atom type."
  (let* ((agg (chem:make-aggregate))
         (mol1 (chem:make-molecule :m1))
         (mol2 (chem:make-molecule :m2))
         (res1 (chem:make-residue :r1))
         (res2 (chem:make-residue :r2))
         (a1 (chem:make-atom :a1 :c))
         (a2 (chem:make-atom :a2 :c)))
    (chem:add-atom res1 a1)
    (chem:add-atom res2 a2)
    (chem:add-matter mol1 res1)
    (chem:add-matter mol2 res2)
    (chem:add-molecule agg mol1)
    (chem:add-molecule agg mol2)
    (chem:setf-force-field-name agg :rosetta)
    (chem:setf-force-field-name mol1 :rosetta)
    (chem:setf-force-field-name mol2 :rosetta)
    ;; Set both :given-atom-type (for nonbond) and :lk_solvation_atom_type (for lksolvation)
    (chem:set-property a1 :given-atom-type atom-type)
    (chem:set-property a2 :given-atom-type atom-type)
    (chem:set-property a1 :lk-solvation-atom-type atom-type)
    (chem:set-property a2 :lk-solvation-atom-type atom-type)
    (values agg a1 a2)))

;;; ============================================================================
;;; The Rosetta LK Gaussian desolvation function
;;;
;;; The base function comes from Lazaridis & Karplus 1999, Appendix eq A3.
;;; Rosetta adopts a sign convention (Alford et al. 2017, eq 11) where
;;; desolvation is a POSITIVE energy penalty:
;;;
;;;   g_desolv_i(r) = -V_j * DG_i^free / (2 * pi^(3/2) * lambda_i * r^2)
;;;                   * exp( -(r - R_i)^2 / lambda_i^2 )
;;;
;;; Note the leading minus sign vs the original LK formula. In the LK paper,
;;; the solvation energy is DG_slv = DG_ref - Sum(f_i), so f_i with DG < 0
;;; (polar atoms) is negative. Rosetta flips the sign so that g_desolv > 0
;;; represents the desolvation penalty.
;;;
;;; Variables:
;;;   DG_i^free = vapor-to-water transfer free energy of atom i  (lk_dgfree)
;;;   V_j       = atomic volume of the desolvating atom j        (lk_volume)
;;;   lambda_i  = correlation length for atom i                  (lk_lambda)
;;;   R_i       = van der Waals radius of atom i                 (lj_radius)
;;;   r         = interatomic distance
;;; ============================================================================

(defun lk-gaussian (r dg-free-i volume-j lambda-i radius-i)
  "Evaluate the Rosetta LK Gaussian desolvation function g_desolv_i(r).
Alford et al. 2017, eq 11; sign convention: positive = desolvation penalty."
  (let* ((prefactor (/ (* (- dg-free-i) volume-j)
                       (* 2.0d0 (expt pi 1.5d0) lambda-i)))
         (delta (- r radius-i))
         (lambda2 (* lambda-i lambda-i)))
    (* prefactor (/ 1.0d0 (* r r)) (exp (/ (- (* delta delta)) lambda2)))))

;;; ============================================================================
;;; Analytical derivative of the LK Gaussian
;;;
;;; From calculus, differentiating f_i(r) = P * exp(-(r-R)^2/L^2) / r^2
;;; where P = DG_i * V_j / (2 * pi^(3/2) * lambda_i):
;;;
;;;   df_i/dr = f_i(r) * [ -2/r  -  2*(r - R_i) / lambda_i^2 ]
;;;           = -2 * f_i(r) * [ 1/r + (r - R_i) / lambda_i^2 ]
;;;
;;; This follows from the product rule on g(r) = exp(...)/r^2:
;;;   d/dr [exp(-u^2)/r^2] = exp(-u^2) * [-2u/r^2 * du/dr + (-2/r^3)]
;;; with u = (r - R_i)/lambda_i, du/dr = 1/lambda_i.
;;; ============================================================================

(defun lk-gaussian-derivative (r dg-free-i volume-j lambda-i radius-i)
  "Evaluate df_i/dr for the LK Gaussian. Derived analytically from lk-gaussian."
  (let* ((f (lk-gaussian r dg-free-i volume-j lambda-i radius-i))
         (delta (- r radius-i))
         (lambda2 (* lambda-i lambda-i)))
    (* -2.0d0 f (+ (/ 1.0d0 r) (/ delta lambda2)))))

;;; ============================================================================
;;; Cubic Hermite interpolation
;;;
;;; Given function values and derivatives at two endpoints [a, b]:
;;;   f(a), f'(a), f(b), f'(b)
;;;
;;; The unique cubic p(t) in the shifted variable t = r - a, with h = b - a, is:
;;;   p(t) = (2*f_a - 2*f_b + df_b*h + df_a*h) * (t/h)^3
;;;        + (-3*f_a + 3*f_b - 2*df_a*h - df_b*h) * (t/h)^2
;;;        + df_a * t
;;;        + f_a
;;;
;;; This is standard cubic Hermite interpolation (see e.g. Burden & Faires,
;;; "Numerical Analysis", Ch. 3.4).
;;;
;;; To evaluate as a polynomial in r (not t), we expand t = r - a and collect
;;; terms: p(r) = A*r^3 + B*r^2 + C*r + D
;;; ============================================================================

(defun cubic-hermite-coefficients (a b f-a df-a f-b df-b)
  "Compute coefficients (A B C D) of the cubic p(r) = A*r^3 + B*r^2 + C*r + D
that interpolates f(a)=f-a, f'(a)=df-a, f(b)=f-b, f'(b)=df-b.
Returns (values A B C D)."
  (let* ((h (- b a))
         (h2 (* h h))
         (h3 (* h2 h))
         ;; Coefficients in shifted basis t = r - a
         (a-shift (/ (+ (* 2.0d0 f-a) (* -2.0d0 f-b) (* df-b h) (* df-a h)) h3))
         (b-shift (/ (+ (* -3.0d0 f-a) (* 3.0d0 f-b) (* -2.0d0 df-a h) (* -1.0d0 df-b h)) h2))
         (c-shift df-a)
         (d-shift f-a)
         ;; Convert to standard polynomial in r by substituting t = r - a
         ;; p(r) = a-shift*(r-a)^3 + b-shift*(r-a)^2 + c-shift*(r-a) + d-shift
         (aa a-shift)
         (bb (- b-shift (* 3.0d0 a-shift a)))
         (cc (+ c-shift (* -2.0d0 b-shift a) (* 3.0d0 a-shift a a)))
         (dd (+ d-shift (- (* c-shift a)) (* b-shift a a) (- (* a-shift a a a)))))
    (values aa bb cc dd)))

(defun eval-cubic (r aa bb cc dd)
  "Evaluate p(r) = dd + cc*r + r^2*(bb + aa*r)."
  (+ dd (* cc r) (* r r (+ bb (* aa r)))))

;;; ============================================================================
;;; Complete piecewise LK solvation energy
;;;
;;; From Alford et al. 2017, eq 12:
;;;
;;;   g_desolv(i,j,d) =
;;;     f_desolv(i,j, sigma_ij)                 if d <= sigma_ij - c0
;;;     f_poly^{solv,low}(i,j, d)               if sigma_ij - c0 < d <= sigma_ij + c1
;;;     f_desolv(i,j, d)                         if sigma_ij + c1 < d <= 4.5 A
;;;     f_poly^{solv,high}(i,j, d)              if 4.5 A < d <= 6.0 A
;;;     0                                        if 6.0 A < d
;;;
;;; The low cubic spline interpolates between:
;;;   - At r = r0_low = sigma_ij - c0: value = f(sigma_ij), derivative = 0
;;;     (contact region is flat, so zero derivative at transition)
;;;   - At r = r1_low = sigma_ij + c1: value = f(r1_low), derivative = f'(r1_low)
;;;     (matches Gaussian value and slope)
;;;
;;; The high cubic spline interpolates between:
;;;   - At r = r2_high = 4.5: value = f(r2_high), derivative = f'(r2_high)
;;;     (matches Gaussian value and slope)
;;;   - At r = r3_high = 6.0: value = 0, derivative = 0
;;;     (smooth transition to zero)
;;;
;;; Total energy for a pair: g_desolv(i,j,d) + g_desolv(j,i,d)
;;; (Alford 2017, eq 13, without w^conn which is 1 for our two-atom test)
;;; ============================================================================

(defun lksolvation--expected-energy (r radius-i dg-free-i lambda-i volume-i
                                       radius-j dg-free-j lambda-j volume-j)
  "Compute the expected LK solvation energy for an atom pair at distance r.
Implements the 5-region piecewise function from Alford 2017 eq 12,
with the Gaussian from Lazaridis & Karplus 1999 eq A3,
and cubic Hermite splines at the transition boundaries."
  (let* ((sigma-ij (+ radius-i radius-j))
         (r0-low (- sigma-ij *lk-c0*))
         (r1-low (+ sigma-ij *lk-c1*))
         (r2-high *lk-r-solv-low*)
         (r3-high *lk-r-solv-high*))
    (cond
      ;; ---- Region 1: Contact (r <= sigma_ij - c0) ----
      ;; Energy is constant: f_desolv evaluated at r = sigma_ij
      ((<= r r0-low)
       (+ (lk-gaussian sigma-ij dg-free-i volume-j lambda-i radius-i)
          (lk-gaussian sigma-ij dg-free-j volume-i lambda-j radius-j)))

      ;; ---- Region 2: Low cubic spline (r0_low < r <= r1_low) ----
      ;; Hermite spline connecting contact plateau to Gaussian
      ((<= r r1-low)
       (let ((fi-contact (lk-gaussian sigma-ij dg-free-i volume-j lambda-i radius-i))
             (fi-r1      (lk-gaussian r1-low   dg-free-i volume-j lambda-i radius-i))
             (dfi-r1     (lk-gaussian-derivative r1-low dg-free-i volume-j lambda-i radius-i))
             (fj-contact (lk-gaussian sigma-ij dg-free-j volume-i lambda-j radius-j))
             (fj-r1      (lk-gaussian r1-low   dg-free-j volume-i lambda-j radius-j))
             (dfj-r1     (lk-gaussian-derivative r1-low dg-free-j volume-i lambda-j radius-j)))
         ;; Boundary conditions:
         ;;   At r0_low: value = f_contact, derivative = 0 (flat contact)
         ;;   At r1_low: value = f(r1_low), derivative = f'(r1_low)
         (multiple-value-bind (ai bi ci di)
             (cubic-hermite-coefficients r0-low r1-low fi-contact 0.0d0 fi-r1 dfi-r1)
           (multiple-value-bind (aj bj cj dj)
               (cubic-hermite-coefficients r0-low r1-low fj-contact 0.0d0 fj-r1 dfj-r1)
             (+ (eval-cubic r ai bi ci di)
                (eval-cubic r aj bj cj dj))))))

      ;; ---- Region 3: Main Gaussian (r1_low < r <= r2_high) ----
      ;; Direct evaluation of LK Gaussian, eq 11/A3
      ((<= r r2-high)
       (+ (lk-gaussian r dg-free-i volume-j lambda-i radius-i)
          (lk-gaussian r dg-free-j volume-i lambda-j radius-j)))

      ;; ---- Region 4: High cubic spline (r2_high < r <= r3_high) ----
      ;; Hermite spline connecting Gaussian to zero
      ((<= r r3-high)
       (let ((fi-r2  (lk-gaussian r2-high dg-free-i volume-j lambda-i radius-i))
             (dfi-r2 (lk-gaussian-derivative r2-high dg-free-i volume-j lambda-i radius-i))
             (fj-r2  (lk-gaussian r2-high dg-free-j volume-i lambda-j radius-j))
             (dfj-r2 (lk-gaussian-derivative r2-high dg-free-j volume-i lambda-j radius-j)))
         ;; Boundary conditions:
         ;;   At r2_high: value = f(r2_high), derivative = f'(r2_high)
         ;;   At r3_high: value = 0, derivative = 0
         (multiple-value-bind (ai bi ci di)
             (cubic-hermite-coefficients r2-high r3-high fi-r2 dfi-r2 0.0d0 0.0d0)
           (multiple-value-bind (aj bj cj dj)
               (cubic-hermite-coefficients r2-high r3-high fj-r2 dfj-r2 0.0d0 0.0d0)
             (+ (eval-cubic r ai bi ci di)
                (eval-cubic r aj bj cj dj))))))

      ;; ---- Region 5: Beyond cutoff (r > r3_high) ----
      (t 0.0d0))))

;;; ============================================================================
;;; Energy function component extraction
;;; ============================================================================

(defun lksolvation--component (energy-function)
  (or (find-if (lambda (comp) (typep comp 'chem:energy-rosetta-lksolvation))
               (chem:all-components energy-function))
      (error "No rosetta lksolvation component found in energy-function")))

;;; ============================================================================
;;; XY pair generation (for plotting/debugging)
;;; ============================================================================

(defun lksolvation-model-xy (&key (start 1.0d0) (end 8.0d0) (step 0.05d0)
                                  (lj-radius *lk-test-lj-radius*)
                                  (lk-dgfree *lk-test-lk-dgfree*)
                                  (lk-lambda *lk-test-lk-lambda*)
                                  (lk-volume *lk-test-lk-volume*)
                                  (pathname nil))
  "Generate distance vs energy pairs from the pure LK solvation model.
No C++ code is invoked -- this is entirely from the paper equations.
When PATHNAME is provided, writes tab-separated output to that file.
Returns a list of (distance energy) pairs."
  (let ((pairs nil))
    (flet ((generate (stream)
             (when stream
               (format stream "# LK solvation model (Lazaridis-Karplus 1999 / Alford 2017)~%")
               (format stream "# lj_radius=~f lk_dgfree=~f lk_lambda=~f lk_volume=~f~%"
                       lj-radius lk-dgfree lk-lambda lk-volume)
               (format stream "# sigma_ij=~f  c0=~f c1=~f r_solv_low=~f r_solv_high=~f~%"
                       (+ lj-radius lj-radius) *lk-c0* *lk-c1* *lk-r-solv-low* *lk-r-solv-high*)
               (format stream "# distance~tenergy~%"))
             (loop for dist from start to end by step
                   for energy = (lksolvation--expected-energy
                                 dist
                                 lj-radius lk-dgfree lk-lambda lk-volume
                                 lj-radius lk-dgfree lk-lambda lk-volume)
                   do (push (list dist energy) pairs)
                      (when stream
                        (format stream "~,6F~t~,10,,,,,'eE~%" dist energy)))))
      (if pathname
          (with-open-file (fout pathname :direction :output :if-exists :supersede)
            (generate fout)
            (format t "Wrote ~d points to ~a~%" (length pairs) pathname))
          (generate nil)))
    (nreverse pairs)))

(defun lksolvation-generate-xy-pairs (&key (start 0.1d0) (end 8.0d0) (step 0.01d0)
                                        (energy-source :actual) pathname)
  "Generate x,y pairs comparing C++ vs model. ENERGY-SOURCE is :actual (C++) or :expected (Lisp model)."
  (unless pathname
    (error "You must provide a pathname"))
  (lksolvation--ensure-force-field)
  (multiple-value-bind (agg a1 a2)
      (lksolvation--make-two-atom-aggregate *lk-test-type*)
    (let* ((setup (list :rosetta (list 'chem:energy-rosetta-lksolvation)))
           (energy-function (chem:make-energy-function :matter agg
                                                       :use-excluded-atoms nil
                                                       :assign-types nil
                                                       :setup setup))
           (component (lksolvation--component energy-function))
           (pos (chem:make-nvector (chem:get-nvector-size energy-function)))
           (pairs nil))
      (with-open-file (stream pathname :direction :output :if-exists :supersede)
        (format stream "# lksolvation energy~%")
        (loop for dist from start to end by step
              do (chem:set-position a1 (geom:vec 0.0d0 0.0d0 0.0d0))
                 (chem:set-position a2 (geom:vec dist 0.0d0 0.0d0))
                 (chem:load-coordinates-into-vector energy-function pos)
                 (let* ((energy (ecase energy-source
                                  (:actual (chem:energy-component-evaluate-energy
                                            energy-function component pos))
                                  (:expected (lksolvation--expected-energy
                                              dist
                                              *lk-test-lj-radius* *lk-test-lk-dgfree*
                                              *lk-test-lk-lambda* *lk-test-lk-volume*
                                              *lk-test-lj-radius* *lk-test-lk-dgfree*
                                              *lk-test-lk-lambda* *lk-test-lk-volume*)))))
                   (push (list dist energy) pairs)
                   (when stream
                     (format stream "~,6f ~,12,,,,,'eE~%" dist energy)))))
      (nreverse pairs))))

;;; ============================================================================
;;; Main test
;;; ============================================================================

(core:set-simd-width 1)

(defun test-lksolvation ()
  (let* ((sigma-ij (+ *lk-test-lj-radius* *lk-test-lj-radius*))
         (r0-low (- sigma-ij *lk-c0*))
         (r1-low (+ sigma-ij *lk-c1*))
         (r2-high *lk-r-solv-low*)
         (r3-high *lk-r-solv-high*)
         ;; Test distances covering all 5 regions
         (distances (list (* 0.5d0 r0-low) ; well inside contact region
                          (* 0.99d0 r0-low) ; just inside contact region
                          (* 0.5d0 (+ r0-low r1-low)) ; middle of low spline
                          (* 0.5d0 (+ r1-low r2-high)) ; middle of Gaussian region
                          3.5d0         ; another Gaussian point
                          (* 0.5d0 (+ r2-high r3-high)) ; middle of high spline
                          5.8d0             ; near end of high spline
                          (+ r3-high 0.5d0))) ; beyond cutoff (zero)
         (tol 1.0d-8)
         (all-ok t))
    (lksolvation--ensure-force-field)
    (format t "LK solvation test (Lazaridis-Karplus 1999 / Alford 2017)~%")
    (format t "  sigma_ij=~,4f r0_low=~,4f r1_low=~,4f r2_high=~,4f r3_high=~,4f~%"
            sigma-ij r0-low r1-low r2-high r3-high)
    (multiple-value-bind (agg a1 a2)
        (lksolvation--make-two-atom-aggregate *lk-test-type*)
      (let* ((setup (list :rosetta (list 'chem:energy-rosetta-lksolvation)))
             (energy-function (chem:make-energy-function :matter agg
                                                         :use-excluded-atoms nil
                                                         :assign-types nil
                                                         :setup setup))
             (component (lksolvation--component energy-function))
             (pos (chem:make-nvector (chem:get-nvector-size energy-function))))
        (loop for dist in distances
              for region in '("contact" "contact-edge" "low-spline" "gaussian"
                              "gaussian" "high-spline" "high-spline-edge" "zero")
              do (chem:set-position a1 (geom:vec 0.0d0 0.0d0 0.0d0))
                 (chem:set-position a2 (geom:vec dist 0.0d0 0.0d0))
                 (chem:load-coordinates-into-vector energy-function pos)
                 (let* ((actual (chem:energy-component-evaluate-energy
                                 energy-function component pos))
                        (expected (lksolvation--expected-energy
                                   dist
                                   *lk-test-lj-radius* *lk-test-lk-dgfree*
                                   *lk-test-lk-lambda* *lk-test-lk-volume*
                                   *lk-test-lj-radius* *lk-test-lk-dgfree*
                                   *lk-test-lk-lambda* *lk-test-lk-volume*))
                        (delta (abs (- actual expected))))
                   (format t "  ~12a dist=~,4f actual=~,12e expected=~,12e delta=~,4e ~a~%"
                           region dist actual expected delta
                           (if (< delta tol) "OK" "FAIL"))
                   (when (>= delta tol)
                     (setf all-ok nil))))))
    #+tests
    (test-true lksolvation-energy-scan all-ok)))

(test-lksolvation)
