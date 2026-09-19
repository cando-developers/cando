(in-package #:clasp-tests)

;;; Public-constructor regressions for :USE-EXCLUDED-ATOMS.
;;; No external force-field files or other regression helpers are required.

(defun aex--ref (vector index)
  (let ((xyz (geom:vec-array vector (* 3 (floor index 3)))))
    (ecase (mod index 3)
      (0 (geom:get-x xyz)) (1 (geom:get-y xyz)) (2 (geom:get-z xyz)))))

(defun aex--vector (values)
  (let ((vector (chem:make-nvector (length values))))
    (loop for i from 0 below (length values) by 3
          do (geom:vec-set (geom:vec (elt values i) (elt values (+ i 1))
                                   (elt values (+ i 2))) vector i))
    vector))

(defun aex--close (actual expected &optional (tolerance 2.0d-7))
  (unless (<= (abs (- actual expected))
              (* tolerance (max 1.0d0 (abs expected))))
    (error "Amber exclusion regression: expected ~s, got ~s" expected actual))
  t)

(defun aex--scale (vdw electrostatic dielectric)
  (let ((scale (chem:make-energy-scale)))
    (chem:set-vdw-scale scale vdw)
    (chem:set-electrostatic-scale scale electrostatic)
    (chem:set-dielectric-constant scale dielectric)
    scale))

(defun aex--fixture (kind callback &key (excluded t))
  "Run CALLBACK with EF, coordinates, nonbond component and aggregate.
PAIR has no bonds; CHAIN has four atoms; RING4/RING5 are closed rings.
Only nonbond energies are enabled, so synthetic bonded parameters are unneeded."
  (let* ((leap.core:*force-fields* (make-hash-table))
         (ff (chem:force-field/make))
         (nb (chem:make-ffnonbond :aex-carbon))
         (agg (chem:make-aggregate))
         (mol (chem:make-molecule :aex))
         (res (chem:make-residue :aex))
         (n (ecase kind (:pair 2) (:chain 4) (:ring4 4) (:ring5 5)))
         (atoms (make-array n)))
    (chem:ffnonbond/set-radius-angstroms nb 1.5d0)
    (chem:ffnonbond/set-epsilon-kcal nb 0.2d0)
    (chem:ffnonbond/set-mass nb 12.0d0)
    (chem:ffnonbond-db-add (chem:get-nonbond-db ff) nb)
    (leap.core:add-force-field-or-modification
     ff :force-field-name :aex :combined-force-field-class-name 'chem:combined-force-field)
    (chem:add-matter mol res)
    (chem:add-molecule agg mol)
    (chem:setf-force-field-name agg :aex)
    (chem:setf-force-field-name mol :aex)
    (dotimes (i n)
      (let ((atom (chem:make-atom (intern (format nil "AEX~d" i) :keyword) :c)))
        (setf (aref atoms i) atom)
        (chem:add-atom res atom)
        (chem:set-property atom :given-atom-type :aex-carbon)
        (chem:set-charge atom (cond ((zerop i) 0.3d0)
                                   ((= i (1- n)) -0.4d0) (t 0.0d0)))
        ;; Oblique, noncoincident positions exercise all Cartesian derivatives.
        (chem:set-position atom
                           (geom:vec (* 1.4d0 i) (* 0.2d0 i i) (* -0.3d0 i)))))
    (unless (eq kind :pair)
      (dotimes (i (1- n))
        (chem:bond-to (aref atoms i) (aref atoms (1+ i)) :single-bond))
      (when (member kind '(:ring4 :ring5))
        (chem:bond-to (aref atoms (1- n)) (aref atoms 0) :single-bond)))
    (when (eq kind :pair)
      (chem:set-position (aref atoms 1) (geom:vec 3.2d0 -1.6d0 3.2d0)))
    (let* ((ef (chem:make-energy-function
                :matter agg :assign-types nil :use-excluded-atoms excluded
                :spec '(:amber :default nil (chem:energy-nonbond t))))
           (position (chem:make-nvector (chem:get-nvector-size ef))))
      (chem:load-coordinates-into-vector ef position)
      (funcall callback ef position (chem:get-nonbond-component ef) agg))))

(defun aex--reference (position kind vdw electrostatic dielectric direction)
  "Independent radial LJ+Coulomb energy, force and H*direction.
The chain has only pair (0,3), scaled by 1/2 and 1/1.2; PAIR is unscaled."
  (let* ((n (if (eq kind :pair) 2 4))
         (last (* 3 (1- n)))
         (delta (loop for k below 3 collect (- (aex--ref position k)
                                              (aex--ref position (+ last k)))))
         (r (sqrt (reduce #'+ delta :key (lambda (x) (* x x)))))
         (lj-scale (* vdw (if (eq kind :chain) 0.5d0 1.0d0)))
         (a (* lj-scale 0.2d0 (expt 3.0d0 12)))
         (c (* lj-scale 0.4d0 (expt 3.0d0 6)))
         (q (* (/ electrostatic dielectric) (expt 18.2223d0 2)
               0.3d0 -0.4d0 (if (eq kind :chain) (/ 1.0d0 1.2d0) 1.0d0)))
         (energy (+ (/ a (expt r 12)) (- (/ c (expt r 6))) (/ q r)))
         (first (+ (/ (* -12 a) (expt r 13)) (/ (* 6 c) (expt r 7))
                   (- (/ q (* r r)))))
         (second (+ (/ (* 156 a) (expt r 14)) (/ (* -42 c) (expt r 8))
                    (/ (* 2 q) (expt r 3))))
         (force (make-array (* 3 n) :initial-element 0.0d0))
         (hd (make-array (* 3 n) :initial-element 0.0d0)))
    (dotimes (k 3)
      (setf (aref force k) (* (- first) (/ (elt delta k) r))
            (aref force (+ last k)) (- (aref force k)))
      (dotimes (j 3)
        (incf (aref hd k)
              (* (+ (* (- second (/ first r)) (/ (* (elt delta k) (elt delta j)) (* r r)))
                    (if (= k j) (/ first r) 0.0d0))
                 (- (aref direction j) (aref direction (+ last j))))))
      (setf (aref hd (+ last k)) (- (aref hd k))))
    (values energy force hd)))

(defun aex--check-potential (kind scales)
  (aex--fixture
   kind
   (lambda (ef position nb agg)
     (declare (ignore agg))
     (let* ((size (chem:get-nvector-size ef))
            (direction (make-array size :element-type 'double-float)))
       (dotimes (i size)
         (setf (aref direction i) (* (if (evenp i) 1.0d0 -1.0d0) (+ 0.1d0 (* 0.03d0 i)))))
       (assert (= (chem:number-of-terms14 nb) (if (eq kind :chain) 1 0)))
       (dolist (settings scales)
         (destructuring-bind (vdw electrostatic dielectric) settings
           (let ((scale (aex--scale vdw electrostatic dielectric))
                 (force (chem:make-nvector size))
                 (hd (chem:make-nvector size)))
             (multiple-value-bind (expected expected-force expected-hd)
                 (aex--reference position kind vdw electrostatic dielectric direction)
               (aex--close (chem:evaluate-energy ef position :energy-scale scale) expected)
               (aex--close (chem:evaluate-energy-force ef position :energy-scale scale
                                                       :calc-force t :force force) expected)
               (dotimes (i size) (aex--close (aex--ref force i) (aref expected-force i)))
               (aex--close
                (chem:scoring-function/evaluate-all
                 ef position :energy-scale scale :calc-force t :force force
                 :calc-diagonal-hessian t :calc-off-diagonal-hessian t
                 :hdvec hd :dvec (aex--vector direction)) expected)
               (dotimes (i size)
                 (aex--close (aex--ref force i) (aref expected-force i))
                 (aex--close (aex--ref hd i) (aref expected-hd i)))))))
       ;; Repeated evaluation and substantial motion must not switch modes or
       ;; build a spatial pair list. The reference follows the changed geometry.
       (geom:vec-set (geom:vec 15.0d0 -3.0d0 2.0d0) position (- size 3))
       (aex--close (chem:evaluate-energy ef position)
                   (aex--reference position kind 1.0d0 1.0d0 1.0d0 direction))
       (assert (chem:nonbond-uses-excluded-atoms-p nb))
       (assert (zerop (chem:nonbond-pair-list-epoch nb)))
       (assert (zerop (chem:number-of-terms nb)))
       t))))

(defun aex--check-ring (kind)
  (aex--fixture kind
                (lambda (ef position nb agg)
                  (declare (ignore agg))
                  (assert (zerop (chem:number-of-terms14 nb)))
                  (aex--close (chem:evaluate-energy ef position) 0.0d0)
                  (assert (chem:nonbond-uses-excluded-atoms-p nb))
                  (assert (zerop (chem:nonbond-pair-list-epoch nb)))
                  t)))

(defun aex--error-containing (thunk text)
  (handler-case (progn (funcall thunk) nil)
    (error (condition) (not (null (search text (princ-to-string condition)))))))

(test-true amber-excluded-component-factory
           (aex--fixture
            :chain
            (lambda (ef position nb agg)
              (declare (ignore position nb agg))
              (let ((table (chem:atom-table ef))
                    (calls 0))
                (multiple-value-bind (expected-counts expected-indices)
                    (chem:calculate-excluded-atom-list table t)
                  (multiple-value-bind (counts indices)
                      (chem:calculate-excluded-atom-list
                       table (lambda (component-class)
                               (incf calls)
                               (assert (eq component-class (find-class 'chem:energy-nonbond)))
                               t))
                    (and (= calls 1) (equalp counts expected-counts)
                         (equalp indices expected-indices))))))
            :excluded nil))

(test-true amber-excluded-reject-filtering
           (aex--fixture
            :chain
            (lambda (ef position nb agg)
              (declare (ignore position nb agg))
              (loop for factory in
                    (list nil
                          (lambda (component-class) (declare (ignore component-class)) nil)
                          (lambda (component-class)
                            (declare (ignore component-class))
                            (lambda (&rest pair) (declare (ignore pair)) t)))
                    always (aex--error-containing
                            (lambda ()
                              (chem:calculate-excluded-atom-list (chem:atom-table ef) factory))
                            "require the EnergyNonbond interaction factory to return T")))
            :excluded nil))

(test-true amber-excluded-chain
           (aex--check-potential :chain '((1.0d0 1.0d0 1.0d0))))
(test-true amber-excluded-unbonded-pair
           (aex--check-potential :pair '((1.0d0 1.0d0 1.0d0))))
(test-true amber-excluded-ring4 (aex--check-ring :ring4))
(test-true amber-excluded-ring5 (aex--check-ring :ring5))
(test-true amber-excluded-runtime-scales
           (loop for kind in '(:pair :chain)
                 always (aex--check-potential
                         kind '((0.0d0 1.0d0 1.0d0) (1.0d0 0.0d0 1.0d0)
                                (0.3d0 0.7d0 4.0d0) (0.0d0 0.0d0 2.0d0)))))

(test-true amber-excluded-invalid-dielectric
           (loop for kind in '(:pair :chain)
                 always
                 (aex--fixture
                  kind
                  (lambda (ef position nb agg)
                    (declare (ignore nb agg))
                    (loop for dielectric in (list 0.0d0 -1.0d0
                                                  ext:double-float-positive-infinity
                                                  (ext:bits-to-double-float #x7ff8000000000000))
                          always
                          (aex--error-containing
                           (lambda ()
                             (chem:evaluate-energy ef position
                                                   :energy-scale (aex--scale 1.0d0 1.0d0 dielectric)))
                           "dielectric must be finite and positive"))))))

(test-true amber-nonbond14-missing-parameters
           ;; Build the fixture with the unguarded pair-list constructor, then
           ;; call the same 1-4 builder with an empty database. This specifically
           ;; tests its false-return handling, not an earlier atom-typing error.
           (aex--fixture
            :chain
            (lambda (ef position nb agg)
              (declare (ignore position))
              (let ((before (chem:number-of-terms14 nb)))
                (and
                 (aex--error-containing
                  (lambda ()
                    (chem:construct-14-interaction-terms
                     nb (chem:atom-table ef) agg
                     (chem:get-nonbond-db (chem:force-field/make)) t (chem:atom-types ef)))
                  "Missing nonbond parameters for 1-4 pair")
                 (= before (chem:number-of-terms14 nb)))))
            :excluded nil))
