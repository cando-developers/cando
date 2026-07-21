#+tests(in-package #:clasp-tests)

;;;; Regression: cached-smirnoff-force-field HIT energy == MISS energy
;;;;
;;;; Parameterize the SAME molecule twice through one cached-smirnoff-force-field:
;;;;   PASS 1 - cache empty -> signature uncovered -> MISS: delegates to real
;;;;            SMIRNOFF and harvests bonded + nonbond params into the cache.
;;;;   PASS 2 - same cache instance, signature now covered -> HIT: applies bonded
;;;;            terms + memoized vdw typing purely from the cache.
;;;; The total energy and every per-component energy must match to tolerance.
;;;;
;;;; The dihedral component holds BOTH proper and improper terms (split by each
;;;; term's PROPER flag), so a per-component total would hide a proper/improper
;;;; divergence in the other.  We therefore also split the dihedral energy into
;;;; proper vs improper by hand-summing the terms (impropers exercise the new
;;;; map-impropers apply path in commit-bonded-from-cache).  Both passes are
;;;; summed the identical way, so the MISS-vs-HIT comparison is valid regardless
;;;; of any kernel-formula subtlety - only self-consistency between passes matters.
;;;;
;;;; To PROVE pass 2 is a real HIT (and not a second miss that would also match),
;;;; we assert the coverage predicate is true after pass 1.  That predicate is the
;;;; exact condition define-for-molecule-using-force-field branches on, so coverage
;;;; guarantees pass 2 takes the HIT branch.  (We can't null out the wrapped real
;;;; SMIRNOFF to force it: the HIT path still needs it for the nonbond merge, whose
;;;; force-fields-as-list delegates to that slot.)

;; real SMIRNOFF, registered under :smirnoff
(leap:load-smirnoff-params
 (probe-file "sys:extensions;cando;src;lisp;regression-tests;data;force-field.offxml"))
(core:set-simd-width 1)

;; a molecule with residues; give every atom a synthetic
;; :given-atom-type = (atom-name . constitution-context).  Using the residue name
;; as the context means all atoms in a residue share one context, which is exactly
;; what molecule-signature keys on.
(defparameter *cc-agg*
  (cando:load-mol2 "sys:extensions;cando;src;lisp;regression-tests;data;hexapeptide-smirnoff.mol2"))
(defparameter *cc-mol* (cando:mol *cc-agg* 0))

(chem:map-residues
 nil
 (lambda (residue)
   (let ((context (chem:get-name residue)))
     (chem:map-atoms
      nil
      (lambda (atom)
        (chem:set-property atom :given-atom-type (cons (chem:get-name atom) context)))
      residue)))
 *cc-mol*)

;; wrap real :smirnoff in a cached FF, register it, point the molecule at it
(defparameter *cc-ff*
  (make-instance 'smirnoff:cached-smirnoff-force-field
                 :smirnoff (chem:find-force-field :smirnoff)
                 :smirnoff-name :smirnoff))
(leap.core:add-combined-force-field *cc-ff* :smirnoff-cache-test)
(chem:setf-force-field-name *cc-mol* :smirnoff-cache-test)

;; Hand-sum the dihedral component's terms, bucketed by the PROPER flag.
;; Per-term kernel energy is  v * (1 + cos(n*phi - phase))  (_Dihedral_termCode.cc).
;; Returns (values proper-energy improper-energy).
(defun dihedral-bucket-energies (ef pos)
  (let ((proper 0.0) (improper 0.0)
        (dih (chem:energy-function/get-dihedral-component ef)))
    (when dih
      (chem:walk-dihedral-terms
       dih
       (lambda (index a1 a2 a3 a4 i3x1 i3x2 i3x3 i3x4 v n phase proper-p)
         (declare (ignore index a1 a2 a3 a4))
         (let* ((phi (geom:calculate-dihedral (geom:vec-array pos i3x1)
                                              (geom:vec-array pos i3x2)
                                              (geom:vec-array pos i3x3)
                                              (geom:vec-array pos i3x4)))
                (e   (* v (+ 1.0 (cos (- (* n phi) phase))))))
           (if proper-p (incf proper e) (incf improper e))))))
    (values proper improper)))

;; Parameterize + evaluate.  Returns four values:
;;   total-energy, component-ht (class-name -> energy), proper-dih, improper-dih
(defun cc-energy ()
  (let* ((ef  (chem:make-energy-function :matter *cc-agg*))
         (pos (chem:make-nvector (chem:get-nvector-size ef)))
         (ht  (make-hash-table)))
    (chem:load-coordinates-into-vector ef pos)
    (let ((total (chem:evaluate-energy ef pos)))
      (loop for comp in (chem:all-components ef)
            for e = (chem:energy-component-evaluate-energy ef comp pos)
            do (setf (gethash (class-name (class-of comp)) ht) e))
      (multiple-value-bind (proper improper) (dihedral-bucket-energies ef pos)
        (values total ht proper improper)))))

;; --- PASS 1: MISS + harvest ---
(multiple-value-bind (total ht proper improper) (cc-energy)
  (defparameter *energy-miss* total)
  (defparameter *components-miss* ht)
  (defparameter *proper-miss* proper)
  (defparameter *improper-miss* improper))

;; the molecule is now covered -> pass 2 is guaranteed to take the HIT branch
(defparameter *covers-after-miss*
  (and (smirnoff::bonded-cache-covers-molecule-p *cc-mol* (smirnoff::cache *cc-ff*)) t))

;; --- PASS 2: HIT (cache only) ---
(multiple-value-bind (total ht proper improper) (cc-energy)
  (defparameter *energy-hit* total)
  (defparameter *components-hit* ht)
  (defparameter *proper-hit* proper)
  (defparameter *improper-hit* improper))

(defparameter *cc-delta* (abs (- *energy-miss* *energy-hit*)))
(format t "cache MISS energy = ~f~%cache HIT  energy = ~f~%total delta = ~g~%"
        *energy-miss* *energy-hit* *cc-delta*)

;; --- per-component MISS-vs-HIT breakdown ------------------------------------
;; Union of component class-names seen in either pass; a term that appears in one
;; pass but not the other (missing/extra) is itself a failure via the |diff| test.
(defparameter *cc-component-names*
  (let ((names '()))
    (maphash (lambda (k v) (declare (ignore v)) (pushnew k names)) *components-miss*)
    (maphash (lambda (k v) (declare (ignore v)) (pushnew k names)) *components-hit*)
    names))

(defparameter *cc-worst-component-delta* 0.0)
(dolist (name *cc-component-names*)
  (let* ((em (gethash name *components-miss* 0.0))
         (eh (gethash name *components-hit* 0.0))
         (d  (abs (- em eh))))
    (setf *cc-worst-component-delta* (max *cc-worst-component-delta* d))
    (format t "~30s  miss=~14,6f  hit=~14,6f  |diff|=~g~%" name em eh d)))

;; --- proper vs improper dihedral split --------------------------------------
(defparameter *proper-delta*   (abs (- *proper-miss*   *proper-hit*)))
(defparameter *improper-delta* (abs (- *improper-miss* *improper-hit*)))
(format t "dihedral PROPER    miss=~14,6f  hit=~14,6f  |diff|=~g~%"
        *proper-miss* *proper-hit* *proper-delta*)
(format t "dihedral IMPROPER  miss=~14,6f  hit=~14,6f  |diff|=~g~%"
        *improper-miss* *improper-hit* *improper-delta*)
;; self-consistency: the two buckets should sum to the dihedral component total
;; (informational - loose tol absorbs the kernel's linear-dihedral guard/clamping)
(let ((split-vs-comp
        (abs (- (+ *proper-miss* *improper-miss*)
                (gethash 'chem:energy-dihedral *components-miss* 0.0)))))
  (format t "dihedral split-vs-component (miss) |diff|=~g~%" split-vs-comp))

;; --- assertions -------------------------------------------------------------
;; pass 1 really harvested (0 signatures -> at least 1), proving it missed
;; ...and the molecule is now covered, so pass 2 provably takes the HIT branch
#+tests(test-true smirnoff-cache-covered-after-miss *covers-after-miss*)
;; the cache-only HIT reproduces the real-SMIRNOFF total energy
#+tests(test-true smirnoff-cache-hit-equals-miss (< *cc-delta* 1.0e-6))
;; ...and reproduces it component by component (localizes any divergence)
#+tests(test-true smirnoff-cache-hit-components (< *cc-worst-component-delta* 1.0e-6))
;; ...and the proper and improper dihedral buckets each match independently,
;; so an improper-only regression (the map-impropers apply path) can't hide
;; inside the combined dihedral component total
#+tests(test-true smirnoff-cache-hit-proper-dihedral   (< *proper-delta*   1.0e-6))
#+tests(test-true smirnoff-cache-hit-improper-dihedral (< *improper-delta* 1.0e-6))


(ql:quickload :path-design)

(defun time-make-energy-function (aggregate
                                  &key (n 50) (warmup 3)
                                    (variants '((:cached   . :smirnoff-spiros)
                                                (:fallback . :smirnoff))))
  "Benchmark chem:make-energy-function on AGGREGATE under each force-field variant.

  AGGREGATE must already have its atoms stamped with :given-atom-type (as the spiro
  assembler does), so the cached path can key on them.  The same aggregate is reused
  for every variant; only its molecules' force-field-name is switched.

  For each (LABEL . FORCE-FIELD-NAME): set every molecule's force-field-name, run
  WARMUP untimed builds (the first cached build harvests the cache, so the timed runs
  are genuine hits - keep WARMUP >= 1 for the cached variant), then time N builds with
  wall-clock (get-internal-real-time, which includes GC).

  Returns an alist of (label n total-seconds seconds-per-call); also prints a table
  and, when both :cached and :fallback are present, the per-call speedup."
  (flet ((build-once (name)
           (let ((force-field-overrides (make-hash-table)))
             (chem:do-molecules (mol aggregate)
               (setf (gethash mol force-field-overrides) name))
             (chem:make-energy-function :matter aggregate :force-field-overrides force-field-overrides))))
    (let ((results '()))
      (dolist (variant variants)
        (let ((label  (car variant))
              (ffname (cdr variant)))
          (dotimes (i warmup) (build-once ffname)) ; harvests cache on 1st cached build
          (let ((start (get-internal-real-time)))
            (dotimes (i n) (build-once ffname))
            (let* ((secs (/ (float (- (get-internal-real-time) start))
                            internal-time-units-per-second))
                   (per-call (/ secs (max n 1))))
              (format t "~12a  ff=~20a  n=~4d  total=~8,4f s  per-call=~10,6f s  (~,1f/s)~%"
                      label ffname n secs per-call (if (zerop secs) 0.0 (/ n secs)))
              (push (list label n secs per-call) results)))))
      (let ((cached (assoc :cached results)) (fallback (assoc :fallback results)))
        (when (and cached fallback (plusp (fourth cached)))
          (format t "speedup (fallback / cached per-call) = ~,2fx~%"
                  (/ (fourth fallback) (fourth cached)))))
      (nreverse results))))

(time 
 (let* ((oligomer (topology:make-oligomer path-design:*spiro8-space*))
        (my-oligomer-shape (topology:make-oligomer-shape oligomer))
        (assembler (topology:make-assembler (list my-oligomer-shape) :receptor-only t))
        (agg       (topology:aggregate assembler))) ; adjust accessor to your API
   (time-make-energy-function agg :n 100 :warmup 3)
   (let* ((ff    (chem:find-force-field :smirnoff-spiros))
          (cache (smirnoff::cache ff))
          (mol   (cando:mol agg 0)))
     (format t "bond-table entries = ~d~%" (hash-table-count (smirnoff::bond-table cache)))
     (format t "nonbond entries    = ~d~%" (hash-table-count (smirnoff::nonbond-table cache)))
     (format t "covered? ~s~%" (smirnoff::bonded-cache-covers-molecule-p mol cache))) 
   ))


(time 
 (let* ((oligomer (topology:make-oligomer path-design:*spiro8-space*))
        (my-oligomer-shape (topology:make-oligomer-shape oligomer))
        (assembler (topology:make-assembler (list my-oligomer-shape) :receptor-only t))
        (agg       (topology:aggregate assembler))) ; adjust accessor to your API
   (time-make-energy-function agg :n 5000 :warmup 3 :variants '((:cache . :smirnoff-spiros)))
   ))

p
