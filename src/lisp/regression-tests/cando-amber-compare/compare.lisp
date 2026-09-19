;;;; Interactive comparison: loading this file does not run the calculation.
(eval-when (:compile-toplevel :load-toplevel :execute)
  (ql:quickload '(:leap :spiros/spiros)))

(defpackage #:cando-amber-compare
  (:use #:cl)
  (:export #:compare-ligand #:compare-receptor #:compare-complex
           #:dump-dihedral-terms #:*data-directory*))
(in-package #:cando-amber-compare)

(defparameter *data-directory*
  (merge-pathnames "data/"
                  #.(uiop:pathname-directory-pathname
                     (or *compile-file-truename* *load-truename*
                         (error "Load or compile this file from disk")))))

(defun dump-dihedral-terms (component &key (stream *standard-output*))
  "Dump each Fourier term of an ENERGY-DIHEDRAL as a readable plist.
Return the number of terms. Atom indices are one-based; offsets are
zero-based. Canonicalize proper terms by reversing all four atoms
when the first index exceeds the last. Sort lexicographically by the four
atom indices, then periodicity. Preserve all duplicates; omit term IDs.
Round the dumped V to three decimal places without modifying the component.
Include the actual cached sine/cosine used by the evaluator as well as the
exported phase and its sine/cosine. This does not evaluate or modify terms."
  (check-type component chem:energy-dihedral)
  (let* ((vectors (chem:extract-vectors-as-alist component))
         (phases (cdr (assoc :phase vectors)))
         (proper (cdr (assoc :proper vectors)))
         (atoms (mapcar (lambda (key) (cdr (assoc key vectors)))
                        '(:atom1 :atom2 :atom3 :atom4)))
         (rows nil)
         (*print-readably* t)
         (*print-pretty* nil)
         (*print-length* nil)
         (*print-level* nil)
         (*print-base* 10)
         (*print-radix* nil)
         (*read-default-float-format* 'single-float))
    (format stream "; Dihedral Fourier terms: ~d; V in kcal/mol; phase in radians.~%"
            (length phases))
    (dotimes (i (length phases))
      ;; This accessor returns cos THEN sin (unlike the C++ constructor).
      (multiple-value-bind (cached-cos cached-sin v n i1 i2 i3 i4)
          (chem::safe-amber-energy-dihedral-term component i)
        (let* ((offsets (list i1 i2 i3 i4))
               (names (mapcar (lambda (vec) (chem:get-name (aref vec i))) atoms))
               (phase (aref phases i))
               (expected-cos (cos phase))
               (expected-sin (sin phase)))
          (dolist (offset offsets)
            (unless (and (integerp offset) (>= offset 0) (zerop (mod offset 3)))
              (error "Invalid dihedral coordinate offset ~s in term ~d" offset i)))
          (when (and (aref proper i) (> i1 i4))
            (setf offsets (reverse offsets)
                  names (reverse names)))
          (push (list :proper (aref proper i)
                       :atoms (mapcar (lambda (offset) (1+ (/ offset 3))) offsets)
                       :names names
                       :offsets offsets :v (/ (round (* v 1000d0)) 1000d0)
                       :periodicity n :phase phase
                       :cached-cos (/ (round (* cached-cos 1000d0)) 1000d0)
                       :cached-sin (/ (round (* cached-sin 1000d0)) 1000d0)
                       :phase-cos (/ (round (* expected-cos 1000d0)) 1000d0)
                       :phase-sin (/ (round (* expected-sin 1000d0)) 1000d0)
                       :cos-delta (/ (round (* (- cached-cos expected-cos) 1000d0)) 1000d0)
                       :sin-delta (/ (round (* (- cached-sin expected-sin) 1000d0)) 1000d0))
                rows))))
    (setf rows
          (stable-sort
           (nreverse rows)
           (lambda (left right)
             (loop for a in (getf left :atoms)
                   for b in (getf right :atoms)
                   when (/= a b) do (return (< a b))
                   finally (return (< (getf left :periodicity)
                                      (getf right :periodicity)))))))
    (dolist (row rows)
      (write row :stream stream)
      (terpri stream))
    (finish-output stream)
    (length phases)))

(defun finite-energy (value)
  (unless (and (realp value)
               (<= (- most-positive-double-float) value most-positive-double-float))
    (error "Nonfinite energy: ~s" value))
  (coerce value 'double-float))

(defun component-sum (components name)
  (loop for (key . value) in components when (eq key name)
        sum (finite-energy value) into total
        finally (return (coerce total 'double-float))))

(defun evaluate-components (energy-function coordinates)
  "Gas-phase MM only; no restraints, GB or SA. Return labelled energies and raw entries."
  (unless (chem:nonbond-uses-excluded-atoms-p
           (chem:get-nonbond-component energy-function))
    (error "Comparison requires excluded-atoms nonbond evaluation"))
  (let ((full (chem:make-energy-components))
        (elec (chem:make-energy-components))
        (scale (chem:make-energy-scale)))
    (chem:set-vdw-scale scale 1.0d0)
    (chem:set-electrostatic-scale scale 1.0d0)
    (chem:set-dielectric-constant scale 1.0d0)
    (let* ((total (finite-energy
                   (chem:evaluate-energy energy-function coordinates
                     :energy-components full :energy-scale scale :disable-restraints t)))
           (raw (chem:energy-components/components full)))
      (chem:set-vdw-scale scale 0.0d0)
      (finite-energy
       (chem:evaluate-energy energy-function coordinates
         :energy-components elec :energy-scale scale :disable-restraints t))
      (let* ((electrostatic (chem:energy-components/components elec))
             (bond (component-sum raw 'chem:energy-stretch))
             (angle (+ (component-sum raw 'chem:energy-angle)
                       (component-sum raw 'chem:energy-linear-angle)))
             (dihed (component-sum raw 'chem:energy-dihedral))
             (nb (component-sum raw 'chem::energy-nonbond-total))
             (nb14 (component-sum raw 'chem::energy-nonbond14))
             (eel (component-sum electrostatic 'chem::energy-nonbond-total))
             (eel14 (component-sum electrostatic 'chem::energy-nonbond14)))
        (dolist (label '(chem::energy-nonbond-total chem::energy-nonbond14))
          (unless (and (assoc label raw) (assoc label electrostatic))
            (error "Missing energy component ~s" label)))
        (values (list (cons :bond bond) (cons :angle angle) (cons :dihed dihed)
                      (cons :vdwaals (- nb eel)) (cons :eel eel)
                      (cons :vdw14 (- nb14 eel14)) (cons :eel14 eel14)
                      (cons :other (- total (+ bond angle dihed nb nb14)))
                      (cons :total total))
                raw)))))

(defun matching-coordinates (source target)
  "Copy the aggregate's saved frame into both EF orders, rejecting order mismatches.
This fixture has short, uncompressed atom names; do not guess a mapping."
  (let* ((a (chem:atom-table source)) (b (chem:atom-table target))
         (n (chem:get-number-of-atoms a))
         (xyz (chem:make-nvector (chem:get-nvector-size source)))
         (other (chem:make-nvector (chem:get-nvector-size target))))
    (unless (and (plusp n) (= n (chem:get-number-of-atoms b))
                 (= (* 3 n) (chem:get-nvector-size source)
                    (chem:get-nvector-size target)))
      (error "Saved aggregate and PRMTOP atom/coordinate counts disagree"))
    ;; Constructed tables use zero-based boundaries (including an end sentinel);
    ;; the PRMTOP reader retains one-based RESIDUE_POINTER entries.
    (unless (equal (remove-duplicates
                    (loop for x across (chem:atom-table-residue-pointers a)
                          when (< x n) collect x))
                   (loop for x across (chem:atom-table-residue-pointers b)
                         when (<= x n) collect (1- x)))
      (error "Saved aggregate and PRMTOP residue boundaries disagree"))
    (chem:load-coordinates-into-vector source xyz)
    (let ((seen-a (make-hash-table)) (seen-b (make-hash-table)))
      (dotimes (i n)
        (unless (and (string= (string (chem:elt-atom-name a i))
                             (string (chem:elt-atom-name b i)))
                     (= (chem:elt-atomic-number a i) (chem:elt-atomic-number b i)))
          (error "Atom order mismatch at atom ~d: ~s / ~s"
                 (1+ i) (chem:elt-atom-name a i) (chem:elt-atom-name b i)))
        (let ((ia (chem:elt-atom-coordinate-index-times3 a i))
              (ib (chem:elt-atom-coordinate-index-times3 b i)))
          (loop for offset in (list ia ib) for seen in (list seen-a seen-b)
                do (unless (and (integerp offset) (<= 0 offset (- (* 3 n) 3))
                                (zerop (mod offset 3)) (not (gethash offset seen)))
                     (error "Duplicate or invalid coordinate offset: ~s" offset))
                   (setf (gethash offset seen) t))
          (let ((position (geom:vec-array xyz ia)))
            (mapc #'finite-energy (list (geom:get-x position) (geom:get-y position)
                                       (geom:get-z position)))
            (geom:vec-set position other ib)))))
    (values xyz other)))

(defun compare-species (species &key (data-directory *data-directory*)
                           (absolute-tolerance 1.0d-4) (relative-tolerance 1.0d-6)
                           dihedral-dump-directory
                           (error-on-mismatch t) (stream *standard-output*))
  "Compare a rebuilt species EF with the EF read from its PRMTOP.
Return rows, rebuilt EF, topology EF, and both raw component alists.
Each row is (:component ... :cando ... :prmtop ... :delta ... :pass ...).
All rows print before a mismatch signals; use :ERROR-ON-MISMATCH NIL to inspect.
When DIHEDRAL-DUMP-DIRECTORY is non-NIL, create that directory and write
<species>-rebuilt-dihedrals.lisp and <species>-prmtop-dihedrals.lisp there before
evaluating energies. Existing dump files are overwritten.
This compares two Cando evaluations, not an invocation of Amber/Sander."
  (unless (and (realp absolute-tolerance) (>= absolute-tolerance 0)
               (realp relative-tolerance) (>= relative-tolerance 0))
    (error "Comparison tolerances must be nonnegative real numbers"))
  (let* ((name (ecase species
                 (:ligand "ligand") (:receptor "receptor") (:complex "complex")))
         (directory (uiop:ensure-directory-pathname data-directory))
         (aggregate (cando.serialize:load-cando
                     (merge-pathnames (format nil "~a.cando" name) directory))))
    (check-type aggregate chem:aggregate)
    ;; Disable construction, not just evaluation: the chiral setup checks neighbors.
    (let* ((rebuilt (chem:make-energy-function :matter aggregate
                                              :spec '(:amber :default t
                                                      :remove (chem:energy-chiral-restraint))
                                              :use-excluded-atoms t))
           (loaded (leap.topology::read-amber-parm-format
                    (merge-pathnames (format nil "~a.prmtop" name) directory))))
      (when (or (chem:bounding-box-bound-p (chem:atom-table rebuilt))
                (chem:bounding-box-bound-p (chem:atom-table loaded)))
        (error "This comparison supports only nonperiodic fixtures"))
      (multiple-value-bind (xyz other) (matching-coordinates rebuilt loaded)
        (when dihedral-dump-directory
          (let ((dump-directory (uiop:ensure-directory-pathname dihedral-dump-directory)))
            (loop for energy-function in (list rebuilt loaded)
                  for source in '("rebuilt" "prmtop")
                  for filename = (format nil "~a-~a-dihedrals.lisp" name source)
                  for pathname = (merge-pathnames filename dump-directory)
                  do (ensure-directories-exist pathname)
                     (with-open-file (out pathname :direction :output
                                                   :if-exists :supersede
                                                   :if-does-not-exist :create)
                       (dump-dihedral-terms (chem:get-dihedral-component energy-function)
                                            :stream out))
                     (format stream "~&Wrote dihedral dump: ~a~%" pathname))))
        (multiple-value-bind (left raw-left) (evaluate-components rebuilt xyz)
          (multiple-value-bind (right raw-right) (evaluate-components loaded other)
            (let ((rows nil))
              (format stream "~&~a: same saved coordinates; kcal/mol; no restraints or GB/SA.~%"
                      (string-capitalize name))
              (format stream "~12a ~18a ~18a ~18a ~a~%"
                      "Component" "Cando rebuilt" "PRMTOP read" "Read - rebuilt" "Match")
              (dolist (entry left)
                (let* ((a (cdr entry)) (b (cdr (assoc (car entry) right)))
                       (delta (- b a))
                       (pass (<= (abs delta) (+ absolute-tolerance
                                               (* relative-tolerance (max (abs a) (abs b)))))))
                  (push (list :component (car entry) :cando a :prmtop b :delta delta :pass pass) rows)
                  (format stream "~12a ~18,8f ~18,8f ~18,8f ~a~%"
                          (car entry) a b delta (if pass "PASS" "FAIL"))))
              (setf rows (nreverse rows))
              (finish-output stream)
              #+(or)
              (when (and error-on-mismatch (some (lambda (row) (not (getf row :pass))) rows))
                (error "~a Cando/PRMTOP energy mismatch: ~{~a~^, ~}" name
                       (loop for row in rows unless (getf row :pass) collect (getf row :component))))
              (values rows rebuilt loaded raw-left raw-right))))))))

(defun compare-ligand (&rest options)
  "Compare ligand.cando with ligand.prmtop. Accepts the keywords of COMPARE-SPECIES:
:DATA-DIRECTORY, :ABSOLUTE-TOLERANCE, :RELATIVE-TOLERANCE,
:DIHEDRAL-DUMP-DIRECTORY, :ERROR-ON-MISMATCH, and :STREAM.
Return rows, rebuilt EF, imported EF, and both raw component alists."
  (apply #'compare-species :ligand options))

(defun compare-receptor (&rest options)
  "Compare receptor.cando with receptor.prmtop, using the same options and
five return values as COMPARE-LIGAND. Preserve saved molecule force fields."
  (apply #'compare-species :receptor options))

(defun compare-complex (&rest options)
  "Compare complex.cando with complex.prmtop, using the same options and
five return values as COMPARE-LIGAND. Preserve each molecule's force field."
  (apply #'compare-species :complex options))
