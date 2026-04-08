#+tests(in-package #:clasp-tests)

;;; ---------------------------------------------------------------------------
;;; Regression tests for chem:parse-energy-spec, chem:copy-energy-function,
;;; and chem:make-energy-function :spec — comparing old API vs new API.
;;; ---------------------------------------------------------------------------

;;; --- Helper: compare per-component energies between two energy functions ---

(defun cef--compare-component-energies (ef1 ef2 coords &key (tolerance 0.0001))
  "Return T if every component in ef1 has a matching component in ef2
with the same energy (within tolerance). Also checks total energy."
  (let ((e1 (chem:evaluate-energy ef1 coords))
        (e2 (chem:evaluate-energy ef2 coords)))
    (unless (< (abs (- e1 e2)) tolerance)
      (format t "  Total energy mismatch: old=~f new=~f delta=~f~%" e1 e2 (- e1 e2))
      (return-from cef--compare-component-energies nil))
    (let ((comps1 (chem:all-components ef1))
          (comps2 (chem:all-components ef2)))
      (unless (= (length comps1) (length comps2))
        (format t "  Component count mismatch: old=~d new=~d~%" (length comps1) (length comps2))
        (return-from cef--compare-component-energies nil))
      (loop for c1 in comps1
            for name = (class-name (class-of c1))
            for ce1 = (chem:energy-component-evaluate-energy ef1 c1 coords)
            for c2 = (chem:find-component-or-nil ef2 name)
            do (if (null c2)
                   (progn (format t "  Missing component ~s in new ef~%" name)
                          (return-from cef--compare-component-energies nil))
                   (let ((ce2 (chem:energy-component-evaluate-energy ef2 c2 coords)))
                     (unless (< (abs (- ce1 ce2)) tolerance)
                       (format t "  Component ~s mismatch: old=~f new=~f~%" name ce1 ce2)
                       (return-from cef--compare-component-energies nil))))))
    t))

;;; =========================================================================
;;; Test 1: make-energy-function — AMBER old way vs :spec :amber
;;; =========================================================================

(let* ((agg (chem:sample-aggregate))
       (coords (chem:matter/extract-coordinates agg))
       ;; Old way
       (old-ef (chem:%make-energy-function :matter agg :use-excluded-atoms nil))
       ;; New way
       (new-ef (chem:make-energy-function :matter agg :use-excluded-atoms nil
                                          :spec :amber)))
  #+tests(test-true make-ef-amber-old-vs-spec
                    (cef--compare-component-energies old-ef new-ef coords)))

;;; =========================================================================
;;; Test 2: make-energy-function — Rosetta old way vs :spec :rosetta
;;;         Uses the rosetta-nonbond test infrastructure.
;;; =========================================================================

(defun cef--ensure-rosetta-force-field ()
  "Set up a minimal rosetta force field for testing."
  (let* ((ff-name :cef-rosetta-test)
         (atom-type :cef-rnb)
         (force-field (chem:force-field/make))
         (nonbond-db (chem:get-nonbond-db force-field))
         (ffnonbond (chem:make-ffnonbond atom-type)))
    (chem:ffnonbond/set-radius-angstroms ffnonbond 1.5d0)
    (chem:ffnonbond/set-epsilon-kcal ffnonbond 0.2d0)
    (chem:ffnonbond-db-add nonbond-db ffnonbond)
    (when (chem:find-force-field ff-name nil)
      (leap:clear-force-field ff-name))
    (leap.core:add-force-field-or-modification
     force-field
     :force-field-name ff-name
     :combined-force-field-class-name 'chem:combined-force-field)
    (values ff-name atom-type)))

(defun cef--make-two-atom-aggregate (ff-name atom-type)
  "Create a two-atom aggregate in separate molecules for rosetta nonbond testing."
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
    (chem:setf-force-field-name agg ff-name)
    (chem:setf-force-field-name mol1 ff-name)
    (chem:setf-force-field-name mol2 ff-name)
    (chem:set-property a1 :given-atom-type atom-type)
    (chem:set-property a2 :given-atom-type atom-type)
    (chem:set-position a1 (geom:vec 0.0d0 0.0d0 0.0d0))
    (chem:set-position a2 (geom:vec 3.0d0 0.0d0 0.0d0))
    (values agg a1 a2)))

(multiple-value-bind (ff-name atom-type)
    (cef--ensure-rosetta-force-field)
  (multiple-value-bind (agg a1 a2)
      (cef--make-two-atom-aggregate ff-name atom-type)
    (declare (ignore a1 a2))
    (let* ((coords (chem:matter/extract-coordinates agg))
           ;; Old way: :setup + :keep-interaction-factory
           (old-ef (chem:%make-energy-function :matter agg
                                               :use-excluded-atoms nil
                                               :assign-types nil
                                               :setup (list :rosetta
                                                            (list 'chem:energy-rosetta-nonbond
                                                                  :rep-weight 0.5))))
           ;; New way: :spec
           (new-ef (chem:make-energy-function :matter agg
                                              :use-excluded-atoms nil
                                              :assign-types nil
                                              :spec '(:rosetta
                                                      (chem:energy-rosetta-nonbond
                                                       :rep-weight 0.5)))))
      #+tests(test-true make-ef-rosetta-old-vs-spec
                        (cef--compare-component-energies old-ef new-ef coords)))))

;;; =========================================================================
;;; Test 3: make-energy-function — Rosetta with :only (component selection)
;;;         Old way uses keep-interaction-factory, new way uses :only.
;;; =========================================================================

(multiple-value-bind (ff-name atom-type)
    (cef--ensure-rosetta-force-field)
  (multiple-value-bind (agg a1 a2)
      (cef--make-two-atom-aggregate ff-name atom-type)
    (declare (ignore a1 a2))
    (let* ((coords (chem:matter/extract-coordinates agg))
           ;; Old way: keep-interaction-factory that only keeps rosetta-nonbond
           (old-ef (chem:%make-energy-function
                    :matter agg
                    :use-excluded-atoms nil
                    :assign-types nil
                    :keep-interaction-factory
                    (lambda (aclass)
                      (if (eq aclass (find-class 'chem:energy-rosetta-nonbond))
                          t
                          nil))
                    :setup '(:rosetta)))
           ;; New way: :spec with :default nil and explicit component
           (new-ef (chem:make-energy-function
                    :matter agg
                    :use-excluded-atoms nil
                    :assign-types nil
                    :spec '(:rosetta
                            :default nil
                            (chem:energy-rosetta-nonbond t))))
           (old-energy (chem:evaluate-energy old-ef coords))
           (new-energy (chem:evaluate-energy new-ef coords)))
      ;; Both should have only rosetta-nonbond component with terms
      #+tests(test-true make-ef-only-same-energy
                        (< (abs (- old-energy new-energy)) 0.0001))
      ;; New ef should not have rosetta-elec or lksolvation
      #+tests(test-true make-ef-only-no-elec
                        (null (chem:find-component-or-nil new-ef 'chem:energy-rosetta-elec)))
      #+tests(test-true make-ef-only-no-lk
                        (null (chem:find-component-or-nil new-ef 'chem:energy-rosetta-lksolvation))))))

;;; =========================================================================
;;; Test 4: copy-energy-function — full copy, old way vs new way
;;; =========================================================================

(let* ((agg (chem:sample-aggregate))
       (coords (chem:matter/extract-coordinates agg))
       (base-ef (chem:%make-energy-function :matter agg :use-excluded-atoms nil))
       ;; Old way: copy-filter with T
       (old-copy (chem:copy-filter base-ef t))
       ;; New way: copy-energy-function with T
       (new-copy (chem:copy-energy-function base-ef t)))
  #+tests(test-true copy-ef-full-old-vs-new
                    (cef--compare-component-energies old-copy new-copy coords)))

;;; =========================================================================
;;; Test 5: copy-energy-function — with setup params (rep-weight),
;;;         old 3-arg copy-filter vs new spec
;;; =========================================================================

(multiple-value-bind (ff-name atom-type)
    (cef--ensure-rosetta-force-field)
  (multiple-value-bind (agg a1 a2)
      (cef--make-two-atom-aggregate ff-name atom-type)
    (declare (ignore a1 a2))
    (let* ((coords (chem:matter/extract-coordinates agg))
           (base-ef (chem:%make-energy-function :matter agg
                                                :use-excluded-atoms nil
                                                :assign-types nil
                                                :setup '(:rosetta)))
           ;; Old way: copy-filter with keep=T and setup list
           (old-copy (chem:copy-filter base-ef t
                                       (list :rosetta
                                             (list 'chem:energy-rosetta-nonbond
                                                   :rep-weight 0.3))))
           ;; New way: copy-energy-function with spec
           (new-copy (chem:copy-energy-function
                      base-ef
                      '(:rosetta
                        (chem:energy-rosetta-nonbond :rep-weight 0.3)))))
      #+tests(test-true copy-ef-rep-weight-old-vs-new
                        (cef--compare-component-energies old-copy new-copy coords)))))

;;; =========================================================================
;;; Test 6: copy-energy-function — with keep-interaction-factory filtering,
;;;         old function-based vs new function passthrough
;;; =========================================================================

(let* ((agg (chem:sample-aggregate))
       (coords (chem:matter/extract-coordinates agg))
       (base-ef (chem:%make-energy-function :matter agg :use-excluded-atoms nil))
       (factory (lambda (aclass)
                  (if (or (eq aclass (find-class 'chem:energy-stretch))
                          (eq aclass (find-class 'chem:energy-angle)))
                      t
                      nil)))
       ;; Old way
       (old-copy (chem:copy-filter base-ef factory))
       ;; New way: function passthrough
       (new-copy (chem:copy-energy-function base-ef factory)))
  #+tests(test-true copy-ef-filter-function-old-vs-new
                    (cef--compare-component-energies old-copy new-copy coords)))

;;; =========================================================================
;;; Test 7: copy-energy-function — :remove vs keep-interaction-factory
;;; =========================================================================

(let* ((agg (chem:sample-aggregate))
       (coords (chem:matter/extract-coordinates agg))
       (base-ef (chem:%make-energy-function :matter agg :use-excluded-atoms nil))
       ;; Old way: factory that skips dihedral
       (old-copy (chem:copy-filter
                  base-ef
                  (lambda (aclass)
                    (if (eq aclass (find-class 'chem:energy-dihedral))
                        nil
                        t))))
       ;; New way: :remove
       (new-copy (chem:copy-energy-function
                  base-ef
                  '(:amber :default t :remove (chem:energy-dihedral)))))
  #+tests(test-true copy-ef-remove-old-vs-new
                    (cef--compare-component-energies old-copy new-copy coords)))

;;; =========================================================================
;;; Test 8: copy-energy-function — :add creates component not in source
;;; =========================================================================

(let* ((agg (chem:sample-aggregate))
       (base-ef (chem:%make-energy-function :matter agg :use-excluded-atoms nil))
       ;; Verify the base doesn't have flat-bottom-anchor-restraint
       (base-has-fbr (chem:find-component-or-nil base-ef 'chem:energy-flat-bottom-anchor-restraint))
       ;; Copy with :add
       (new-copy (chem:copy-energy-function
                  base-ef
                  '(:amber :add (chem:energy-flat-bottom-anchor-restraint))))
       (copy-has-fbr (chem:find-component-or-nil new-copy 'chem:energy-flat-bottom-anchor-restraint)))
  #+tests(test-true copy-ef-add-base-lacks-component (null base-has-fbr))
  #+tests(test-true copy-ef-add-new-has-component (not (null copy-has-fbr))))

;;; =========================================================================
;;; Test 9: make-energy-function — :spec conflicts with old args signal errors
;;; =========================================================================

(let ((agg (chem:sample-aggregate)))
  #+tests(test-true make-ef-spec-and-factory-errors
                    (handler-case
                        (progn (chem:make-energy-function :matter agg
                                                          :spec :amber
                                                          :keep-interaction-factory t)
                               nil)
                      (error () t)))
  #+tests(test-true make-ef-spec-and-setup-errors
                    (handler-case
                        (progn (chem:make-energy-function :matter agg
                                                          :spec :amber
                                                          :setup '(:amber))
                               nil)
                      (error () t))))

(format t "copy-energy-function regression tests complete.~%")
