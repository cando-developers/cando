#+tests(in-package #:clasp-tests)

;;;; cached-types.lisp
;;;;
;;;; Regression test for the per-atom nonbond parameter cache
;;;; (EnergyRosettaNonbond_O::ensureParameterCache / tryAddTermCached and the
;;;; pairList.h hookup).  The cache resolves each atom's nonbond type -> radius/epsilon
;;;; ONCE instead of calling getType + find-type per pair; it must not change any energy.
;;;;
;;;; This test loads a structure, builds a Rosetta energy function, and checks the
;;;; rosetta-nonbond / rosetta-elec / rosetta-lksolvation component energies against
;;;; saved reference values.
;;;;
;;;; WORKFLOW
;;;;   1. Build cando with the cached-types code DISABLED (or before enabling it).
;;;;   2. Load this file.  It prints an "==== REFERENCE energies ====" block.
;;;;   3. Paste those numbers into the *ref-...* parameters below.
;;;;   4. Re-enable the cached-types code, rebuild, reload this file.  The (test-true ...)
;;;;      forms then verify the cached path reproduces the same energies to *tolerance*.
;;;;
;;;; NOTE on components: a plain load-mol2 + :setup (:rosetta) types nonbond (given atom
;;;; types) and elec (charges from the mol2), but lksolvation needs a :lk-solvation-atom-type
;;;; property that a generic mol2 does not carry.  If the lksolvation line prints <absent>
;;;; or "0 terms", the structure/setup needs lk typing (or pick a structure that has it);
;;;; nonbond+elec still fully exercise the Level-2 cache.

(eval-when (:compile-toplevel :load-toplevel :execute)
  ;; These packages must exist while COMPILE-FILE reads their qualified symbols.
  (ql:quickload :open-force-field/load)
  (ql:quickload :spiros/spiros))

(defun cached-types--make-force-field-registry ()
  "Return a private registry containing canonical SMIRNOFF and Rosetta data.

Other regression files install small synthetic force fields under the standard
:SMIRNOFF and :ROSETTA names.  QUICKLOAD does not undo those changes when its
systems are already loaded, so this test must establish its own parameter state."
  (let ((leap.core:*force-fields* (make-hash-table))
        ;; Force the public loader to populate this private registry even when
        ;; OPEN-FORCE-FIELD/LOAD has already run in the surrounding image.
        (open-force-field::*loaded-already* nil))
    ;; Aggregates default to :DEFAULT even when their molecules name a more
    ;; specific force field.  Reproduce LEaP's normal empty aggregate-level
    ;; force field so nonbond merging has the same inputs as a fresh image.
    (leap.core:add-force-field-or-modification
     (chem:force-field/make)
     :force-field-name :default
     :combined-force-field-class-name 'chem:combined-force-field)
    (open-force-field:load-smirnoff)
    (rosetta.params:load-rosetta-params)
    leap.core:*force-fields*))

(defparameter *cached-types-force-fields*
  (cached-types--make-force-field-registry))

(defparameter *cached-types-structure*
  (merge-pathnames "data/struct-0000.mol2"
                   #.(or *compile-file-truename* *load-truename*)))

;; --- Reference energies: paste from a cache-DISABLED run.  NIL = not captured yet. ---
(defparameter *ref-rosetta-nonbond*     1.0360468883d0)
(defparameter *ref-rosetta-elec*        -155.7431664998d0)
(defparameter *ref-rosetta-lksolvation* 5.9383679081d0)

(defparameter *cached-types-tolerance* 1d-6)

(defparameter *cached-types-soak-iterations* 2
  "Number of repeated cache evaluations before the reference checks.

Two iterations exercise cache reuse without turning the regression suite into
a benchmark.  Set this to a larger value for an interactive performance soak.")

(defun cached-types--find-component (energy-function class-name)
  "The rosetta component of type CLASS-NAME, or NIL if none is present."
  (find-if (lambda (c) (typep c class-name)) (chem:all-components energy-function)))

(defun cached-types--component-energy (energy-function pos class-name)
  "Return (values energy number-of-terms) for the CLASS-NAME component, or (values NIL NIL)."
  (let ((comp (cached-types--find-component energy-function class-name)))
    (if comp
        (values (chem:energy-component-evaluate-energy energy-function comp pos)
                (chem:number-of-terms comp))
        (values nil nil))))

(defun cached-types--energies ()
  "Load the structure, build a Rosetta energy function, and return
   (values nb-energy nb-terms el-energy el-terms lk-energy lk-terms)."
  (let ((leap.core:*force-fields* *cached-types-force-fields*))
    (let* ((agg (let ((ag (cando:load-mol2 *cached-types-structure*)))
                (chem:setf-force-field-name (cando:mol ag) :smirnoff)
                ag))
         (rings (chem:identify-rings agg))
         (chem:*current-rings* rings)
         (ff-name (or (chem:force-field-name agg) :mdl))
         (aromaticity (chem:identify-aromatic-rings agg ff-name))
         (chem:*current-aromaticity-information* aromaticity)
         (assigned (rosetta.params:assign-lk-types
                    agg (spiros::load-spiro-lk-type-rules) :error-on-missing t)))
    (maphash (lambda (atom lk-type)
               (when lk-type
                 (chem:set-property atom :lk-solvation-atom-type lk-type)))
             assigned)
    ;; --- now build the rosetta EF; lksolvation will have terms ---
      (let* ((ef (chem:make-energy-function :matter agg :use-excluded-atoms nil
                                            :assign-types t :setup '(:rosetta)))
             (pos (chem:make-nvector (chem:get-nvector-size ef))))
        (chem:load-coordinates-into-vector ef pos)
        (dotimes (ii *cached-types-soak-iterations*)
          (setf (aref pos 0) (+ (aref pos 0) 1d-4))
          (multiple-value-bind (nb nbt) (cached-types--component-energy ef pos 'chem:energy-rosetta-nonbond)
            (multiple-value-bind (el elt) (cached-types--component-energy ef pos 'chem:energy-rosetta-elec)
              (multiple-value-bind (lk lkt) (cached-types--component-energy ef pos 'chem:energy-rosetta-lksolvation)
                (values nb nbt el elt lk lkt))))))
    (let* ((ef (chem:make-energy-function :matter agg :use-excluded-atoms nil
                                          :assign-types t :setup '(:rosetta)))
           (pos (chem:make-nvector (chem:get-nvector-size ef))))
      (chem:load-coordinates-into-vector ef pos)
      (multiple-value-bind (nb nbt) (cached-types--component-energy ef pos 'chem:energy-rosetta-nonbond)
        (multiple-value-bind (el elt) (cached-types--component-energy ef pos 'chem:energy-rosetta-elec)
          (multiple-value-bind (lk lkt) (cached-types--component-energy ef pos 'chem:energy-rosetta-lksolvation)
            (values nb nbt el elt lk lkt))))))))
    

(defun cached-types--fmt (v) (if v (format nil "~,10f" v) "<absent>"))

(multiple-value-bind (nb nbt el elt lk lkt) (cached-types--energies)
  (format t "~&==== cached-types REFERENCE energies (paste into the *ref-...* params) ====~%")
  (format t "  rosetta-nonbond      = ~a   (~a terms)~%" (cached-types--fmt nb) (or nbt 0))
  (format t "  rosetta-elec         = ~a   (~a terms)~%" (cached-types--fmt el) (or elt 0))
  (format t "  rosetta-lksolvation  = ~a   (~a terms)~%" (cached-types--fmt lk) (or lkt 0))
  (format t "============================================================================~%")
  (finish-output)

  (flet ((approx= (a b)
           (and (numberp a) (numberp b) (< (abs (- a b)) *cached-types-tolerance*))))
    ;; Each test runs only once its reference has been filled in.
    #+tests
    (progn
      (when *ref-rosetta-nonbond*
        (test-true cached-types-rosetta-nonbond
                   (approx= nb *ref-rosetta-nonbond*)))
      (when *ref-rosetta-elec*
        (test-true cached-types-rosetta-elec
                   (approx= el *ref-rosetta-elec*)))
      (when *ref-rosetta-lksolvation*
        (test-true cached-types-rosetta-lksolvation
                   (approx= lk *ref-rosetta-lksolvation*))))))



;; --- between-matters references (rebuildPairListBetweenMattersImpl path) ---
(defparameter *ref-bm-rosetta-nonbond*     nil)
(defparameter *ref-bm-rosetta-elec*        nil)
(defparameter *ref-bm-rosetta-lksolvation* nil)

(defun cached-types--between-matters-energies ()
  "Build the rosetta EF, then restrict each component to pairs between two residues via
     chem:set-matters — this routes evaluation through rebuildPairListBetweenMattersImpl
     (the path the design's single-scan / pair-scan population uses).
     Returns (values nb nbt el elt lk lkt)."
  (let ((leap.core:*force-fields* *cached-types-force-fields*))
    (let* ((agg (let ((ag (cando:load-mol2 *cached-types-structure*)))
                (chem:setf-force-field-name (cando:mol ag) :smirnoff)
                ag)) 
         (rings (chem:identify-rings agg))
         (chem:*current-rings* rings)
         (ff-name (or (chem:force-field-name agg) :mdl))
         (aromaticity (chem:identify-aromatic-rings agg ff-name))
         (chem:*current-aromaticity-information* aromaticity)
         (assigned (rosetta.params:assign-lk-types
                    agg (spiros::load-spiro-lk-type-rules) :error-on-missing t)))
    (maphash (lambda (atom lk-type)
               (when lk-type (chem:set-property atom :lk-solvation-atom-type lk-type)))
             assigned)
    (let* ((ef  (chem:make-energy-function :matter agg :use-excluded-atoms nil
                                           :assign-types t :setup '(:rosetta)))
           (pos (chem:make-nvector (chem:get-nvector-size ef)))
           (mol  (cando:mol agg 0))
           (mat1 (cando:res mol 0)) ; two residues -> between-matters pairs
           (mat2 (cando:res mol 1)))
      (chem:load-coordinates-into-vector ef pos)
      (dolist (class-name '(chem:energy-rosetta-nonbond
                            chem:energy-rosetta-elec 
                            chem:energy-rosetta-lksolvation))
        (let ((comp (cached-types--find-component ef class-name)))
          (when comp
            (chem:set-matters comp mat1 mat2)))) ; flips comp to between-matters rebuild
      (multiple-value-bind (nb nbt) (cached-types--component-energy ef pos 'chem:energy-rosetta-nonbond)
        (multiple-value-bind (el elt) (cached-types--component-energy ef pos 'chem:energy-rosetta-elec)
          (multiple-value-bind (lk lkt) (cached-types--component-energy ef pos 'chem:energy-rosetta-lksolvation)
            (values nb nbt el elt lk lkt))))))))

(multiple-value-bind (nb nbt el elt lk lkt) (cached-types--between-matters-energies)
  (format t "~&==== cached-types BETWEEN-MATTERS energies (rebuildPairListBetweenMattersImpl) ====~%")
  (format t "  bm-rosetta-nonbond      = ~a   (~a terms)~%" (cached-types--fmt nb) (or nbt 0))
  (format t "  bm-rosetta-elec         = ~a   (~a terms)~%" (cached-types--fmt el) (or elt 0))
  (format t "  bm-rosetta-lksolvation  = ~a   (~a terms)~%" (cached-types--fmt lk) (or lkt 0))
  (format t "===============================================================================~%")
  (finish-output)
  (flet ((approx= (a b) (and (numberp a) (numberp b) (< (abs (- a b)) *cached-types-tolerance*))))
    #+tests
    (progn
      (test-true cached-types-bm-rosetta-nonbond-present
                 (and (numberp nb) (integerp nbt) (plusp nbt)))
      (test-true cached-types-bm-rosetta-elec-present
                 (and (numberp el) (integerp elt) (plusp elt)))
      (test-true cached-types-bm-rosetta-lksolvation-present
                 (and (numberp lk) (integerp lkt) (plusp lkt)))
      (when *ref-bm-rosetta-nonbond*
        (test-true cached-types-bm-rosetta-nonbond (approx= nb *ref-bm-rosetta-nonbond*)))
      (when *ref-bm-rosetta-elec*
        (test-true cached-types-bm-rosetta-elec (approx= el *ref-bm-rosetta-elec*)))
      (when *ref-bm-rosetta-lksolvation*
        (test-true cached-types-bm-rosetta-lksolvation (approx= lk *ref-bm-rosetta-lksolvation*))))))
