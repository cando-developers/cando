(in-package :chem.extras)

(export 'chem::coordinates :chem)

(defparameter chem:*save-positions* nil)

(defgeneric chem:compute-merged-nonbond-force-field-for-aggregate (aggregate atom-types))
(defgeneric chem:compute-merged-lksolvation-force-field-for-aggregate (aggregate))

;;(defmethod chem:lookup-atom-properties-radius (atom &optional (force-field-name :default))
  ;; "Return (values atom-radius radius-is-default-t-nil)"
;;  (values 1.5 t))

(defun chem:lookup-atom-properties-radius (atom nonbond-db &optional (force-field-name :default))
  ;; "Return (values atom-radius radius-is-default-t-nil)"
  (if (chem:has-type nonbond-db (chem:get-type atom))
      (let* ((type-index (chem:find-type-index nonbond-db (chem:get-type atom)))
             (ffnonbond (chem:get-ffnonbond-using-type-index nonbond-db type-index))
             (atom-radius (chem:ffnonbond/get-radius-angstroms ffnonbond)))
        (values atom-radius nil))
      (progn
        (values 1.5 t))))


;;; ------------------------------------------------------------
;;;
;;; Force field generic functions
;;;
;;; Find nonbond atom type position in non-bond force-field
;;;

(defgeneric chem:find-atom-type-position (nonbond-force-field type))

(defmethod chem:find-atom-type-position ((nonbond-force-field chem:ffnonbond-db) type)
  (chem:ffnonbond-find-atom-type-position nonbond-force-field type))

(defgeneric chem:find-type (nonbond-force-field-component type)
  (:documentation "This plays the same role as FFNonbondDb_O::findType"))

(defmethod chem:find-type (nonbond-force-field-component type)
  (error "Provide a method for chem:find-type ~s ~s" nonbond-force-field-component type))

(defmethod chem:find-type ((nbforce-field chem:ffnonbond-db) type)
  (chem:ffnonbond-find-type nbforce-field type))

#|
(defgeneric chem:nonbond-force-field-database (force-field)
  (:documentation "Return an object that represents the nonbond component of the force-field"))

(defmethod chem:nonbond-force-field-database ((force-field chem:force-field))
  (chem:force-field-get-nonbond-db force-field))
|#


;;; ------------------------------------------------------------
;;;
;;;  add-shadowing-force-field
;;;
;;;  Adds a force-field to a combined-force-field
;;;



(defgeneric chem:clear-combined-force-field (dest)
  (:documentation "Clear everything out of the combined-force-field"))

(defmethod chem:clear-combined-force-field ((dest chem:combined-force-field))
  (chem:combined-force-field/clear dest))

(defgeneric chem:force-field-component-merge (dest source)
  (:documentation "Merge two force-field components of the specified kind"))

(defmethod chem:force-field-component-merge ((dest chem:ffnonbond-db) (source chem:ffnonbond-db))
  (chem:ffnonbond-db-merge dest source))

(defgeneric chem:add-shadowing-force-field (combined-force-field force-field force-field-info)
  (:documentation "Add a force-field to a combined-force-field.  combined-force-fields contain
multiple force-fields and know how a more recently added force-field shadows a later addition."))

(defmethod chem:add-shadowing-force-field ((combined-force-field chem:combined-force-field)
                                           (force-field chem:force-field)
                                           force-field-info)
  (chem:combined-force-field/add-shadowing-force-field combined-force-field force-field force-field-info))


(defgeneric chem:force-fields-as-list (combined-force-field)
  (:documentation "Return the internal force-fields as a list with more recently added ones before older ones.
 The caller needs to know how they shadow each other."))

(defmethod chem:force-fields-as-list ((combined-force-field chem:combined-force-field))
  (chem:combined-force-field/force-fields-as-list combined-force-field))


(defgeneric chem:nonbond-component (force-field)
  (:documentation "Return the nonbond component of the force-field"))

(defmethod chem:nonbond-component ((force-field chem:force-field))
  (chem:get-nonbond-db force-field))



(defgeneric chem:assign-force-field-types (combined-force-field molecule atom-types)
  (:documentation  "Assign force-field types"))

(defmethod chem:assign-force-field-types ((combined-force-field chem:combined-force-field) molecule atom-types)
  (chem:combined-force-field/assign-force-field-types combined-force-field molecule atom-types))



#+(or)
(defgeneric chem:assign-molecular-force-field-parameters (energy-function first-force-field-part all-force-field-parts molecule)
  (:documentation "Assign the bond, angle and dihedral parameters"))
#+(or)
(defmethod chem:assign-molecular-force-field-parameters (energy-function (first-force-field-part chem:force-field) all-force-field-parts molecule)
  (chem:generate-standard-assign-molecular-force-field-parameters energy-function force-field molecule))


(defgeneric chem:generate-molecule-energy-function-tables (energy-function molecule combined-force-field keep-interaction-factory)
  (:documentation "Generate the molecule energy-function tables"))

(defmethod chem:generate-molecule-energy-function-tables (energy-function molecule (combined-force-field chem:combined-force-field) keep-interaction-factory)
  "Combine AMBER and GAFF force field and frcmods and run parmchk2"
  (let ((merged-force-field (chem:force-field/make)))
    (loop for partial-force-field in (chem:combined-force-field/force-fields-as-list combined-force-field)
          do (chem:force-field/force-field-merge merged-force-field partial-force-field))
    (let ((ffstretch-db (chem:get-stretch-db merged-force-field))
          (ffangle-db (chem:get-angle-db merged-force-field))
          (ffptor-db (chem:get-ptor-db merged-force-field))
          (ffitor-db (chem:get-itor-db merged-force-field)))
      #+(or)(warn "At this point we should run parmchk2 on the stretch/angle/ptor/itor components of the merged-force-field - any missing parameters should be provided by parmchk2")
      (chem:energy-function/generate-standard-energy-function-tables energy-function molecule ffstretch-db ffangle-db ffptor-db ffitor-db :keep-interaction-factory keep-interaction-factory :atom-types (chem:atom-types energy-function)))))


;;; ------------------------------------------------------------
;;;
;;; node-table protocol
;;;
;;;


(defgeneric chem:make-node-table-from-graph (graph nonbond-force-field atom-types)
  (:documentation "Make a node-table from the graph. Return the node-table object."))

(defgeneric chem:node-table-size (node-table)
  (:documentation "Return the number of nodes in the node-table"))

(defgeneric chem:node-table-coordinate-index-times3 (node-table index)
  (:documentation "Return a coordinate index (an index into a coordinate vector, a multiple of 3)
for the node in the table."))

(defgeneric chem:node-table-node-at-index (node-table index)
  (:documentation "Return the node in the node-table at the index."))

(defgeneric chem:node-get-position (node)
  (:documentation "Return the vector position of the node."))

(defgeneric chem:node-set-position (node vector3)
  (:documentation "Set the vector position of the node."))


(defmethod chem:make-node-table-from-graph ((molecule chem:molecule) nonbond-force-field atom-types)
  "Construct a chem:atom-table from a molecule."
  (let ((atom-table (chem:make-atom-table)))
    (chem:construct-from-molecule atom-table molecule nonbond-force-field nil atom-types)
    atom-table))

(defmethod chem:node-table-size ((atom-table chem:atom-table))
  (chem:get-number-of-atoms atom-table))

(defmethod chem:node-table-coordinate-index-times3 ((atom-table chem:atom-table) index)
  (chem:get-coordinate-index-times3-for-atom-at-index atom-table index))

(defmethod chem:node-table-node-at-index ((atom-table chem:atom-table) index)
  (chem:elt-atom atom-table index))

(defmethod chem:node-get-position ((atom chem:atom))
  (chem:get-position atom))

(defmethod chem:node-set-position ((atom chem:atom) vector3)
  (chem:set-position atom vector3))



;;; ------------------------------------------------------------
;;;
;;; Generic functions for bounding-box for aggregates and atom-tables
;;;

(defgeneric chem:bounding-box (object)
  (:documentation "Return the bounding-box for the object")
  (:method (object)
    (declare (ignore object))
    nil))

(defgeneric chem:bounding-box-bound-p (object)
  (:documentation "Return T if the bounding-box is bound to a value."))

(defgeneric (setf chem:bounding-box) (bounding-box object)
  (:documentation "Set the bounding box to the value."))

(defgeneric chem:mak-unbound-bounding-box (object)
  (:documentation "Make the bounding-box unbound."))

(defmethod chem:bounding-box ((atom-table chem:atom-table))
  (chem:atom-table-bounding-box atom-table))

(defmethod chem:bounding-box-bound-p ((atom-table chem:atom-table))
  (chem:atom-table-bounding-box-bound-p atom-table))

(defmethod (setf chem:bounding-box) (bounding-box (atom-table chem:atom-table))
  (chem:atom-table-set-bounding-box atom-table bounding-box)
  bounding-box)

(defmethod chem:mak-unbound-bounding-box ((atom-table chem:atom-table))
  (chem:atom-table-mak-unbound-bounding-box atom-table))


(defmethod chem:bounding-box ((aggregate chem:aggregate))
  (chem:aggregate-bounding-box aggregate))

(defmethod chem:bounding-box-bound-p ((aggregate chem:aggregate))
  (chem:aggregate-bounding-box-bound-p aggregate))

(defmethod (setf chem:bounding-box) (bounding-box (aggregate chem:aggregate))
  (chem:aggregate-set-bounding-box aggregate bounding-box)
  bounding-box)

(defmethod chem:mak-unbound-bounding-box ((aggregate chem:aggregate))
  (chem:aggregate-mak-unbound-bounding-box aggregate))


(defmethod (setf chem:bounding-box) (bounding-box (octtree chem:add-ion-octree))
  (chem:aggregate-set-bounding-box octtree bounding-box)
  bounding-box)


(defmethod chem:bounding-box ((rigid-body-energy-function chem:rigid-body-energy-function))
  (chem:rigid-body-energy-function-bounding-box rigid-body-energy-function))

(defmethod chem:bounding-box-bound-p ((rigid-body-energy-function chem:rigid-body-energy-function))
  (chem:rigid-body-energy-function-bounding-box-bound-p rigid-body-energy-function))

(defmethod (setf chem:bounding-box) (bounding-box (rigid-body-energy-function chem:rigid-body-energy-function))
  (chem:rigid-body-energy-function-set-bounding-box rigid-body-energy-function bounding-box)
  bounding-box)

(defmethod chem:mak-unbound-bounding-box ((rigid-body-energy-function chem:rigid-body-energy-function))
  (chem:rigid-body-energy-function-mak-unbound-bounding-box rigid-body-energy-function))


(defmethod chem:bounding-box ((energy-function chem:energy-function))
  (chem:energy-function-bounding-box energy-function))

(defmethod chem:bounding-box-bound-p ((energy-function chem:energy-function))
  (chem:energy-function-bounding-box-bound-p energy-function))

(defmethod (setf chem:bounding-box) (bounding-box (energy-function chem:energy-function))
  (chem:energy-function-set-bounding-box energy-function bounding-box)
  bounding-box)

(defmethod chem:mak-unbound-bounding-box ((energy-function chem:energy-function))
  (chem:energy-function-mak-unbound-bounding-box energy-function))




(defparameter chem:*cando-database* (chem:make-cando-database :default)
  "This is where topology's are stored")


(defun chem::register-topology (topology &optional name)
  "Register a topology with cando "
  (unless name
    (setf name (chem:get-name topology)))
  #+(or)(unless (keywordp name)
          (setf name (intern (string-upcase (string name)) :keyword)))
  (export name (symbol-package name))
  (chem:setf-find-topology chem:*cando-database* name topology))

(export 'chem::register-topology :chem)


(defun chem::register-foldamer-monomer-context-info (foldamer-monomer-context-info)
  "Register a foldamer-monomer-context-info object with cando "
  (let ((foldamer (chem:foldamer-monomer-context-info-foldamer foldamer-monomer-context-info))
        (monomer-context (chem:foldamer-monomer-context-info-monomer-context foldamer-monomer-context-info)))
    (chem:setf-find-foldamer-monomer-context chem:*cando-database* foldamer monomer-context foldamer-monomer-context-info)))

(export 'chem::register-foldamer-monomer-context-info :chem)



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Loops for molecules, residues, atoms
;;;


(defmacro chem::do-molecules ((molecule-var matter) &body body)
  `(flet ((do-one-molecule (,molecule-var)
            ,@body))
     (if (typep ,matter 'chem:molecule)
         (funcall #'do-one-molecule ,matter) 
         (chem:map-molecules
          nil
          #'do-one-molecule
          ,matter))))

(defmacro chem::do-residues ((residue-var matter) &body body)
  `(flet ((do-one-residue (,residue-var)
            ,@body))
     (if (typep ,matter 'chem:residue)
         (funcall #'do-one-residue ,matter) 
         (chem:map-residues
          nil
          #'do-one-residue
          ,matter))))

(defmacro chem::do-atoms ((atom-var matter) &body body)
  `(flet ((do-one-atom (,atom-var)
            ,@body))
     (if (typep ,matter 'chem:atom)
         (funcall #'do-one-atom ,matter) 
         (chem:map-atoms
          nil
          #'do-one-atom
          ,matter))))


(export '(chem::do-atoms chem::do-residues chem::do-molecules) :chem)



;;
;; Identify all rings, isolate atoms in rings and apply all of the aromaticity rules
;; from Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
;; Return the rings for bond type assignment
;; Make this work using the hash-table returned by chem:identify-aromatic-rings


(defmacro chem::with-aromaticity-information ((var matter force-field-name) &body body)
  "Provide a dynamic environment where the aromaticity:is-aromatic function works.
This macro first calculates the aromaticity information of the matter and then
evaluates the body in that dynamic environment."
  `(let* ((,var (chem:identify-aromatic-rings ,matter ,force-field-name))
          (chem:*current-aromaticity-information* ,var))
     ,@body))


(export 'chem::with-aromaticity-information :chem)


(defgeneric chem:summarize-energy-function (energy-function pos active-atom-mask))

(defparameter *save-pos* nil)
(defmethod chem:summarize-energy-function ((energy-function chem:energy-function) pos active-atom-mask)
  (setf *save-pos* pos)
  (with-output-to-string (sout)
    (let* ((energy-components (chem:make-energy-components))
           (total-energy (chem:evaluate-all energy-function pos :energy-components energy-components)))
      (format sout "Total E. ~10,5f  mask ~a~%" total-energy (not (null active-atom-mask)))
      (format sout "~s~%" (chem:energy-components/components energy-components)))))
#||
      (loop for component in components
            for comp-energy = (chem:energy-component-evaluate-energy energy-function component pos active-atom-mask)
            do (format sout "  ~30a ~a ~10,5f~%" (class-name (class-of component))
                       (chem:energy-component/is-enabled component)
                       comp-energy)
            when (typep component 'chem:energy-nonbond)
              do (format sout "      ee ~10,5f    vdw ~10,5f~%"
                         (chem:energy-nonbond/get-vdw-energy component)
                         (chem:energy-nonbond/get-electrostatic-energy component)
                         )
            )
      )))
||#


(defmethod print-object ((energy-function chem:energy-function) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (energy-function stream :type t :identity t))))

(defmethod print-object ((closure core::closure) stream)
  (print-unreadable-object (closure stream :type t :identity t)
    (format stream "(LAMBDA ~s)" (core:function-lambda-list closure))))


(defun amber-energy-components (energy-components)
  (let* ((comps (chem:energy-components/components energy-components))
         (bond (cdr (assoc 'chem:energy-stretch comps)))
         (angle (cdr (assoc 'chem:energy-angle comps)))
         (dihed (cdr (assoc 'chem:energy-dihedral comps)))
         (vdwaals (+ (or (cdr (assoc 'chem:energy-vdw-excluded-atoms comps)) 0.0)
                     (or (cdr (assoc 'chem:energy-vdw comps)) 0.0)))
         (eel (+ (or (cdr (assoc 'chem:energy-electrostatic-excluded-atoms comps)) 0.0)
                 (or (cdr (assoc 'chem:energy-electrostatic comps)) 0.0)))
         (vdw1-4 (cdr (assoc 'chem:energy-vdw14 comps)))
         (eel1-4 (cdr (assoc 'chem:energy-electrostatic14 comps)))
         (restraint (+ (or (cdr (assoc 'chem:energy-anchor-restraint comps)) 0.0)
                       (or (cdr (assoc 'chem:energy-chiral-restraint comps)) 0.0)
                       (or (cdr (assoc 'chem:energy-dihedral-restraint comps)) 0.0)
                       ))
         )
    (with-output-to-string (sout)
      (format sout "BOND     = ~12,4f  ANGLE   = ~12,4f   DIHED     = ~12,4f~%" bond angle dihed)
      (format sout "VDWAALS  = ~12,4f  EEL     = ~12,4f   HBOND     = ~12,4f~%" vdwaals eel 0.0)
      (format sout "1-4 VDW  = ~12,4f  1-4 EEL = ~12,4f   RESTRAINT = ~12,4f~%" vdw1-4 eel1-4 restraint)
    )))

(define-condition chem:parameterization-error (error) ()
  (:documentation "Base condition for parameterization errors"))

(define-condition chem:missing-type (chem:parameterization-error)
  ((molecule :initarg :molecule :accessor molecule)
   (a1-name :initarg :a1-name :accessor a1-name)
   (a1-element :initarg :a1-element :accessor a1-element))
  (:report (lambda (obj stream)
             (let* ((name (chem:get-name (molecule obj)))
                    (properties (chem:matter/properties (molecule obj)))
                    (maybe-description (getf properties :description)))
               (format stream "Molecule name: ~s cannot identify OpenFF nonbond parameter for atom(element)s: (~a(~a)); ~a"
                       (chem:get-name (molecule obj))
                       (a1-name obj)
                       (a1-element obj)
                       (if maybe-description
                           (format nil "description: ~a" maybe-description)
                           "")
                       )))))

(defun chem:missing-type-error (a1 molecule)
  (let* ((a1-name (chem:get-name a1))
         (a1-element (chem:get-element a1)))
    (error 'chem:missing-type
           :molecule molecule
           :a1-name a1-name
           :a1-element a1-element
           )))

(define-condition chem:missing-stretch (chem:parameterization-error)
  ((molecule :initarg :molecule :accessor molecule)
   (a1-name :initarg :a1-name :accessor a1-name)
   (a2-name :initarg :a2-name :accessor a2-name)
   (a1-element :initarg :a1-element :accessor a1-element)
   (a2-element :initarg :a2-element :accessor a2-element)
   (a1-a2-bond-order :initarg :a1-a2-bond-order :accessor a1-a2-bond-order))
  (:report (lambda (obj stream)
             (let* ((name (chem:get-name (molecule obj)))
                    (properties (chem:matter/properties (molecule obj)))
                    (maybe-description (getf properties :description)))
               (format stream "Molecule name: ~s cannot identify OpenFF stretch parameter for list of two atom(element)s: (~a(~a) ~a(~a)); ~a"
                       (chem:get-name (molecule obj))
                       (a1-name obj)
                       (a1-element obj)
                       (a2-name obj)
                       (a2-element obj)
                       (if maybe-description
                           (format nil "description: ~a" maybe-description)
                           "")
                       )))))

(defun chem:missing-stretch-error (key molecule)
  (let* ((a1 (first key))
         (a2 (second key))
         (a1-name (chem:get-name a1))
         (a2-name (chem:get-name a2))
         (a1-element (chem:get-element a1))
         (a2-element (chem:get-element a2))
         (a1-a2-bond-order (chem:bond-order-to a1 a2))
         )
    (error 'chem:missing-stretch
           :molecule molecule
           :a1-name a1-name
           :a2-name a2-name
           :a1-element a1-element
           :a2-element a2-element
           :a1-a2-bond-order a1-a2-bond-order
           )))

(define-condition chem:missing-angle (chem:parameterization-error)
  ((molecule :initarg :molecule :accessor molecule)
   (a1-name :initarg :a1-name :accessor a1-name)
   (a2-name :initarg :a2-name :accessor a2-name)
   (a3-name :initarg :a3-name :accessor a3-name)
   (a1-element :initarg :a1-element :accessor a1-element)
   (a2-element :initarg :a2-element :accessor a2-element)
   (a3-element :initarg :a3-element :accessor a3-element)
   (a1-a2-bond-order :initarg :a1-a2-bond-order :accessor a1-a2-bond-order)
   (a2-a3-bond-order :initarg :a2-a3-bond-order :accessor a2-a3-bond-order))
  (:report (lambda (obj stream)
             (let* ((name (chem:get-name (molecule obj)))
                    (properties (chem:matter/properties (molecule obj)))
                    (maybe-description (getf properties :description)))
               (format stream "Molecule name: ~s cannot identify OpenFF angle parameter for list of three atom(element)s: (~a(~a) ~a(~a) ~a(~a)); ~a"
                       (chem:get-name (molecule obj))
                       (a1-name obj)
                       (a1-element obj)
                       (a2-name obj)
                       (a2-element obj)
                       (a3-name obj)
                       (a3-element obj)
                       (if maybe-description
                           (format nil "description: ~a" maybe-description)
                           "")
                       )))))

(defun chem:missing-angle-error (key molecule)
  (let* ((a1 (first key))
         (a2 (second key))
         (a3 (third key))
         (a1-name (chem:get-name a1))
         (a2-name (chem:get-name a2))
         (a3-name (chem:get-name a3))
         (a1-element (chem:get-element a1))
         (a2-element (chem:get-element a2))
         (a3-element (chem:get-element a3))
         (a1-a2-bond-order (chem:bond-order-to a1 a2))
         (a2-a3-bond-order (chem:bond-order-to a2 a3))
         )
    (error 'chem:missing-angle
           :molecule molecule
           :a1-name a1-name
           :a2-name a2-name
           :a3-name a3-name
           :a1-element a1-element
           :a2-element a2-element
           :a3-element a3-element
           :a1-a2-bond-order a1-a2-bond-order
           :a2-a3-bond-order a2-a3-bond-order
           )))

(define-condition chem:missing-dihedral (chem:parameterization-error)
  ((molecule :initarg :molecule :accessor molecule)
   (a1-name :initarg :a1-name :accessor a1-name)
   (a2-name :initarg :a2-name :accessor a2-name)
   (a3-name :initarg :a3-name :accessor a3-name)
   (a4-name :initarg :a4-name :accessor a4-name)
   (a1-element :initarg :a1-element :accessor a1-element)
   (a2-element :initarg :a2-element :accessor a2-element)
   (a3-element :initarg :a3-element :accessor a3-element)
   (a4-element :initarg :a4-element :accessor a4-element)
   (a1-a2-bond-order :initarg :a1-a2-bond-order :accessor a1-a2-bond-order)
   (a2-a3-bond-order :initarg :a2-a3-bond-order :accessor a2-a3-bond-order)
   (a3-a4-bond-order :initarg :a3-a4-bond-order :accessor a3-a4-bond-order))
  (:report (lambda (obj stream)
             (let* ((name (chem:get-name (molecule obj)))
                    (properties (chem:matter/properties (molecule obj)))
                    (maybe-description (getf properties :description)))
               (format stream "Molecule name: ~s cannot identify OpenFF dihedral parameter for list of four atom(element)s: (~a(~a) ~a(~a) ~a(~a) ~a(~a)); ~a"
                       (chem:get-name (molecule obj))
                       (a1-name obj)
                       (a1-element obj)
                       (a2-name obj)
                       (a2-element obj)
                       (a3-name obj)
                       (a3-element obj)
                       (a4-name obj)
                       (a4-element obj)
                       (if maybe-description
                           (format nil "description: ~a" maybe-description)
                           "")
                       )))))

(defun chem:missing-dihedral-error (a1 a2 a3 a4 molecule)
  (let ((a1-name (chem:get-name a1))
        (a2-name (chem:get-name a2))
        (a3-name (chem:get-name a3))
        (a4-name (chem:get-name a4))
        (a1-element (chem:get-element a1))
        (a2-element (chem:get-element a2))
        (a3-element (chem:get-element a3))
        (a4-element (chem:get-element a4))
        (a1-a2-bond-order (chem:bond-order-to a1 a2))
        (a2-a3-bond-order (chem:bond-order-to a2 a3))
        (a3-a4-bond-order (chem:bond-order-to a3 a4))
        )
    (error 'chem:missing-dihedral
           :molecule molecule
           :a1-name a1-name
           :a2-name a2-name
           :a3-name a3-name
           :a4-name a4-name
           :a1-element a1-element
           :a2-element a2-element
           :a3-element a3-element
           :a4-element a4-element
           :a1-a2-bond-order a1-a2-bond-order
           :a2-a3-bond-order a2-a3-bond-order
           :a3-a4-bond-order a3-a4-bond-order
           )))


;;; ---------------------------------------------------------------------------
;;; Energy spec: a unified s-expression for configuring energy functions.
;;;
;;; A spec can be:
;;;   :amber or :rosetta     — PRESET, all components, all interactions
;;;   T                      — keep all (for copy-filter)
;;;   NIL                    — keep none (for copy-filter)
;;;   (PRESET . clauses)     — preset with customization
;;;
;;; Clauses in a spec list:
;;;   :default T/NIL/fn      — behavior for components not explicitly listed
;;;   :add (class ...)       — add components beyond the preset
;;;   :remove (class ...)    — remove components from the preset
;;;   (class-name ...)       — per-component entry: T, :filter fn, :param val
;;;
;;; A per-component entry is (class-name . options) where options can be:
;;;   T                      — keep all interactions for this component
;;;   :filter fn             — filter interactions with this closure
;;;   :key value ...         — set parameters (passed through to setup)
;;;
;;; chem:parse-energy-spec returns three values:
;;;   1. keep-interaction-factory — T, NIL, or a function of (aclass)
;;;   2. setup — list suitable for the C++ setup parameter (:rosetta ...)
;;;   3. add-components — list of class names to add (for copy-filter)
;;; ---------------------------------------------------------------------------

(defun chem::%%parse-spec-clauses (clauses)
  "Parse the clauses of an energy spec list.
Returns (values default add remove component-entries).
component-entries is a list of (class-name . options) entries."
  (let ((default t)
        (add nil)
        (remove nil)
        (component-entries nil))
    (loop with rest = clauses
          while rest
          do (let ((item (car rest)))
               (cond
                 ((eq item :default)
                  (setf default (cadr rest))
                  (setf rest (cddr rest)))
                 ((eq item :add)
                  (setf add (cadr rest))
                  (setf rest (cddr rest)))
                 ((eq item :remove)
                  (setf remove (cadr rest))
                  (setf rest (cddr rest)))
                 ((consp item)
                  ;; Per-component entry: (class-name . options)
                  (push item component-entries)
                  (setf rest (cdr rest)))
                 (t
                  (error "Unexpected item in energy spec: ~s" item)))))
    (values default add remove (nreverse component-entries))))

(defun chem::%%build-component-filter-map (component-entries)
  "Build a hash table mapping component class to its filter value.
The filter value is T (keep all) or a closure.
Also returns a list of (class-name . param-plist) for setup."
  (let ((filter-map (make-hash-table :test 'eq))
        (setup-params nil))
    (dolist (entry component-entries)
      (let* ((class-name (car entry))
             (options (cdr entry))
             (class (if (symbolp class-name)
                        (find-class class-name)
                        class-name)))
        ;; Parse the options for this component
        (cond
          ;; (class-name t) or (class-name) — keep all
          ((or (null options) (eq (car options) t))
           (setf (gethash class filter-map) t))
          ;; (class-name :filter fn :param val ...)
          (t
           (let ((filter t)
                 (params nil))
             (loop with rest = options
                   while rest
                   do (cond
                        ((eq (car rest) :filter)
                         (setf filter (cadr rest))
                         (setf rest (cddr rest)))
                        (t
                         ;; It's a parameter key-value pair
                         (push (cadr rest) params)
                         (push (car rest) params)
                         (setf rest (cddr rest)))))
             (setf (gethash class filter-map) filter)
             (when params
               (push (cons class-name (nreverse params)) setup-params)))))))
    (values filter-map (nreverse setup-params))))

(defun chem::%%build-keep-interaction-factory (default add remove filter-map)
  "Build a keep-interaction-factory function from parsed spec components.
Returns T, NIL, or a closure of (aclass).
The returned function, when called with a component class, returns:
  NIL      — skip this component entirely
  T        — keep all interactions
  closure  — filter interactions with this closure"
  (let ((remove-set (when remove
                      (let ((ht (make-hash-table :test 'eq)))
                        (dolist (class-name remove)
                          (setf (gethash (find-class class-name) ht) t))
                        ht)))
        (add-set (when add
                   (let ((ht (make-hash-table :test 'eq)))
                     (dolist (class-name add)
                       (setf (gethash (find-class class-name) ht) t))
                     ht)))
        (has-entries (plusp (hash-table-count filter-map))))
    ;; If no customization at all, just return the default
    (when (and (eq default t) (null remove) (null add) (not has-entries))
      (return-from chem::%%build-keep-interaction-factory t))
    (when (and (null default) (null add) (not has-entries))
      (return-from chem::%%build-keep-interaction-factory nil))
    ;; Build a factory function
    (lambda (aclass)
      (block result
        ;; Check removal first
        (when (and remove-set (gethash aclass remove-set))
          (return-from result nil))
        ;; Check explicit component entry
        (multiple-value-bind (filter found)
            (gethash aclass filter-map)
          (when found
            (return-from result filter)))
        ;; Check if it's an added component (not in filter-map)
        (when (and add-set (gethash aclass add-set))
          (return-from result default))
        ;; Fall through to default
        default))))

(defun chem:parse-energy-spec (spec)
  "Parse an energy spec into values for the existing C++ API.
Returns (values keep-interaction-factory setup add-components).

SPEC can be:
  :amber            — default AMBER energy function
  :rosetta          — Rosetta energy function, all components
  T                 — keep all (for copy-filter)
  NIL               — keep none (for copy-filter)
  (:rosetta ...)    — Rosetta with customization
  (:amber ...)      — AMBER with customization"
  (cond
    ;; T and NIL pass through directly
    ((eq spec t)
     (values t nil nil))
    ((null spec)
     (values nil nil nil))
    ;; Bare keyword preset
    ((keywordp spec)
     (values t (list spec) nil))
    ;; Function — treat as a raw keep-interaction-factory (backward compat)
    ((functionp spec)
     (values spec nil nil))
    ;; Full spec list
    ((consp spec)
     (let ((preset (car spec))
           (clauses (cdr spec)))
       (unless (keywordp preset)
         (error "Energy spec must start with a preset keyword (:amber or :rosetta), got ~s" preset))
       (multiple-value-bind (default add remove component-entries)
           (chem::%%parse-spec-clauses clauses)
         (multiple-value-bind (filter-map setup-params)
             (chem::%%build-component-filter-map component-entries)
           (let ((keep-interaction-factory
                   (chem::%%build-keep-interaction-factory default add remove filter-map))
                 (setup (cons preset setup-params)))
             (values keep-interaction-factory setup add))))))
    (t
     (error "Invalid energy spec: ~s — expected T, NIL, a keyword, a function, or a list" spec))))

;;; ---------------------------------------------------------------------------
;;; Wrappers that accept :spec and delegate to the C++ API
;;; ---------------------------------------------------------------------------

(defun chem:make-energy-function (&rest args &key matter spec
                                    use-excluded-atoms keep-interaction-factory
                                    assign-types setup
                                    disable-components enable-components)
  "Create an energy function for MATTER.
If :spec is provided, it is parsed to derive :keep-interaction-factory and :setup.
Otherwise, :keep-interaction-factory and :setup are used directly (backward compatible).
It is an error to provide :spec together with :keep-interaction-factory or :setup.

See chem:parse-energy-spec for the spec syntax.

Examples:
  (chem:make-energy-function :matter agg)                      ; default AMBER
  (chem:make-energy-function :matter agg :spec :rosetta)       ; Rosetta preset
  (chem:make-energy-function :matter agg
    :spec '(:rosetta
            (chem:energy-rosetta-nonbond :rep-weight 0.5)))"
  (declare (ignore matter use-excluded-atoms assign-types
                   disable-components enable-components))
  (if spec
      (progn
        (when keep-interaction-factory
          (error "Cannot provide both :spec and :keep-interaction-factory to make-energy-function"))
        (when setup
          (error "Cannot provide both :spec and :setup to make-energy-function"))
        (multiple-value-bind (parsed-factory parsed-setup)
            (chem:parse-energy-spec spec)
          (let ((cleaned-args (loop for (key val) on args by #'cddr
                                    unless (member key '(:spec :keep-interaction-factory :setup))
                                      nconc (list key val))))
            (apply #'chem:%make-energy-function
                   :keep-interaction-factory parsed-factory
                   :setup parsed-setup
                   cleaned-args))))
      ;; Legacy path: pass through directly
      (apply #'chem:%make-energy-function
             (loop for (key val) on args by #'cddr
                   unless (eq key :spec)
                     nconc (list key val)))))

(defun chem:copy-energy-function (energy-function spec)
  "Copy an energy function, applying the given energy spec.
SPEC controls which components are included, how interactions are filtered,
and what parameters are set. See chem:parse-energy-spec for spec syntax.

Examples:
  (chem:copy-energy-function ef t)              ; full copy
  (chem:copy-energy-function ef nil)            ; empty copy
  (chem:copy-energy-function ef :rosetta)       ; copy all rosetta components
  (chem:copy-energy-function ef
    '(:rosetta
      :add (chem:energy-anchor-restraint)
      (chem:energy-rosetta-nonbond :rep-weight 0.2
                                   :filter #'my-filter)))"
  (multiple-value-bind (keep-interaction-factory setup add-components)
      (chem:parse-energy-spec spec)
    (chem:copy-filter energy-function keep-interaction-factory setup add-components)))

;;; ---------------------------------------------------------------------------

(defun chem:find-lksolvation-type (lksolvation-db type &key (errorp t))
  "Used by energyRosettaLKSolvation.cc"
  (let ((ff (chem:fflksolvation-find-type lksolvation-db type)))
    (if ff
        ff
        (when errorp
          (error "Could not find LKSolvation parameters for :lk-solvation-atom-type ~s" type)))))


