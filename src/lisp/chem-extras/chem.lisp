(in-package :chem.extras)


(defgeneric chem:compute-merged-nonbond-force-field-for-aggregate (aggregate atom-types))

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
           (total-energy (chem:evaluate-all energy-function pos :component-energy energy-components)))
      (format sout "Total E. ~10,5f  mask ~a~%" total-energy (not (null active-atom-mask)))
      (format sout "~s~%" (chem:energy-components/component-energies energy-components)))))
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


(defun amber-energy-components (energy-components)
  (let* ((comps (chem:energy-components/component-energies energy-components))
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
