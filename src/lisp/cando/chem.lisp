(in-package :cando)


(defgeneric chem:compute-merged-nonbond-force-field-for-aggregate (aggregate))

;;(defmethod chem:lookup-atom-properties-radius (atom &optional (force-field-name :default))
  ;; "Return (values atom-radius radius-is-default-t-nil)"
;;  (values 1.5 t))

(defun chem:lookup-atom-properties-radius (atom nonbond-db &optional (force-field-name :default))
  ;; "Return (values atom-radius radius-is-default-t-nil)"
  (if (chem:has-type nonbond-db (chem:get-type atom))
      (let* ((type-index (chem:find-type-index nonbond-db (chem:get-type atom)))
             (ffnonbond (chem:get-ffnonbond-using-type-index nonbond-db type-index))
             (atom-radius (chem:get-radius-angstroms ffnonbond)))
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
  (chem:combined-force-field-clear dest))

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
  (chem:combined-force-field-add-shadowing-force-field combined-force-field force-field force-field-info))


(defgeneric chem:force-fields-as-list (combined-force-field)
  (:documentation "Return the internal force-fields as a list with more recently added ones before older ones.
 The caller needs to know how they shadow each other."))

(defmethod chem:force-fields-as-list ((combined-force-field chem:combined-force-field))
  (chem:combined-force-field-force-fields-as-list combined-force-field))


(defgeneric chem:nonbond-component (force-field)
  (:documentation "Return the nonbond component of the force-field"))

(defmethod chem:nonbond-component ((force-field chem:force-field))
  (chem:get-nonbond-db force-field))



(defgeneric chem:assign-force-field-types (combined-force-field molecule)
  (:documentation  "Assign force-field types"))

(defmethod chem:assign-force-field-types ((combined-force-field chem:combined-force-field) molecule)
  (chem:combined-force-field-assign-force-field-types combined-force-field molecule))



#+(or)
(defgeneric chem:assign-molecular-force-field-parameters (energy-function first-force-field-part all-force-field-parts molecule)
  (:documentation "Assign the bond, angle and dihedral parameters"))
#+(or)
(defmethod chem:assign-molecular-force-field-parameters (energy-function (first-force-field-part chem:force-field) all-force-field-parts molecule)
  (chem:generate-standard-assign-molecular-force-field-parameters energy-function force-field molecule))


(defgeneric chem:generate-molecule-energy-function-tables (energy-function molecule combined-force-field active-atoms)
  (:documentation "Generate the molecule energy-function tables"))

(defmethod chem:generate-molecule-energy-function-tables (energy-function molecule (combined-force-field chem:combined-force-field) active-atoms)
  (chem:generate-standard-energy-function-tables energy-function molecule combined-force-field active-atoms))


;;; ------------------------------------------------------------
;;;
;;; node-table protocol
;;;
;;;


(defgeneric chem:make-node-table-from-graph (graph nonbond-force-field)
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


(defmethod chem:make-node-table-from-graph ((molecule chem:molecule) nonbond-force-field)
  "Construct a chem:atom-table from a molecule."
  (let ((atom-table (core:make-cxx-object 'chem:atom-table)))
    (chem:construct-from-molecule atom-table molecule nonbond-force-field nil)
    atom-table))

(defmethod chem:node-table-size ((atom-table chem:atom-table))
  (chem:get-number-of-atoms atom-table))

(defmethod chem:node-table-coordinate-index-times3 ((atom-table chem:atom-table) index)
  (chem:get-coordinate-index-for-atom-at-index atom-table index))

(defmethod chem:node-table-node-at-index ((atom-table chem:atom-table) index)
  (chem:elt-atom atom-table index))

(defmethod chem:node-get-position ((atom chem:atom))
  (chem:get-position atom))

(defmethod chem:node-set-position ((atom chem:atom) vector3)
  (chem:set-position atom vector3))

