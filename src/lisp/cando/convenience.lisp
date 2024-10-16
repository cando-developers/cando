(in-package :cando)

(defgeneric agg (object)
  (:documentation "Return the aggregate that the object represents"))

(defgeneric mol (object &optional id)
  (:documentation "Return the molecule that the object/id represents"))

(defgeneric res (object id)
  (:documentation "Return the residue that the object/id represents"))

(defgeneric atm (object id)
  (:documentation "Return the atom that the object/id represents"))

(defmethod agg ((object chem:aggregate))
  object)

(defmethod mol ((object chem:aggregate) &optional id)
  (etypecase id
    (integer (chem:content-at object id))
    (null (if (= (chem:content-size object) 1)
              (chem:content-at object 0)
              (error "If no index is provided to mol there can only be one molecule and there are ~a molecules" (chem:content-size object))))))

(defmethod res ((object chem:molecule) (id integer))
  (chem:content-at object id))

(defmethod atm ((object chem:residue) (id integer))
  (chem:content-at object id))

(defmethod atm ((object chem:residue) (id symbol))
  (chem:content-with-name object id))

(defmethod molid ((object chem:aggregate) (id integer))
  (chem:content-with-id object id))

(defmethod resid ((object chem:molecule) (id integer))
  (chem:content-with-id object id))

(defmethod atmid ((object chem:residue) (id integer))
  (chem:content-with-id object id))




