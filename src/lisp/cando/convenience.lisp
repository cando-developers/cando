(in-package :cando)

(defgeneric agg (object)
  (:documentation "Return the aggregate that the object represents"))

(defgeneric mol (object id)
  (:documentation "Return the molecule that the object/id represents"))

(defgeneric res (object id)
  (:documentation "Return the residue that the object/id represents"))

(defgeneric atm (object id)
  (:documentation "Return the atom that the object/id represents"))


(defmethod mol ((object chem:aggregate) (id integer))
  (chem:content-at object id))

(defmethod res ((object chem:molecule) (id integer))
  (chem:content-at object id))

(defmethod atm ((object chem:residue) (id integer))
  (chem:content-at object id))

(defmethod molid ((object chem:aggregate) (id integer))
  (chem:content-with-id object id))

(defmethod resid ((object chem:molecule) (id integer))
  (chem:content-with-id object id))

(defmethod atmid ((object chem:residue) (id integer))
  (chem:content-with-id object id))




