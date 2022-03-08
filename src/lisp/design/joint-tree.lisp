
(in-package :kin)

(defgeneric children (object)
  (:documentation "Return list of children of object"))

(defmethod children ((joint joint))
  (kin:joint-children joint))

(defclass joint-tree ()
  ((root :initarg :root :accessor root)
   (atom-id-map :initarg :atom-id-map
                :initform (chem:make-atom-id-map)
                :accessor atom-id-map)))

(defun make-joint-tree ()
  (make-instance 'joint-tree))

(eval-when (:compile-toplevel :execute :load-toplevel)
  (export '(joint-tree children) :kin)
  )


