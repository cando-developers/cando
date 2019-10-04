(in-package :leap.atom-types)


(defclass type-element-info ()
  ((atom-type :initarg :atom-type :accessor atom-type)
   (element :initarg :element :accessor element)
   (hybridization :initarg :hybridization :accessor hybridization)))


(defvar *type-element-info* (make-hash-table))

(defun add-one-atom-type (atom-type element hybridization)
  (let ((info (make-instance 'type-element-info
                             :atom-type atom-type
                             :element element
                             :hybridization hybridization)))
    (setf (gethash atom-type *type-element-info*) info)))

