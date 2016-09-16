(in-package :amber)

(defparameter *objects* (make-hash-table :test #'equal))

(define-condition object-not-found ()
  ((name :initarg :name :reader name)))

(defun register-variable (name object)
  "* Arguments
- name :: Symbol
- object :: an object
* Description
Associate the name with the object"
  (when (gethash name *objects*)
    (warn "Overwriting variable ~a" name))
  (setf (gethash name *objects*) object))

(defun lookup-variable (name &optional (errorp t) error-value)
  "* Arguments
- name : Symbol
* Description
Lookup the object in the *objects*."
  (multiple-value-bind (val found)
      (gethash name *objects*)
    (if found
        val
        (if errorp
            (error 'object-not-found :name name)
            error-value))))

