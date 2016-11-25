(in-package :amber)

(defparameter *path* nil)

(defun add-path (path)
  "Add a path to the *path* dynamic variable"
  (declare (pathname path))
  (push path *path*))

(defun search-path (filename)
  "Search the *path* list for a directory that contains filename and return it.  Otherwise return NIL."
  (declare (pathname filename))
  (loop for path in (reverse *path*)
     for full-path = (merge-pathnames filename path)
     when (probe-file full-path)
     return it))
       
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

