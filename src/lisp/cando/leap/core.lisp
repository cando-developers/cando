(in-package :leap.core)

;;; ------------------------------------------------------------
;;;
;;; Maintain a list of search paths for files
;;;

(defvar *path* nil)
(defvar *default-force-field* nil)

(defun add-path (path)
  "Add a path to the *path* dynamic variable"
  (declare (pathname path))
  (push path *path*))

(defun search-path (filename)
  "Search the *path* list for a directory that contains filename and return it.
  Otherwise return NIL."
  (declare (pathname filename))
  (loop for path in (reverse *path*)
     for full-path = (merge-pathnames filename path)
     when (probe-file full-path)
     return it))

;;; ------------------------------------------------------------
;;;
;;; Maintain a hash-table that maps object names to objects
;;; 
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


;;; ------------------------------------------------------------
;;;
;;; Maintain a hash-table of named force-field-lists that are merged when
;;;    they are needed.
;;;
;;; The environment supports a collection of named force-fields that are
;;; assembled by merging a list of a force-field with any number of force-field modifications (frcmods).
;;; The force-fields and modifications are maintained as a list and
;;; only merged when the merged force-field is needed.
;;; There is a default merged force-field called :default.

(defvar *force-field-lists* (make-hash-table)
  "Every time a force field is loaded, it is put into a new list and
associated with a keyword symbol name.  Every frcmod is loaded it is pushed onto
a named force-field list.  When parameters are needed a name must be provided and
the list of frcmods and force-field are used to create a merged force-field that is
used to provide parameters.  There is one default force-field called :default.")

(defun add-force-field-or-modification (force-field-or-frcmod &optional (force-field-name :default))
  (push force-field-or-frcmod (gethash force-field-name *force-field-lists*)))

(defun merged-force-field (&optional (force-field-name :default))
  "Merge the force-field-list with _force-field-name_ and return it."
  (let ((force-field-list (gethash force-field-name *force-field-lists*))
        (reversed-force-field-list (reverse *force-field-list*))
        (merged-force-field (chem:make-force-field)))
    (dolist (ff reversed-force-field-list)
      (chem:force-field-merge merged-force-field ff))
    merged-force-field))

(defun coerce-to-force-field (force-field-designator)
  "force-field-designator can be a force-field or a list of force-fields. In the case of a list of force-fields they will be merged and the merged force-field will be returned."
  (cond
    ((typeof force-field-designator 'chem:force-field)
     force-field-designator)
    ((and (consp force-field-designator)
          (every (lambda (x) (typeof x 'chem:force-field))))
     (let ((merged-force-field (chem:make-force-field))
           (reversed-force-field-list (reverse force-field-designator)))
       (dolist (ff reversed-force-field-list)
         (chem:force-field-merge merged-force-field ff))
       merged-force-field))
    (error "Unsupported force-field designator ~a" force-field-designator)))
      
