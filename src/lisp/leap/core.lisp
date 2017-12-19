(in-package :leap.core)

;;; ------------------------------------------------------------
;;;
;;; Maintain a list of search paths for files
;;;

(defvar *path* nil)
(defvar *default-force-field* nil)

(defun clear-path ()
  (setf *path* nil))

(defun add-path (path)
  "Add a path to the *path* dynamic variable"
  (declare (pathname path))
  (push (translate-logical-pathname path) *path*))

(defun search-path (filename)
  "Search the *path* list for a directory that contains filename and return it.
  Otherwise return NIL."
  (declare (pathname filename))
  (loop for path in (reverse *path*)
     for full-path = (merge-pathnames filename path)
     when (probe-file full-path)
     return it))

;;; --
;;;
;;; Leap environment


(defclass leap-environment ()
  ((%functions :initarg  :functions
               :type     hash-table
               :reader   %functions
               :initform (make-hash-table :test #'equalp))
   (%variables :initarg  :variables
               :type     hash-table
               :reader   %variables
               :initform (make-hash-table :test #'equalp))))


(defparameter *leap-env* (make-instance 'leap-environment))

(defun function-lookup (name environment)
  (or (gethash name (%functions environment))
      (error "The function ~S is not defined" name)))

(defun (setf function-lookup) (new-value name environment)
  (setf (gethash name (%functions environment)) new-value))

(defun variable-lookup (name environment &optional errorp error-value)
  (or (gethash name (%variables environment))
   (if errorp
       (error "The variable ~S is not defined" name)
       error-value)))

(defun (setf variable-lookup) (new-value name environment)
  (setf (gethash name (%variables environment)) new-value))

(defun evaluate (builder ast environment)
  (architecture.builder-protocol:walk-nodes
   builder
   (lambda (recurse relation relation-args node kind relations
            &key name value &allow-other-keys)
     (declare (ignore relation relation-args node relations))
     (case kind
       (:literal
        value)
       (:list
        (first (funcall recurse :relations '(:element))))
       (:function
        (apply (function-lookup name environment)
               (first (funcall recurse :relations '(:argument)))))
       (:assignment
        (let ((value (first (first (funcall recurse :relations '(:value))))))
          (if (symbolp value)
              (setf (variable-lookup name environment) (lookup-variable value))
              (setf (variable-lookup name environment) value))))
       (t
        (funcall recurse))))

   ast))







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
  (setf (variable-lookup name *leap-env*) object))

(defun lookup-variable (name &optional (errorp t) error-value)
  "* Arguments
- name : Symbol
* Description
Lookup the object in the *objects*."
  (variable-lookup name *leap-env* errorp error-value))


(defun leap-lookup-variable-reader-macro (stream char)
  (let* ((*package* (find-package :keyword))
         (var-name (read stream t nil t)))
    `(lookup-variable ',var-name)))


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

(defvar *force-fields* (make-hash-table)
  "Every time a force field is loaded, it is put into a new list and
associated with a keyword symbol name.  As every frcmod is loaded it is pushed onto
a named force-field list.  When parameters are needed a name must be provided and
the list of frcmods and force-field are used to create a merged force-field that is
used to provide parameters.  There is one default force-field called :default.")

(defun add-force-field-or-modification (force-field-or-frcmod &optional (force-field-name :default))
  (push force-field-or-frcmod (gethash force-field-name *force-fields*)))

(defun merged-force-field (&optional (force-field-name :default))
  "Merge the force-field-list with _force-field-name_ and return it."
  (let* ((force-field-list (gethash force-field-name *force-fields*))
         (reversed-force-field-list (reverse force-field-list)))
    (if (= (length force-field-list) 1)
        (car force-field-list)
        (let ((merged-force-field (chem:make-force-field)))
          (dolist (ff reversed-force-field-list)
            (chem:force-field-merge merged-force-field ff))
          merged-force-field))))

