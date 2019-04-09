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

(defun lookup-variable* (name environment &optional errorp error-value)
  "Lookup a variable in the given environment"
  (let ((topology (cando:lookup-topology name nil)))
    (if topology
        topology
        (or (gethash name (%variables environment))
            (if errorp
                (error "The variable ~S is not defined" name)
                error-value)))))

(defun all-variables ()
  (let (vars)
    (maphash (lambda (k v) (push k vars)) (%variables *leap-env*))
    (flet ((accumulate-topologys (top)
             (push (chem:get-name top) vars)))
      (cando:walk-topologys #'accumulate-topologys))
    (let ((sorted-vars (sort vars #'string<)))
      sorted-vars)))

(defun (setf lookup-variable*) (new-value name environment)
  (if (typep new-value 'chem:topology)
      (cando:register-topology name new-value)
      (setf (gethash name (%variables environment)) new-value)))

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
              (setf (lookup-variable* name environment) (lookup-variable value))
              (setf (lookup-variable* name environment) (if (typep value 'number)
                                                            (float value 1d0)
                                                            value)))))
       (t
        (funcall recurse))))
   ast))







;;; ------------------------------------------------------------
;;;
;;; Maintain a hash-table that maps object names to objects
;;; 
(define-condition object-not-found ()
  ((name :initarg :name :reader name)))

(defun register-variable (name object)
  "* Arguments
- name :: Symbol
- object :: an object
* Description
Associate the name with the object"
  (setf (lookup-variable* name *leap-env*) object))

(defun lookup-variable (name &optional (errorp t) error-value)
  "* Arguments
- name : Symbol
* Description
Lookup the object in the variable space."
  (lookup-variable* name *leap-env* errorp error-value))


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

(defun clear-force-field (&optional (force-field-name :default))
  (chem:clear-combined-force-field (gethash force-field-name *force-fields*))
  ;;; The default conbined-force-field always has one empty force-field for global nonbond info
  (when (eq force-field-name :default)
    (add-force-field-or-modification (chem:make-force-field))))

(defun add-combined-force-field (combined-force-field force-field-name)
  (setf (gethash force-field-name *force-fields*) combined-force-field))

(defun chem:find-force-field (force-field-name &optional (errorp t))
  (let ((ff (gethash force-field-name *force-fields*)))
    (when (and errorp (null ff))
      (error "Could not find force-field with name ~s" force-field-name))
    ff))

(defun add-force-field-or-modification (force-field-or-frcmod &key (force-field-name :default) force-field-info
                                                                combined-force-field-class-name)
  (let ((cff (chem:find-force-field force-field-name nil)))
    (unless cff
      (add-combined-force-field (if (eq combined-force-field-class-name 'chem:combined-force-field)
                                              (sys:make-cxx-object combined-force-field-class-name)
                                              (make-instance combined-force-field-class-name))
                                          force-field-name)))
  (let ((combined-force-field (chem:find-force-field force-field-name)))
    (chem:add-shadowing-force-field combined-force-field force-field-or-frcmod force-field-info))
  force-field-or-frcmod)

(eval-when (:load-toplevel :execute)
  (add-force-field-or-modification (chem:make-force-field) :force-field-name :default
                                   :combined-force-field-class-name 'chem:combined-force-field))

(defun merged-force-field (&optional (force-field-name :default))
  "Merge the force-field-list with _force-field-name_ and return it."
  (error "Deprecated merged-force-field - use chem:find-force-field and work with the result")
  (let* ((force-field-list (gethash force-field-name *force-fields*))
         (reversed-force-field-list (reverse force-field-list)))
    (if (= (length force-field-list) 1)
        (car force-field-list)
        (let ((merged-force-field (chem:make-force-field)))
          (dolist (ff reversed-force-field-list)
            (chem:force-field-merge merged-force-field ff))
          merged-force-field))))


(defun nonbond-force-field-component (force-field-name)
  (let* ((combined-force-field (chem:find-force-field force-field-name))
         (force-field-parts (chem:force-fields-as-list combined-force-field))
         (nbmerged (chem:make-ffnonbond-db)))
    (mapc (lambda (force-field-part)
            (let ((nonbond-part (chem:nonbond-component force-field-part)))
              (chem:force-field-component-merge nbmerged nonbond-part)))
          force-field-parts)
    nbmerged))

(defun force-fields ()
  (maphash (lambda (name parts)
             (format t "~a ~a~%" name parts))
           *force-fields*))

(defun ensure-path (obj)
  "Convert obj into a pathname and search the amber path list for the file and return it"
  (declare (type (or symbol pathname string) obj))
  (let ((pn (cond
              ((symbolp obj) (pathname (string obj)))
              ((pathname obj) obj)
              ((stringp obj) (pathname obj))
              (t (pathname (string obj))))))
    (leap.core:search-path pn)))

