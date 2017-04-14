(cl:in-package #:leap-parser)

(defclass leap-environment ()
  ((%functions :initarg  :functions
               :type     hash-table
               :reader   %functions
               :initform (make-hash-table :test #'equal))
   (%variables :initarg  :variables
               :type     hash-table
               :reader   %variables
               :initform (make-hash-table :test #'equal))))

(defmethod lookup ((name        string)
                   (kind        (eql :function))
                   (environment leap-environment))
  (or (gethash name (%functions environment))
      (error "The function ~S is not defined" name)))

(defmethod (setf lookup) ((new-value   function)
                          (name        string)
                          (kind        (eql :function))
                          (environment leap-environment))
  (setf (gethash name (%functions environment)) new-value))

(defmethod lookup ((name        string)
                   (kind        (eql :variable))
                   (environment leap-environment))
  (or (gethash name (%variables environment))
      (error "The variable ~S is not defined" name)))

(defmethod (setf lookup) ((new-value   t)
                          (name        string)
                          (kind        (eql :variable))
                          (environment leap-environment))
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
        (apply (lookup name :function environment)
               (first (funcall recurse :relations '(:argument)))))
       (:assignment
        (setf (lookup name :variable environment)
              (first (first (funcall recurse :relations '(:value))))))
       (t
        (funcall recurse))))

   ast))

(defparameter *input*
  (alexandria:read-file-into-string
   (asdf:system-relative-pathname :leap-parser "examples/test.leap")))

(defparameter *ast*
  (architecture.builder-protocol:with-builder ('list)
    (esrap:parse 'leap *input*)))

(let ((env (make-instance 'leap-environment)))
  (setf (lookup "logFile" :function env)
        (lambda (filename)
          (print (list "logFile" filename)))
        (lookup "addAtomTypes" :function env)
        (lambda (types)
          (print (list "addAtomTypes" types))))
  (evaluate 'list *ast* env))
