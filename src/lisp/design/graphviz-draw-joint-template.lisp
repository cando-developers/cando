(in-package :design.graphviz-draw-joint-template)

(defparameter *table* nil)

(defun id (ast)
  (symbol-name (gethash ast *table*)))

(defun draw-joint-template (joint-template filename)
  (with-open-file (stream filename :direction :output :if-exists :supersede)
    (format stream "digraph G {~%   ordering = out; ~%")
    (let ((*table* (make-hash-table :test #'eq)))
      (stream-draw-joint-template joint-template stream))
    (format stream "}~%")))

(defmethod stream-draw-joint-template :around (joint-template stream)
  (when (null (gethash joint-template *table*))
    (setf (gethash joint-template *table*) (gensym))
    (format stream "  ~a [shape = box];~%"
            (id joint-template))
    (call-next-method)
    (loop for child in (kin:children joint-template)
          for i from 1
          do (stream-draw-joint-template child stream)
             (format stream "   ~a -> ~a [label = \"~d\"];~%"
                     (id joint-template) (id child) i))))

(defmethod label (joint-template)
  (let ((name (string (class-name (class-of joint-template)))))
    (format nil "~a/~a[~a]" name (kin:name joint-template) (kin:id joint-template))))

(defmethod stream-draw-joint-template ((joint-template kin:joint-template) stream)
  (format stream "   ~a [label = \"~a\"];~%" (id joint-template) (label joint-template)))
