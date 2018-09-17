

(in-package :design.graphviz-draw-joint-tree)

(defparameter *table* nil)

(defun id (ast)
  (symbol-name (gethash ast *table*)))

(defun draw-joint-tree (joint filename)
  (with-open-file (stream filename :direction :output :if-exists :supersede)
    (format stream "digraph G {~%   ordering = out; ~%")
    (let ((*table* (make-hash-table :test #'eq)))
      (stream-draw-joint (kin:atom-tree-root joint) stream))
    (format stream "}~%")))

(defmethod stream-draw-joint :around (joint stream)
  (when (null (gethash joint *table*))
    (setf (gethash joint *table*) (gensym))
    (format stream "  ~a [shape = box];~%"
            (id joint))
    (call-next-method)
    (loop for child in (kin:children joint)
          for i from 1
          do (stream-draw-joint child stream)
             (format stream "   ~a -> ~a [label = \"~d\"];~%"
                     (id child) (id joint) i))))

(defmethod label ((joint kin:joint))
  (let ((name (string (class-name (class-of joint)))))
    (format nil "~a/~a" name (kin:name joint))))

(defmethod stream-draw-joint (joint stream)
  (format stream "   ~a [label = \"~a\"];~%" (id joint) (label joint)))

(defmethod stream-draw-joint ((joint kin:joint) stream)
  (format stream "   ~a [label = \"~a\"];~%" (id joint) (label joint)))
