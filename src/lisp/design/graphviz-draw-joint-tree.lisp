

(in-package :design.graphviz-draw-joint-tree)

(defparameter *table* nil)

(defun id (ast)
  (symbol-name (gethash ast *table*)))

(defun draw-joint-tree (atom-tree filename)
  (with-open-file (stream filename :direction :output :if-exists :supersede)
    (format stream "digraph G {~%   ordering = out; ~%")
    (let ((*table* (make-hash-table :test #'eq)))
      (stream-draw-joint (kin:atom-tree-root atom-tree) atom-tree stream))
    (format stream "}~%")))

(defmethod stream-draw-joint :around (joint atom-tree stream)
  (when (null (gethash joint *table*))
    (setf (gethash joint *table*) (gensym))
    (format stream "  ~a [shape = box];~%"
            (id joint))
    (call-next-method)
    (loop for child in (kin:children joint)
          for i from 1
          do (stream-draw-joint child atom-tree stream)
             (format stream "   ~a -> ~a [label = \"~d\"];~%"
                     (id child) (id joint) i))))

(defmethod label ((joint kin:joint))
  (let ((name (string (class-name (class-of joint)))))
    (format nil "~a/~a" name (kin:name joint))))

(defmethod stream-draw-joint (joint atom-tree stream)
  (format stream "   ~a [label = \"~a\"];~%" (id joint) (label joint)))

(defmethod stream-draw-joint ((joint kin:joint) atom-tree stream)
  (format stream "   ~a [label = \"~a\"];~%" (id joint) (label joint)))

(defmethod stream-draw-joint ((joint kin:bonded-atom) atom-tree stream)
  (let* ((pos (kin:get-position joint))
         (start-pos (kin:get-property-or-default joint :start-pos nil))
         (highlight (if start-pos
                        (if (> (geom:vlength (geom:v- start-pos pos)) 0.01)
                            (format nil "style=filled;color=red;")
                            ""))))
    (format stream "   ~a [~alabel = \"~a~%~a~%#(~4,2f ~4,0f ~4,0f)~%~s~%~s~%~s\"];~%"
            (id joint)
            highlight
            (label joint)
            (kin:atom-id joint)
            (kin:get-distance joint)
            (/ (kin:get-theta joint) 0.0174533)
            (/ (kin:get-phi joint) 0.0174533)
            (kin:get-stub joint)
            (kin:get-position joint)
            (kin:properties joint))))

(defmethod stream-draw-joint ((joint kin:jump-atom) atom-tree stream)
  (format stream "   ~a [label = \"~a~%~a~%relative-transform~%~s~%lab-frame~%~s~%~s\"];~%"
          (id joint)
          (label joint)
          (kin:atom-id joint)
          (kin:get-parent-relative-transform joint)
          (kin:get-lab-frame joint)
          (kin:properties joint)))
