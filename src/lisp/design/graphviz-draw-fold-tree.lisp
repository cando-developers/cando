
(in-package :design.graphviz-draw-fold-tree)

(defparameter *table* nil)

(defun graphid (ast)
  (symbol-name (gethash ast *table*)))

(defun draw-fold-tree (fold-tree filename)
  (with-open-file (stream filename :direction :output :if-exists :supersede)
    (format stream "digraph G {~%   ordering = out; ~%")
    (let ((*table* (make-hash-table :test #'eq)))
      (stream-draw-fold (kin:root fold-tree) fold-tree stream))
    (format stream "}~%")))

(defmethod stream-draw-fold :around (fold fold-tree stream)
  (when (null (gethash fold *table*))
    (setf (gethash fold *table*) (gensym))
    (format stream "  ~a [shape = box];~%"
            (graphid fold))
    (call-next-method)
    (loop for child in (kin:children fold)
          for i from 1
          do (stream-draw-fold child fold-tree stream)
             (format stream "   ~a -> ~a [label = \"~d\"];~%"
                     (graphid child) (graphid fold) i))))

(defmethod label ((fold kin:monomer-base-node))
  (let ((name (string (class-name (class-of fold)))))
    (format nil "~a" name)))

(defmethod stream-draw-fold (fold fold-tree stream)
  (format stream "   ~a [label = \"~a\"];~%" (graphid fold) (label fold)))

(defmethod stream-draw-fold ((fold kin:monomer-node) fold-tree stream)
  (format stream "   ~a [label = \"~a~%~a\"];~%" (graphid fold) (label fold) (kin:stereoisomer-name fold)))

