(in-package :graphviz-molecule-graph)

(defmethod draw-node (id index (node mg:atom-cluster) stream)
  (let* ((width 0)
         (name (with-output-to-string (sout)
                 (maphash (lambda (atom dummy)
                            (declare (ignore dummy))
                            (let ((name (string (chem:get-name atom))))
                              (format sout "~a " name)
                              (if (> width 10)
                                  (progn
                                    (terpri sout)
                                    (setf width 0))
                                  (incf width (1+ (length name))))))
                          (mg:atom-cluster-atoms node)))))
    (format stream "~a [label = \"~a: ~a\";shape=box];~%" id index name)))

(defmethod draw-node (id index (node chem:atom) stream)
  (format stream "~a [label = \"~a: ~a\"];~%" id index (chem:get-name node)))

(defmethod draw-edge (source-id target-id label stream)
  (format stream "~a -- ~a [label = \"~a\"];~%" source-id target-id (string label)))


(defun draw-graph-to-stream (graph stream)
  (let ((ids (make-hash-table))
        (index 0))
    (chem:walk-vertices
     graph
     (lambda (node)
       (let ((id (gensym)))
         (setf (gethash node ids) id)
         (draw-node id index node stream)
         (incf index))))
    (chem:walk-edges
     graph
     (lambda (v1 v2 bond-order)
       (let ((inode1 (gethash v1 ids))
             (inode2 (gethash v2 ids)))
         (draw-edge inode1 inode2 bond-order stream))))))


(defun draw-graph-to-file (graph filename)
  (with-open-file (fout filename :direction :output :if-exists :supersede)
    (format fout "graph {~%")
    (draw-graph-to-stream graph fout)
    (format fout "}~%")))


