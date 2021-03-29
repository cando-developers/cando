
(in-package :smarts-jupyter)


(defun graph (smarts)
  (let ((nodes nil)
        (edges nil)
        (root-id nil)
        (chem-info-to-node (make-hash-table))
        (chem-info-to-node-id (make-hash-table)))
    (chem:walk-chem-info-with-parent
     smarts
     (lambda (parent child)
       (let* ((parent-node (gethash parent chem-info-to-node))
              (parent-node-id (gethash parent chem-info-to-node-id))
              (node-id (gensym))
              (label (cond
                       ((typep child 'chem:atom-test)
                        (format nil "~s" child))
                       (t (format nil "~s" child))))
              (node (make-instance 'cytoscape:element
                                   :group "nodes"
                                   :data (j:make-object "id" (string node-id)
                                                        "label" label))))
         (unless root-id
           (setf root-id (string node-id)))
         (setf (gethash child chem-info-to-node) node
               (gethash child chem-info-to-node-id) node-id)
         (push node nodes)
         (when parent-node-id
           (let ((edge (make-instance 'cytoscape:element
                                      :group "edges"
                                      :data (j:make-object "source" (string parent-node-id)
                                                           "target" (string node-id)))))
             (push edge edges))))))
    (make-instance 'cytoscape:cytoscape-widget
                   :graph-layouts (list (make-instance 'cytoscape:dagre-layout
                                                       :padding -30
                                                       :roots (list root-id)))
                   :graph-style "* { label: data(label); font-size: 10; }"
                   :elements (append (nreverse nodes) (nreverse edges)))))

