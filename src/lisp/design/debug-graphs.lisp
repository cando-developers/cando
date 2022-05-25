(in-package :design.debug)


(defun property-label (atm)
  (with-output-to-string (sout)
    (let* ((props (chem:properties atm))
           (weight (getf props :weight))
           (bond1-atom (getf props :bond1-atom))
           (out-plug (getf props :out-plug))
           (root-atom (getf props :root-atom))
           (origin (getf props :origin))
           (high-priority (getf props :high-priority))
           (highest-priority (getf props :highest-priority))
           (build-before (getf props :build-before))
           (entity-to-delay-children-for (getf props :entity-to-delay-children-for))
           )
      (when (> weight 1)
        (format sout ":weight ~d~%" weight))
      (when highest-priority
        (format sout ":highest-priority t~%"))
      (when high-priority
        (format sout ":high-priority t~%"))
      (when out-plug
        (format sout ":out-plug ~a~%" out-plug))
      (when root-atom
        (format sout ":root-atom ~a~%" root-atom))
      (when origin
        (format sout ":origin ~a~%" origin))
      (when bond1-atom
        (format sout ":bond1-atom ~a~%" bond1-atom))
      (when build-before
        (format sout ":build-before ~a~%" build-before))
      (when entity-to-delay-children-for
        (format sout ":entity-to-delay-children-for ~a~%" entity-to-delay-children-for))
      )))


(defun residue-graph-nodes-edges (residue)
  (let ((nodes nil)
        (root-id nil)
        (atom-to-node (make-hash-table))
        (atom-to-node-id (make-hash-table)))
    (cando:do-atoms (atm residue)
      (let* ((node-id (gensym))
             (label (chem:get-name atm))
             (node (cytoscape:make-element
                    :group "nodes"
                    :data (j:make-object "id" (string node-id) "label" (format nil "~a~%~a" (string label) (property-label atm))))))
        (push node nodes)
        (unless root-id
          (setf root-id (string node-id)))
        (setf (gethash atm atom-to-node) node
              (gethash atm atom-to-node-id) node-id)))
    (let ((edges nil))
      (chem:map-bonds 'nil
                      (lambda (atm1 atm2 bond-order)
                        (declare (ignore bond-order))
                        (let* ((node1 (gethash atm1 atom-to-node-id))
                               (children1 (getf (chem:properties atm1) :children))
                               (node2 (gethash atm2 atom-to-node-id))
                               (children2 (getf (chem:properties atm2) :children))
                               (node2-child-of-node1 (member atm2 children1))
                               (edge (when (and node1 node2)
                                       (if node2-child-of-node1
                                           (let ((child-index (position atm2 children1)))
                                             (cytoscape:make-element
                                              :group "edges"
                                              :data (j:make-object "source" (string node2)
                                                                   "target" (string node1)
                                                                   "label" (format nil " ~a" child-index))))
                                           (let ((child-index (position atm1 children2)))
                                             (cytoscape:make-element
                                              :group "edges"
                                              :data (j:make-object "source" (string node1)
                                                                   "target" (string node2)
                                                                   "label" (format nil " ~a" child-index))))))))
                          (when edge (push edge edges))))
                      residue)
      #+(or)(progn
        (format t "Number of nodes: ~a  edges: ~s~%" (length nodes) (length edges))
        (loop for node in nodes
              do (format t "node: ~s~%" (alexandria:hash-table-plist (cytoscape:data node))))
        (loop for edge in edges
              do (format t "edge: ~s~%" (alexandria:hash-table-plist (cytoscape:data edge)))))
      (values nodes edges root-id))))

(defun residue-graph (residue)
  #+(or)(format t "Drawing graph ~a~%" residue)
  (multiple-value-bind (nodes edges root-id)
      (residue-graph-nodes-edges residue)
    (make-instance 'cytoscape:cytoscape-widget
                   :graph-layouts (list (make-instance 'cytoscape:cose-layout
                                                       :roots (list root-id)))
                   :graph-style "* { label: data(label); font-size: 10;text-wrap: wrap; } edge { curve-style: straight; target-arrow-shape: triangle; }"
                   :elements (append (nreverse nodes) (nreverse edges))
                   :layout (make-instance 'jupyter-widgets:layout :width "auto" :height "1000px"))))






(defun matter-graph-nodes-edges (matter)
  (let ((nodes nil)
        (root-id nil)
        (atom-to-node (make-hash-table))
        (atom-to-node-id (make-hash-table)))
    (cando:do-atoms (atm matter)
      (let* ((node-id (gensym))
             (label (chem:get-name atm))
             (node (cytoscape:make-element
                    :group "nodes"
                    :data (j:make-object "id" (string node-id) "label" (format nil "~a~%~a" (string label) (chem:get-type atm))))))
        (push node nodes)
        (unless root-id
          (setf root-id (string node-id)))
        (setf (gethash atm atom-to-node) node
              (gethash atm atom-to-node-id) node-id)))
    (let ((edges nil))
      (chem:map-bonds 'nil
                      (lambda (atm1 atm2 bond-order)
                        (declare (ignore bond-order))
                        (let* ((node1 (gethash atm1 atom-to-node-id))
                               (children1 (getf (chem:properties atm1) :children))
                               (node2 (gethash atm2 atom-to-node-id))
                               (children2 (getf (chem:properties atm2) :children))
                               (node2-child-of-node1 (member atm2 children1))
                               (edge (when (and node1 node2)
                                       (if node2-child-of-node1
                                           (let ((child-index (position atm2 children1)))
                                             (cytoscape:make-element
                                              :group "edges"
                                              :data (j:make-object "source" (string node2)
                                                                   "target" (string node1)
                                                                   "label" (format nil " ~a" child-index))))
                                           (let ((child-index (position atm1 children2)))
                                             (cytoscape:make-element
                                              :group "edges"
                                              :data (j:make-object "source" (string node1)
                                                                   "target" (string node2)
                                                                   "label" (format nil " ~a" child-index))))))))
                          (when edge (push edge edges))))
                      matter)
      #+(or)(progn
        (format t "Number of nodes: ~a  edges: ~s~%" (length nodes) (length edges))
        (loop for node in nodes
              do (format t "node: ~s~%" (alexandria:hash-table-plist (cytoscape:data node))))
        (loop for edge in edges
              do (format t "edge: ~s~%" (alexandria:hash-table-plist (cytoscape:data edge)))))
      (values nodes edges root-id))))

(defun matter-graph (matter)
  #+(or)(format t "Drawing graph ~a~%" matter)
  (multiple-value-bind (nodes edges root-id)
      (matter-graph-nodes-edges matter)
    (make-instance 'cytoscape:cytoscape-widget
                   :graph-layouts (list (make-instance 'cytoscape:cose-layout
                                                       :roots (list root-id)))
                   :graph-style "* { label: data(label); font-size: 10;text-wrap: wrap; } edge { curve-style: straight; target-arrow-shape: triangle; }"
                   :elements (append (nreverse nodes) (nreverse edges))
                   :layout (make-instance 'jupyter-widgets:layout :width "auto" :height "1000px"))))



(defgeneric label (joint-template node-id))

(defmethod label (joint-template node-id)
  (format nil "~a~%~a" (class-name (class-of joint-template)) (design.joint-tree:atom-name joint-template)))

(defmethod label ((joint-template design.joint-tree:bonded-joint-template) node-id)
  (format nil "BONDJT~%~a" (design.joint-tree:atom-name joint-template)))

(defmethod label ((joint-template design.joint-tree:jump-joint-template) node-id)
  (format nil "JUMP~%~a" (design.joint-tree:atom-name joint-template)))

(defmethod label ((joint-template design.joint-tree:in-plug-bonded-joint-template) node-id)
  (format nil "IN~%~a~%~a" (design.joint-tree:atom-name joint-template) (design.joint-tree:in-plug joint-template)))

(defmethod label ((joint-template design.joint-tree:complex-bonded-joint-template) node-id)
  (let* ((ajoint (design.joint-tree:parent joint-template))
         (ajoint-name (when ajoint (string (design.joint-tree:atom-name ajoint))))
         (input-stub-joints (design.joint-tree:input-stub-joints joint-template))
         (bjoint (aref input-stub-joints 0))
         (bjoint-name (when bjoint (string (design.joint-tree:atom-name bjoint))))
         (cjoint (aref input-stub-joints 1))
         (cjoint-name (when cjoint (string (design.joint-tree:atom-name cjoint)))))
  (format nil "CBJT~%~a~%~a" (design.joint-tree:atom-name joint-template)
          (list ajoint-name bjoint-name cjoint-name))))



(defun render-joint-template (parent joint-template child-index nodes-edges joint-template-to-node-id)
  (when (null (gethash joint-template joint-template-to-node-id))
    (let* ((node-id (gensym "AA"))
           (node (cytoscape:make-element
                  :group "nodes"
                  :data (j:make-object "id" (string node-id) "label" (label joint-template node-id)))))
      (setf (gethash joint-template joint-template-to-node-id) node-id)
      (push node (car nodes-edges))
      (let ((edge (when parent
                    (let ((node1 (gethash parent joint-template-to-node-id))
                          (node2 (gethash joint-template joint-template-to-node-id))
                          )
                      (cytoscape:make-element :group "edges"
                                              :data (j:make-object "source" (string node2)
                                                                   "target" (string node1)
                                                                   "label" (format nil "~d" child-index)))))))
        (when edge
          (push edge (cdr nodes-edges)))
        (loop for child in (design.joint-tree:children joint-template)
              for sub-child-index from 0
              do (render-joint-template joint-template child sub-child-index nodes-edges joint-template-to-node-id))
        node-id))))

(defun do-render-joint-template (joint-template)
  (let* ((nodes-edges (cons nil nil))
         (joint-template-to-node-id (make-hash-table))
         (root-id (render-joint-template nil joint-template 0 nodes-edges joint-template-to-node-id))
         (nodes (car nodes-edges))
         (edges (cdr nodes-edges)))
    #+(or)(progn
            (format t "root-id: ~a~%" root-id)
            (loop for node in nodes
                  do (format t "node: ~s~%" (alexandria:hash-table-plist (cytoscape:data node))))
            (loop for edge in edges
                  do (format t "edge: ~s~%" (alexandria:hash-table-plist (cytoscape:data edge)))))
    (values nodes edges root-id)))

(defun joint-template-graph (joint-template)
  (multiple-value-bind (nodes edges root-id)
      (do-render-joint-template joint-template)
    #+(or)(format t "Number of nodes: ~a  edges: ~s root-id: ~s~%" (length nodes) (length edges) (string root-id))
    (make-instance 'cytoscape:cytoscape-widget
                   :graph-layouts (list (make-instance 'cytoscape:cose-layout
                                                       :roots (list (string root-id))))
                   :graph-style "* { label: data(label); font-size: 10;text-wrap: wrap; } edge { curve-style: straight; target-arrow-shape: triangle; }"
                   :elements (append (reverse nodes) (reverse edges))
                   :layout (make-instance 'jupyter-widgets:layout :width "auto" :height "1000px"))))
