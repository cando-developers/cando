(in-package :lomap)



(defun edge-in-spanning-tree-p (edge back-span-info)
  ;; Check the first vertex of edge
  (multiple-value-bind (vertex1 vertex2)
      (edge-vertices edge)
    (let* ((vertex vertex1)
           (bsi (gethash vertex back-span-info)))
      (when bsi
        (let ((back-vertex (spanning:back-vertex bsi)))
          (when (eq back-vertex vertex2)
            (return-from edge-in-spanning-tree-p t)))))
    (let* ((vertex vertex2)
           (bsi (gethash vertex back-span-info)))
      (when bsi
        (let ((back-vertex (spanning:back-vertex bsi)))
          (when (eq back-vertex vertex1)
            (return-from edge-in-spanning-tree-p t)))))
    nil))

(defun edges-outside-of-spanning-tree (graph back-span-info)
  (loop for edge in (edges graph)
        when (not (edge-in-spanning-tree-p edge back-span-info))
          collect edge))


(defun generate-backspan-bitvec (vertex back-span-info graph)
  (let* ((num-vertices (length (vertices graph)))
        (bitvec (make-array (bitvec-length num-vertices) :element-type 'bit :adjustable nil :initial-element 0)))
    (loop named walk
          for vertex-index = (index vertex)
          for bsi = (gethash vertex back-span-info)
          do (if bsi
               (let* ((back-vertex (spanning:back-vertex bsi))
                      (back-vertex-index (index back-vertex))
                      (bit-index (index-to-bit vertex-index back-vertex-index num-vertices)))
                 (setf (elt bitvec bit-index) 1)
                 (setf vertex back-vertex))
               (return-from generate-backspan-bitvec bitvec)))))

(defun all-nodes-in-fundamental-cycles-p (graph spanning-tree)
  ;;question: to do check might need to geberate spanning tree here instead of pass one
  ;;becuase this will be a new "chopped graph" when checking constraints.
  (let* ((outedges (edges-outside-of-spanning-tree graph spanning-tree))
         (num-vertices (length (vertices graph)))
         (all-cycles-bitvec (make-array (bitvec-length num-vertices)
                                        :element-type 'bit
                                        :adjustable nil
                                        :initial-element 0)))
    (loop for outedge-i in outedges
          do (let* ((s1 (generate-backspan-bitvec (vertex1 outedge-i) spanning-tree graph))
                    (s2 (generate-backspan-bitvec (vertex2 outedge-i) spanning-tree graph))
                    (s12 (bit-xor s1 s2)))
                    (setf (elt s12 (index-to-bit (index (vertex1 outedge-i)) (index (vertex2 outedge-i)) num-vertices)) 1)
               (setf all-cycles-bitvec (bit-ior s12 all-cycles-bitvec))))
    #+debug-lomap(format t "all-cycles-bitvec: ~a~%" all-cycles-bitvec)
  (let ((nodes-in-cycles (make-array num-vertices :element-type 't :adjustable nil :initial-element nil)))
    (loop for bit-index from 0 below (length all-cycles-bitvec)
          for bit-val = (elt all-cycles-bitvec bit-index)
          do (when (= bit-val 1)
               (multiple-value-bind (node1-index node2-index)
                   (bit-to-index bit-index num-vertices)
                 (let ((vertex1 (elt (vertices graph) node1-index))
                       (vertex2 (elt (vertices graph) node2-index)))
                   #+debug-lomap
                   (progn
                     (format t "edge ~a - ~a~%" vertex1 vertex2)
                     (format t "Node-in-cycle ~a ~a~%" node1-index vertex1)
                     (format t "Node-in-cycle ~a ~a~%" node2-index vertex2)))
                 (setf (elt nodes-in-cycles node1-index) t
                       (elt nodes-in-cycles node2-index) t))))
    (let ((all-in (every #'identity nodes-in-cycles)))
      #+debug-lomap(format t "nodes-in-cycles -> ~a    all-in -> ~a~%" nodes-in-cycles all-in)
      all-in))))

(defun graph-wider-than-p (graph maxdist)
  "Return T if the graph-diameter is less than or equal to maxdist"
  (loop for vertex in (vertices graph)
        do (multiple-value-bind (spanning-tree width)
               (spanning:calculate-spanning-tree graph vertex)
             (when (> width maxdist)
               (return-from graph-wider-than-p T))))
  nil)


(defun connected-components (graph)
  "Return a list of lists of vertices.  Each list of vertices is a connected component in the graph"
  (let ((vertices (copy-list (vertices graph)))
        all-components)
    (loop
      do (let* ((one-vertex (first vertices))
                (spanning-tree (spanning:calculate-spanning-tree graph one-vertex))
                component
                untouched-vertices)
           (loop for vertex in vertices
                 do (multiple-value-bind (back-span foundp)
                        (gethash vertex spanning-tree)
                      (if foundp
                          (push vertex component)
                          (push vertex untouched-vertices))))
           (push component all-components)
           (when (null untouched-vertices)
             (return-from connected-components all-components))
           (setf vertices untouched-vertices)))))

(defun number-connected-components (graph)
  (length (connected-components graph)))
