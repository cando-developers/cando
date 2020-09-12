(in-package :spanning)

(defclass back-span ()
  ((distance-to-root :initarg :distance-to-root :accessor distance-to-root)
   (back-vertex :initarg :back-vertex :accessor back-vertex)))

(defparameter *debug-spanning-tree* nil)
(defparameter *largest-distance* nil)


(defgeneric vertex-edges (vertex graph)
  (:documentation "Return the list of edges that connect the vertex to other vertices in the graph"))

(defgeneric other-vertex (edge vertex)
  (:documentation "Given an edge and a vertex, return the other vertex in the edge"))

(defun span-from-node (vertex distance back-span-info graph)
  (when *debug-spanning-tree*
    (format t "On vertex: ~a~%" vertex))
  (setf distance (1+ distance))
  (setf *largest-distance* (max *largest-distance* distance))
  ;; First visit all the other edges and create a back-span to the current vertex
  (let ((other-vertices nil))
    (loop for edge in (vertex-edges vertex graph)
          for other-vertex = (other-vertex edge vertex)
          do (multiple-value-bind (value foundp)
                 (gethash other-vertex back-span-info)
               (unless foundp
                 (when *debug-spanning-tree*
                   (format t "   child: ~a~%" other-vertex))
                 (let ((back-span (make-instance 'back-span
                                                 :distance-to-root distance
                                                 :back-vertex vertex)))
                   (setf (gethash other-vertex back-span-info) back-span))
                 (push other-vertex other-vertices))))
    ;; then loop over the other vertices and visit them
    (loop for other-vertex in other-vertices
          do (span-from-node other-vertex distance back-span-info graph))))

(defun calculate-spanning-tree (graph root-vertex &key debug)
  "Return a hash-table of back-span records that describe a spanning tree from the root.
Also return a second value that is the longest path in the spanning tree from the root."
  (let ((*debug-spanning-tree* debug)
        (*largest-distance* 0))
    (let ((back-span-info (make-hash-table)))
      (setf (gethash root-vertex back-span-info) nil)
      (span-from-node root-vertex 0 back-span-info graph)
      (values back-span-info *largest-distance*))))

