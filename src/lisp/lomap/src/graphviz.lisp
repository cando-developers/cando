(in-package :lomap)

(defparameter *neato-scale* 0.04)

(defun draw-graph-to-stream (graph stream &key edge-bitvec back-span-info extra-edges)
  ;;maybe change back-span-info to spanning-tree to make it more readable.
  (format stream "digraph goofy {~%")
  (let ((use-neato (some #'xypos (vertices graph))))
    (format stream "layout=\"neato\"~%"))
  (let ((ids (make-hash-table))
        (num-vertices (length (vertices graph))))
    (loop for vertex in (vertices graph)
          do (let ((id (gensym)))
               (setf (gethash vertex ids) id)
               (if (xypos vertex)
                 (let* ((pos (xypos vertex))
                        (xpos (* *neato-scale* (first pos)))
                        (ypos (* *neato-scale* (second pos))))
                   (format stream "~a [label = \"~a\", pos=\"~8,4f,~8,4f!\"];~%" id (chem:get-name (molecule vertex)) xpos ypos))
                 (format stream "~a [label = \"~a\"];~%" id (chem:get-name (molecule vertex))))
               ))
    (flet ((draw-edge (edge color)
             (let* ((id1 (gethash (vertex1 edge) ids))
                    (id2 (gethash (vertex2 edge) ids))
                    (weight (if edge-bitvec
                                (let* ((index1 (index (vertex1 edge)))
                                       (index2 (index (vertex2 edge)))
                                       (bit-index (index-to-bit index1 index2 num-vertices))
                                       (bit (elt edge-bitvec bit-index)))
                                  (if (= bit 0) 1 3))
                                1)))
               ;; Below here id1 and id2 may be swapped
               (let* ((edge-symbol (if back-span-info
                                       (cond
                                         ((let ((back-span (gethash (vertex1 edge) back-span-info)))
                                            (if (and back-span (eq (back-vertex back-span) (vertex2 edge)))
                                                "->"
                                                nil)))
                                         ((let ((back-span (gethash (vertex2 edge) back-span-info)))
                                            (if (and back-span (eq (back-vertex back-span) (vertex1 edge)))
                                                (progn
                                                  (rotatef id1 id2)
                                                  "->")
                                                nil)))
                                         (t nil))
                                       nil))
                      (dir (if edge-symbol
                               "forward"
                               "none"))
                      (edge-symbol (if edge-symbol
                                       edge-symbol
                                       "->")))
                 (format stream "~a ~a ~a [dir=~a,label = \"~,3f\",penwidth=~a,color=~s];~%"
                         id1 edge-symbol id2 dir (sim-score edge) weight color)))))
      (loop for edge in (edges graph)
            do (draw-edge edge "#000000"))
      (loop for edge in extra-edges
            do (draw-edge edge "#ff0000"))))
  (format stream "}~%"))

(defun draw-graph-to-file (graph filename &key edge-bitvec back-span-info extra-edges)
  (with-open-file (fout filename :direction :output :if-exists :supersede)
    (draw-graph-to-stream graph fout :edge-bitvec edge-bitvec :back-span-info back-span-info :extra-edges extra-edges)))
