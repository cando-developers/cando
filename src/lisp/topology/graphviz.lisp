(in-package :topology)

(defun generate-graph (graph stream)
  (format stream "digraph mygraph {~%")
  (maphash (lambda (key node)
             (format stream "  ~s~%" (symbol-name (name node))))
           (nodes graph))
  (loop for edge in (edges graph)
        for from-node = (from-node edge)
        for to-node = (to-node edge)
        for edge-type = (edge-type edge)
        do (format stream "~s -> ~s [label=~s]~%"
                   (symbol-name (name from-node))
                   (symbol-name (name to-node))
                   (symbol-name edge-type)))
  (format stream "}~%"))


(defun draw-graph (graph filename)
  (with-open-file (fout filename :direction :output)
    (generate-graph graph fout)))
