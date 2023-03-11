(in-package :foldamer)

(defun render-dag (dag stream)
  (format stream "digraph mygraph {~%")
  (maphash (lambda (key node)
             (declare (ignore key))
             (cond
               ((typep node 'cap-node)
                (format stream "  ~s~%" (symbol-name (name node))))
               ((typep node 'node)
                (format stream "  ~s~%" (symbol-name (name node))))
               ))
           (nodes dag))
  (loop for edge in (edges dag)
        for from-node = (from-node edge)
        for to-node = (to-node edge)
        for edge-type = (raw-name edge)
        do (format stream "~s -> ~s [label=~s]~%"
                   (symbol-name (name from-node))
                   (symbol-name (name to-node))
                   (symbol-name edge-type)))
  (format stream "}~%"))


(defun draw-dag (dag filename)
  (with-open-file (fout filename :direction :output)
    (render-dag dag fout)))


(defun node-name (node)
  (format nil "~a" (id node)))

(defun render-multiple-dags (dags stream)
  (format stream "digraph mygraph {~%")
  (loop for dag in dags
        for index from 0
        do (format stream "subgraph {~%")
        do (mapc (lambda (node)
                   (cond
                     ((eq node (root dag))
                      (format stream "  ~s[label=~s,shape=box,fillcolor=lightblue,style=filled]~%" (node-name node) (string (name node))))
                     ((typep node 'cap-node)
                      (format stream "  ~s[label=~s,shape=diamond]~%" (node-name node) (string (name node))))
                     ((typep node 'node)
                      (format stream "  ~s[label=~s]~%" (node-name node) (string (name node))))
                     ))
                 (nodes dag))
        do (loop for edge in (edges dag)
                 for from-node = (from-node edge)
                 for to-node = (to-node edge)
                 for edge-type = (raw-name edge)
                 do (format stream "  ~s -> ~s [label=~s]~%"
                            (node-name from-node)
                            (node-name to-node)
                            (symbol-name edge-type)))
        do (format stream "}~%")
        )
  (format stream "}~%"))

(defun draw-foldamer (foldamer filename)
  (with-open-file (fout filename :direction :output)
    (let ((dags (loop for training-oligomer-space in (foldamer::training-oligomer-spaces foldamer)
                      collect (foldamer::expression-dag training-oligomer-space))))
      (render-multiple-dags dags fout))))


;;;
;;; Render topology
;;;

(in-package topology.graphviz)
(defgeneric make-graph (object))
(defgeneric graph-label (graph))
(defgeneric graph-name (graph))
(defgeneric nodes (graph))
(defgeneric node-label (object graph))
(defgeneric node-id (object graph))
(defgeneric directed-edges (object))
(defgeneric directed-edge-from (object))
(defgeneric directed-edge-to (object))
(defgeneric directed-edge-label (object))
(defgeneric undirected-edges (object))
(defgeneric undirected-edge-from (object graph))
(defgeneric undirected-edge-to (object graph))
(defgeneric undirected-edge-label (object graph))

(defun render-subgraph (dag stream)
  (format stream "subgraph cluster_~a {~%" (graph-name dag))
  (format stream "style=filled; color=lightgrey;~%")
  (format stream "label=~s~%" (graph-label dag))
  (loop for node in (nodes dag)
        do (format stream "  ~a [label=~s]~%" (node-id node dag) (node-label node dag)))
  (let ((directed-edges (directed-edges dag)))
    (when directed-edges
      (format stream "  subgraph Rel1 {~%")
      (format stream "     edge [color=red]~%")
      (loop for edge in directed-edges
            for from-node = (directed-edge-from edge)
            for to-node = (directed-edge-to edge)
            for edge-type = (directed-edge-label edge)
            do (format stream "     ~a -> ~a [label=~s]~%"
                       (node-id from-node dag)
                       (node-id to-node dag)
                       (if edge-type edge-type "")))
      (format stream "  }~%")))
  (let ((undirected-edges (undirected-edges dag)))
    (when undirected-edges
      (format stream "  subgraph Rel2 {~%")
      (format stream "     edge [dir=none, color=blue]~%")
      (loop for edge in undirected-edges
            for from-node = (undirected-edge-from edge dag)
            for to-node = (undirected-edge-to edge dag)
            for edge-type = (undirected-edge-label edge dag)
            do (format stream "     ~a -> ~a [label=~s]~%"
                       (node-id from-node dag)
                       (node-id to-node dag)
                       (if edge-type edge-type "")))
      (format stream "  }~%")))
  (format stream "}~%"))




;;; Render different objects


(defun dot-svg-foldamer-joint-trees (foldamer filename)
  (let ((path (probe-file filename)))
    (if path
        (progn
          (with-open-file (fout filename :direction :output)
            (format fout "digraph mine {~%")
            (loop for topology in (foldamer:topologys foldamer)
                  for graph = (make-graph topology)
                  do (render-subgraph graph fout))
            (format fout "}~%"))
          (let* ((stream (ext:run-program "/bin/bash" (list "-c" (format nil "dot -Tsvg -O ~a" (namestring path)))))
                 (result (read-line stream nil)))
            (when result
              (format t "~a~%" result))))
        (error "Could not find file ~a~%" filename))))
