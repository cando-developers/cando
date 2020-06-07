

(ql:quickload :lomap)


(in-package :lomap)

(defparameter *mols* (sdf:load-sdf-as-list-of-molecules (format nil "~a/Development/fep-benchmark/eg5/ligands.sdf" (ext:getenv "HOME"))))

;;; Smaller graph
(defparameter *mols4* (subseq *mols* 6 10))
(defparameter *mat4* (similarity-matrix *mols4*))
(defparameter *graph* (similarity-graph *mols4* *mat4*))

;;; Larger graph
(defparameter *mols11* (subseq *mols* 6 16))
(defparameter *mat11* (similarity-matrix *mols11*))
(defparameter *graph* (similarity-graph *mols11* *mat11*))

(graph-set-positions *graph*
                     '(("CHEMBL1096002" 213.29 104.6)
                       ("CHEMBL1084143" 104.4 18)
                       ("CHEMBL1085692" 161.95 47.491)
                       ("CHEMBL1085666" 313.98 194.15)
                       ("CHEMBL1083836" 270.01 156.43)
                       ("CHEMBL1077204" 77.91 73.688)
                       ("CHEMBL1096003" 137.04 100.11)
                       ("CHEMBL1083517" 351.84 163.16)
                       ("CHEMBL1086410" 289.73 110.38)
                       ("CHEMBL1077227" 347.21 114.61)
                       ))
 *graph*

(progn
  (defparameter *new-graph* (lomap-graph *graph* :debug t))
  (format t "Done~%"))



(draw-graph-to-file *graph* "/tmp/tmp.dot")

(defparameter *spanning-tree*
  (calculate-spanning-tree
   *graph*
   (first (vertices *graph*))
   :debug t))


(second (vertices *graph*))
(defparameter *span-bitvec* (generate-backspan-bitvec
                             (second (vertices *graph*))
                             *spanning-tree*
                             *graph*))
(print *span-bitvec*)

(draw-graph-to-file *graph* "/tmp/edge.dot" :edge-bitvec *span-bitvec*)

(draw-span-graph-to-file *graph* "/tmp/edge-span.dot" :edge-bitvec *span-bitvec* :back-span-info *spanning-tree*)

(defparameter *outedges* (edges-outside-of-spanning-tree *graph* *spanning-tree*))

(defparameter *e* (first *outedges*))
(defparameter *s1* (generate-backspan-bitvec
                    (vertex1 *e*)
                    *spanning-tree*
                    *graph*))
(defparameter *s2* (generate-backspan-bitvec
                    (vertex2 *e*)
                    *spanning-tree*
                    *graph*))
(defparameter *s12* (bit-xor *s1* *s2*))
(setf (elt *s12* (index-to-bit (index (vertex1 *e*))
                                      (index (vertex2 *e*))
                                      (length (vertices *graph*)))) 1)
;; Now *s12* has a fundamental cycle
(draw-span-graph-to-file *graph* "/tmp/s12.dot" :edge-bitvec *s12* :back-span-info *spanning-tree*)

(nodes-in-fundamental-cycles-p *graph* *spanning-tree*)
(graph-diameter-p *graph*)

(defparameter *a* (make-array 100 :element-type 'bit))
(setf (elt *a* 32) 1)
(setf (elt *a* 31) 1)
(defparameter *b* (make-array 100 :element-type 'bit))

;;; Try running lomap
(progn
  (defparameter *new-graph* 1)
  (lomap-graph *graph*)
  (print "Hello")
  )
(progn
  (defparameter *new-graph* (lomap-graph *graph* :debug t))
  (format t "Done~%"))
*new-graph*

(edge-bitvec-from-edges *graph*)
(print (length (vertices *graph*)))

*mols11*

(print (index-to-bit 9 8 10))
  
