(in-package :molecule-graph)

(defun atoms-match-p (v1 v2)
  (and (eq (chem:get-element v1)
           (chem:get-element v2))
       (chem:atom-within-angstroms v1 v2 0.2)))

(defun vertex-match-callback (v1 v2)
  (cond
    ((and (typep v1 'chem:atom)
          (typep v2 'chem:atom))
     (atoms-match-p v1 v2))
    ((and (atom-cluster-p v1)
          (atom-cluster-p v2))
     (eq (atom-cluster-equivalent v1) v2))
    (t nil)))


(defparameter *equivalent* nil)

(defun match-callback (results graph1 graph2)
  (setq *equivalent* results)
  (< (length results) 6))


(defstruct atom-cluster equivalent atoms)

(defun merge-clusters (target-cluster source-node)
  "Merge the atoms in source-hash into target-hash"
  (cond
    ((atom-cluster-p source-node)
     (maphash (lambda (atom dummy)
                (declare (ignore dummy))
                (setf (gethash atom (atom-cluster-atoms target-cluster)) atom))
              (atom-cluster-atoms source-node)))
    (t (setf (gethash source-node (atom-cluster-atoms target-cluster)) source-node))))

(defun make-atom-cluster-from-matches (matches graph1 graph2 equivalences)
  (let ((cluster1 (make-atom-cluster :atoms (make-hash-table)))
        (cluster2 (make-atom-cluster :atoms (make-hash-table))))
    (loop for count from 0 below (/ (length matches) 2)
          for index = 0 then (+ index 2)
          for v1 = (aref matches index)
          for v2 = (aref matches (1+ index))
          for node1 = (chem:get-vertex graph1 v1)
          for node2 = (chem:get-vertex graph2 v2)
          do (merge-clusters cluster1 node1)
             (merge-clusters cluster2 node2)
             (when (and (typep node1 'chem:atom) (typep node2 'chem:atom))
               (setf (gethash (cons node1 node2) equivalences) t)))
    (setf (atom-cluster-equivalent cluster1) cluster2)
    (setf (atom-cluster-equivalent cluster2) cluster1)
    (values cluster1 cluster2)))

(defun find-equivalent (graph1 graph2)
  (chem:boost-graph-mcgregor-common-subgraphs
   graph1 graph2
   t
   #'vertex-match-callback
   #'match-callback)
  *equivalent*)


(defun add-remaining-vertices (new-graph graph cluster atoms-to-node-indices)
  (chem:walk-vertices
   graph
   (lambda (node)
     (cond
       ((atom-cluster-p node)
        ;; Are the atoms in the node in the new cluster?
        (let ((subset t))
          (maphash (lambda (atom dummy)
                     (declare (ignore dummy))
                     (unless (gethash atom atoms-to-node-indices)
                       (setf subset nil)))
                   (atom-cluster-atoms node))
          (unless subset
            (let ((cluster-index (chem:add-vertex new-graph node)))
              (maphash (lambda (atom dummy)
                         (declare (ignore dummy))
                         (setf (gethash atom atoms-to-node-indices) cluster-index))
                       (atom-cluster-atoms node))))))
       (t (unless (gethash node (atom-cluster-atoms cluster))
            (let ((node-index (chem:add-vertex new-graph node)))
              (setf (gethash node atoms-to-node-indices) node-index))))))))

(defun add-edges (new-graph atom-to-node-indices original-graph)
  (chem:walk-edges
   original-graph
   (lambda (node1 node2 bond-order)
     (let ((index1 (gethash node1 atom-to-node-indices))
           (index2 (gethash node2 atom-to-node-indices)))
       (unless index1
         (error "Could not find ~a in atom-to-node-indices" node1))
       (unless index2
         (error "Could not find ~a in atom-to-node-indices" node2))
       (unless (= index1 index2)
         (chem:add-edge new-graph index1 index2 bond-order))))))

(defun collapse-graph (graph cluster original-graph)
  "Construct a new graph from the graph and the cluster"
  (let ((new-graph (chem:make-molecule-graph))
        (atom-to-node-indices (make-hash-table)))
    (let ((index (chem:add-vertex new-graph cluster)))
      (maphash (lambda (atom dummy)
                 (declare (ignore dummy))
                 (setf (gethash atom atom-to-node-indices) index))
               (atom-cluster-atoms cluster)))
    (add-remaining-vertices new-graph graph cluster atom-to-node-indices)
    (add-edges new-graph atom-to-node-indices original-graph)
    new-graph))
      
(defun collapse-graphs-once (&key graph1 graph2 original-graph1 original-graph2 equivs)
  "Collapse the equivalent parts of the graph1 and graph2"
  (let ((equiv (find-equivalent graph1 graph2)))
    (multiple-value-bind (cluster1 cluster2)
        (make-atom-cluster-from-matches equiv graph1 graph2 equivs)
      (let ((new-graph1 (collapse-graph graph1 cluster1 original-graph1))
            (new-graph2 (collapse-graph graph2 cluster2 original-graph2)))
        (values new-graph1 new-graph2 (= (length equiv) 2))))))

(defun compare-graphs (original-graph1 original-graph2)
  "Repeatedly collapse equivalent parts of the graph until nothing is left to collapse"
  (let ((equivs (make-hash-table :test #'equal))
        (graph1 original-graph1)
        (graph2 original-graph2)
        done)
    (loop repeat 100
          for iteration from 0
          do (multiple-value-setq (graph1 graph2 done)
               (collapse-graphs-once :graph1 graph1
                                     :graph2 graph2
                                     :original-graph1 original-graph1
                                     :original-graph2 original-graph2
                                     :equivs equivs))
          until done)
    (let ((equiv-atoms (make-hash-table))
          diff1 diff2)
      (maphash (lambda (pair dummy)
                 (declare (ignore dummy))
                 (let ((a1 (car pair))
                       (a2 (cdr pair)))
                   (setf (gethash a1 equiv-atoms) t)
                   (setf (gethash a2 equiv-atoms) t)))
               equivs)
      (chem:walk-vertices
       original-graph1
       (lambda (atom)
         (unless (gethash atom equiv-atoms)
           (push atom diff1))))
      (chem:walk-vertices
       original-graph2
       (lambda (atom)
         (unless (gethash atom equiv-atoms)
           (push atom diff2))))
      (values equivs diff1 diff2 graph1 graph2))))



(defun compare-molecules (molecule1 molecule2)
  "Compare the two molecules.
Return: (values equivalences diff1 diff2)
equivalences - a hash table of cons cells. Each cons cell is a pair of equivalent atoms in molecule1 and molecule2.
diff1 - The atoms of molecule1 that do not have equivalences to anything in molecule 2.
diff2 - The atoms of molecule2 that do not have equivalences to anything in molecule 1."
  (let ((graph1 (chem:make-molecule-graph-from-molecule molecule1))
        (graph2 (chem:make-molecule-graph-from-molecule molecule2)))
    (compare-graphs graph1 graph2)))
