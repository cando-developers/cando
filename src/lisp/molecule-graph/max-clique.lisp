
(in-package :molecule-graph.max-clique)

(defstruct molecule-graph atoms edges)

(defun build-molecule-graph (molecule &key (exclude-hydrogens t))
  (let ((atoms (make-array 256 :adjustable t :fill-pointer 0)))
    (cando:do-atoms (atm molecule)
      (unless (and exclude-hydrogens (= (chem:get-atomic-number atm) 1))
        (vector-push-extend atm atoms)))
    (let ((edges (make-array (list (length atoms) (length atoms)) :element-type 'ext:byte8 :initial-element 0)))
      (loop for xi below (length atoms)
            for xatom = (elt atoms xi)
            do (loop for yi below (length atoms)
                     for yatom = (elt atoms yi)
                     do (if (eq xatom yatom)
                            (setf (aref edges xi yi) 1)
                            (when (chem:is-bonded-to xatom yatom)
                              (setf (aref edges xi yi) 1)))))
      (make-molecule-graph :atoms (copy-seq atoms) :edges edges))))

(defun topological-distance-matrix (graph)
  "Return the topological distance matrix for a molecule graph and a hashtable of
atoms to indices into the matrix."
  ;; Save the index of each atom
  (let ((atom-indices (make-hash-table)))
    (loop for xi below (length (molecule-graph-atoms graph))
          for xatom = (elt (molecule-graph-atoms graph) xi)
          do (setf (gethash xatom atom-indices) xi))
    ;; Run a spanning loop from every atom and build a matrix of distances back
    ;; to the root atom.
    (let ((topological-distance-matrix (make-array (list (hash-table-count atom-indices) (hash-table-count atom-indices)) :initial-element 0)))
      (loop for xi below (hash-table-count atom-indices)
            for xatom = (elt (molecule-graph-atoms graph) xi)
            do (let ((spanning-loop (chem:make-spanning-loop xatom)))
                 (loop for next-atom = (chem:advance-loop-and-process spanning-loop)
                       if next-atom
                         do (let* ((atm (chem:get-atom spanning-loop))
                                   (atm-index (gethash atm atom-indices)))
                              (when atm-index
                                (let ((back-count (chem:get-back-count atm)))
                                  (setf (aref topological-distance-matrix xi atm-index) back-count))))
                       else
                         do (progn
                              (return)))))
      (values topological-distance-matrix atom-indices))))

(defstruct cross-product graph1 graph2 nodes edges)

(defun modular-cross-product (graph1 graph2 &key atom-match-callback)
  "Return a cross-product for the two molecules.
If atom-match-callback is NIL then all atoms are considered equivalent."
  (let ((pairs (make-array 1024 :adjustable t :fill-pointer 0)))
    (if atom-match-callback
        (loop for xi below (length (molecule-graph-atoms graph1))
              for xatom = (elt (molecule-graph-atoms graph1) xi)
              do (loop for yi below (length (molecule-graph-atoms graph2))
                       for yatom = (elt (molecule-graph-atoms graph2) yi)
                       do (when (funcall atom-match-callback xatom yatom)
                            (vector-push-extend (cons xi yi) pairs))))
        (loop for xi below (length (molecule-graph-atoms graph1))
              do (loop for yi below (length (molecule-graph-atoms graph2))
                       do (vector-push-extend (cons xi yi) pairs))))
    (let* ((pairs-length (length pairs))
           (edges (make-array (list pairs-length pairs-length)
                              :element-type 'ext:byte8 :initial-element 0)))
      (loop for xii below pairs-length
            do (let* ((xpair (elt pairs xii))
                      (u (car xpair))
                      (v (cdr xpair)))
                 (loop for yii below pairs-length
                       do (let* ((ypair (elt pairs yii))
                                 (up (car ypair))
                                 (vp (cdr ypair))
                                 (u-up (aref (molecule-graph-edges graph1) u up))
                                 (v-vp (aref (molecule-graph-edges graph2) v vp))
                                 (adjacent (if (and (/= u up) (/= v vp) (= u-up v-vp)) 1 0)))
                            #+(or)
                            (progn
                              (format t "-------xii yii -> ~a ~a~%" xii yii)
                              (format t "u(~a) up(~a) u-up -> ~a~%" u up u-up)
                              (format t "v(~a) vp(~a) v-vp -> ~a~%" v vp v-vp)
                              (format t "xii(~a) yii(~a)  adjacent -> ~a~%" xii yii adjacent))
                            (setf (aref edges xii yii) adjacent)))))
      (make-cross-product :graph1 graph1 :graph2 graph2 :nodes pairs :edges edges))))

(defun modular-cross-product-with-topological-constraint (graph1 graph2 theta &key atom-match-callback)
  "Return a modular-cross-product for the two graphs with a topological constraint.
The topological constraint is best described here: dx.doi.org/10.1021/ci2001023 J. Chem. Inf. Model.2011, 51, 1775–1787.
You have two graphs graph1 and graph2.
The modular-cross-product is a new graph with vertices (v1,v2) where (v1 element of graph1) and (v2 element-of graph2).
A constraint to created vertices can be defined using atom-match-callback.
If atom-match-callback is NIL then all atoms are considered equivalent and all (v1,v2) vertices are created.
The edges ((v1,v2) - (v1',v2')) in the new graph are created in two situations.
1. Connected edges: If bond(v1,v1') and bond(v2,v2') then create the edge.
2. Disconnected edges: if !bond(v1,v1') and !bond(v2,v2') and abs(T(v1,v1')-T(v2,v2'))<=theta.
The last constraint abs(T(v1,v1')-T(v2,v2'))<=theta is the topological constraint. 
"
  (let ((pairs (make-array 1024 :adjustable t :fill-pointer 0)))
    (if atom-match-callback
              (loop for xi below (length (molecule-graph-atoms graph1))
                    for xatom = (elt (molecule-graph-atoms graph1) xi)
                    do (loop for yi below (length (molecule-graph-atoms graph2))
                             for yatom = (elt (molecule-graph-atoms graph2) yi)
                             do (when (funcall atom-match-callback xatom yatom)
                                  (vector-push-extend (cons xi yi) pairs))))
              (loop for xi below (length (molecule-graph-atoms graph1))
                    do (loop for yi below (length (molecule-graph-atoms graph2))
                             do (vector-push-extend (cons xi yi) pairs))))
    (let* ((pairs-length (length pairs))
           (edges (make-array (list pairs-length pairs-length)
                              :element-type 'ext:byte8 :initial-element 0)))
      (loop for xii below pairs-length
            do (let* ((xpair (elt pairs xii))
                      (u (car xpair))
                      (v (cdr xpair)))
                 (loop for yii below pairs-length
                       do (let* ((ypair (elt pairs yii))
                                 (up (car ypair))
                                 (vp (cdr ypair))
                                 (u-up (aref (molecule-graph-edges graph1) u up))
                                 (v-vp (aref (molecule-graph-edges graph2) v vp))
                                 (connected-edge (and (/= u up) (/= v vp) (= 1 u-up v-vp))))
                            (when connected-edge
                              (setf (aref edges xii yii) 1))))))
      ;; Created Disconnected edges with topological constraint
      (multiple-value-bind (topological-matrix-graph1 atom-index-map-graph1)
          (topological-distance-matrix graph1)
        (multiple-value-bind (topological-matrix-graph2 atom-index-map-graph2)
            (topological-distance-matrix graph2)
          (loop for xii below pairs-length
                do (let* ((xpair (elt pairs xii))
                          (u (car xpair))
                          (v (cdr xpair)))
                     (loop for yii below pairs-length
                           do (let* ((ypair (elt pairs yii))
                                     (up (car ypair))
                                     (vp (cdr ypair))
                                     (u-up (aref (molecule-graph-edges graph1) u up))
                                     (v-vp (aref (molecule-graph-edges graph2) v vp))
                                     (disconnected-edge (and (/= u up) (/= v vp) (= 0 u-up v-vp))))
                                (when disconnected-edge
                                  (let* ((top-dist1 (aref topological-matrix-graph1 u up))
                                         (top-dist2 (aref topological-matrix-graph2 v vp))
                                         (top-dist-delta (abs (- top-dist1 top-dist2)))
                                         (satisfy-top-dist (<= top-dist-delta theta)))
                                    #+(or)(format t "u,v (~a,~a)  up,vp (~a,~a) top-dist1(~a) top-dist2(~a) top-dist-delta(~a) satisfy-top-dist(~a)~%" u v up vp top-dist1 top-dist2 top-dist-delta satisfy-top-dist)
                                    (when satisfy-top-dist
                                      (setf (aref edges xii yii) 1))))))))))
      (make-cross-product :graph1 graph1 :graph2 graph2 :nodes pairs :edges edges))))

(defun molecular-modular-cross-product (molecule1 molecule2 &key (exclude-hydrogens t) atom-match-callback)
  "Return a cross-product for the two molecules"
  (let ((graph1 (build-molecule-graph molecule1 :exclude-hydrogens exclude-hydrogens))
        (graph2 (build-molecule-graph molecule2 :exclude-hydrogens exclude-hydrogens)))
    (modular-cross-product graph1 graph2 :atom-match-callback atom-match-callback)))


(defun test-cross-product ()
  (flet ((make-bonds (atoms edges pairs)
           (loop for pair in pairs
                 for a1 = (elt atoms (first pair))
                 for a2 = (elt atoms (second pair))
                 do (chem:bond-to a1 a2 :single-bond)
                 do (setf (aref edges (first pair) (second pair)) 1
                          (aref edges (second pair) (first pair)) 1))))
    (let ((g1-atoms (make-array 3 :initial-contents (list (chem:make-atom :ac1 :c)
                                                          (chem:make-atom :ac2 :c)
                                                          (chem:make-atom :ac3 :c))))
          (g2-atoms (make-array 3 :initial-contents (list (chem:make-atom :bc1 :c)
                                                          (chem:make-atom :bc2 :c)
                                                          (chem:make-atom :bc3 :c))))
          (g1-edges (make-array (list 3 3) :element-type 'ext:byte8 :initial-element 0))
          (g2-edges (make-array (list 3 3) :element-type 'ext:byte8 :initial-element 0)))
      (make-bonds g1-atoms g1-edges '((0 1) (1 2)))
      (make-bonds g2-atoms g2-edges '((0 1) (1 2)))
      (let ((g1 (make-molecule-graph :atoms g1-atoms :edges g1-edges))
            (g2 (make-molecule-graph :atoms g2-atoms :edges g2-edges)))
        (values (modular-cross-product g1 g2) g1 g2)))))

(defun build-dimacs (cross-product)
  (let* ((dims (array-dimensions (cross-product-edges cross-product)))
         (xdim (first dims))
         (ydim (second dims)))
    (let ((num-vertices (length (cross-product-nodes cross-product)))
          (num-edges (let ((count 0))
                       (loop for xi below xdim
                             do (loop for yi below ydim
                                      do (when (> (aref (cross-product-edges cross-product) xi yi) 0)
                                           (incf count))))
                       count)))
      (let ((dimacs (chem:make-dimacs num-vertices num-edges)))
        (loop for xi below xdim
              do (loop for yi below ydim
                       do (if (> (aref (cross-product-edges cross-product) xi yi) 0)
                              (chem:dimacs-add-edge dimacs xi yi))))
        dimacs))))

(defun write-dimacs (cross-product stream)
  (let* ((dims (array-dimensions (cross-product-edges cross-product)))
         (xdim (first dims))
         (ydim (second dims)))
    (let ((num-vertices (length (cross-product-nodes cross-product)))
          (num-edges (let ((count 0))
                       (loop for xi below xdim
                             do (loop for yi below ydim
                                      do (when (> (aref (cross-product-edges cross-product) xi yi) 0)
                                           (incf count))))
                       count)))
      (format stream "p edge ~a ~a~%" num-vertices num-edges)
      (loop for xi below xdim
            do (loop for yi below ydim
                     do (if (> (aref (cross-product-edges cross-product) xi yi) 0)
                            (format stream "e ~a ~a~%" (1+ xi) (1+ yi))))))))

(defun extract-atom-pairs (cross-product indices)
  (loop for index across indices
        for real-index = index
        for pair = (elt (cross-product-nodes cross-product) real-index)
        for xii = (car pair)
        for yii = (cdr pair)
        for xatom = (elt (molecule-graph-atoms (cross-product-graph1 cross-product)) xii)
        for yatom = (elt (molecule-graph-atoms (cross-product-graph2 cross-product)) yii)
        collect (cons xatom yatom)))

(defun rename-atoms (mol prefix)
  (cando:do-atoms (atm mol)
                  (let* ((name (chem:get-name atm))
                         (new-name (intern (format nil "~a~a" prefix (string name)))))
                    (chem:set-name atm new-name))))


(defun dothing (dimacs)
  (chem:find-maximum-clique-search dimacs (core:num-logical-processors) 1))


(defun compare-molecules (molecule1 molecule2 &key topological-constraint-theta (exclude-hydrogens t) atom-match-callback)
  "Return three values:
a list of pairs of matching atoms 
a list of mismatching atoms from molecule1
a list of mismatching atoms from molecule2.
If atom-match-callback is NIL then all atoms are considered equivalent.
Otherwise pass a function that takes two atoms and returns T if they are matchable."
  (let* ((graph1 (build-molecule-graph molecule1 :exclude-hydrogens exclude-hydrogens))
         (graph2 (build-molecule-graph molecule2 :exclude-hydrogens exclude-hydrogens))
         (cross-product (if topological-constraint-theta
                            (modular-cross-product-with-topological-constraint graph1 graph2 topological-constraint-theta :atom-match-callback atom-match-callback)
                            (modular-cross-product graph1 graph2 :atom-match-callback atom-match-callback)))
         (dimacs (build-dimacs cross-product))
         (maximum-clique (chem:find-maximum-clique-search dimacs (core:num-logical-processors) 1))
         (matches (extract-atom-pairs cross-product maximum-clique))
         (mol1-matches (make-hash-table))
         (mol2-matches (make-hash-table)))
    (loop for pair in matches
          for atom1 = (car pair)
          for atom2 = (cdr pair)
          do (setf (gethash atom1 mol1-matches) t
                   (gethash atom2 mol2-matches) t))
    (let ((mol1-mismatches (loop for atm across (molecule-graph-atoms graph1)
                                 unless (gethash atm mol1-matches)
                                   collect atm))
          (mol2-mismatches (loop for atm across (molecule-graph-atoms graph2)
                                 unless (gethash atm mol2-matches)
                                   collect atm)))
      (values matches mol1-mismatches mol2-mismatches))))
          
         
    
