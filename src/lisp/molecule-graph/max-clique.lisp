
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

(defun molecular-modular-cross-product (molecule1 molecule2 &key (exclude-hydrogens t) atom-match-callback)
  "Return a cross-product for the two molecules"
  (let ((graph1 (build-molecule-graph molecule1 :exclude-hydrogens exclude-hydrogens))
        (graph2 (build-molecule-graph molecule2 :exclude-hydrogens exclude-hydrogens)))
    (modular-cross-product graph1 graph2 :atom-match-callback atom-match-callback)))


(defun test-cross-product ()
  (flet ((make-bonds (edges pairs)
           (loop for pair in pairs
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
      (make-bonds g1-edges '((0 1) (1 2)))
      (make-bonds g2-edges '((0 1) (1 2)))
      (let ((g1 (make-molecule-graph :atoms g1-atoms :edges g1-edges))
            (g2 (make-molecule-graph :atoms g2-atoms :edges g2-edges)))
        (modular-cross-product g1 g2)))))

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


(defun compare-molecules (molecule1 molecule2 &key (exclude-hydrogens t) atom-match-callback)
  "Return three values:
a list of pairs of matching atoms 
a list of mismatching atoms from molecule1
a list of mismatching atoms from molecule2.
If atom-match-callback is NIL then all atoms are considered equivalent.
Otherwise pass a function that takes two atoms and returns T if they are matchable."
  (let* ((graph1 (build-molecule-graph molecule1 :exclude-hydrogens exclude-hydrogens))
         (graph2 (build-molecule-graph molecule2 :exclude-hydrogens exclude-hydrogens))
         (cross-product (modular-cross-product graph1 graph2 :atom-match-callback atom-match-callback))
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
          
         
    
