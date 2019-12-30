
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

(defstruct clique-graph graph1 graph2 nodes edges)

(defun build-clique-graph-from-graphs (graph1 graph2 &key match-elements)
  "Return a clique-graph for the two molecules"
  (let ((pairs (make-array 1024 :adjustable t :fill-pointer 0)))
    (if match-elements
        (loop for xi below (length (molecule-graph-atoms graph1))
              for xatom = (elt (molecule-graph-atoms graph1) xi)
              do (loop for yi below (length (molecule-graph-atoms graph2))
                       for yatom = (elt (molecule-graph-atoms graph2) yi)
                       do (when (eq (chem:get-element xatom) (chem:get-element yatom))
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
      (format t "Ending~%")
      (make-clique-graph :graph1 graph1 :graph2 graph2 :nodes pairs :edges edges))))

(defun build-clique-graph (molecule1 molecule2 &key (exclude-hydrogens t))
  "Return a clique-graph for the two molecules"
  (let ((graph1 (build-molecule-graph molecule1 :exclude-hydrogens exclude-hydrogens))
        (graph2 (build-molecule-graph molecule2 :exclude-hydrogens exclude-hydrogens)))
    (build-clique-graph-from-graphs graph1 graph2)))


(defun test-clique-graph ()
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
        (build-clique-graph-from-graphs g1 g2)))))


(defun write-dimacs (clique-graph stream)
  (let* ((dims (array-dimensions (clique-graph-edges clique-graph)))
         (xdim (first dims))
         (ydim (second dims)))
    (let ((num-vertices (length (clique-graph-nodes clique-graph)))
          (num-edges (let ((count 0))
                       (loop for xi below xdim
                             do (loop for yi below ydim
                                      do (when (> (aref (clique-graph-edges clique-graph) xi yi) 0)
                                           (incf count))))
                       count)))
      (format stream "p edge ~a ~a~%" num-vertices num-edges)
      (loop for xi below xdim
            do (loop for yi below ydim
                     do (if (> (aref (clique-graph-edges clique-graph) xi yi) 0)
                            (format stream "e ~a ~a~%" (1+ xi) (1+ yi))))))))

(defun translate-result (clique-graph indices)
  (loop for index in indices
        for real-index = index
        for pair = (elt (clique-graph-nodes clique-graph) real-index)
        for xii = (car pair)
        for yii = (cdr pair)
        for xatom = (elt (molecule-graph-atoms (clique-graph-graph1 clique-graph)) xii)
        for yatom = (elt (molecule-graph-atoms (clique-graph-graph2 clique-graph)) yii)
        do (format t "~a <-> ~a~%" xatom yatom)))



(defun rename-atoms (mol prefix)
  (cando:do-atoms (atm mol)
                  (let* ((name (chem:get-name atm))
                         (new-name (intern (format nil "~a~a" prefix (string name)))))
                    (chem:set-name atm new-name))))
