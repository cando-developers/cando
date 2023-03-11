(in-package :topology)


(defclass conformation-selector ()
  ((all-conformations :initarg :all-conformations :accessor all-conformations)
   (consistent-conformations :initform nil :initarg :consistent-conformations :accessor consistent-conformations)
   (consistent-index :initform 0 :initarg :consistent-index :accessor consistent-index)))



(defun classify-dihedral-angle (angle-degrees)
  (when (< angle-degrees 0.0)
    (setf angle-degrees (+ 360.0 angle-degrees)))
  (let* ((rounded (round (/ angle-degrees 30.0)))
         (result (if (= rounded 12) 0 rounded)))
    result))


(defun classify-dihedral-angles (angles)
  (loop with accumulator = 0
        for angle in angles
        for classified-angle = (classify-dihedral-angle angle)
        do (setf accumulator (+ (* accumulator 12) classified-angle))
        finally (return-from classify-dihedral-angles accumulator)))

(defparameter *saved-fi* nil)
(defun fill-internals-from-fragments-debug (conf fragments &optional index)
  "Fill internal coordinates from the fragments"
  (let ((monomer-to-conformation-selector (make-hash-table))
        (all-classified (make-hash-table :test 'eql)))
    (loop with atagg = (ataggregate conf)
          with atmol = (elt (atmolecules atagg) 0)
          for monomer across (monomers (oligomer conf))
          for monomer-context = (gethash monomer (monomer-contexts conf))
          for _mc = (format t "monomer-context = ~a~%" monomer-context)
          for frags = (gethash monomer-context (monomer-context-to-fragment-conformations fragments))
          for monomer-position = (gethash monomer (monomer-positions conf))
          for monomer-index = (monomer-position-residue-index monomer-position)
          for atres = (elt (atresidues atmol) monomer-index)
          for fragment-conformations = (let ((fp (gethash monomer-context (monomer-context-to-fragment-conformations fragments))))
                                (unless fp
                                  (error "Cannot find context ~a in conformations" monomer-context))
                                fp)
          for conformations = (fragments fragment-conformations)
          for conformation-selector = (make-instance 'conformation-selector
                                                     :all-conformations conformations)
          do (setf (gethash monomer monomer-to-conformation-selector) conformation-selector))
    (let ((ordered-monomers (ordered-monomers (oligomer conf))))
      (loop with previous-monomer = nil
            with previous-conformation-selector = nil
            for monomer in ordered-monomers
            for conformation-selector = (gethash monomer monomer-to-conformation-selector)
            do (progn
                 (if (null previous-monomer)
                     (setf (consistent-index conformation-selector) 0
                           (consistent-conformations conformation-selector) (all-conformations conformation-selector))
                     (let* ((previous-conformation (elt (consistent-conformations previous-conformation-selector)
                                                        (consistent-index previous-conformation-selector)))
                            (out-of-focus (out-of-focus-internals previous-conformation)))
                     
                       (format t "Look for consistent-conformations for ~a : ~a~%" monomer (length (all-conformations conformation-selector)))
                       (format t "out-of-focus ~s~%" out-of-focus)
                       (loop for conf in (all-conformations conformation-selector)
                             for internals = (internals conf)
                             for internal0 = (first internals)
                             for internal1 = (second internals)
                             for internal2 = (third internals)
                             for dihedral0 = (/ (topology:dihedral internal0) 0.0174533)
                             for dihedral1 = (/ (topology:dihedral internal1) 0.0174533)
                             for dihedral2 = (/ (topology:dihedral internal2) 0.0174533)
                             for classify0 = (classify-dihedral-angle dihedral0)
                             for classify1 = (classify-dihedral-angle dihedral1)
                             for classify2 = (classify-dihedral-angle dihedral2)
                             for classified = (classify-dihedral-angles (list dihedral0 dihedral1 dihedral2))
                             
                             do (format t "internals = ~a, ~a, ~a = ~a~%"
                                        (topology:name internal0)
                                        (topology:name internal1)
                                        (topology:name internal2)
                                        classified)
                             do (incf (gethash classified all-classified 0))))))
               (setf previous-monomer monomer
                     previous-conformation-selector conformation-selector)))
    monomer-to-conformation-selector))

#|
  (loop with atagg = (ataggregate conf)
        with atmol = (elt (atmolecules atagg) 0)
        for monomer across (monomers (oligomer conf))
        for monomer-context = (gethash monomer (monomer-contexts conf))
        for _mc = (format t "monomer-context = ~a~%" monomer-context)
        for frags = (gethash monomer-context (monomer-context-to-fragment-conformations fragments))
        for monomer-index = (gethash monomer (monomer-positions conf))
        for atres = (elt (atresidues atmol) monomer-index)
        for fragment-conformations = (let ((fp (gethash monomer-context (monomer-context-to-fragment-conformations fragments))))
                              (unless fp
                                (error "Cannot find context ~a in conformations" monomer-context))
                              fp)
        for _fp = (format t "fragment-conformations = ~a~%" fragment-conformations)
        for rand-limit = (length (fragments fragment-conformations))
        for rand-index = (if index index (random rand-limit))
        for fragment-internals = (let ((fragments (fragments fragment-conformations)))
                                   (unless fragments
                                     (error "fragments is NIL for context ~a" monomer-context))
                                   (elt (fragments fragment-conformations) rand-index))
        for _fi = (format t "fragment-internals = ~a~%" fragment-internals)
        do (push fragment-internals *saved-fi*)
        do (loop for joint across (joints atres)
                 for internal in (internals fragment-internals)
                 do (fill-joint-internals joint internal))
        ))
|#
