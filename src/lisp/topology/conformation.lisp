(in-package :topology)

(defclass conformation ()
  ((monomer-positions :initarg :monomer-positions :accessor monomer-positions)
   (monomer-contexts :type hash-table :initarg :monomer-contexts :accessor monomer-contexts)
   (oligomer :initarg :oligomer :accessor oligomer)
   (aggregate :initarg :aggregate :accessor aggregate)
   (ataggregate :initarg :ataggregate :accessor ataggregate)
   (joint-tree :initarg :joint-tree :accessor joint-tree)))

(defun make-conformation (oligomer &key focus-monomer)
  "Build a conformation for the oligomer.   If focus-monomer is set (to a monomer) then
   only worry about that monomers monomer-context."
  (let* ((oligomer-space (oligomer-space oligomer))
         (foldamer (foldamer oligomer-space))
         (monomer-contexts (make-hash-table)))
    (if focus-monomer
        (let ((monomer-context (foldamer:foldamer-monomer-context focus-monomer oligomer foldamer)))
          (setf (gethash focus-monomer monomer-contexts) monomer-context))
        (loop for monomer across (monomers oligomer-space)
              for monomer-context = (foldamer:foldamer-monomer-context monomer oligomer foldamer)
              do (setf (gethash monomer monomer-contexts) monomer-context)))
    (multiple-value-bind (aggregate molecule monomer-positions)
        (topology:build-molecule oligomer)
      (let* ((ataggregate (make-instance 'ataggregate :aggregate aggregate))
             (joint-tree (make-joint-tree)))
        (resize-atmolecules ataggregate 1)
        (let ((molecule-index 0))
          ;; This is where I would invoke Conformation_O::buildMoleculeUsingOligomer
          ;; Use the monomers-to-topologys
          (let* ((atmolecule (build-atmolecule-using-oligomer oligomer molecule molecule-index monomer-positions joint-tree)))
            (put-atmolecule ataggregate atmolecule molecule-index)
            (let ((conf (make-instance 'conformation
                                       :monomer-positions monomer-positions
                                       :monomer-contexts monomer-contexts
                                       :oligomer oligomer
                                       :aggregate aggregate
                                       :ataggregate ataggregate
                                       :joint-tree joint-tree)))
              conf
              )))))))

(defun walk-atoms-joints (conformation callback)
  (let ((aggregate (aggregate conformation))
        (ataggregate (ataggregate conformation)))
    (loop for molecule-index below (chem:content-size aggregate)
          for molecule = (chem:content-at aggregate molecule-index)
          for atmolecule = (aref (atmolecules ataggregate) molecule-index)
          do (loop for residue-index below (chem:content-size molecule)
                   for residue = (chem:content-at molecule residue-index)
                   for atresidue = (aref (atresidues atmolecule) residue-index)
                   do (loop for atom-index below (chem:content-size residue)
                           for atom = (chem:content-at residue atom-index)
                           for joint = (aref (joints atresidue) atom-index)
                           do (funcall callback atom joint (list molecule-index residue-index atom-index)))))))

(defun copy-atom-positions-into-joints (conformation)
  (walk-atoms-joints conformation
                     (lambda (atm jnt atomid)
                       (declare (ignore atomid))
                       (kin:set-position jnt (chem:get-position atm)))))

(defun copy-joint-positions-into-atoms (conformation)
  (walk-atoms-joints conformation
                     (lambda (atm jnt atomid)
                       (declare (ignore atomid))
                       (chem:set-position atm (kin:get-position jnt)))))

(defun update-joint-tree-internal-coordinates (conformation)
  (let ((ataggregate (ataggregate conformation)))
    (walk-ataggregate-joints ataggregate
                             (lambda (joint atom-id)
                               (declare (ignore atom-id))
                               (kin:update-internal-coord joint)))))

(defun build-all-atom-tree-external-coordinates (conformation)
  (let ((joint (root (joint-tree conformation))))
    (kin:update-xyz-coords joint)))


(defun zero-all-atom-tree-external-coordinates (conf)
  "Set the external coordinates for each joint to the origin"
  (let ((ataggregate (ataggregate conf)))
    (walk-ataggregate-joints ataggregate
                             (lambda (joint atom-id)
                               (declare (ignore atom-id))
                               (kin:set-position joint (geom:vec 0.0 0.0 0.0))))))


(defun fill-internals-from-fragments (conf fragments &optional (index 0))
  "Fill internal coordinates from the fragments"
  (loop with atagg = (ataggregate conf)
        with atmol = (elt (atmolecules atagg) 0)
        for monomer across (monomers (oligomer conf))
        for monomer-context = (gethash monomer (monomer-contexts conf))
        for frags = (gethash monomer-context (monomer-context-to-fragment-pool fragments))
        for monomer-index = (gethash monomer (monomer-positions conf))
        for atres = (elt (atresidues atmol) monomer-index)
        for fragment-pool = (let ((fp (gethash monomer-context (monomer-context-to-fragment-pool fragments))))
                              (unless fp
                                (error "Cannot find context ~a in conformations" monomer-context))
                              fp)
        for rand-limit = (length (fragments fragment-pool))
        for rand-index = (if index index (random rand-limit))
        for fragment-internals = (let ((fragments (fragments fragment-pool)))
                                   (unless fragments
                                     (error "fragments is NIL for context ~a" monomer-context))
                                   (elt (fragments fragment-pool) rand-index))
        do (loop for joint across (joints atres)
                 for internal in (internals fragment-internals)
                 do (fill-joint-internals joint internal))
        ))

(defun build-conformation (oligomer-space fragment-conformations)
  (format t "Entered build-conformation~%")
  (let* ((oligomer (topology:make-oligomer oligomer-space 0))
         (conf (topology:make-conformation oligomer)))
    (topology::fill-internals-from-fragments conf fragment-conformations)
    (topology:zero-all-atom-tree-external-coordinates conf)
    (topology:build-all-atom-tree-external-coordinates conf)
    (topology:copy-joint-positions-into-atoms conf)
    (topology:aggregate conf)))

(defun search-conformations (oligomer-space fragment-conformations monomer-names sdf-filename &optional (number-struct 50) (number-conf 100))
  (with-open-file (fout sdf-filename :direction :output)
    (loop for struct-count below number-struct
          for rand-sequence = (random (topology:number-of-sequences oligomer-space))
          for _a = (format t "rand-sequence ~a~%" rand-sequence)
          for oligomer = (topology:make-oligomer oligomer-space rand-sequence)
          for conf = (topology:make-conformation oligomer)
          do (topology::fill-internals-from-fragments conf fragment-conformations 0)
             (loop for count below number-conf
                   do (topology::fill-internals-from-fragments-for-monomers-named conf fragment-conformations monomer-names)
                   do (topology:zero-all-atom-tree-external-coordinates conf)
                   do (topology:build-all-atom-tree-external-coordinates conf)
                   do (topology:copy-joint-positions-into-atoms conf)
                   do (sdf:write-sdf-stream (topology:aggregate conf) fout)))))

(defun fill-internals-from-fragments-for-monomers-named (conf fragments monomer-names)
  "Fill internal coordinates from the fragments"
  (loop with atagg = (ataggregate conf)
        with atmol = (elt (atmolecules atagg) 0)
        for monomer across (monomers (oligomer conf))
        when (member (topology:current-stereoisomer-name monomer (oligomer conf)) monomer-names)
          do (let* ((monomer-context (gethash monomer (monomer-contexts conf)))
                    (frags (gethash monomer-context (monomer-context-to-fragment-pool fragments)))
                    (monomer-index (gethash monomer (monomer-positions conf)))
                    (atres (elt (atresidues atmol) monomer-index))
                    (fragment-pool (gethash monomer-context (monomer-context-to-fragment-pool fragments)))
                    (rand-limit (length (fragments fragment-pool)))
                    (rand-index (random rand-limit))
                    (fragment-internals (elt (fragments fragment-pool) rand-index)))
               (loop for joint across (joints atres)
                     for internal in (internals fragment-internals)
                     do (fill-joint-internals joint internal))
               )))

