(in-package :topology)

(defstruct monomer-position
  molecule-index
  residue-index)

(defclass conformation ()
  ((monomer-positions :initarg :monomer-positions :accessor monomer-positions)
   (monomer-contexts :type hash-table :initarg :monomer-contexts :accessor monomer-contexts)
   (oligomers :initarg :oligomers :accessor oligomers)
   (aggregate :initarg :aggregate :accessor aggregate)
   (energy-function :initarg :energy-function :accessor energy-function)
   (ataggregate :initarg :ataggregate :accessor ataggregate)
   (joint-tree :initarg :joint-tree :accessor joint-tree)))

(cando.serialize:make-class-save-load conformation)

(defgeneric foldamer-monomer-context (focus-monomer oligomer foldamer)
  (:documentation "Return a monomer-context for a monomer in the oligomer using the foldamer.
Specialize the foldamer argument to provide methods"))

(defun make-conformation (oligomers &key monomer-order)
  "Build a conformation for the oligomers."
  (let* ((aggregate (chem:make-aggregate :all))
         (monomer-positions (make-hash-table))
         (oligomer-molecules (loop for oligomer in oligomers
                                   for molecule-index from 0
                                   for molecule = (topology:build-molecule oligomer aggregate molecule-index monomer-positions)
                                   collect (cons oligomer molecule)))
         (monomer-contexts (make-hash-table))
         (ataggregate (let ((atagg (make-instance 'ataggregate :aggregate aggregate)))
                        (resize-atmolecules atagg (length oligomers))
                        atagg))
         (joint-tree (make-joint-tree))
         (energy-function (chem:make-energy-function :matter aggregate))
         )
    (loop for oligomer-molecule in oligomer-molecules
          for oligomer = (car oligomer-molecule)
          for molecule = (cdr oligomer-molecule)
          for oligomer-space = (oligomer-space oligomer)
          for foldamer = (foldamer oligomer-space)
          for molecule-index from 0
          do (loop for monomer across (monomers oligomer-space)
                   for monomer-context = (foldamer-monomer-context monomer oligomer foldamer)
                   do (setf (gethash monomer monomer-contexts) monomer-context))
             ;; This is where I would invoke Conformation_O::buildMoleculeUsingOligomer
             ;; Use the monomers-to-topologys
          do (let ((atmolecule (build-atmolecule-using-oligomer oligomer molecule molecule-index monomer-positions joint-tree (chem:atom-table energy-function))))
               (put-atmolecule ataggregate atmolecule molecule-index))
          finally (return (make-instance 'conformation
                                         :monomer-positions monomer-positions
                                         :monomer-contexts monomer-contexts
                                         :oligomers oligomers
                                         :aggregate aggregate
                                         :energy-function energy-function
                                         :ataggregate ataggregate
                                         :joint-tree joint-tree)))
    ))

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

#+(or)
(defun copy-atom-positions-into-joints (conformation)
  (walk-atoms-joints conformation
                     (lambda (atm jnt atomid)
                       (declare (ignore atomid))
                       (kin:set-position jnt (chem:get-position atm)))))

#+(or)(defun copy-joint-positions-into-atoms (conformation)
  (walk-atoms-joints conformation
                     (lambda (atm jnt atomid)
                       (declare (ignore atomid))
                       (chem:set-position atm (kin:position coords jnt)))))

(defun update-joint-tree-internal-coordinates (conformation coordinates)
  (let ((ataggregate (ataggregate conformation)))
    (walk-ataggregate-joints ataggregate
                             (lambda (joint atom-id)
                               (declare (ignore atom-id))
                               (kin:update-internal-coord joint coordinates)))))

(defun build-all-atom-tree-external-coordinates (conformation oligomer coords)
  (loop for one-oligomer in (oligomers conformation)
        for joint in (roots (joint-tree conformation))
        when (eq oligomer one-oligomer)
          do (kin:update-xyz-coords joint coords)))

(defun build-atresidue-atom-tree-external-coordinates (atresidue coords)
  (loop for atom-index below (length (joints atresidue))
        for joint = (aref (joints atresidue) atom-index)
        do (kin:update-xyz-coords joint coords)))

(defun zero-all-atom-tree-external-coordinates (conf)
  "Set the external coordinates for each joint to the origin"
  (error "Implement me")
  (let ((ataggregate (ataggregate conf)))
    (walk-ataggregate-joints ataggregate
                             (lambda (joint atom-id)
                               (declare (ignore atom-id))
                               (kin:set-position joint (geom:vec 0.0 0.0 0.0))))))

(defun fill-internals-from-oligomer-shape (conf oligomer-shape)
  "Fill internal coordinates from the fragments"
  (let ((fragments (matched-fragment-conformations-map oligomer-shape)))
    (loop for oligomer in (oligomers conf)
          when (eq oligomer (oligomer oligomer-shape))
            do (loop with atagg = (ataggregate conf)
                     for monomer in (ordered-monomers oligomer)
                     for monomer-context = (gethash monomer (monomer-contexts conf))
                     for fragment-conformations = (let ((frag (gethash monomer-context (monomer-context-to-fragment-conformations fragments))))
                                                    (unless frag
                                                      (error "Could not find monomer-context ~a" monomer-context))
                                                    frag)
                     for monomer-position = (gethash monomer (monomer-positions conf))
                     for molecule-index = (monomer-position-molecule-index monomer-position)
                     for monomer-index = (monomer-position-residue-index monomer-position)
                     for atmol = (elt (atmolecules atagg) molecule-index)
                     for atres = (elt (atresidues atmol) monomer-index)
                     for monomer-shape = (let ((ms (gethash monomer (monomer-shape-map oligomer-shape))))
                                           (unless ms (error "Could not get monomer-shape for monomer ~a" monomer))
                                           ms)
                     for fragment-conformation-index = (let ((fi (fragment-conformation-index monomer-shape)))
                                                         (unless fi (break "The fragment-conformation-index is nil - you must call random-fragment-conformation-index on the oligomer-shape before you do anything else: ~a" monomer-shape))
                                                         fi)
                     for fragment-internals = (let ((fragments (fragments fragment-conformations)))
                                                (unless fragments
                                                  (error "fragments is NIL for context ~a" monomer-context))
                                                (elt (fragments fragment-conformations) fragment-conformation-index))
                     do (loop for joint across (joints atres)
                              for internal across (internals fragment-internals)
                              do (fill-joint-internals joint internal))
                     ))))



(defun search-conformations (oligomer-space fragment-conformations monomer-names sdf-filename &optional (number-struct 50) (number-conf 100))
  (error "Implement search-conformations")
  #+(or)(with-open-file (fout sdf-filename :direction :output)
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
  (loop for oligomer in (oligomers conf)
        do (loop with atagg = (ataggregate conf)
                 for monomer across (monomers oligomer)
                 when (member (topology:current-stereoisomer-name monomer oligomer) monomer-names)
                   do (let* ((monomer-context (gethash monomer (monomer-contexts conf)))
                             (frags (gethash monomer-context (monomer-context-to-fragment-conformations fragments)))
                             (monomer-position (gethash monomer (monomer-positions conf)))
                             (monomer-index (monomer-position-residue-index monomer-position))
                             (molecule-index (monomer-position-molecule-index monomer-position))
                             (atmol (elt (atmolecules atagg) molecule-index))
                             (atres (elt (atresidues atmol) monomer-index))
                             (fragment-conformations (gethash monomer-context (monomer-context-to-fragment-conformations fragments)))
                             (rand-limit (length (fragments fragment-conformations)))
                             (rand-index (random rand-limit))
                             (fragment-internals (elt (fragments fragment-conformations) rand-index)))
                        (loop for joint across (joints atres)
                              for internal in (internals fragment-internals)
                              do (fill-joint-internals joint internal))
                        ))))


(defun copy-externals (conf monomers-to-residues)
  (maphash (lambda (monomer residue)
             (let* ((monomer-position (gethash monomer (topology:monomer-positions conf)))
                    (molecule-index (topology:monomer-position-molecule-index monomer-position))
                    (residue-index (topology:monomer-position-residue-index monomer-position))
                    (conf-molecule (aref (topology:atmolecules (topology:ataggregate conf)) molecule-index))
                    (conf-atresidue (aref (topology:atresidues conf-molecule) residue-index))
                    (conf-residue (topology:residue conf-atresidue)))
               (chem:do-atoms (atm conf-residue)
                 (let* ((source-atm (chem:atom-with-name residue (chem:get-name atm)))
                        (source-pos (chem:get-position source-atm)))
                   (chem:set-position atm source-pos)))))
           monomers-to-residues))
                   
                   
