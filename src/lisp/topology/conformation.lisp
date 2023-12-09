(in-package :topology)


(defclass orientation ()
  ((to-origin :initarg :to-origin :accessor to-origin
              :documentation "This is the transform that takes the object to the origin")
   (from-origin :initarg :from-origin :accessor from-origin
                              :documentation "This is the second transform that takes the object after to-origin is applied")))

(defun make-orientation (&key (from-origin  (geom:make-matrix-identity))
                           (to-origin (geom:make-matrix-identity)))
  (make-instance 'orientation
                 :from-origin from-origin
                 :to-origin to-origin))

(defmethod print-object ((obj orientation) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (format stream ":from-origin~%~s~%:to-origin~%~s~%" (from-origin obj) (to-origin obj)))))

(defgeneric kin:orientation-transform (orientation))

(defmethod kin:orientation-transform ((orientation orientation))
  (let* ((from-origin (topology:from-origin orientation))
         (to-origin (to-origin orientation))
         (transform (geom:m*m from-origin to-origin)))
    transform))

(defclass orientations ()
  ((orientations :initform (make-hash-table) :initarg :orientations :accessor orientations)))

(defclass monomer-position ()
  ((molecule-index :initarg :molecule-index :reader molecule-index)
   (residue-index :initarg :residue-index :reader residue-index)))

(defmethod print-object ((obj monomer-position) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a" (molecule-index obj) (residue-index obj))))

(defgeneric at-position (obj pos)
  (:documentation "Return the object by following the monomer-position"))

(defmethod at-position ((obj chem:aggregate) pos)
  (let* ((mol (chem:content-at obj (molecule-index pos)))
         (res (chem:content-at mol (residue-index pos))))
    res))

(defmethod at-position ((ataggregate ataggregate) pos)
  (let* ((atmol (aref (atmolecules ataggregate) (molecule-index pos)))
         (atres (aref (atresidues atmol) (residue-index pos))))
  atres))


(defclass assembler-base (cando.serialize:serializable)
  ((monomer-positions :initarg :monomer-positions :accessor monomer-positions)
   (monomer-contexts :type hash-table :initarg :monomer-contexts :accessor monomer-contexts)
   (aggregate :initarg :aggregate :accessor aggregate)
   (energy-function :initarg :energy-function :accessor energy-function)
   (ataggregate :initarg :ataggregate :accessor ataggregate)
   (joint-tree :initarg :joint-tree :accessor joint-tree)))

(defclass assembler (assembler-base)
  ((oligomer-shapes :initarg :oligomer-shapes :accessor oligomer-shapes)
   (adjustments :initarg :adjustments :accessor adjustments)))


(defun aggregate* (assembler &optional coordinates)
  "Return a copy of the aggregate with coordinates if provided"
  (let ((agg (chem:matter-copy (topology:aggregate assembler))))
    (when coordinates (chem:matter-apply-coordinates agg coordinates))
    agg))

(defclass training-assembler (assembler-base)
  ((oligomers :initarg :oligomers :accessor oligomers)))

(defclass focused-training-assembler (training-assembler)
  ((focus-monomer :initarg :focus-monomer :accessor focus-monomer))
  (:documentation "An assembler with a focus monomer"))

(defun make-focused-training-assembler (oligomers focus-monomer)
  (let ((assembler (make-training-assembler oligomers :focus-monomer focus-monomer)))
    (change-class assembler 'focused-training-assembler :focus-monomer focus-monomer)
    assembler))

(defgeneric oligomer-containing-monomer (assembler monomer &optional errorp))

(defmethod oligomer-containing-monomer ((assembler training-assembler) monomer &optional errorp)
  "Return the oligomer that contains the monomer"
  (loop for oligomer in (oligomers assembler)
        do (loop for mon across (monomers oligomer)
                 when (eq mon monomer)
                   do (return-from oligomer-containing-monomer oligomer)))
  (if errorp
      (error "Could not find ~s in ~s" monomer assembler)
      nil))

(defmethod oligomer-containing-monomer ((assembler assembler) monomer &optional errorp)
  "Return the oligomer that contains the monomer"
  (loop for oligomer-shape in (oligomer-shapes assembler)
        for oligomer = (oligomer oligomer-shape)
        do (loop for mon across (monomers oligomer)
                 when (eq mon monomer)
                   do (return-from oligomer-containing-monomer oligomer)))
  (if errorp
      (error "Could not find ~s in ~s" monomer assembler)
      nil))
(defun make-coordinates-for-assembler (assembler)
  (let ((number-of-atoms (chem:number-of-atoms (topology:aggregate assembler))))
    (make-array (* 3 number-of-atoms) :element-type (geom:vecreal-type)
                :initial-element (geom:vecreal 0.0))))

(defun assembler-atresidue (assembler monomer)
  "Return the atresidue corresponding to the monomer"
  (let* ((pos (let ((p (gethash monomer (monomer-positions assembler))))
                (unless p
                  (error "The monomer ~s is not in the assembler" monomer))
                p))
         (mol-index (molecule-index pos))
         (res-index (residue-index pos))
         (atagg (ataggregate assembler))
         (atmol (aref (atmolecules atagg) mol-index))
         (atres (aref (atresidues atmol) res-index)))
    atres))

(defun focused-assembler-atresidue (focused-assembler)
  "Return the atresidue corresponding to the focus-monomer"
  (assembler-atresidue focused-assembler (focus-monomer focused-assembler)))

(defgeneric foldamer-monomer-context (focus-monomer oligomer foldamer)
  (:documentation "Return a monomer-context for a monomer in the oligomer using the foldamer.
Specialize the foldamer argument to provide methods"))

(defun make-assembler (oligomer-shapes &key monomer-order tune-energy-function)
  "Build a assembler for the oligomers.
energy-function-factory - provide a function that takes an aggregate and returns an energy-function."
  (unless (every (lambda (os) (typep os 'oligomer-shape)) oligomer-shapes)
    (error "You must provide a list of oligomer-shapes"))
  (let* ((full-orientations oligomer-shapes)
         (aggregate (chem:make-aggregate :all))
         (monomer-positions (make-hash-table))
         (oligomers (mapcar (lambda (oligomer-shape)
                              (oligomer oligomer-shape))
                            oligomer-shapes))
         (foldamers (mapcar (lambda (oligomer)
                              (let* ((oligomer-space (oligomer-space oligomer))
                                     (foldamer (foldamer oligomer-space)))
                                foldamer))
                            oligomers))
         (oligomer-molecules (loop for oligomer in oligomers
                                   for oligomer-space = (oligomer-space oligomer)
                                   for foldamer = (foldamer oligomer-space)
                                   for molecule-index from 0
                                   for molecule = (topology:build-molecule oligomer
                                                                           :aggregate aggregate
                                                                           :molecule-index molecule-index
                                                                           :monomer-positions-accumulator monomer-positions)
                                   do (chem:setf-force-field-name molecule (oligomer-force-field-name foldamer))
                                   collect (cons oligomer molecule)))
         (monomer-contexts (make-hash-table))
         (ataggregate (let ((atagg (make-instance 'ataggregate :aggregate aggregate)))
                        (resize-atmolecules atagg (length oligomers))
                        atagg))
         (joint-tree (make-joint-tree))
         (energy-function (let ((ef (chem:make-energy-function :keep-interaction keep-interaction :matter aggregate)))
                            (when tune-energy-function (funcall tune-energy-function ef))
                            ef))
         (adjustments (make-instance 'adjustments))
         (assembler (loop for oligomer-molecule in oligomer-molecules
                          for oligomer = (car oligomer-molecule)
                          for one-orientation in full-orientations
                          for molecule = (cdr oligomer-molecule)
                          for oligomer-space = (oligomer-space oligomer)
                          for foldamer = (foldamer oligomer-space)
                          for molecule-index from 0
                          do (loop for monomer across (monomers oligomer-space)
                                   for monomer-context = (foldamer-monomer-context monomer oligomer foldamer)
                                   do (setf (gethash monomer monomer-contexts) monomer-context))
                             ;; This is where I would invoke Conformation_O::buildMoleculeUsingOligomer
                             ;; Use the monomers-to-topologys
                          do (let* ((atmolecule (build-atmolecule-using-oligomer
                                                oligomer
                                                molecule
                                                molecule-index
                                                monomer-positions
                                                joint-tree
                                                (chem:atom-table energy-function)
                                                adjustments
                                                one-orientation)))
                               (put-atmolecule ataggregate atmolecule molecule-index))
                          finally (return (make-instance 'assembler
                                                         :monomer-positions monomer-positions
                                                         :monomer-contexts monomer-contexts
                                                         :oligomer-shapes oligomer-shapes
                                                         :aggregate aggregate
                                                         :energy-function energy-function
                                                         :ataggregate ataggregate
                                                         :joint-tree joint-tree
                                                         :adjustments adjustments)))))
    ;; The assembler is built - initialize the adjustments
    (loop for atmol across (atmolecules ataggregate)
          for atmol-index from 0
          do (loop for atres across (atresidues atmol)
                   for atres-index from 0
                   for internal-adjusts = (gethash atres (internal-adjustments adjustments))
                   for external-adjusts = (gethash atres (external-adjustments adjustments))
                   when internal-adjusts
                     do (loop for adjust in internal-adjusts
                              do (initialize-adjustment adjust assembler))
                   when external-adjusts
                     do (loop for adjust in external-adjusts
                              do (initialize-adjustment adjust assembler))))
    assembler))


(defun make-training-assembler (oligomers &key focus-monomer)
  "Build a assembler for the oligomers."
  (unless (every (lambda (os) (typep os 'oligomer)) oligomers)
    (error "You must provide a list of oligomers"))
  (let* ((aggregate (chem:make-aggregate :all))
         (monomer-positions (make-hash-table))
         (oligomer-molecules (loop for oligomer in oligomers
                                   for molecule-index from 0
                                   for molecule = (topology:build-molecule oligomer
                                                                           :aggregate aggregate
                                                                           :molecule-index molecule-index
                                                                           :monomer-positions-accumulator monomer-positions)
                                   collect (cons oligomer molecule)))
         (monomer-contexts (make-hash-table))
         (ataggregate (let ((atagg (make-instance 'ataggregate :aggregate aggregate)))
                        (resize-atmolecules atagg (length oligomers))
                        atagg))
         (joint-tree (make-joint-tree))
         (energy-function (chem:make-energy-function :matter aggregate))
         (assembler (loop for oligomer-molecule in oligomer-molecules
                          for oligomer = (car oligomer-molecule)
                          for molecule = (cdr oligomer-molecule)
                          for oligomer-space = (oligomer-space oligomer)
                          for foldamer = (foldamer oligomer-space)
                          for molecule-index from 0
                          do (if focus-monomer
                                 (progn
                                   (let ((focus-monomer-context (foldamer-monomer-context focus-monomer oligomer foldamer)))
                                     (setf (gethash focus-monomer monomer-contexts) focus-monomer-context)
                                     ))
                                 (loop for monomer across (monomers oligomer-space)
                                       for monomer-context = (foldamer-monomer-context monomer oligomer foldamer)
                                       do (setf (gethash monomer monomer-contexts) monomer-context)))
                             ;; This is where I would invoke Conformation_O::buildMoleculeUsingOligomer
                             ;; Use the monomers-to-topologys
                          do (let* ((atmolecule (build-atmolecule-using-oligomer
                                                 oligomer
                                                 molecule
                                                 molecule-index
                                                 monomer-positions
                                                 joint-tree
                                                 (chem:atom-table energy-function)
                                                 nil
                                                 nil)))
                               (put-atmolecule ataggregate atmolecule molecule-index))
                          finally (return (make-instance 'training-assembler
                                                         :monomer-positions monomer-positions
                                                         :monomer-contexts monomer-contexts
                                                         :oligomers oligomers
                                                         :aggregate aggregate
                                                         :energy-function energy-function
                                                         :ataggregate ataggregate
                                                         :joint-tree joint-tree)))))
    assembler))

(defun walk-atoms-joints (assembler molecule-index callback)
  (let* ((aggregate (aggregate assembler))
        (ataggregate (ataggregate assembler))
        (molecule (chem:content-at aggregate molecule-index))
         (atmolecule (aref (atmolecules ataggregate) molecule-index)))
    (loop for residue-index below (chem:content-size molecule)
          for residue = (chem:content-at molecule residue-index)
          for atresidue = (aref (atresidues atmolecule) residue-index)
          do (loop for atom-index below (chem:content-size residue)
                   for atom = (chem:content-at residue atom-index)
                   for joint = (aref (joints atresidue) atom-index)
                   do (funcall callback atom joint (list molecule-index residue-index atom-index))))))

(defun old-walk-atoms-joints (assembler molecule-index callback)
  (let ((aggregate (aggregate assembler))
        (ataggregate (ataggregate assembler)))
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
(defun copy-atom-positions-into-joints (assembler)
  (walk-atoms-joints assembler
                     (lambda (atm jnt atomid)
                       (declare (ignore atomid))
                       (kin:set-position jnt (chem:get-position atm)))))

(defun find-joint-for-atom (assembler atom)
  (let ((ataggregate (ataggregate assembler))
        (aggregate (aggregate assembler)))
    (loop for atmol across (atmolecules ataggregate)
          for imol from 0
          for mol = (chem:content-at aggregate imol)
          do (loop for atres across (atresidues atmol)
                   for ires from 0
                   for res = (chem:content-at mol ires)
                   do (loop for joint across (joints atres)
                            for iatm from 0
                            for atm = (chem:content-at res iatm)
                            when (eq atm atom)
                              do (return-from find-joint-for-atom joint))))
    (error "Could not find atom ~s" atom)
    ))

(defun copy-joint-positions-into-atoms (assembler coords oligomer-shape)
  (let ((pos (position oligomer-shape (oligomer-shapes assembler))))
    (walk-atoms-joints assembler pos
                       (lambda (atm jnt atomid)
                         (declare (ignore atomid))
                         (chem:set-position atm (kin:joint/position jnt coords))))))

(defun copy-all-joint-positions-into-atoms (assembler coords)
  (loop for index below (length (oligomer-shapes assembler))
        do (walk-atoms-joints assembler index
                     (lambda (atm jnt atomid)
                       (declare (ignore atomid))
                       (chem:set-position atm (kin:joint/position jnt coords))))))

(defun update-atresidue-joint-tree-internal-coordinates (assembler atresidue coordinates)
  (walk-atresidue-joints atresidue
                         (lambda (joint atom-id)
                           (declare (ignore atom-id))
                           (kin:update-internal-coord joint coordinates)
                           )))

(defun update-atmolecule-joint-tree-internal-coordinates (assembler atmolecule coordinates)
  (walk-atmolecule-joints atmolecule
                          (lambda (joint atom-id)
                            (declare (ignore atom-id))
                            (kin:update-internal-coord joint coordinates))))

(defun update-ataggregate-joint-tree-internal-coordinates (assembler coordinates)
  (let ((ataggregate (ataggregate assembler)))
    (walk-ataggregate-joints ataggregate
                             (lambda (joint atom-id)
                               (declare (ignore atom-id))
                               (kin:update-internal-coord joint coordinates)))))

(defun build-external-coordinates (assembler &key oligomer-shape (coords (topology:make-coordinates-for-assembler assembler)))
  (if oligomer-shape
      (progn
        (build-atom-tree-external-coordinates assembler coords oligomer-shape)
        (adjust-atom-tree-external-coordinates assembler coords oligomer-shape))
      (loop for oligomer-shape in (oligomer-shapes assembler)
            do (build-external-coordinates assembler :oligomer-shape oligomer-shape :coords coords))))

(defun build-atom-tree-external-coordinates (assembler coords oligomer-shape)
  (let* ((one-oligomer (oligomer oligomer-shape))
         (joint (gethash one-oligomer (root-map (joint-tree assembler)))))
    (unless joint
      (error "Could not find oligomer ~s in root-map ~s" one-oligomer (root-map (joint-tree assembler))))
    (kin:update-xyz-coords joint coords)))

(defun build-all-external-coordinates (assembler &key (coords (topology:make-coordinates-for-assembler assembler)))
  (loop for oligomer-shape in (oligomer-shapes assembler)
        do (build-atom-tree-external-coordinates assembler coords oligomer-shape)))

(defun adjust-atom-tree-external-coordinates (assembler coords oligomer-shape)
  (let* ((pos (position oligomer-shape (oligomer-shapes assembler)))
         (atmol (elt (atmolecules (ataggregate assembler)) pos)))
    (loop for atres across (atresidues atmol)
          for adjustments = (gethash atres (external-adjustments (adjustments assembler)))
          do (loop for adjustment in adjustments
                   do (external-adjust adjustment assembler coords)))))

(defun adjust-all-atom-tree-external-coordinates (assembler coords)
  (loop for oligomer-shape in (oligomer-shapes assembler)
        do (adjust-atom-tree-external-coordinates assembler coords oligomer-shape)))

(defun build-all-atom-tree-external-coordinates-and-adjust (assembler coords)
  (build-all-atom-tree-external-coordinates assembler coords)
  (adjust-all-atom-tree-external-coordinates assembler coords))

(defun build-atom-tree-external-coordinates-and-adjust (assembler coords oligomer-shape)
  (build-atom-tree-external-coordinates assembler coords oligomer-shape)
  (adjust-atom-tree-external-coordinates assembler coords oligomer-shape))

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

(defun adjust-internals (assembler oligomer-shape)
  "Do the adjustment of internal coordinates"
  (let* ((pos (position oligomer-shape (oligomer-shapes assembler)))
         (atmol (aref (atmolecules (ataggregate assembler)) pos)))
    (loop for atres across (atresidues atmol)
          for adjustments = (gethash atres (internal-adjustments (adjustments assembler)))
          do (loop for adjustment in adjustments
                   do (internal-adjust adjustment assembler)))))

(defun adjust-all-internals (assembler)
  "Do the adjustment of internal coordinates"
  (loop for oligomer-shape in (oligomer-shapes assembler)
        do (adjust-internals assembler oligomer-shape)))

(defun fill-internals-from-oligomer-shape (assembler oligomer-shape &optional monomers)
  "Fill internal coordinates from the fragments"
  (let ((fragments (rotamers-map oligomer-shape)))
    (loop for ass-oligomer-shape in (oligomer-shapes assembler)
          for oligomer = (oligomer ass-oligomer-shape)
          when (eq ass-oligomer-shape oligomer-shape)
            do (loop with atagg = (ataggregate assembler)
                     for monomer in (ordered-monomers oligomer)
                     when (or (null monomers)
                              (member monomer monomers))
                       do (let* ((monomer-context (gethash monomer (monomer-contexts assembler)))
                                 (context-rotamers (let ((frag (gethash monomer-context (monomer-context-to-context-rotamers fragments))))
                                                           (unless frag
                                                             (error "Could not find monomer-context ~a" monomer-context))
                                                           frag))
                                 (monomer-position (gethash monomer (monomer-positions assembler)))
                                 (molecule-index (molecule-index monomer-position))
                                 (monomer-index (residue-index monomer-position))
                                 (atmol (elt (atmolecules atagg) molecule-index))
                                 (atres (elt (atresidues atmol) monomer-index))
                                 (monomer-shape (let ((ms (gethash monomer (monomer-shape-map oligomer-shape))))
                                                  (unless ms (error "Could not get monomer-shape for monomer ~a" monomer))
                                                  ms))
                                 (fragment-conformation-index (let ((fi (fragment-conformation-index monomer-shape)))
                                                                (unless fi (break "The fragment-conformation-index is nil - you must call random-fragment-conformation-index on the oligomer-shape before you do anything else: ~a" monomer-shape))
                                                                fi))
                                 (fragment-internals (let ((fragments (rotamers context-rotamers)))
                                                       (unless (> (length fragments) 0)
                                                         (error "fragments is empty for context ~a" monomer-context))
                                                       (unless (< fragment-conformation-index (length fragments))
                                                         (break "Check fragments ~a because fragment-conformation-index ~a is out of bounds"
                                                                fragments fragment-conformation-index))
                                                       (elt fragments fragment-conformation-index)))
                                 )
                            (apply-fragment-internals-to-atresidue fragment-internals atres)
                            )))))

#|
;;;Idea to make monomer-shape subclasses control how internal coordinates get generated.

(defmethod apply-monomer-shape-to-atresidue ((monomer-shape monomer-shape) atresidue context-rotamers)
  (let* ((fragments (rotamers context-rotamers))
         (fragment-conformation-index (fragment-conformation-index monomer-shape))
         (fragment-internals (elt fragments fragment-conformation-index)))
    (apply-fragment-internals-to-atresidue fragment-internals atresidue)))
|#


(defun fill-internals-from-oligomer-shape-and-adjust (assembler oligomer-shape &optional monomers)
  (fill-internals-from-oligomer-shape assembler oligomer-shape monomers)
  (adjust-internals assembler oligomer-shape))

(defun search-conformations (oligomer-space context-rotamers monomer-names sdf-filename &optional (number-struct 50) (number-conf 100))
  (error "Implement search-conformations")
  #+(or)(with-open-file (fout sdf-filename :direction :output)
    (loop for struct-count below number-struct
          for rand-sequence = (random (topology:number-of-sequences oligomer-space))
          for _a = (format t "rand-sequence ~a~%" rand-sequence)
          for oligomer = (topology:make-oligomer oligomer-space rand-sequence)
          for conf = (topology:make-assembler oligomer)
          do (topology::fill-internals-from-fragments conf context-rotamers 0)
             (loop for count below number-conf
                   do (topology::fill-internals-from-fragments-for-monomers-named conf context-rotamers monomer-names)
                   do (topology:zero-all-atom-tree-external-coordinates conf)
                   do (topology:build-all-atom-tree-external-coordinates conf)
                   do (topology:copy-joint-positions-into-atoms conf)
                   do (sdf:write-sdf-stream (topology:aggregate conf) fout)))))

#+(or)
(defun fill-internals-from-fragments-for-monomers-named (assembler fragments monomer-names)
  "Fill internal coordinates from the fragments"
  (loop for ass-oligomer-shape in (oligomer-shapes assembler)
        for oligomer = (oligomer ass-oligomer-shape)
        do (loop with atagg = (ataggregate assembler)
                 for monomer across (monomers oligomer)
                 when (member (topology:current-stereoisomer-name monomer oligomer) monomer-names)
                   do (let* ((monomer-context (gethash monomer (monomer-contexts assembler)))
                             (frags (gethash monomer-context (monomer-context-to-context-rotamers fragments)))
                             (monomer-position (gethash monomer (monomer-positions assembler)))
                             (monomer-index (residue-index monomer-position))
                             (molecule-index (molecule-index monomer-position))
                             (atmol (elt (atmolecules atagg) molecule-index))
                             (atres (elt (atresidues atmol) monomer-index))
                             (context-rotamers (gethash monomer-context (monomer-context-to-context-rotamers fragments)))
                             (rand-limit (length (fragments context-rotamers)))
                             (rand-index (random rand-limit))
                             (fragment-internals (elt (fragments context-rotamers) rand-index)))
                        (loop for joint across (joints atres)
                              for internal in (internals fragment-internals)
                              do (fill-joint-internals joint internal))
                        ))))


(defun copy-externals (assembler monomers-to-residues)
  (maphash (lambda (monomer residue)
             (let* ((monomer-position (let ((mp (gethash monomer (topology:monomer-positions assembler))))
                                        (unless mp (error "Could not find monomer ~s in monomer-positions of assembler ~s" monomer (alexandria:hash-table-keys (topology:monomer-positions assembler))))
                                        mp))
                    (molecule-index (topology:molecule-index monomer-position))
                    (residue-index (topology:residue-index monomer-position))
                    (assembler-molecule (aref (topology:atmolecules (topology:ataggregate assembler)) molecule-index))
                    (assembler-atresidue (aref (topology:atresidues assembler-molecule) residue-index))
                    (assembler-residue (topology:residue assembler-atresidue)))
               (chem:do-atoms (atm assembler-residue)
                 (let* ((source-atm (chem:atom-with-name residue (chem:get-name atm)))
                        (source-pos (chem:get-position source-atm)))
                   (chem:set-position atm source-pos)))))
           monomers-to-residues))
