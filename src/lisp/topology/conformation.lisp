(in-package :topology)


(defclass orientation ()
  ((to-origin :initarg :to-origin :accessor to-origin
              :documentation "This is the transform that takes the object to the origin")
   (from-origin :initarg :from-origin :accessor from-origin
                              :documentation "This is the second transform that takes the object after to-origin is applied")))

(defun make-orientation (&key (from-origin (geom:make-matrix-identity))
                           (to-origin (geom:make-matrix-identity)))
  (make-instance 'orientation
                 :from-origin from-origin
                 :to-origin to-origin))

(defun copy-orientation (orientation)
  (make-instance 'orientation
                 :from-origin (geom:copy-matrix (from-origin orientation))
                 :to-origin (geom:copy-matrix (to-origin orientation))))

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
   (monomer-subset :initform nil :initarg :monomer-subset :accessor monomer-subset)
   (adjustments :initarg :adjustments :accessor adjustments)))

(defclass subset-assembler (assembler)
  ()
  (:documentation "This class is an assembler that builds a subset of another assembler"))

(defun aggregate* (assembler coordinates &key (name :all))
  "Return a copy of the aggregate with coordinates if provided"
  (let ((agg (chem:matter-copy (topology:aggregate assembler))))
    (when coordinates (chem:matter-apply-coordinates agg coordinates))
    (chem:set-name agg name)
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

(defun make-coordinates-for-number-of-atoms (number-of-atoms)
    (make-array (* 3 number-of-atoms) :element-type (geom:vecreal-type)
                                      :initial-element (geom:vecreal 0.0)))

(defun make-coordinates-for-assembler (assembler)
  (let ((number-of-atoms (chem:number-of-atoms (topology:aggregate assembler))))
    (make-coordinates-for-number-of-atoms number-of-atoms)))

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

(defgeneric foldamer-rotamers-database (foldamer)
  (:documentation "Return the rotamers-database for the foldamer"))

(defclass monomer-subset ()
  ((monomers :initarg :monomers :reader monomers)))

(defclass monomer-subset-with-assembler (monomer-subset)
  ((backbone-assembler :initarg :backbone-assembler :reader backbone-assembler)
   (coordinates :initarg :coordinates :reader coordinates)))

(defun make-monomer-subset-with-assembler (backbone-assembler coordinates list-of-monomers)
  "Make a monomer-subset, currently a hash-table with keys of monomers and values of T"
  (let ((ht (make-hash-table)))
    (loop for monomer in list-of-monomers
          for backbone-pos = (gethash monomer (monomer-positions backbone-assembler))
          for backbone-atresidue = (topology:at-position (topology:ataggregate backbone-assembler) backbone-pos)
          for first-joint = (aref (topology:joints backbone-atresidue) 0)
          for backbone-parent = (kin:parent first-joint)
          do (setf (gethash monomer ht) backbone-parent))
    (make-instance 'monomer-subset-with-assembler
                   :backbone-assembler backbone-assembler
                   :coordinates coordinates
                   :monomers ht)))

(defun in-monomer-subset (monomer-subset &rest monomers)
  (cond
    ((null monomer-subset)
     ;; If monomer-subset is NULL then everything is in the subset
     t)
    ((typep monomer-subset 'monomer-subset)
     (every (lambda (mon) (gethash mon (monomers monomer-subset))) monomers))
    (t (error "Illegal value for monomer-subset ~s - must be NIL or a hash-table"))))

(defun make-assembler (oligomer-shapes &key monomer-subset tune-energy-function (keep-interaction t))
  "Build a assembler for the oligomers.
energy-function-factory - provide a function that takes an aggregate and returns an energy-function.
tune-energy-function - A function that takes the energy-function and an assembler and modifies the energy-function."
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
         (monomers-to-residues (make-hash-table))
         (oligomer-molecules (loop for oligomer in oligomers
                                   for oligomer-space = (oligomer-space oligomer)
                                   for foldamer = (foldamer oligomer-space)
                                   for molecule-index from 0
                                   for molecule = (topology:build-molecule oligomer
                                                                           :monomer-subset monomer-subset
                                                                           :aggregate aggregate
                                                                           :molecule-index molecule-index
                                                                           :monomers-to-residues monomers-to-residues
                                                                           :monomer-positions-accumulator monomer-positions)
                                   do (chem:setf-force-field-name molecule (oligomer-force-field-name foldamer))
                                   collect (cons oligomer molecule))))
    (let* ((monomer-contexts (make-hash-table))
           (ataggregate (let ((atagg (make-instance 'ataggregate :aggregate aggregate)))
                          (resize-atmolecules atagg (length oligomers))
                          atagg))
           (joint-tree (make-joint-tree))
           (adjustments (make-instance 'adjustments))
           (energy-function (chem:make-energy-function :keep-interaction keep-interaction :matter aggregate))
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
                            do (let* ((atmolecule (build-atmolecule-using-oligomer oligomer
                                                                                   monomer-subset
                                                                                   molecule
                                                                                   molecule-index
                                                                                   monomer-positions
                                                                                   joint-tree
                                                                                   (chem:atom-table energy-function)
                                                                                   adjustments
                                                                                   one-orientation)))
                                 (put-atmolecule ataggregate atmolecule molecule-index))
                            finally (return (make-instance (if monomer-subset
                                                               'subset-assembler
                                                               'assembler)
                                                           :monomer-positions monomer-positions
                                                           :monomer-contexts monomer-contexts
                                                           :oligomer-shapes oligomer-shapes
                                                           :aggregate aggregate
                                                           :ataggregate ataggregate
                                                           :monomer-subset monomer-subset
                                                           :joint-tree joint-tree
                                                           :adjustments adjustments)))))
      ;; The energy-function may be adjusted for the assembler
      (when tune-energy-function
        (funcall tune-energy-function energy-function assembler))
      (setf (energy-function assembler) energy-function)
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
      assembler)))


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
                                                 nil
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
         (joints (gethash one-oligomer (root-map (joint-tree assembler)))))
    (when (null joints)
      (error "Could not find oligomer ~s in root-map ~s" one-oligomer (root-map (joint-tree assembler))))
    (loop for joint in joints
          do (kin:update-xyz-coords joint coords))))

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
  (loop for oligomer-shape in (oligomer-shapes assembler)
        do (build-atom-tree-external-coordinates assembler coords oligomer-shape)
        do (adjust-atom-tree-external-coordinates assembler coords oligomer-shape)))

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

(defgeneric apply-monomer-shape-to-atresidue-internals (assembler oligomer-shape monomer-shape monomer-context atresidue coordinates &key verbose)
  (:documentation "Specialize this on different monomer-shape classes.
Fill in internal coordinates into the atresidue for the monomer-shape.
Some specialized methods will need coordinates for the assembler"))



(defun fill-internals-from-oligomer-shape (assembler oligomer-shape &optional verbose)
  "Fill internal coordinates from the fragments"
  (when verbose (let ((*print-pretty* nil)) (format t "fill-internals-from-oligomer-shape ~s~%" oligomer-shape)))
  (let ((coordinates (topology:make-coordinates-for-assembler assembler)))
    (loop for ass-oligomer-shape in (oligomer-shapes assembler)
          when (eq ass-oligomer-shape oligomer-shape)
            do (let* ((oligomer (oligomer ass-oligomer-shape))
                      (ordered-monomers (ordered-monomers oligomer)))
                 (when verbose
                   (let ((*print-pretty* nil))
                     (format t "fill-internals-from-oligomer-shape ~s ordered-monomers: ~s~%" ass-oligomer-shape ordered-monomers)))
                 (loop with atagg = (ataggregate assembler)
                       ;; It's really important that we use the ordered-monomers so that the monomer-shapes
                       ;;  will install internal coordinates in the order from the root outwards.
                       ;;  so that any preceeding monomer-shapes are built before any following ones.
                       for monomer in ordered-monomers
                       when (in-monomer-subset (monomer-subset assembler) monomer)
                         do (let* ((monomer-context (gethash monomer (monomer-contexts assembler)))
                                   (monomer-position (gethash monomer (monomer-positions assembler)))
                                   (molecule-index (molecule-index monomer-position))
                                   (residue-index (residue-index monomer-position))
                                   (atmol (elt (atmolecules atagg) molecule-index))
                                   (atres (elt (atresidues atmol) residue-index))
                                   (monomer-shape (let ((ms (gethash monomer (monomer-shape-map oligomer-shape))))
                                                    (unless ms (error "Could not get monomer-shape for monomer ~a" monomer))
                                                    ms)))
                              (when verbose (format t "applying internals for monomer: ~s~%" monomer))
                              (apply-monomer-shape-to-atresidue-internals assembler oligomer-shape monomer-shape monomer-context atres coordinates :verbose verbose)))))))

#|
;;;Idea to make monomer-shape subclasses control how internal coordinates get generated.

(defmethod apply-monomer-shape-to-atresidue ((monomer-shape monomer-shape) atresidue context-rotamers)
  (let* ((fragments (rotamers context-rotamers))
         (rotamer-index (rotamer-index monomer-shape))
         (fragment-internals (elt fragments rotamer-index)))
    (apply-fragment-internals-to-atresidue fragment-internals atresidue)))
|#


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


(defmacro do-atresidue-residue ((atresidue residue assembler) &body body)
  "Iterate over each atresidue and residue pair in the assembler"
  (let ((oligomer-shape (gensym "oligomer-shape"))
        (monomer-shape-info (gensym "monomer-shape-info"))
        (monomer (gensym "monomer"))
        (monomer-pos (gensym "monomer-pos"))
        )
    `(loop for ,oligomer-shape in (topology:oligomer-shapes ,assembler)
           do (loop for ,monomer-shape-info across (topology:monomer-shape-info-vector ,oligomer-shape)
                    for ,monomer = (topology:monomer ,monomer-shape-info)
                    for ,monomer-pos = (gethash ,monomer (topology:monomer-positions ,assembler))
                    for ,atresidue = (topology:at-position (topology:ataggregate ,assembler) ,monomer-pos)
                    for ,residue = (topology:at-position (topology:aggregate ,assembler) ,monomer-pos)
                    do (progn
                         ,@body)))))

(defmacro do-joint-atom ((joint atm assembler) &body body)
  "Iterate over each atresidue and residue pair in the assembler"
  (let ((atresidue (gensym "atresidue"))
        (residue (gensym "residue"))
        (name (gensym "name")))
    (declare (ignorable joint atm))
    `(do-atresidue-residue (,atresidue ,residue ,assembler)
       (loop for ,joint across (topology:joints ,atresidue)
             for ,name = (kin:joint/name ,joint)
             for ,atm = (chem:atom-with-name ,residue ,name)
             do (progn
                  ,@body)))))

(defun copy-atom-positions-to-xyz-joint-pos (assembler)
  (do-joint-atom (joint atm assembler)
    (when (typep joint 'kin:xyz-joint)
      (let ((pos (chem:get-position atm)))
        (when (geom:is-defined pos)
          (kin:xyz-joint/set-pos joint pos))))))

(defun copy-xyz-joint-pos-to-atom-positions (assembler)
  (do-joint-atom (joint atm assembler)
    (when (typep joint 'kin:xyz-joint)
      (let ((pos (kin:xyz-joint/get-pos joint)))
        (when (geom:is-defined pos)
          (chem:set-position atm pos))))))



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

(defmethod analyze-assembler ((assembler assembler) &optional name)
  (loop for oligomer-shape in (oligomer-shapes assembler)
        for oligomer-shape-index from 0
        do (analyze-oligomer-shape oligomer-shape))
  (let (
        (backbone-internals-count (make-array (length (oligomer-shapes assembler)) :initial-element 0))
        (backbone-internals-defined (make-array (length (oligomer-shapes assembler)) :initial-element 0))
        (sidechain-internals-count (make-array (length (oligomer-shapes assembler)) :initial-element 0))
        (sidechain-internals-defined (make-array (length (oligomer-shapes assembler)) :initial-element 0))
        )
  (loop for oligomer-shape in (oligomer-shapes assembler)
        for oligomer-shape-index from 0
        for backbone-internals-count = 0
        for sidechain-internals-count = 0
        for backbone-internals-defined = 0
        for sidechain-internals-defined = 0
        do (loop for monomer-shape across (monomer-shape-vector oligomer-shape)
                 for monomer-shape-info across (monomer-shape-info-vector oligomer-shape)
                 for monomer = (monomer monomer-shape-info)
                 for monomer-shape-kind = (monomer-shape-kind monomer-shape-info)
                 for monomer-pos = (gethash monomer (monomer-positions assembler))
                 for atresidue = (topology:at-position (topology:ataggregate assembler) monomer-pos)
                 if (eq :backbone monomer-shape-kind)
                   do (loop for joint across (joints atresidue)
                            do (incf backbone-internals-count)
                                when (kin:joint/definedp joint)
                                  do (incf backbone-internals-defined))
                 else
                   do (loop for joint across (joints atresidue)
                            do (incf sidechain-internals-count)
                            when (kin:joint/definedp joint)
                                 do (incf sidechain-internals-defined)))
        do (analyze-oligomer-shape oligomer-shape)
        do (format t "  backbone-internals-count ~d~%" backbone-internals-count)
        do (format t "  backbone-internals-defined ~d~%" backbone-internals-defined)
        do (format t "  sidechain-internals-count ~d~%" sidechain-internals-count)
        do (format t "  sidechain-internals-defined ~d~%" sidechain-internals-defined)
           )))

(defun assembler-check-joints (assembler)
  (loop for atmol across (atmolecules (ataggregate assembler))
        do (loop for atres across (atresidues atmol)
                 for undefs = (loop for joint across (joints atres)
                                    when (not (kin:joint/definedp joint))
                                      collect joint)
                 do (when undefs
                      (format t "atres: ~s~%" atres)
                      (loop for joint in undefs
                            do (format t "   joint: ~s undefined~%" joint))))))
