(in-package :topology)


(defclass orientation ()
  ((to-origin :initarg :to-origin :reader to-origin
              :documentation "This is the transform that takes the object to the origin")
   (adjust-rotation :initarg :adjust-rotation :reader adjust-rotation
           :documentation "Apply this to the coordinates after the to-origin transform was applied")
   (adjust-translation :initarg :adjust-translation :reader adjust-translation
           :documentation "Apply this to the coordinates after the adjust-rotation transform was applied")
   (from-origin :initarg :from-origin :reader from-origin
                :documentation "This is the second transform that takes the object after adjusting to-origin is applied")
   (transform-cache :initarg :transform-cache :reader transform-cache))
  (:documentation "Represent the orientation of a molecule in a design calculation.
There are three transforms that define the orientation.
The `to-origin` transform moves something on the molecule onto the origin.
The `adjust` transform is for small adjustments of the molecule on the origin.
The `from-origin` transform moves the molecule to wherever it needs to be in the laboratory
frame to say, position it within a protein receptor."))

(defun make-orientation (&key (from-origin (geom:make-matrix-identity))
                           (adjust-rotation (geom:make-matrix-identity))
                           (adjust-translation (geom:make-matrix-identity))
                           (to-origin (geom:make-matrix-identity)))
  "Make an ORIENTATION object, if any of FROM-ORIGIN, ADJUST, or TO-ORIGIN are not specified then they
default to the identity matrix." 
  (let ((transform (calculate-orientation-transform from-origin adjust-translation adjust-rotation to-origin)))
    (make-instance 'orientation
                   :from-origin from-origin
                   :adjust-rotation adjust-rotation
                   :adjust-translation adjust-translation
                   :to-origin to-origin
                   :transform-cache transform)))

(defun copy-orientation (orientation)
  "Copy an ORIENTATION."
  (make-instance 'orientation
                 :from-origin (geom:copy-matrix (from-origin orientation))
                 :adjust-rotation (geom:copy-matrix (adjust-rotation orientation))
                 :adjust-translation (geom:copy-matrix (adjust-translation orientation))
                 :to-origin (geom:copy-matrix (to-origin orientation))
                 :transform-cache (geom:copy-matrix (transform-cache orientation))
                 ))

(defmethod print-object ((obj orientation) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (format stream ":from-origin~%~s~%:adjust-rot~%~s~%:adjust-trans~%~s~%:to-origin~%~s~%" (from-origin obj) (adjust-rotation obj) (adjust-translation obj) (to-origin obj)))))

(defgeneric kin:orientation-transform (orientation))

(defun calculate-orientation-transform (from-origin adjust-translation adjust-rotation to-origin)
  (let* ((transform0 (geom:m*m adjust-rotation to-origin))
         (transform1 (geom:m*m adjust-translation transform0))
         (transform (geom:m*m from-origin transform1)))
    transform))

(defmethod kin:orientation-transform ((orientation orientation))
  "Combine the components of the ORIENTATION into a 4x4 homogeneous matrix."
  (transform-cache orientation))

(defclass orientations ()
  ((orientations :initform (make-hash-table) :initarg :orientations :accessor orientations))
  (:documentation "A hash-table that maps OLIGOMER-SHAPEs to ORIENTATIONs"))

(defun ensure-complete-orientations (orientations oligomer-shapes)
  (loop for oligomer-shape in oligomer-shapes
        unless (gethash oligomer-shape (orientations orientations))
          do (error "The orientation for oligomer-shape ~s is missing" oligomer-shape)))

(defun make-orientations (pairs)
  (let ((ht (make-hash-table)))
    (loop for cur = pairs then (cddr cur)
          for oligomer-shape = (car cur)
          for maybe-orientation = (cadr cur)
          for orientation = (or maybe-orientation (make-orientation))
          unless (typep orientation 'orientation)
          do (error "~s must be an orientation or nil" orientation)
          when (null cur)
          do (return nil)
          do (setf (gethash oligomer-shape ht) orientation))
    (make-instance 'orientations :orientations ht)))

(defun oligomer-space-orientations (orientations)
  "Convert the orientations hash-table to be keyed on oligomer-space, rather than oligomer or oligomer-shape"
  (let ((ht (make-hash-table)))
    (maphash (lambda (key value)
               (setf (gethash (oligomer-space key) ht) value))
             (orientations orientations))
    (make-instance 'orientations :orientations ht)))

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
   (orientations :initarg :orientations :reader orientations)
   (monomer-subset :initform nil :initarg :monomer-subset :accessor monomer-subset)
   (adjustments :initarg :adjustments :accessor adjustments))
  (:documentation "The assembler class maintains a list of OLIGOMER-SHAPEs and a hash-table of
OLIGOMER-SHAPE to ORIENTATIONs.

The most important functions are UPDATE-INTERNALS and UPDATE-EXTERNALS."
 ))

(defun ligand-oligomer-shape (assembler)
  "Return the ligand oligomer-shape for the ASSEMBLER."
  (first (oligomer-shapes assembler)))

(defun ligand-oligomer-shape-orientation (assembler)
  "Return the ligand oligomer-shape and orientation for the ASSEMBLER."
  (let ((ligand-oligomer-shape (first (oligomer-shapes assembler))))
    (values ligand-oligomer-shape (gethash ligand-oligomer-shape (orientations (orientations assembler))))))

(defun receptor-oligomer-shape (assembler)
  "Return the receptor oligomer-shape for the ASSEMBLER."
  (second (oligomer-shapes assembler)))

(defun receptor-oligomer-shape-orientation (assembler)
  "Return the receptor oligomer-shape and orientation for the ASSEMBLER."
  (let ((receptor-oligomer-shape (second (oligomer-shapes assembler))))
    (values receptor-oligomer-shape (if receptor-oligomer-shape
                                        (gethash receptor-oligomer-shape (orientations (orientations assembler)))
                                        nil))))

(defun complex-oligomer-shapes-orientations (assembler)
  (multiple-value-bind (ligand-oligomer-shape ligand-orientation)
      (ligand-oligomer-shape-orientation assembler)
    (multiple-value-bind (receptor-oligomer-shape receptor-orientation)
        (receptor-oligomer-shape-orientation assembler)
      (values ligand-oligomer-shape ligand-orientation receptor-oligomer-shape receptor-orientation))))

(defgeneric lookup-orientation (assembler-or-orientations oligomer-thing)
  (:documentation "Return the orientation for the OLIGOMER-THING in the ASSEMBLER-OR-ORIENTATION"))

(defmethod lookup-orientation ((assembler assembler) oligomer-thing)
  "Return the orientation for the OLIGOMER-THING (oligomer-space, oligomer, or oligomer-shape) in the ASSEMBLER."
  (or (gethash oligomer-thing (orientations (orientations assembler)))
      (error "Could not find ~s as a key in ~s" oligomer-thing assembler)))


(defmethod lookup-orientation ((orientations orientations) oligomer-thing)
  (or (gethash oligomer-thing (orientations orientations))
      (error "Could not find ~s as a key in ~s" oligomer-thing orientations)))

(defclass subset-assembler (assembler)
  ()
  (:documentation "This class is an assembler that builds a subset of another assembler"))

(defgeneric orientation (oligomer-shape holder))

(defmethod orientation ((orientation orientation) (orientations t))
  orientation)

(defmethod orientation ((identity (eql :identity)) (dummy t))
  (declare (ignore dummy))
  (make-orientation))

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

(defgeneric orientation-for-oligomer-shape (thing oligomer-shape))

(defmethod orientation-for-oligomer-shape ((assembler assembler) oligomer-shape)
  "Return the orientation for the OLIGOMER-SHAPE in the ASSEMBLER"
  (or (gethash oligomer-shape (orientations (orientations assembler)))
      (error "Could not find the orientation for ~s in ~s" oligomer-shape assembler)))

(defmethod orientation-for-oligomer-shape ((orientations orientations) oligomer-shape)
  "Return the orientation for the OLIGOMER-SHAPE in the ORIENTATIONS"
  (or (gethash oligomer-shape (orientations orientations))
      (error "Could not find the orientation for ~s in ~s" oligomer-shape orientations)))

(defun make-coordinates-for-number-of-atoms (number-of-atoms)
    (make-array (* 3 number-of-atoms) :element-type (geom:vecreal-type)
                                      :initial-element (geom:vecreal 0.0)))

(defun make-coordinates-for-assembler (assembler)
  "Make a vector of coordinates that can store the external coordinates of the structure in the assembler."
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

(defun in-monomer-subset (monomer-subset monomer)
  (cond
    ((null monomer-subset)
     ;; If monomer-subset is NULL then everything is in the subset
     t)
    ((typep monomer-subset 'monomer-subset)
     (gethash monomer (monomers monomer-subset)))
    (t (error "Illegal value for monomer-subset ~s - must be NIL or a hash-table"))))

(defun maybe-update-backbone-dihedral-cache (assembler)
  (loop for oligomer-shape in (oligomer-shapes assembler)
        do (topology:with-orientation (gethash oligomer-shape (topology:orientations (topology:orientations assembler)))
             (do-oligomer-shape (monomer-shape monomer monomer-context monomer-shape-kind) oligomer-shape
               (when (typep monomer-shape 'backbone-residue-shape)
                 (let* ((monomer-pos (let ((mp (gethash monomer (monomer-positions assembler))))
                                       (unless mp
                                         (error "Could not find monomer ~s in monomer-positions with keys ~s" monomer (alexandria:hash-table-keys (monomer-positions assembler))))
                                       mp))
                        (atresidue (at-position (ataggregate assembler) monomer-pos))
                        (topology (topology atresidue))
                        (constitution (constitution topology))
                        (residue-properties (residue-properties constitution))
                        (dihedrals (getf residue-properties :dihedrals))
                        (dihedrals-to-cache (loop for dihedral in dihedrals
                                                  when (and (typep dihedral 'dihedral-info-atom)
                                                            (member (name dihedral) '(:phi :psi :phi-1 :psi-1)))
                                                    collect dihedral))
                        (dihedral-cache (make-backbone-dihedral-cache)))
                   (loop for dihedral in dihedrals-to-cache
                         for dihedral-name = (name dihedral)
                         for atom-name = (atom-name dihedral)
                         for joint = (joint-with-name atresidue atom-name)
                         for parent = (kin:parent joint)
                         for grand-parent = (kin:parent parent)
                         for great-grand-parent = (kin:parent grand-parent)
                         for jpos = (kin:xyz-joint/transformed-pos joint)
                         for ppos = (kin:xyz-joint/transformed-pos parent)
                         for gppos = (kin:xyz-joint/transformed-pos grand-parent)
                         for ggppos = (kin:xyz-joint/transformed-pos great-grand-parent)
                         for dih = (geom:calculate-dihedral jpos ppos gppos ggppos)
                         for dih-deg = (bin-dihedral-deg (rad-to-deg dih))
                         do (add-to-cache dihedral-cache dihedral-name dih-deg))
                   (setf (backbone-dihedral-cache-deg monomer-shape) dihedral-cache)
                   ))))))

(defun make-assembler (oligomer-shapes &key (orientations nil orientations-p) monomer-subset energy-function-factory (monomer-contexts nil monomer-contexts-p))
  "Build a assembler for the OLIGOMER-SHAPES.
OLIGOMER-SHAPES - A list of OLIGOMER-SHAPEs that the ASSEMBLER will build.
MONOMER-CONTEXTS - A map of monomers to monomer-contexts copied from another assembler (avoids recalculating them).
ORIENTATIONS - An ORIENTATIONS object that maps OLIGOMER-SHAPES to ORIENTATIONs.
USE-EXCLUDED-ATOMS - A parameter passed to make-energy-function.
ENERGY-FUNCTION-FACTORY - If defined, call this with the aggregate to make the energy-function."
  (cond
    ((not (or (= (length oligomer-shapes) 1) orientations))
     (error "You must provide orientations when there is more than one oligomer-shape"))
    ((and (= (length oligomer-shapes) 1) (null orientations))
     (setf orientations (make-orientations (list (first oligomer-shapes) (make-orientation))))))
  (ensure-complete-orientations orientations oligomer-shapes)
  (unless (every (lambda (os) (typep os 'oligomer-shape)) oligomer-shapes)
    (error "You must provide a list of oligomer-shapes"))
  (let* ((aggregate (chem:make-aggregate :all))
         (monomer-positions (make-hash-table))
         (foldamers (mapcar (lambda (oligomer-shape)
                              (let* ((oligomer (oligomer oligomer-shape))
                                     (oligomer-space (oligomer-space oligomer))
                                     (foldamer (foldamer oligomer-space)))
                                foldamer))
                            oligomer-shapes))
         (monomers-to-residues (make-hash-table))
         (oligomer-shapes-molecules (loop for oligomer-shape in oligomer-shapes
                                          for oligomer = (oligomer oligomer-shape)
                                          for oligomer-space = (oligomer-space oligomer)
                                          for foldamer = (foldamer oligomer-space)
                                          for molecule-index from 0
                                          for molecule = (topology:build-molecule
                                                          oligomer
                                                          :monomer-subset monomer-subset
                                                          :aggregate aggregate
                                                          :molecule-index molecule-index
                                                          :monomers-to-residues monomers-to-residues
                                                          :monomer-positions-accumulator monomer-positions)
                                          do (chem:setf-force-field-name molecule (oligomer-force-field-name foldamer))
                                          collect (cons oligomer-shape molecule))))
    (let* ((ataggregate (let ((atagg (make-instance 'ataggregate :aggregate aggregate)))
                          (resize-atmolecules atagg (length oligomer-shapes))
                          atagg))
           (new-monomer-contexts (if monomer-contexts-p
                                     monomer-contexts
                                     (make-hash-table)))
           (joint-tree (make-joint-tree))
           (adjustments (make-instance 'adjustments))
           (energy-function (if energy-function-factory
                                (funcall energy-function-factory aggregate)
                                (chem:make-energy-function :matter aggregate)))
           (assembler (loop for oligomer-shape-molecule in oligomer-shapes-molecules
                            for oligomer-shape = (car oligomer-shape-molecule)
                            for oligomer = (oligomer oligomer-shape)
                            for molecule = (cdr oligomer-shape-molecule)
                            for oligomer-space = (oligomer-space oligomer)
                            for foldamer = (foldamer oligomer-space)
                            for molecule-index from 0
                            for monomer-to-monomer-shape-map = (monomer-shape-map oligomer-shape)
                            do (cond
                                 (monomer-contexts-p nil)
                                 (monomer-subset
                                  (let ((monomer-contexts (make-hash-table)))
                                    (loop for monomer in (alexandria:hash-table-keys (topology:monomers monomer-subset))
                                          when (topology:monomer-position monomer (topology:oligomer-space oligomer))
                                            do (let ((monomer-context (foldamer-monomer-context monomer oligomer foldamer)))
                                                 (setf (gethash monomer new-monomer-contexts) monomer-context)))))
                                 (t
                                  (loop for monomer across (monomers oligomer-space)
                                        for monomer-context = (foldamer-monomer-context monomer oligomer foldamer)
                                        do (setf (gethash monomer new-monomer-contexts) monomer-context))))
                               ;. This is where I would invoke Conformation_O::buildMoleculeUsingOligomer
                               ;; Use the monomers-to-topologys
                            do (let* ((atmolecule (build-atmolecule-using-oligomer oligomer
                                                                                   monomer-to-monomer-shape-map
                                                                                   monomer-subset
                                                                                   molecule
                                                                                   molecule-index
                                                                                   monomer-positions
                                                                                   joint-tree
                                                                                   (chem:atom-table energy-function)
                                                                                   adjustments)))
                                 (put-atmolecule ataggregate atmolecule molecule-index))
                            finally (return (make-instance (if monomer-subset
                                                               'subset-assembler
                                                               'assembler)
                                                           :monomer-positions monomer-positions
                                                           :monomer-contexts new-monomer-contexts
                                                           :oligomer-shapes oligomer-shapes
                                                           :orientations orientations
                                                           :aggregate aggregate
                                                           :ataggregate ataggregate
                                                           :monomer-subset monomer-subset
                                                           :joint-tree joint-tree
                                                           :energy-function energy-function
                                                           :adjustments adjustments)))))
      ;; The energy-function may be adjusted for the assembler
#|      (when tune-energy-function
        (funcall tune-energy-function energy-function assembler))
      (setf (energy-function assembler) energy-function)|#
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
      (unless monomer-subset
        ;; update the internals so we can build dihedral caches
        (loop for oligomer-shape in (oligomer-shapes assembler)
              do (update-internals assembler oligomer-shape))
        ;; update the dihedral caches
        (maybe-update-backbone-dihedral-cache assembler))
      assembler)))


(defun make-ligand-receptor-assembler (&key ligand receptor orientations monomer-subset energy-function-factory #|tune-energy-function (keep-interaction-factory-factory t)|#)
  "Define an assembler for a complex between a LIGAND and a RECEPTOR."
  ;; The ligand is the first oligomer-shape and the receptor is the second one
  (make-assembler (list ligand receptor)
                  :orientations orientations
                  :monomer-subset monomer-subset
                  :energy-function-factory energy-function-factory))
#|                  :tune-energy-function tune-energy-function
                  :keep-interaction-factory-factory keep-interaction-factory-factory))|#

(defun make-training-assembler (oligomers &key focus-monomer)
  "Build a assembler for the oligomers. This is used for building training molecules."
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
                          do (let* ((atmolecule (build-atmolecule-using-oligomer oligomer
                                                                                 nil ; no monomer-to-monomer-shape-map
                                                                                 nil
                                                                                 molecule
                                                                                 molecule-index
                                                                                 monomer-positions
                                                                                 joint-tree
                                                                                 (chem:atom-table energy-function)
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

(defun find-atom-for-joint (assembler joint)
  "Return the atom corresponding to the JOINT in the aggregate of the ASSEMBLER."
  (destructuring-bind (mol-idx res-idx atm-idx) (kin:id joint)
    (let ((aggregate (aggregate assembler)))
      (let* ((mol (chem:content-at aggregate mol-idx))
             (res (chem:content-at mol res-idx))
             (atm (chem:content-at res atm-idx)))
        atm))))

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


(defun build-atom-tree-external-coordinates* (assembler coords oligomer-shape maybe-orientation)
  (let* ((orientation (orientation maybe-orientation assembler))
         (one-oligomer (oligomer oligomer-shape))
         (joints (gethash one-oligomer (root-map (joint-tree assembler)))))
    (when (null joints)
      (error "Could not find oligomer ~s in root-map ~s" one-oligomer (root-map (joint-tree assembler))))
    (with-orientation orientation
      (loop for joint in joints
            do (kin:update-xyz-coords joint coords)))))


(defun build-atom-tree-for-monomer-shape-external-coordinates* (assembler coords oligomer-shape monomer-shape maybe-orientation)
  (let* ((orientation (orientation maybe-orientation assembler))
         (one-oligomer (oligomer oligomer-shape))
         (joints (gethash one-oligomer (root-map (joint-tree assembler)))))
    (when (null joints)
      (error "Could not find oligomer ~s in root-map ~s" one-oligomer (root-map (joint-tree assembler))))
    (with-orientation orientation
      (loop for joint in joints
            do (kin:update-xyz-coords joint coords)))))


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
        do (build-atom-tree-external-coordinates-and-adjust assembler coords oligomer-shape oligomer-shape)))

(defun build-atom-tree-external-coordinates-and-adjust (assembler coords oligomer-shape maybe-orientation)
  (build-atom-tree-external-coordinates* assembler coords oligomer-shape maybe-orientation)
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

(defgeneric apply-monomer-shape-to-atresidue-internals (assembler oligomer-shape monomer-shape monomer-context atresidue &key verbose)
  (:documentation "Specialize this on different monomer-shape classes.
Fill in internal coordinates into the atresidue for the monomer-shape.
Some specialized methods will need coordinates for the assembler"))


#+(or)
(defun fill-internals-from-oligomer-shape (assembler oligomer-shape
                                           &key (root-monomer (topology:root-monomer (topology:oligomer oligomer-shape)) root-monomer-p)
                                             verbose)
  "Fill internal coordinates from the monomer-shapes in OLIGOMER-SHAPE of the ASSEMBLER.
If ROOT-MONOMER is provided, then start from that."
  (when verbose (let ((*print-pretty* nil)) (format t "fill-internals-from-oligomer-shape ~s~%" oligomer-shape)))
  (let ((coordinates (topology:make-coordinates-for-assembler assembler)))
    (loop for ass-oligomer-shape in (oligomer-shapes assembler)
          when (eq ass-oligomer-shape oligomer-shape)
            do (let* ((oligomer (oligomer ass-oligomer-shape))
                      (ordered-monomers (ordered-monomers oligomer :root-monomer root-monomer)))
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

(defun fill-internals-from-oligomer-shape (assembler oligomer-shape &key verbose)
  "Fill internal coordinates from the monomer-shapes in OLIGOMER-SHAPE of the ASSEMBLER."
  (when verbose (let ((*print-pretty* nil)) (format t "fill-internals-from-oligomer-shape ~s~%" oligomer-shape)))
  (loop for ass-oligomer-shape in (oligomer-shapes assembler)
        when (eq ass-oligomer-shape oligomer-shape)
          do (let* ((oligomer (oligomer ass-oligomer-shape))
                    (monomers (monomers (oligomer-space oligomer))))
               (loop with atagg = (ataggregate assembler)
                     ;;  IGNORING THE FOLLOWING
                     ;; It's really important that we use the ordered-monomers so that the monomer-shapes
                     ;;  will install internal coordinates in the order from the root outwards.
                     ;;  so that any preceeding monomer-shapes are built before any following ones.
                     for monomer across monomers
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
                            (apply-monomer-shape-to-atresidue-internals assembler oligomer-shape monomer-shape monomer-context atres :verbose verbose))))))

(defun update-internals (assembler oligomer-shape &key verbose)
  "Update the internal coordinates of the assembler for the oligomer-shape
or whatever type the second argument is.
The MONOMER-SHAPEs of the OLIGOMER-SHAPE are used to update the internal coordinates in
the atom-tree.

ASSEMBLER - the assembler to update.
OLIGOMER-SHAPE - An oligomer-shape (or permissible-rotamers - I think this is wrong) in the assembler."
  (fill-internals-from-oligomer-shape assembler oligomer-shape :verbose verbose)
  (topology:with-orientation (topology:lookup-orientation assembler oligomer-shape)
    (adjust-internals assembler oligomer-shape)))

(defun update-internals-for-monomer-shape (assembler oligomer-shape monomer-shape &key verbose)
  "Fill internal coordinates from the MONOMER-SHAPE in OLIGOMER-SHAPE of the ASSEMBLER."
  (let* ((atagg (ataggregate assembler))
         (monomer-shape-pos (or (gethash monomer-shape (monomer-shape-to-index oligomer-shape))
                                (error "Could not find monomer-shape ~s in oligomer-shape ~s" monomer-shape oligomer-shape)))
         (monomer-info (aref (monomer-shape-info-vector oligomer-shape) monomer-shape-pos))
         (monomer (monomer monomer-info))
         (monomer-context (gethash monomer (monomer-contexts assembler)))
         (monomer-position (gethash monomer (monomer-positions assembler)))
         (molecule-index (molecule-index monomer-position))
         (residue-index (residue-index monomer-position))
         (atmol (elt (atmolecules atagg) molecule-index))
         (atres (elt (atresidues atmol) residue-index)))
    (when verbose (format t "applying internals for monomer: ~s~%" monomer))
    (apply-monomer-shape-to-atresidue-internals assembler oligomer-shape monomer-shape monomer-context atres :verbose verbose))
  (topology:with-orientation (topology:lookup-orientation assembler oligomer-shape)
    (adjust-internals assembler oligomer-shape)))

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
          (kin:xyz-joint/set-atom joint atm))))

(defun copy-xyz-joint-pos-to-atom-positions (assembler)
  (do-joint-atom (joint atm assembler)
    (when (typep joint 'kin:xyz-joint)
      (let* ((joint-atom (kin:xyz-joint/get-atom joint))
             (pos (chem:atom/get-position joint-atom)))
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
  (let ((backbone-internals-count (make-array (length (oligomer-shapes assembler)) :initial-element 0))
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

(defun assembler-dump-internals (assembler)
  (loop for oligomer-shapes in (oligomer-shapes assembler)
        for index from 0
        do (walk-atoms-joints assembler index
                              (lambda (atom joint atomid)
                                (format t "~s ~s~%" atom joint)))))


(defun assembler-dump-monomer-shapes (assembler)
  (loop for oligomer-shape in (oligomer-shapes assembler)
        for index from 0
        do (format t "oligomer-shape ~s~%" oligomer-shape)
        do (loop for monomer-shape across (topology:monomer-shape-vector oligomer-shape)
                 do (format t "monomer-shape ~s~%" monomer-shape))))

(defun monomer-center (assembler monomer coords)
  "Return the geometric center of the MONOMER in the ASSEMBLER with COORDS."
  (let* ((pos (alexandria:ensure-gethash monomer (monomer-positions assembler)
                                         (error "Could not find ~s in ~s"
                                                monomer
                                                (alexandria:hash-table-keys (monomer-positions assembler)))))
         (atres (at-position (ataggregate assembler) pos))
         (centerx 0.0s0)
         (centery 0.0s0)
         (centerz 0.0s0)
         (num 0))
    (loop for atatm across (joints atres)
          for index3 = (kin:position-index-x3 atatm)
          for x = (aref coords index3)
          for y = (aref coords (+ 1 index3))
          for z = (aref coords (+ 2 index3))
          do (incf centerx x)
          do (incf centery y)
          do (incf centerz z)
          do (incf num))
    (geom:vec (/ centerx num) (/ centery num) (/ centerz num))))

(defun update-externals (assembler &key oligomer-shape
                                     (orientation :identity orientationp)
                                     (coords (topology:make-coordinates-for-assembler assembler)))
  "Update the external coordinates in COORDS using the ASSEMBLER and ORIENTATION.
IF OLIGOMER-SHAPE is provided then just build externals for that OLIGOMER-SHAPE.
If OLIGOMER-SHAPE is not provided then build them all and use the OLIGOMER-SHAPE as the ORIENTATION key.
Return the COORDS."
  (if oligomer-shape
      (progn
        (when (and oligomer-shape (not orientationp))
          (error "You must provide orientation when you provide oligomer-shape"))
        (build-atom-tree-external-coordinates* assembler coords oligomer-shape orientation)
        (adjust-atom-tree-external-coordinates assembler coords oligomer-shape))
      (loop for oligomer-shape in (oligomer-shapes assembler)
            do (update-externals assembler :oligomer-shape oligomer-shape
                                           :orientation oligomer-shape
                                           :coords coords)))
  coords)


(defun update-externals-for-monomer-shape (assembler oligomer-shape monomer-shape coords orientation)
  "Update the external coordinates for one MONOMER-SHAPE in the OLIGOMER-SHAPE in the ASSEMBLER using the ORIENTATION."
  (let* ((atagg (ataggregate assembler))
         (monomer-shape-pos (or (gethash monomer-shape (monomer-shape-to-index oligomer-shape))
                                (error "Could not find monomer-shape ~s in oligomer-shape ~s" monomer-shape oligomer-shape)))
         (monomer-info (aref (monomer-shape-info-vector oligomer-shape) monomer-shape-pos))
         (monomer (monomer monomer-info))
         (monomer-context (gethash monomer (monomer-contexts assembler)))
         (monomer-position (gethash monomer (monomer-positions assembler)))
         (molecule-index (molecule-index monomer-position))
         (residue-index (residue-index monomer-position))
         (atmol (elt (atmolecules atagg) molecule-index))
         (atres (elt (atresidues atmol) residue-index))
         (joints (joints atres))
         (joint0 (elt joints 0)))
    (with-orientation orientation
      (kin:update-xyz-coords joint0 coords))
    (adjust-atom-tree-external-coordinates assembler coords oligomer-shape)))
