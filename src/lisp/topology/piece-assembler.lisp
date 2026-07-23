(in-package :topology)

;;;; piece-assembler.lisp
;;;;
;;;; Build a PIECE-ASSEMBLER from a complex-assembler.  A piece-assembler holds
;;;; one throwaway single-monomer OLIGOMER-SHAPE per (locus monomer-index rotamer-index)
;;;; - every chemical sidechain that could sit at every mobile sidechain locus, on both
;;;; the ligand and the receptor.  Each such shape holds one PIECE-SHAPE that
;;;; carries the allowed rotamers and the anchor points it builds off of.
;;;;
;;;; A locus is (monomer-index-in-its-oligomer-space . owner) where owner is
;;;; :LIGAND or :RECEPTOR.  LOCUS-TO-OLIGOMER-SHAPE is keyed by
;;;; (locus monomer-index rotamer-index) because a ligand locus has many choices.
;;;;
;;;; The ligand and receptor halves run through the same COLLECT-PIECES loop.
;;;; They differ only in (a) which oligomer-shape they walk, (b) the owner
;;;; keyword, and (c) whether a monomer-shape must already be a ROTAMER-SHAPE
;;;; to be considered mobile - which is how the receptor's ~10 moveable
;;;; sidechains are distinguished from its frozen RESIDUE-SHAPEs.  Receptor
;;;; monomers carry exactly one topology name (:monomers (list residue-name)),
;;;; so their inner monomer-index loop simply runs once.


(defclass piece-assembler (assembler-base)
  ((oligomer-shapes :initarg :oligomer-shapes :reader oligomer-shapes)
   (locus-to-oligomer-shape :initarg :locus-to-oligomer-shape :reader locus-to-oligomer-shape)
   )
  )


(defstruct piece-id
  "Identifies one piece of a PIECE-ASSEMBLER.

LOCUS         - (locus-position-in-its-oligomer-space . owner), owner :LIGAND or :RECEPTOR.
MONOMER-NAME  - the topology name of the chemical species at that locus.
MONOMER-INDEX - index into (MONOMERS monomer), i.e. which name was chosen.  This is
                the same quantity ADD-SINGLE-SCAN-ENERGY and MAKE-MRKEY call
                monomer-index, so it needs no translation.
ROTAMER-INDEX - index into the monomer-context's ROTAMER-VECTOR, likewise the same
                quantity the mrkey uses.

Stamped on each piece's throwaway monomer as its ID, so a piece is self-describing
in backtraces and in the LOCUS-TO-OLIGOMER-SHAPE map."
  locus monomer-name monomer-index rotamer-index)

(defun make-one-monomer-oligomer (template-space one-monomer-name &key id)
    "Build a throwaway oligomer-space/oligomer holding exactly ONE monomer.
  Bypasses make-oligomer-space (whose verify-oligomer-space would reject a
  context-less lone sidechain) the same way amber-protein builds receptors."
    (let* ((space   (make-instance 'oligomer-space
                                   :foldamer-name (foldamer-name template-space)
                                   :name :one-monomer))
           (monomer (make-instance 'monomer
                                   :id (or id one-monomer-name)
                                   :monomers (list one-monomer-name))))
      (vector-push-extend monomer (monomers space))
      (values (make-instance 'oligomer
                             :monomer-indexes (make-array 1 :element-type 'ext:byte32
                                                            :initial-element 0)
                             :oligomer-space space)
              monomer
              space)))

(defclass piece-shape (monomer-shape cando.serialize:serializable)
  ((monomer-index :initform 0 :initarg :monomer-index :accessor monomer-index)
   (rotamer-index :initform 0 :initarg :rotamer-index :accessor rotamer-index)
   (rotamer-internals :initarg :rotamer-internals :reader rotamer-internals)
   (anchor-i3xs :initarg :anchor-i3xs :reader anchor-i3xs)
   (piece-constitution-context :initform nil
                               :initarg :piece-constitution-context
                               :reader piece-constitution-context)
   (piece-topology    :initarg :piece-topology  :reader piece-topology)
   (piece-atom-info   :initarg :piece-atom-info :reader piece-atom-info)
   ;; The monomer this piece stands in for, in the COMPLEX-ASSEMBLER's oligomer-space -
   ;; NOT the fresh one-monomer monomer that make-one-monomer-oligomer created.
   ;; ADD-SINGLE-SCAN-ENERGY keys MONOMER-TO-LOCUS by EQ on this object, so passing it
   ;; is what makes the energies object assign the same loci/mrkindexes as the current
   ;; pair-scan machinery.
   (original-monomer  :initarg :original-monomer :reader original-monomer)))


(defmethod make-appropriate-joint ((joint-template bonded-joint-template) (piece-shape piece-shape) (parent-joint (eql nil)) atomid atom-name atom-table)
  (declare (ignore joint-template parent-joint))
  (let ((anchor-i3xs (anchor-i3xs piece-shape)))
    (kin:make-anchored-bonded-joint atomid atom-name atom-table
                                    (parent anchor-i3xs) (grand-parent anchor-i3xs) (great-grand-parent anchor-i3xs))))

(defmethod make-appropriate-joint ((joint-template bonded-joint-template) (nil-or-monomer-shape piece-shape) (parent-joint t) atomid atom-name atom-table)
  (declare (ignore joint-template parent-joint))
  (kin:make-bonded-joint atomid atom-name atom-table))

(defun allowed-rotamers-for-piece (monomer-context shape-key rotamers-db one-monomer-name)
  "Vector of (original-index . rotamer) for MONOMER-CONTEXT at SHAPE-KEY.
  The car is the index INTO the context's rotamer-vector - the mrkeys that
  ADD-PAIR-SCAN-ENERGY builds must use that same index space."
  (let* ((rotamers        (gethash monomer-context (context-to-rotamers rotamers-db)))
         (allowed-indexes (and rotamers (gethash shape-key (shape-key-to-index rotamers))))
         (rotamer-vector  (and rotamers (rotamer-vector rotamers))))
    (cond               
      ((null rotamers)
       (warn "No rotamers for monomer-context ~s" monomer-context) #())
      ((or (null allowed-indexes) (= 0 (length allowed-indexes)))
       (warn "No allowed rotamers for ~s / context ~s / shape-key ~s"
             one-monomer-name monomer-context shape-key) #())
      (t (map 'vector (lambda (i) (cons i (aref rotamer-vector i))) allowed-indexes)))))


(defclass anchor-i3xs ()
  ((parent :initarg :parent :reader parent)
   (grand-parent :initarg :grand-parent :reader grand-parent)
   (great-grand-parent :initarg :great-grand-parent :reader great-grand-parent)))

(defun extract-anchor-i3xs (complex-assembler complex-coords monomer)
  (let* ((pos    (gethash monomer (monomer-positions complex-assembler)))
         (atres  (at-position (ataggregate complex-assembler) pos))
         (joint0 (aref (joints atres) 0))
         (p      (kin:parent joint0))
         (gp     (kin:parent p))
         (ggp    (kin:parent gp))
         (i3xp (kin:joint/position-index-x3 p))
         (i3xgp (kin:joint/position-index-x3 gp))
         (i3xggp (kin:joint/position-index-x3 ggp)))
    (make-instance 'anchor-i3xs
                   :parent i3xp
                   :grand-parent i3xgp
                   :great-grand-parent i3xggp)))

(defun refresh-piece-anchors (piece-assembler complex-coords)
  "Cache the backbone anchor positions on every piece's root AnchoredBondedJoint
  from COMPLEX-COORDS.  The anchor positions are the ONLY thing in a piece-assembler
  that depends on complex-coords - everything else (residues, atom types, joints,
  rotamer sets) is invariant as long as the backbone ROTAMERS are unchanged.  So
  this is the cheap retarget path: a rigid backbone-frame change needs only a call
  to this + a coordinate rebuild, NOT a full rebuild of the piece-assembler.  (A
  backbone ROTAMER change alters shape-keys and thus the piece set, and does need a
  full rebuild.)
  Must run after BUILD-PIECE-ATAGGREGATE (the anchored joints must exist)."
  (loop for roots being the hash-values of (root-map (joint-tree piece-assembler))
        do (dolist (root roots)
             (when (typep root 'kin:anchored-bonded-joint)
               (kin:update-anchor-positions-from-coords root complex-coords))))
  piece-assembler)

(defun build-piece-ataggregate (piece-assembler complex-coords)
    "Build the ataggregate + joint-tree for PIECE-ASSEMBLER and register them.
  Must run AFTER the aggregate and energy-function exist (joints resolve their
  coordinate index through the energy function's atom-table).  Each piece is its
  own one-residue molecule built free (monomer-subset nil), so its root atom gets
  an AnchoredBondedJoint that reads the fixed backbone frame from COMPLEX-COORDS."
    (let* ((all-shapes  (oligomer-shapes piece-assembler))
           (aggregate   (aggregate piece-assembler))
           (ef          (energy-function piece-assembler))
           (atom-table  (chem:atom-table ef))
           (mon-pos     (monomer-positions piece-assembler))
           (ataggregate (make-instance 'ataggregate :aggregate aggregate))
           (joint-tree  (make-joint-tree))
           (adjustments (make-instance 'adjustments)))
      (resize-atmolecules ataggregate (length all-shapes))
      (loop for piece-olig-shape across all-shapes
            for molecule-index from 0
            for one-oligomer = (oligomer piece-olig-shape)
            for molecule     = (chem:content-at aggregate molecule-index)
            for atmolecule   = (build-atmolecule-using-oligomer
                                one-oligomer
                                (monomer-shape-map piece-olig-shape)   ; -> the piece-shape
                                nil                                    ; monomer-subset: build free
                                molecule
                                molecule-index
                                mon-pos
                                joint-tree
                                atom-table
                                adjustments)
            do (put-atmolecule ataggregate atmolecule molecule-index))
      (setf (ataggregate piece-assembler) ataggregate
            (joint-tree  piece-assembler) joint-tree)
      ;; Anchor refresh: each piece root is an AnchoredBondedJoint whose three
      ;; anchor indices point into COMPLEX-COORDS.  Cache the positions now, while
      ;; we still hold that vector; afterwards the pieces are self-sufficient.
      (refresh-piece-anchors piece-assembler complex-coords)
      piece-assembler))



;;; ---------------------------------------------------------------------------
;;; One (locus monomer-index rotamer-index) -> a single-monomer oligomer-shape
;;; ---------------------------------------------------------------------------

(defun make-piece-oligomer-shape (template-space locus one-monomer-name
                                  monomer-context monomer-index rotamer-index rotamer-internals
                                  piece-topology piece-atom-info
                                  piece-constitution-context anchor-i3xs
                                  original-monomer)
  "Build a throwaway one-monomer OLIGOMER-SHAPE wrapping a single PIECE-SHAPE
  for ONE rotamer.  ROTAMER-INDEX is the index into the context's rotamer-vector.
  ORIGINAL-MONOMER is the monomer in the complex-assembler's oligomer-space that
  this piece stands in for - the bridge to the energies object's bookkeeping.
  Returns (values oligomer-shape one-monomer piece-shape)."
  (let ((piece-shape (make-instance 'piece-shape
                                    :anchor-i3xs     anchor-i3xs
                                    :monomer-index     monomer-index
                                    :rotamer-index     rotamer-index
                                    :rotamer-internals rotamer-internals
                                    :piece-constitution-context piece-constitution-context
                                    :piece-topology    piece-topology
                                    :piece-atom-info   piece-atom-info
                                    :original-monomer  original-monomer)))
    (multiple-value-bind (one-oligomer one-monomer)
        (make-one-monomer-oligomer template-space one-monomer-name
                                   :id (make-piece-id :locus locus
                                                      :monomer-name one-monomer-name
                                                      :monomer-index monomer-index
                                                      :rotamer-index rotamer-index))
      (let ((shape-to-index (make-hash-table))
            (shape-map      (make-hash-table))
            (out-mons       (make-hash-table)))
        (setf (gethash piece-shape shape-to-index) 0
              (gethash one-monomer shape-map)      piece-shape
              ;; a lone monomer has no out-couplings
              (gethash one-monomer out-mons)       nil)
        (values (make-instance 'oligomer-shape
                               :oligomer                  one-oligomer
                               :name                      (format nil "~s-~s#~d" locus one-monomer-name rotamer-index)
                               :monomer-shape-vector      (vector piece-shape)
                               :monomer-shape-info-vector (vector (make-instance 'monomer-shape-info
                                                                                 :monomer one-monomer
                                                                                 :monomer-context monomer-context
                                                                                 :monomer-shape-kind :sidechain))
                               :monomer-shape-to-index    shape-to-index
                               :monomer-shape-map         shape-map
                               ;; no in-coupling, so the monomer is its own root
                               :the-root-monomer          one-monomer
                               :in-monomers               (make-hash-table)
                               :out-monomers              out-mons
                               ;; we set the rotamer ourselves and there is no backbone
                               ;; in this shape to be incomplete about
                               :rotamers-state (ensure-valid-rotamers-state
                                                :complete-sidechain-and-backbone-rotamers))
                one-monomer
                piece-shape)))))

;;; ---------------------------------------------------------------------------
;;; The shared collector - drives both the ligand and receptor halves
;;; ---------------------------------------------------------------------------

(defun collect-pieces (complex-assembler complex-coords oligomer-shape owner
                       all-oligomer-shapes locus-to-oligomer-shape monomer-contexts
                       &key require-rotamer-shape verbose)
  "Walk the mobile sidechain monomers of OLIGOMER-SHAPE and build one piece per
(locus, choice), accumulating into ALL-OLIGOMER-SHAPES, LOCUS-TO-OLIGOMER-SHAPE
and MONOMER-CONTEXTS.

OWNER is :LIGAND or :RECEPTOR and becomes the cdr of each locus.
REQUIRE-ROTAMER-SHAPE restricts to monomers whose monomer-shape is already a
ROTAMER-SHAPE - that is how MAKE-ASSEMBLER-FOR-COMPLEX-SCANS marks the receptor
sidechains near the ligand as mobile, leaving the rest frozen RESIDUE-SHAPEs.

Returns the number of pieces collected."
  (let* ((base-oligomer (oligomer oligomer-shape))
         (space         (oligomer-space base-oligomer))
         (foldamer      (find-foldamer (foldamer-name space)))
         (rotamers-db   (rotamers-database oligomer-shape))
         (count         0))
    (loop for locus-pos from 0
          for monomer-shape across (monomer-shape-vector oligomer-shape)
          for info across (monomer-shape-info-vector oligomer-shape)
          for monomer = (monomer info)
          when (and (sidechain-shape-kind-p (monomer-shape-kind info))
                    (or (not require-rotamer-shape)
                        (typep monomer-shape 'rotamer-shape)))
            do (let* ((locus (cons locus-pos owner))
                      ;; the shape-key comes from the flanking BACKBONE rotamers, and the
                      ;; backbone is fixed for the whole scan - so compute it once per locus
                      (original-context (foldamer-monomer-context monomer base-oligomer foldamer))
                      (shape-key (shape-key-for-sidechain-monomer
                                  foldamer original-context oligomer-shape monomer))
                      (anchor-i3xs (extract-anchor-i3xs complex-assembler complex-coords monomer)))
                 (loop for one-monomer-name in (monomers monomer)
                       for monomer-index from 0
                       for mutant-oligomer = (copy-oligomer base-oligomer)
                       do (setf (aref (monomer-indexes mutant-oligomer) locus-pos) monomer-index)
                       do (multiple-value-bind (piece-monomer-context piece-constitution-context)
                              (foldamer-monomer-context monomer mutant-oligomer foldamer)
                            (let* ((context-info (monomer-context-info monomer mutant-oligomer))
                                   (piece-topology (if context-info
                                                       (foldamer-monomer-context-info-topology context-info)
                                                       (monomer-topology monomer mutant-oligomer)))
                                   (piece-atom-info (and context-info
                                                         (foldamer-monomer-context-info-atom-info context-info)))
                                   (allowed (allowed-rotamers-for-piece piece-monomer-context shape-key
                                                                        rotamers-db one-monomer-name)))
                              (when verbose
                                (format t "~s ~s: ~d rotamers ~%" locus one-monomer-name (length allowed)))
                              (loop for entry across allowed
                                    for rotamer-index     = (car entry)
                                    for rotamer-internals = (cdr entry)
                                    do (multiple-value-bind (piece-oligomer-shape one-monomer)
                                           (make-piece-oligomer-shape space locus one-monomer-name
                                                                      piece-monomer-context
                                                                      monomer-index rotamer-index
                                                                      rotamer-internals
                                                                      piece-topology piece-atom-info
                                                                      piece-constitution-context anchor-i3xs
                                                                      monomer)
                                         (vector-push-extend piece-oligomer-shape all-oligomer-shapes)
                                         ;; Key on the PIECE-ID stamped on the monomer, so the map and
                                         ;; the id are one representation.  The table's test must be
                                         ;; EQUALP - EQUAL compares structs by EQ, so a freshly built
                                         ;; piece-id would never match a stored one.
                                         (setf (gethash (id one-monomer) locus-to-oligomer-shape)
                                               piece-oligomer-shape)
                                         (setf (gethash one-monomer monomer-contexts) piece-monomer-context))
                                    do (incf count)))))))
    (when verbose
      (format t "collect-pieces ~s: ~d pieces~%" owner count))
    count))

 (defun build-piece-aggregate (all-oligomer-shapes)
    "One molecule per piece, one residue per molecule. Returns (values aggregate monomer-positions)."
    (let ((aggregate         (chem:make-aggregate :all))
          (monomer-positions (make-hash-table))
          (monomers-to-residues (make-hash-table)))
      (loop for piece-shape-olig across all-oligomer-shapes
            for molecule-index from 0
            for piece-shape = (aref (monomer-shape-vector piece-shape-olig) 0)
            for info        = (aref (monomer-shape-info-vector piece-shape-olig) 0)
            for one-monomer = (monomer info)
            for one-name    = (first (monomers one-monomer))
            for foldamer    = (foldamer (oligomer-space (oligomer piece-shape-olig)))
            do (let ((residue  (build-residue-for-monomer-name (piece-topology piece-shape) one-name
                                                               :atom-info (piece-atom-info piece-shape)))
                     (molecule (chem:make-molecule (intern (name piece-shape-olig) :keyword))))
                 (chem:setf-force-field-name molecule (oligomer-force-field-name foldamer))
                 (chem:add-matter molecule residue)
                 (chem:add-matter aggregate molecule)
                 (setf (gethash one-monomer monomers-to-residues) residue
                       (gethash one-monomer monomer-positions)
                       (make-instance 'monomer-position
                                      :molecule-index molecule-index
                                      :residue-index 0))))
      (values aggregate monomer-positions monomers-to-residues)))


(defun assign-piece-given-atom-types (all-oligomer-shapes monomers-to-residues)
  "Stamp :given-atom-type on every atom of every piece residue.

  Replaces TOPOLOGY:ASSIGN-GIVEN-ATOM-TYPES, which cannot be used here: it
  re-derives the context with FOLDAMER-MONOMER-CONTEXT from the oligomer it is
  handed, and a piece's one-monomer oligomer has no backbone to match against."
  (loop for piece-olig-shape across all-oligomer-shapes
        for piece-oligomer = (oligomer piece-olig-shape)
        for piece-oligomer-space = (oligomer-space piece-oligomer)
        for foldamer    = (foldamer piece-oligomer-space)
        for piece-shape = (aref (monomer-shape-vector piece-olig-shape) 0)
        for one-monomer = (monomer (aref (monomer-shape-info-vector piece-olig-shape) 0))
        for cc          = (piece-constitution-context piece-shape)
        for residue     = (gethash one-monomer monomers-to-residues)
        when (and cc residue)
          do (topology:assign-residue-given-atom-types foldamer residue cc)))

(defun piece-mrkindex (energies original-monomer monomer-index rotamer-index)
  "Deterministic mrkindex for a piece, resolved through the energies object itself."
  (let ((locus (gethash original-monomer (monomer-to-locus energies))))
    (when locus
      (gethash (make-mrkey :locus locus
                           :monomer-index monomer-index
                           :rotamer-index rotamer-index)
               (mrkey-to-mrkindex energies)))))

(defmethod apply-monomer-shape-to-atresidue-internals
    (assembler assembler-internals oligomer-shape (piece-shape piece-shape) monomer-context atresidue &key verbose)
  (declare (ignore oligomer-shape monomer-context verbose))
  ;; piece-shape carries the rotamer object itself - no database lookup by context+index
  (let ((rotamer       (rotamer-internals piece-shape))
        (rotamer-index (rotamer-index piece-shape)))
    (setf (rotamer-index atresidue) rotamer-index)
    (apply-fragment-internals-to-atresidue assembler assembler-internals rotamer rotamer-index atresidue)))

(defmethod adjust-internals ((assembler piece-assembler) internals oligomer-shape)
    (declare (ignore internals oligomer-shape))
  nil)

#+(or)
(defmethod update-externals ((assembler piece-assembler) assembler-internals &key coords)
  (loop for piece-olig-shape across (oligomer-shapes assembler)
        for monomer-shape = (aref (monomer-shape-vector piece-olig-shape) 0)
        for atagg   = (ataggregate assembler)
        for mon     = (monomer (aref (monomer-shape-info-vector piece-olig-shape) 0))
        for mon-pos = (gethash mon (monomer-positions assembler))
        for atres   = (elt (atresidues (elt (atmolecules atagg) (molecule-index mon-pos)))
                           (residue-index mon-pos))
        for joint0  = (elt (joints atres) 0)
        do (update-xyz-coords assembler assembler-internals joint0 coords)))

(defmethod update-externals ((assembler piece-assembler) assembler-internals &key coords)
    ;; Anchored roots carry their own frame; *orientation* stays identity.
    (loop for piece-olig-shape across (oligomer-shapes assembler)
          for mon     = (monomer (aref (monomer-shape-info-vector piece-olig-shape) 0))
          for pos     = (gethash mon (monomer-positions assembler))
          for atres   = (elt (atresidues (elt (atmolecules (ataggregate assembler)) (molecule-index pos)))
                            (residue-index pos))
          for joint0  = (elt (joints atres) 0)
          do (update-xyz-coords assembler assembler-internals joint0 coords)))


#+(or)
(defun check-piece-geometry (piece-assembler piece-id complex-assembler complex-coords)
  "Build coords for one piece and compare its atom positions to the same
  sidechain/rotamer in the complex-assembler.  Returns the max per-atom deviation (Å).
  Near 0 means the AnchoredBondedJoint scheme reproduces the real geometry."
  (let* ((piece-shape-olig (gethash piece-id (locus-to-oligomer-shape piece-assembler)))
         (internals (make-internals-for-assembler piece-assembler))
         (coords    (make-coordinates-for-assembler piece-assembler)))
    (update-internals piece-assembler internals)
    (update-externals piece-assembler internals :coords coords)
    (let* ((ps       (aref (monomer-shape-vector piece-shape-olig) 0))
           (orig-mon (original-monomer ps))
           (mon      (monomer (aref (monomer-shape-info-vector piece-shape-olig) 0)))
           (pos      (gethash mon (monomer-positions piece-assembler)))
           (residue  (chem:content-at (chem:content-at (aggregate piece-assembler)
                                                       (molecule-index pos))
                                      (residue-index pos)))
           ;; the reference: original monomer's residue in the complex-assembler
           (ref-pos  (gethash orig-mon (monomer-positions complex-assembler)))
           (ref-res  (at-position (ataggregate complex-assembler) ref-pos)))
      (copy-all-joint-positions-into-atoms piece-assembler coords)
      (format t "piece ~s: ~d atoms~%" piece-id (chem:number-of-atoms residue))
      ;; NOTE: the reference is the ORIGINAL rotamer state, which may differ from
      ;; this piece's rotamer - so compare only when the rotamer indices match, or
      ;; just eyeball a mol2.  This returns the residue for inspection.
      (values residue ref-res))))

(defun check-piece-geometry (piece-assembler piece-id &key (mol2 "/tmp/piece.mol2"))
    "Build coordinates for ONE piece and write its residue to MOL2 for inspection.
  Runs with *orientation* = identity (the default), so the piece lands in the same
  frame as the complex-coords its anchors were cached from - overlay MOL2 on the
  complex-assembler's sidechain at the same locus to check the anchoring."
    (let* ((piece-olig-shape (or (gethash piece-id (locus-to-oligomer-shape piece-assembler))
                                 (error "No piece for id ~s" piece-id)))
           (piece-shape (aref (monomer-shape-vector piece-olig-shape) 0))
           (internals   (make-internals-for-assembler piece-assembler))
           (coords      (make-coordinates-for-assembler piece-assembler)))
      (update-internals-for-monomer-shape piece-assembler internals piece-olig-shape piece-shape)
      (let* ((mon    (monomer (aref (monomer-shape-info-vector piece-olig-shape) 0)))
             (pos    (gethash mon (monomer-positions piece-assembler)))
             (atres  (elt (atresidues (elt (atmolecules (ataggregate piece-assembler))
                                          (molecule-index pos)))
                         (residue-index pos)))
             (joint0 (elt (joints atres) 0)))
        (update-xyz-coords piece-assembler internals joint0 coords)
        (copy-all-joint-positions-into-atoms piece-assembler coords)
        (let ((residue (chem:content-at (chem:content-at (aggregate piece-assembler)
                                                        (molecule-index pos))
                                       (residue-index pos))))
          (chem:save-mol2 residue mol2 t)
          (format t "wrote ~a (~d atoms), root joint ~a~%"
                  mol2 (chem:number-of-atoms residue) (type-of joint0))
          residue))))




;;; ---------------------------------------------------------------------------
;;; Top level
;;; ---------------------------------------------------------------------------

(defun make-piece-assembler (complex-assembler complex-coords &key verbose)
  "Build a PIECE-ASSEMBLER from COMPLEX-ASSEMBLER and its COMPLEX-COORDS.

Collects every ligand sidechain choice at every ligand sidechain locus, plus
every mobile receptor sidechain (those promoted to ROTAMER-SHAPEs by
MAKE-ASSEMBLER-FOR-COMPLEX-SCANS).  Receptor monomers carry a single topology
name fixed by the PDB, so they contribute one piece each."
  (let ((all-oligomer-shapes    (make-array 16 :adjustable t :fill-pointer 0))
        ;; keyed by PIECE-ID structs - EQUALP descends structs, EQUAL does not
        (locus-to-oligomer-shape (make-hash-table :test 'equalp))
        (monomer-contexts        (make-hash-table)))
    (collect-pieces complex-assembler complex-coords
                    (ligand-oligomer-shape complex-assembler) :ligand
                    all-oligomer-shapes locus-to-oligomer-shape monomer-contexts
                    :require-rotamer-shape nil
                    :verbose verbose)
    (let ((receptor-shape (receptor-oligomer-shape complex-assembler)))
      (when receptor-shape
        (collect-pieces complex-assembler complex-coords
                        receptor-shape :receptor
                        all-oligomer-shapes locus-to-oligomer-shape monomer-contexts
                        :require-rotamer-shape t
                        :verbose verbose)))
    (multiple-value-bind (aggregate monomer-positions monomers-to-residues)
        (build-piece-aggregate all-oligomer-shapes)
      ;; stamp :given-atom-type BEFORE the energy function reads atom types
      (assign-piece-given-atom-types all-oligomer-shapes monomers-to-residues)
      
      (let ((energy-function (chem:make-energy-function :matter aggregate
                                                        :keep-interaction-factory #'rosetta-components-keep ; was nil
                                                        :assign-types t
                                                        :setup (list :rosetta))))
        (let ((pa (make-instance 'piece-assembler
                                 :oligomer-shapes         all-oligomer-shapes
                                 :locus-to-oligomer-shape locus-to-oligomer-shape
                                 :monomer-contexts        monomer-contexts
                                 :aggregate               aggregate
                                 :monomer-positions       monomer-positions
                                 :energy-function         energy-function
                                 )))
          (build-piece-ataggregate pa complex-coords)
          pa
          )))))

(defun count-piece-pairs (piece-assembler &key verbose)
  "Count the cross-locus piece pairs - the number of components a
  one-component-per-pair scheme would need - and assign each piece an mrkindex.
                
  Pieces at the SAME locus are never paired: every choice and rotamer at a locus
  is anchored at the same joint and physically interpenetrates.
        
  Walks OLIGOMER-SHAPES in vector order, which COLLECT-PIECES fills
  deterministically (ligand then receptor; within each, locus -> monomer-index ->
  rotamer-index).  That is the same order the current machinery registers singles
  in, so piece index i here is the mrkindex ADD-SINGLE-SCAN-ENERGY will assign.

  Returns (values n-pairs n-pieces locus-counts)."
  (let ((by-locus (make-hash-table :test 'equal))   ; locus is (locus-pos . owner), a cons
        (n-pieces 0))
    (loop for piece-olig-shape across (oligomer-shapes piece-assembler)
          for one-monomer = (monomer (aref (monomer-shape-info-vector piece-olig-shape) 0))
          for locus       = (piece-id-locus (id one-monomer))
          do (incf (gethash locus by-locus 0))
             (incf n-pieces))
    (let* ((sum-sq  (loop for c being the hash-values of by-locus sum (* c c)))
           (n-pairs (floor (- (* n-pieces n-pieces) sum-sq) 2))
           (counts  (sort (loop for locus being the hash-keys of by-locus
                                  using (hash-value c)
                                collect (cons locus c))
                          (lambda (a b)
                            (let ((la (car a)) (lb (car b)))
                              (if (eq (cdr la) (cdr lb))
                                  (< (car la) (car lb))
                                  (string< (symbol-name (cdr la)) (symbol-name (cdr lb)))))))))
      (when verbose
        (loop for (locus . c) in counts
              do (format t "locus ~s: ~5d pieces~%" locus c))
        (format t "~&~d pieces in ~d loci~%" n-pieces (hash-table-count by-locus))
        (format t "~d cross-locus pairs (all-pairs would be ~d)~%"
                n-pairs (floor (* n-pieces (1- n-pieces)) 2))
        ;; ~200 terms/pair, ~64 bytes/term, three component types
        (format t "one-component-per-pair would resident ~,1f GB of terms~%"
                (/ (* n-pairs 200 64 3) 1d9)))
      (values n-pairs n-pieces counts))))

(defun generate-all-sidechains-mol2 (piece-assembler &key (mol2 "/tmp/all-sidechains.mol2") verbose)
    "Place every piece in PIECE-ASSEMBLER and write the whole aggregate to MOL2.
  Runs with *orientation* = identity (default), so pieces land in the same frame
  their anchors were cached from.  All choices/rotamers at one locus are anchored
  to the same joint, so they overlap in space - expect a fan of sidechains at each
  locus.  Returns the aggregate."
    (let ((internals (make-internals-for-assembler piece-assembler))
          (coords    (make-coordinates-for-assembler piece-assembler))
          (atagg     (ataggregate piece-assembler))
          (mon-pos   (monomer-positions piece-assembler))
          (n 0) (anchored 0))
      (loop for piece-olig-shape across (oligomer-shapes piece-assembler)
            for piece-shape = (aref (monomer-shape-vector piece-olig-shape) 0)
            for mon    = (monomer (aref (monomer-shape-info-vector piece-olig-shape) 0))
            for pos    = (gethash mon mon-pos)
            for atres  = (elt (atresidues (elt (atmolecules atagg) (molecule-index pos)))
                             (residue-index pos))
            for joint0 = (elt (joints atres) 0)
            do (update-internals-for-monomer-shape piece-assembler internals piece-olig-shape piece-shape)
               (update-xyz-coords piece-assembler internals joint0 coords)
               (incf n)      
               (when (typep joint0 'kin:anchored-bonded-joint) (incf anchored))
               (when (and verbose (zerop (mod n 200)))
                 (format t "  placed ~d pieces~%" n)))
      (copy-all-joint-positions-into-atoms piece-assembler coords)
      (chem:save-mol2 (aggregate piece-assembler) mol2 t)
      (format t "wrote ~a: ~d pieces (~d anchored roots), ~d atoms~%"
              mol2 n anchored (chem:number-of-atoms (aggregate piece-assembler)))
      (aggregate piece-assembler)))


(defun generate-sidechains-mol2 (piece-assembler
                                 &key (mol2 "/tmp/sidechains.mol2")
                                   (filter (constantly t))
                                   verbose)
  "Place the pieces of PIECE-ASSEMBLER that satisfy FILTER and write those
  residues to MOL2.  FILTER is called with the PIECE-ID and returns non-nil to
  include the piece.  Runs with *orientation* = identity.  Returns the count placed.
                        
  Only the selected residues are added to the output aggregate, so unlike
  GENERATE-ALL-SIDECHAINS-MOL2 there is no same-locus overlap unless FILTER admits
  more than one piece per locus."            
  (let ((internals (make-internals-for-assembler piece-assembler))
        (coords    (make-coordinates-for-assembler piece-assembler))
        (atagg     (ataggregate piece-assembler))
        (mon-pos   (monomer-positions piece-assembler))
        (out       (chem:make-aggregate :selected-sidechains))
        (n 0) (anchored 0))
    (loop for piece-olig-shape across (oligomer-shapes piece-assembler)
          for piece-shape = (aref (monomer-shape-vector piece-olig-shape) 0)
          for mon    = (monomer (aref (monomer-shape-info-vector piece-olig-shape) 0))
          for id     = (id mon)         ; the piece-id
          when (funcall filter id)
            do (let* ((pos    (gethash mon mon-pos))
                      (atres  (elt (atresidues (elt (atmolecules atagg) (molecule-index pos)))
                                   (residue-index pos)))
                      (joint0 (elt (joints atres) 0)))
                 (update-internals-for-monomer-shape piece-assembler internals piece-olig-shape piece-shape)
                 (update-xyz-coords piece-assembler internals joint0 coords)
                 (incf n)   
                 (when (typep joint0 'kin:anchored-bonded-joint) (incf anchored))
                 (when verbose (format t "  ~s~%" id))))
    ;; push coords into atoms, then copy the selected residues into the output aggregate
    (copy-all-joint-positions-into-atoms piece-assembler coords)
    (loop for piece-olig-shape across (oligomer-shapes piece-assembler)
          for mon = (monomer (aref (monomer-shape-info-vector piece-olig-shape) 0))
          for id  = (id mon)
          when (funcall filter id)
            do (let* ((pos (gethash mon mon-pos))
                      (res (chem:content-at (chem:content-at (aggregate piece-assembler)
                                                             (molecule-index pos))
                                            (residue-index pos)))
                      (mol (chem:make-molecule (chem:get-name res))))
                 (chem:add-matter mol (chem:matter-copy res))
                 (chem:add-matter out mol)))
    (chem:save-mol2 out mol2 t)            
    (format t "wrote ~a: ~d pieces (~d anchored), ~d atoms~%"
            mol2 n anchored (chem:number-of-atoms out))
    n)) 



(defun place-all-pieces (piece-assembler)
  "Position every piece into one coordinate vector.  Returns coords."
  (let ((internals (make-internals-for-assembler piece-assembler))
        (coords    (make-coordinates-for-assembler piece-assembler))
        (atagg     (ataggregate piece-assembler)))
    (loop for os across (oligomer-shapes piece-assembler)
          for ps  = (aref (monomer-shape-vector os) 0)
          for mon = (monomer (aref (monomer-shape-info-vector os) 0))
          for pos = (gethash mon (monomer-positions piece-assembler))
          for atres = (elt (atresidues (elt (atmolecules atagg) (molecule-index pos)))
                           (residue-index pos))
          do (update-internals-for-monomer-shape piece-assembler internals os ps)
             (update-xyz-coords piece-assembler internals (elt (joints atres) 0) coords))
    coords))


(defun piece-atoms (piece-assembler os)
    (let* ((mon (monomer (aref (monomer-shape-info-vector os) 0))) 
           (pos (gethash mon (monomer-positions piece-assembler)))
           (res (chem:content-at (chem:content-at (aggregate piece-assembler) (molecule-index pos))
                                (residue-index pos)))
           (v   (make-array 16 :adjustable t :fill-pointer 0)))
      (chem:do-atoms (a res) (vector-push-extend a v))
      (copy-seq v)))





(defun rosetta-components-keep (aclass)
    "Create exactly the rosetta nonbond/elec/LK components (no bonded, no 14)."
    (or (eq aclass (load-time-value (find-class 'chem:energy-rosetta-nonbond)))
        (eq aclass (load-time-value (find-class 'chem:energy-rosetta-elec)))
        (eq aclass (load-time-value (find-class 'chem:energy-rosetta-lksolvation)))))

(defun test-piece-pair-energy (piece-assembler)
  "Prove the between-matters component path: components exist empty, a pair sets
  matters and evaluates to a finite number, terms appear after evaluation."
  (let* ((ef     (energy-function piece-assembler))
         (nb     (chem:find-component-or-nil ef 'chem:energy-rosetta-nonbond))
         (elec   (chem:find-component-or-nil ef 'chem:energy-rosetta-elec))
         (lk     (chem:find-component-or-nil ef 'chem:energy-rosetta-lksolvation))
         (shapes (oligomer-shapes piece-assembler)))
    ;; (1) the 3 components exist and start empty
    (format t "components: nb=~a elec=~a lk=~a~%" (and nb t) (and elec t) (and lk t))
    (format t "initial terms: nb=~d elec=~d lk=~d~%"
            (and nb (chem:number-of-terms nb)) (and elec (chem:number-of-terms elec))
            (and lk (chem:number-of-terms lk)))
    ;; (2) place all pieces
    (let* ((coords (place-all-pieces piece-assembler))
           (os1    (aref shapes 0))
           (id1    (id (monomer (aref (monomer-shape-info-vector os1) 0))))
           ;; first piece at a DIFFERENT locus
           (os2    (loop for k from 1 below (length shapes)
                         for os = (aref shapes k)
                         for id = (id (monomer (aref (monomer-shape-info-vector os) 0)))
                         unless (equal (piece-id-locus id) (piece-id-locus id1))
                           do (return os)))
           (atoms1 (piece-atoms piece-assembler os1))
           (atoms2 (piece-atoms piece-assembler os2)))
      (format t "pair: ~s  x  ~s~%   (~d vs ~d atoms)~%"
              (id (monomer (aref (monomer-shape-info-vector os1) 0)))
              (id (monomer (aref (monomer-shape-info-vector os2) 0)))
              (length atoms1) (length atoms2))
      ;; (3) set matters on the 3 components
      (dolist (c (list nb elec lk))
        (when c (chem:set-matters c atoms1 atoms2)))
      ;; (4) evaluate - triggers the between-matters pair-list build
      (let ((energy (chem:evaluate-energy ef coords)))
        (format t "energy = ~f~%" energy)
        (format t "terms after eval: nb=~d elec=~d lk=~d~%"
                (and nb (chem:number-of-terms nb)) (and elec (chem:number-of-terms elec))
                (and lk (chem:number-of-terms lk)))
        energy))))


 (defun compare-raw-matrices (e-old e-new &key (tolerance 1d-4))
    "Element-wise diff of two pair matrices.  Valid only if both energies registered
  their singles in the same order (identical mrkindex layout)."
    (let ((ma (pair-scan-energy-lower-triangular-matrix e-old))
          (mb (pair-scan-energy-lower-triangular-matrix e-new)))
      (unless (= (length ma) (length mb))
        (error "matrix sizes differ (~d vs ~d) - single registration differs, cannot compare"
               (length ma) (length mb)))
      (loop with mism = 0 with maxd = 0d0 with worst = nil
            for i below (length ma)
            for a = (aref ma i) for b = (aref mb i)
            when (> (abs (- a b)) tolerance)
              do (incf mism)
                 (when (> (abs (- a b)) maxd) (setf maxd (abs (- a b)) worst i))
                 (when (< mism 20) (format t "cell ~d: old=~,6f new=~,6f  diff=~,2e~%" i a b (abs (- a b))))
            finally (format t "~%~d cells, ~d mismatch (tol ~g), max-diff ~,3e at cell ~a~%"
                            (length ma) mism tolerance maxd worst)
                    (return (values mism maxd)))))

