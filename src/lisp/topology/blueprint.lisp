
(in-package :topology)

(export '() :topology)


;;; ------------------------------------------------------------------
;;; Structural index: (locus, monomer-name, slot) -> flat index
;;;
;;; The flat index IS the mrkindex.  protein-binder.lisp:172 defines mrkindex as "an integer
;;; [0...N) that uniquely maps to a mrkey [where] the range [0...i] represents all of the mrkey
;;; at locus 0" -- already dense and locus-grouped.  Today MRKEY-TO-MRKINDEX is a hash because
;;; the numbering is discovered during fill by VECTOR-PUSH-EXTEND.  Here the layout is fixed at
;;; construction, so the same number is arithmetic:
;;;
;;;     mrkindex = (base locus) + (base monomer) + slot
;;;
;;; and MONOMER-LOCUS-MAX-MRKINDEX becomes (+ (base locus) (extent locus)).
;;;
;;; Nothing below changes when the backbone moves.  That is the point: a geometry change
;;; rewrites SLOT-TO-ROTAMER and invalidates energies, but every index stays valid.

(defclass rotamer-scan ()
  ((mrkindex :initarg :mrkindex :reader mrkindex
             :documentation "This scan's index into every vector keyed by mrkindex - the
interleaved (fa_rest, fa_rep) energies, the mcstate vector, and via LOWER-TRIANGULAR-INDEX%% the
pair table.

Assigned once at construction by BLUEPRINT-MRKINDEX and NEVER changes for the lifetime of the
blueprint.  That immutability is the property the whole design rests on: a backbone change
rewrites conformations and invalidates energies, but every cached energy stays addressable
because no index moves.  READER only - nothing should ever write it.

Derivable from (base locus) + (base monomer) + slot, but stored because a ROTAMER-SCAN has no
back-pointer to its locus or monomer and energy code holding a scan needs this constantly.")
   (residue :initarg :residue :reader residue
            :documentation "The residue in the blueprint's aggregate for THIS slot.  Unique per
slot - unlike SINGLE-SCAN, which owns one residue per (locus, monomer-name) and swaps rotamers
through it, the blueprint materializes every rotamer simultaneously.  That is what lets the
atom table stay fixed while the backbone moves.")
   (atresidue :initarg :atresidue :reader atresidue
              :documentation "The atresidue in the blueprint's ataggregate for THIS slot - the
joints that place this rotamer's atoms.  Unique per slot, for the same reason as RESIDUE.")
   (rotamer-catalogue-index :initform :unfolded :accessor rotamer-catalogue-index
                            :documentation "WHICH rotamer of the monomer-context's ROTAMER-VECTOR
this slot was last folded with - the same number a MONOMER-SHAPE's ROTAMER-INDEX holds, so the two
are directly comparable.

Three index spaces meet at a slot and only this one is shared with the pose:

  SLOT             position in the monomer's ACTIVE subset, which is what MRKINDEX is built from
  shortlist entry  SHAPE-KEY-TO-INDEX[shape-key][slot] - which catalogue entries this backbone allows
  CATALOGUE INDEX  position in ROTAMER-VECTOR.  This slot.  What a monomer-shape records.

Set by UPDATE-BLUEPRINT-SLOT-INTERNALS, at the moment the rotamer is actually applied - not by
FILL-PERSONA, because what matters is what was FOLDED, and a persona can be filled without folding.

  :UNFOLDED            never folded
  :UNUSED              folded with a stand-in - nothing was loaded here
  :GEOMETRY-OPTIMIZED  folded with a rotamer that is not in the catalogue at all

PERSONA-RATE, NOT BLUEPRINT-RATE.  The scan is otherwise permanent, and this is not: a backbone move
changes the shortlist and re-folds, and this must be rewritten with it.  It lives here because this
is where the fold records what it did, and it is a DIAGNOSTIC - nothing should compute from it.")
   (intramolecular-energy-components :initarg :intramolecular-energy-components
                                     :accessor intramolecular-energy-components
                                     :documentation "Energy components evaluating this rotamer's
INTRAmolecular terms - what SINGLE-SCAN's INTRAMOLECULAR-ENERGY-FUNCTION evaluates today.")
   (intermolecular-energy-components :initarg :intermolecular-energy-components
                                     :accessor intermolecular-energy-components
                                     :documentation "Energy components evaluating this rotamer's
INTERmolecular terms.

The intra/inter split is kept explicit rather than flattened into one list because everything
downstream depends on it: energies are stored interleaved as (fa_rest, fa_rep) per mrkindex, the
INTERMOLECULAR-P flag matrix is single-stride over the same indexing, and the bucket
decomposition (backbone/single/pair x intra/inter) can only be computed if both channels exist
separately for every mrkindex."))
  (:documentation "Per-slot materialization record - one for each entry of a BLUEPRINT-MONOMER's
SLOT-TO-ROTAMER.

Components are held here by direct reference, accumulated as they are created.  They are NOT to
be reached through the shared energy function by class: FINDCOMPONENTORNIL returns the FIRST
component of a class (energyFunction.cc:1229-1237), and with one energy function holding every
slot's components that would silently return some other slot's."))

;;; SHARED MRKINDEX INDEXING
;;;
;;; There are exactly two per-slot vectors, both FLAT and both indexed by MRKINDEX over
;;; [0, TOTAL-SLOTS).  Entry i of one describes the same slot as entry i of the other, and
;;; neither may change length after construction:
;;;
;;;   BLUEPRINT's ROTAMER-SCANS   what slot i IS - residue, atresidue, energy components.
;;;                               Structural, invariant for the blueprint's lifetime.
;;;   PERSONA's SLOT-TO-ROTAMER   what is IN slot i right now.  Conformational, rewritten on a
;;;                               backbone change or refinement.
;;;
;;; Flat rather than nested per BLUEPRINT-MONOMER because the pair table iterates mrkindex PAIRS
;;; through LOWER-TRIANGULAR-INDEX%%, so "the scan at mrkindex i" must be one AREF.  Nesting
;;; would force decomposing an mrkindex back into (locus, monomer, slot) at every lookup - the
;;; indirection BLUEPRINT-MRKINDEX exists to remove.
;;;
;;; Splitting them by LIFETIME rather than gathering both into one object means "is this stale?"
;;; has one answer per object instead of per entry.  The usual parallel-vector hazard - two
;;; vectors mutated independently until their indexing drifts - cannot arise: one is never
;;; mutated at all, the other is overwritten in place at fixed length, and both are addressed by
;;; a number neither of them owns.
;;;
;;; Each half reaches the other by mrkindex alone, so neither needs a pointer to the other.  Do
;;; not move a rotamer into a ROTAMER-SCAN to save an AREF: that would put a field with
;;; conformational lifetime inside a structural object and reintroduce exactly the ambiguity
;;; this split exists to remove.

(defclass blueprint-monomer ()
  ((monomer-name :initarg :monomer-name :reader monomer-name
                 :documentation "Structural - it determines which atoms exist.  GLU, ARG and LEU
at one locus are different residues with different atom counts, so this cannot be state.  Which
monomer-name is SELECTED at a locus is already on the pose: MCSTATE holds an mrkindex per locus
and an mrkindex identifies (locus, monomer, slot).")
   (stereoisomers :initform nil :initarg :stereoisomers :reader stereoisomers
                  :documentation "Monomer-names sharing this monomer's constitution - the set
SWITCH-STEREOISOMER may move between.  Same atoms, same bonding, same atom types; they differ in
configuration and possibly in partial charges.

A separate axis from MONOMER-NAME because it varies differently: a sidechain's name says which
residue it is, while a backbone stereoisomer says which mirror of the same residue.  Backbone
needs its own home for this precisely because backbone is NOT in the mrkindex space - putting it
there would have made the pair table 47037 wide instead of a few hundred.

ROTAMER-SLOT-COUNT is maxed over this whole set, not just MONOMER-NAME, so a switch to a richer
stereoisomer fits in slots that already exist and no mrkindex moves.  NIL, or a single-element
list, for a monomer with no alternatives.")
   (topology :initarg :topology :reader topology)
   (matcher :initarg :matcher :reader matcher
            :documentation "The MONOMER-CONTEXT matcher that generates this monomer's context - a
node tree (MONOMER-CONTEXT:CHAIN-NODE / BRANCH-NODE / PLUG-TO-MONOMER-NODE, context.lisp:52-88)
parsed from a training-oligomer-space's pattern.

Structural, and frozen here, because it does not vary with the names at this locus.  A
MONOMER-MATCH-NODE holds NAMES as a SET and tests MEMBER, so one node matches every stereoisomer
of a constitution; the name only affects what ADD-MATCH pushes into PARTS.  The rotamers database
shows this directly - every backbone context is the same three-part shape
(:PRO4RR_-DKP_CGLY, :PRO4SS_-DKP_PROR) and every PME context the same five-part shape, across all
stereoisomer substitutions.  Different winning matchers would mean different pattern trees and so
different shapes.  That is why SWITCH-STEREOISOMER re-runs THIS matcher rather than re-selecting
one.

FOLDAMER-MONOMER-CONTEXT normally rediscovers the matcher on every call by scoring all training
spaces (foldamer.lisp:1196-1204).  Holding it lets a persona refresh skip that search.")
   ;; Deliberately NOT called NMAX -- that name reads as either "highest valid rotamer index"
   ;; or "number of rotamers", and the two differ by one in every expression below.  The count
   ;; is the primitive; the two derived quantities have their own names.
   (rotamer-slot-count :initarg :rotamer-slot-count :reader rotamer-slot-count
                       :documentation "Maximum number of database rotamers+1, maximised over ALL shape-keys for
this (locus, monomer-name).  Sizing on the maximum rather than the current shape-key's count is
what keeps the atom table fixed when the backbone moves.")
   ;; SLOT-TO-ROTAMER, ACTIVE-ROTAMER-COUNT and REFINED-ROTAMER-INDEX-OR-NIL are deliberately
   ;; NOT here - they live on the PERSONA.  They are conformational state, rewritten whenever
   ;; the backbone moves or a refinement runs, while everything in this class is structural and
   ;; fixed for the blueprint's lifetime.  Keeping them out means N personas can share one
   ;; blueprint - possible because a joint stores no coordinates, so each folds with its own
   ;; internals vector - and it means "is this stale?" has one answer per object rather than
   ;; per slot of an object.
   (base :initarg :base :reader base
         :documentation "Flat index of this monomer's slot 0, relative to its locus.")
   (monomer-ordinal :initarg :monomer-ordinal :reader monomer-ordinal
                    :documentation "Dense index of this BLUEPRINT-MONOMER across the whole
blueprint, [0, number-of-blueprint-monomers).  Assigned at construction and never changed.

The persona's ACTIVE-ROTAMER-COUNT and REFINED-ROTAMER-INDEX-OR-NIL are per-monomer, so they need a
compact index to be flat vectors; (+ (base locus) (base monomer)) is unique but sparse, which
would waste a slot-sized vector to hold one number per monomer.")
   )
  (:documentation "One allowed monomer-name at one locus, with its fixed slot range."))

(defmethod print-object ((obj blueprint-monomer) stream)
  ;; BASE is relative to the locus and there is no back-pointer to it, so the ABSOLUTE mrkindex
  ;; range cannot be shown here - BLUEPRINT-LOCUS prints that.  MONOMER-ORDINAL is shown instead
  ;; because it is the index into the persona's per-monomer vectors, i.e. the debugging handle.
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (format stream "~a~@[ (~d stereo)~] ~d slots ord ~d"
                (monomer-name obj)
                (let ((s (stereoisomers obj))) (when (cdr s) (length s)))
                (rotamer-slot-count obj)
                (monomer-ordinal obj)))))

(defclass blueprint-locus ()
  ((locus :initarg :locus :reader locus)
   (kind :initarg :kind :reader kind
         :documentation "One of :FIXED, :BACKBONE or :SIDECHAIN.

The axis that matters is NOT 'does it move' - it is 'is it in the pair table'.  Only :SIDECHAIN
is, and only :SIDECHAIN owns mrkindexes.

  :FIXED      non-moveable receptor backbone and immobile sidechains.  Materialized once as
              RESIDUE-SHAPEs, whose joints are XyzJoints holding positions outright
              (shape.lisp:121-123).  No monomers, no rotamers, EXTENT 0.

  :BACKBONE   the ligand backbone.  It DOES change conformation, rarely, but only ONE backbone
              exists at a time - a backbone rotamer moves everything downstream of it, so twenty
              alternatives at a locus are twenty different chains, not twenty branches off a
              shared trunk.  So it gets BLUEPRINT-MONOMERs, so the persona can track its rotamers
              and stereoisomers, but EXTENT 0: its energy goes to the backbone bucket, and which
              rotamer is installed lives in the persona's BACKBONE-ROTAMER-INDEX, one integer per
              monomer rather than a range of slots.

  :SIDECHAIN  every rotamer materialized simultaneously off a shared backbone.  Owns a
              contiguous mrkindex range and is what the pair table indexes.

Putting backbone in the mrkindex space would have cost 56 slots on the myc system and, worse,
claimed a pair-table entry for every backbone-rotamer/sidechain-rotamer combination as though
they were independent choices evaluable in one structure.  They are not.")
   (original-monomer :initarg :original-monomer :reader original-monomer
                     :documentation "The monomer in the design oligomer-space this locus stands
for.  Keyed by EQ, the way ADD-SINGLE-SCAN-ENERGY keys MONOMER-TO-LOCUS.")
   (oligomer-shape :initarg :oligomer-shape :reader oligomer-shape
                   :documentation "The oligomer-shape this locus came from.  Pass 2 needs it to
reach the monomer-shape-map and the foldamer's rotamers-database.")
   (residue-shape :initarg :residue-shape :initform nil :reader residue-shape
                  :documentation "For a FIXED locus, the RESIDUE-SHAPE that already holds this
monomer's ATOM-COORDINATES (shape.lisp:115-119).  This is the 'conformer' - no new class needed.
NIL for a moveable locus.")
   (residue :initform nil :accessor residue
            :documentation "For a FIXED locus, the ONE residue pass 2 materializes here.  A
moveable locus has no residue of its own - it has one per SLOT, held by the ROTAMER-SCANs - so
this stays NIL there.

Separate from RESIDUE-SHAPE, which is the coordinates rather than the built residue.  Fixed
residues cannot live in the blueprint's ROTAMER-SCANS vector because that is indexed by mrkindex
and a fixed locus owns none.

Written once by pass 2 and never again; the ACCESSOR is for construction, not for search.")
   (atresidue :initform nil :accessor atresidue
              :documentation "For a FIXED locus, the atresidue pass 2 builds for RESIDUE - the
XyzJoints holding this monomer's positions outright.  NIL for a moveable locus, for the same
reason as RESIDUE.")
   (monomers :initarg :monomers :reader monomers
             :documentation "Vector of BLUEPRINT-MONOMER, one per allowed monomer-name.")
   (base :initarg :base :reader base
         :documentation "Flat index (mrkindex) of this locus's first slot.")
   (extent :initarg :extent :reader extent
           :documentation "Number of slots at this locus = sum over monomers of rotamer-slot-count.")
   ;; GENERATION and REFINED-CONFORMATION used to live here.  Both moved to the PERSONA: they
   ;; change when the backbone or a refinement does, which is the persona's rate, and leaving
   ;; them here would have made the blueprint writable by every backbone move.
   )
  (:documentation "One locus: the monomer-names allowed there and their slot ranges."))

(defun moveable-p (blueprint-locus)
  "True when this locus is not fixed in space - :BACKBONE or :SIDECHAIN.
Says nothing about whether it owns mrkindexes; ask OWNS-SLOTS-P for that."
  (not (eq (kind blueprint-locus) :fixed)))

(defun owns-slots-p (blueprint-locus)
  "True when this locus owns a range of mrkindexes, i.e. it is a :SIDECHAIN.

The predicate that guards every index into SLOT-TO-ROTAMER and ROTAMER-SCANS.  A :FIXED or
:BACKBONE locus has EXTENT 0 and its BASE equals the next slot-owning locus's base, so calling
BLUEPRINT-MRKINDEX on one would silently return an index belonging to a different locus."
  (eq (kind blueprint-locus) :sidechain))

(defun locus-fixed-name (blueprint-locus)
  "The monomer-name selected at a FIXED locus, or NIL if it cannot be determined.

Read from the OLIGOMER via CURRENT-STEREOISOMER-NAME (topology-classes.lisp:326-330), which
indexes the monomer's allowed names by the oligomer's MONOMER-INDEXES.  For a fixed locus that
selection is the whole truth - the locus exists precisely because nothing there will change.

Not from RESIDUE-SHAPE's NAME, which looks like the obvious source and is not: it has no
:INITFORM (shape.lisp:116) and nothing populates it.  RESIDUE-SHAPE's own PRINT-OBJECT gives up
and prints the literal \"A-Residue\" with its attempt at a real name commented out
(shape.lisp:133-137).

Not from (MONOMERS (ORIGINAL-MONOMER ...)) either - that is the space's full list of ALLOWED
names, and the fixed tier is exactly where the loci carrying thirty-six of them ended up.

Returns NIL rather than signalling: this feeds PRINT-OBJECT, which must not error on a
half-built object."
  (ignore-errors
   (let ((shape (oligomer-shape blueprint-locus)))
     (when shape
       (current-stereoisomer-name (original-monomer blueprint-locus)
                                  (oligomer shape))))))

(defmethod print-object ((obj blueprint-locus) stream)
  ;; Names first - "locus 2" is a number, "locus 2 METHYL/BNZ" is a thing you recognize.  Then
  ;; the mrkindex range, because everything downstream is keyed by mrkindex and "who owns index
  ;; 137?" is a question the locus number cannot answer.
  ;;
  ;; A fixed locus prints its ONE name - see LOCUS-FIXED-NAME - and no range: its extent is 0,
  ;; and showing 95..94 for ten consecutive loci reads as a bug rather than as an empty extent.
  ;;
  ;; The monomer's ID is whatever the oligomer-space DAG named it - the most recognizable label
  ;; there is when it exists, and :INITFORM NIL (topology-classes.lisp:306) when it does not, so
  ;; ~@[ drops it silently.
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        ;; Keyed on OWNS-SLOTS-P, not MOVEABLE-P: a :BACKBONE locus is moveable but has EXTENT 0,
        ;; so a range would come out as base..base-1.
        (if (owns-slots-p obj)
            (format stream "~d~@[ ~a~] ~{~a~^/~} ~d..~d"
                    (locus obj)
                    (id (original-monomer obj))
                    (map 'list #'monomer-name (monomers obj))
                    (base obj)
                    (+ (base obj) (extent obj) -1))
            (format stream "~d~@[ ~a~] ~a ~a"
                    (locus obj)
                    (id (original-monomer obj))
                    (if (eq (kind obj) :backbone)
                        (format nil "~{~a~^/~}" (map 'list #'monomer-name (monomers obj)))
                        (or (locus-fixed-name obj) ""))
                    (kind obj))))))

(defclass blueprint (assembler)
  ;; ASSEMBLER, not ASSEMBLER-BASE.  OLIGOMER-SHAPES, %LIGAND-OLIGOMER-SHAPE,
  ;; %RECEPTOR-OLIGOMER-SHAPE and ADJUSTMENTS used to be redeclared here, verbatim copies of
  ;; ASSEMBLER's; inheriting them instead is what gives the blueprint UPDATE-EXTERNALS, which is a
  ;; DEFGENERIC with a method on ASSEMBLER.  That method is not a convenience - it is the only
  ;; thing that applies the ligand ORIENTATION TRANSFORM and both ADJUST steps.  Hand-rolling the
  ;; fold without them left the ligand folded in its own local frame, nowhere near the receptor.
  ;;
  ;; The inherited slots RECEPTOR-ONLY and MONOMER-SUBSET default to NIL, which is right:
  ;; IN-MONOMER-SUBSET reads NIL as "everything is in the subset" (assembler.lisp:623-627).
  ;;
  ;; UPDATE-INTERNALS is NOT inherited - it is a plain DEFUN, so there is nothing to specialize,
  ;; and its FILL-INTERNALS-FROM-OLIGOMER-SHAPE dereferences a MONOMER-POSITION that a declined
  ;; sidechain monomer does not have.  UPDATE-BLUEPRINT-INTERNALS stays the blueprint's own path.
  ((loci :initarg :loci :reader loci
         :documentation "Vector of BLUEPRINT-LOCUS indexed by locus number, in mrkindex order -
(base locus) is non-decreasing across the vector.")
   (locus-map :initarg :locus-map :reader locus-map
              :documentation "EQ hash from the design oligomer-space's monomer to its
BLUEPRINT-LOCUS.")
   (total-slots :initarg :total-slots :reader total-slots
                :documentation "Sum of (extent locus) - the number of mrkindexes, and so the
length of any vector indexed by mrkindex.")
   (monomer-count :initarg :monomer-count :reader monomer-count
                  :documentation "Number of BLUEPRINT-MONOMERs across the whole blueprint, and so
the length of any vector indexed by MONOMER-ORDINAL.  The per-monomer counterpart to
TOTAL-SLOTS: a persona needs both, because some of its state is per-slot and some per-monomer.")
   (rotamer-scans :initarg :rotamer-scans :reader rotamer-scans
                  :documentation "ROTAMER-SCAN per slot, indexed by MRKINDEX, length TOTAL-SLOTS.
The structural half of every slot - residue, atresidue, energy components - and the exact mirror
of the persona's SLOT-TO-ROTAMER, which holds the conformational half at the same indices.  See
SHARED MRKINDEX INDEXING above for why both are flat rather than nested per monomer.

Allocated full-length by MAKE-BLUEPRINT with NIL entries, because the LENGTH is a pass 1 fact
even though the CONTENTS are pass 2's: a vector that can never disagree with TOTAL-SLOTS is
worth more than a slot left unbound.  A NIL entry means not yet materialized.

Covers the MOVEABLE tier only.  Fixed loci own no mrkindex, so their single residue and
atresidue hang off the BLUEPRINT-LOCUS instead."))
  (:documentation "An entire design search space, laid out once with a permanent numbering.

Where an assembler describes ONE structure, a blueprint describes every structure reachable
without rebuilding: at each locus, the monomer-names allowed there and every rotamer available to
each.  A search step selects within the blueprint; it does not reshape it.

That is bought with fixed-width sizing.  Each monomer gets NMAX+1 slots - the most rotamers it
can ever need, plus one reserved for a refinement-derived conformation - so slot counts do not
follow the backbone's shape-key.  MRKINDEX therefore reduces to arithmetic,

    (+ (base locus) (base monomer) slot)      ; see BLUEPRINT-MRKINDEX

instead of a hash lookup, and every cached pair energy stays addressable across a backbone move.
A backbone change becomes a REFRESH - rewrite the persona's SLOT-TO-ROTAMER and
ACTIVE-ROTAMER-COUNT - rather than a rebuild.  It costs slots that are materialized but currently
unselectable; at a few hundred total slots that is not a cost worth avoiding.

Loci divide into two tiers.  MOVEABLE-LOCI are the design choices and own all the slots.
FIXED-LOCI - non-moveable receptor backbone and immobile sidechains - are recorded and
materialized, but hold no monomers and no mrkindexes; their coordinates already live in a
RESIDUE-SHAPE and their atoms reach the energy through the backbone bucket.  Skipping them is
what takes the myc system from 47037 slots to a few hundred.

Everything here is IMMUTABLE for the blueprint's lifetime, because the two things that change do
so at their own rates and have their own objects.  A PERSONA holds what changes on a backbone
move, refinement or stereoisomer switch - which rotamer is in each slot, how many are live,
which stereoisomer is instantiated.  A POSE holds what changes every Monte Carlo step - the
coordinates and the selection that produced them.  Several personas can be live against one
blueprint at once, because a joint stores no coordinates and each persona folds with its own
internals vector.

Built in passes.  MAKE-BLUEPRINT runs pass 1 and fills LOCI, LOCUS-MAP and TOTAL-SLOTS only; the
inherited ASSEMBLER-BASE slots stay unbound until pass 2 materializes the residues, aggregate,
ataggregate and energy function, and pass 3 the ROTAMER-SCANs."))

(declaim (inline blueprint-mrkindex))
(defun blueprint-mrkindex (blueprint-locus blueprint-monomer slot)
  "The mrkindex for (locus, monomer-name, slot).

Replaces the MRKEY-TO-MRKINDEX hash lookup with arithmetic, which is only possible because the
slot count is fixed at NMAX+1 rather than varying with the backbone's shape-key."
  (+ (base blueprint-locus) (base blueprint-monomer) slot))

#|
(declaim (inline slot-count refinement-slot))
(defun slot-count (blueprint-monomer)
  "Total slots for this monomer: ROTAMER-COUNT rotamers plus one reserved refinement slot.
This is the length of SLOT-TO-ROTAMER and the amount this monomer contributes to its locus's
EXTENT."
  (1+ (rotamer-count blueprint-monomer)))

(defun refinement-slot (blueprint-monomer)
  "The reserved slot index, always the last one."
  (rotamer-count blueprint-monomer))

(defun slot-selectable-p (blueprint-monomer slot)
  "True when SLOT can be chosen for the CURRENT backbone shape-key.

Selectable slots are [0, VALID-COUNT) plus the refinement slot.  Everything between is
materialized but unaddressed by SLOT-TO-ROTAMER for this shape-key."
  (or (< slot (valid-count blueprint-monomer))
      (= slot (refinement-slot blueprint-monomer))))

|#


;;; ------------------------------------------------------------------
;;; Pass 1 - enumerate and size
;;;
;;; Pure: walks the design assembler and computes the complete layout without materializing
;;; anything.  Everything it produces is structural and fixed for the blueprint's lifetime -
;;; loci, allowed monomer-names, slot counts, bases, ordinals, total-slots.
;;;
;;; Testable on its own.  Run it against the myc system and check the per-locus boundaries
;;; against MONOMER-LOCUS-MAX-MRKINDEX from a current run: if those agree, the addressing is
;;; right before a single residue has been built.
;;;
;;; NOTE on vocabulary.  (MONOMERS oligomer-space) are the loci; (MONOMERS monomer) are the
;;; monomer-NAMES allowed at one locus.  Two different meanings one step apart, so the code
;;; below never uses the raw accessor twice in a row - LOCUS-MONOMER and ALLOWED-NAMES say
;;; which is meant.

(defgeneric max-rotamer-count (rotamers)
  (:documentation "The largest number of rotamers this container can offer under ANY backbone.

Sizing on the maximum rather than on the count for the current backbone is what lets the atom
table stay fixed when the backbone moves - see the SHARED SLOT INDEXING comment above.

Two methods because the two kinds of rotamer container are keyed differently, and the
difference is structural: SIDECHAIN rotamers are backbone-DEPENDENT, so they are keyed by
shape-key and the count varies per key.  BACKBONE rotamers are not keyed by anything - they ARE
the backbone - so the count is just the list length.

Nothing here interprets a shape-key, and nothing should.  MAKE-SHAPE-KEY is (&rest pparts) ->
pparts (internals.lisp:26-27) - a shape-key is whatever list of parts the FOLDAMER decides on.
Amino acids use backbone dihedrals; spiroligomers use combinations of Cremer-Pople ring pucker
parameters.  The count below is a MAPHASH over however many keys exist, which stays correct for
any foldamer.  (Note SIDECHAIN-ROTAMERS' own slot documentation at internals.lisp:453-454 says
\"This is a cons of (phi.psi)\" - that is the amino-acid case, not the general one.)"))

(defmethod max-rotamer-count ((rotamers sidechain-rotamers))
  (let ((most 0))
    (maphash (lambda (shape-key allowed)
               (declare (ignore shape-key))
               (setf most (max most (length allowed))))
             (shape-key-to-index rotamers))
    most))

(defmethod max-rotamer-count ((rotamers backbone-rotamers-base))
  ;; No shape-key indirection - SHAPE-KEY-TO-INDEX is a slot of SIDECHAIN-ROTAMERS only
  ;; (internals.lisp:452-454).  A plain AA backbone carries 36 psi1 x 36 phi = 1296 of these.
  (length (rotamer-vector rotamers)))

(defun context-max-rotamer-count (rotamers-database monomer-context)
  "MAX-ROTAMER-COUNT for MONOMER-CONTEXT, or 0 when the context has no entry in the database -
which gives that monomer a single slot: the reserved refinement slot and nothing else."
  (let ((container (gethash monomer-context (context-to-rotamers rotamers-database))))
    (if (null container)
        0
        (max-rotamer-count container))))

(defun group-names-by-constitution (names)
  "Partition NAMES into groups that share a CONSTITUTION, preserving order.

Two monomer-names sharing a constitution have the same atoms and the same bonding and differ
only in configuration and possibly partial charges - they are stereoisomers, and the blueprint
materializes ONE residue for the group with the pose selecting among them.  Names with distinct
constitutions are different residues (GLU vs ARG vs LEU) and each needs its own materialization
and its own slot range.

Order is preserved because BASE and MONOMER-ORDINAL are assigned by walking the result, and
those must be reproducible for the blueprint's lifetime."
  (let ((groups nil))                       ; alist of (constitution . reversed-names)
    (dolist (name names)
      (let* ((constitution (topology:constitution (chem:find-topology name t)))
             (entry (assoc constitution groups :test #'eq)))
        (if entry
            (push name (cdr entry))
            (setf groups (append groups (list (cons constitution (list name))))))))
    (mapcar (lambda (entry) (reverse (cdr entry))) groups)))

(defun locus-moveable-p (monomer-shape)
  "Is this locus a design choice, or is it fixed in space?

Do not invent a policy here - the oligomer-shape already carries the answer.  A monomer that is
fixed was given a RESIDUE-SHAPE, and MAKE-APPROPRIATE-JOINT specializes on that class to build
XyzJoints instead of BondedJoints (shape.lisp:121-123).  So 'is a residue-shape' and 'is fixed in
space' are the same statement, decided by whoever built the oligomer-shape.  A moveable monomer
got a plain MONOMER-SHAPE driven by a rotamer index."
  (not (typep monomer-shape 'residue-shape)))

;;; ------------------------------------------------------------------
;;; The bound
;;;
;;; ROTAMER-SLOT-COUNT for a (locus, name-group) is
;;;
;;;     1 + max over EVERY monomer-context reachable at that locus with focus in the group
;;;
;;; not the count for the context that happens to apply now.  This is the same medicine the
;;; blueprint already takes for shape-keys: when a quantity varies with something the search
;;; mutates, size on its maximum so no index ever moves.  Over-allocating costs materialized
;;; but unselectable slots, which the design already accepts.
;;;
;;; The point is that it needs no pose.  An earlier attempt computed the context per name with
;;; COPY-OLIGOMER, mirroring piece-assembler.lisp:337-341, and the layout came out dependent on
;;; which name the oligomer happened to have selected - two runs of one system disagreed, 479
;;; slots against 485.  Matching against the SPACE removes the question: MATCHES-IMPL for an
;;; oligomer-space pushes the monomer's WHOLE name list into PARTS (context.lisp:253-260), and
;;; MATCH-ITERATOR is an odometer over exactly those multi-name positions (context.lisp:94-124),
;;; so it enumerates every combination the space permits.  Nothing is selected, so nothing can
;;; be selected wrongly.
;;;
;;; The reach is small.  A spiro sidechain context is :PME_-SIDE_PRO4SR_+DKP_PRO4RR - its own
;;; backbone plus one neighbour across +DKP - so the enumeration is on the order of sixteen
;;; contexts, and computing the bound exactly beats any conservative over-estimate.
;;; ------------------------------------------------------------------

(defgeneric best-monomer-context-matcher (monomer oligomer-or-space foldamer)
  (:documentation "Return (values MATCHER MATCH): the training-oligomer-space matcher that wins for MONOMER, and
its match.  Same selection rule as FOLDAMER-MONOMER-CONTEXT (foldamer.lisp:1196-1204) - highest
NODE-COUNT among those that match - but returns the matcher itself, which that function discards.

NIL NIL when nothing matches.  Matching against a SPACE is stricter than against an oligomer: an
oligomer needs the one current name to be in the pattern, a space needs EVERY allowed name to be
(context.lisp:256-258).  A locus whose design alternatives outrun what the foldamer was trained
on therefore matches nothing here, and gets one slot."))

(defun reachable-monomer-contexts (match)
  "Every (FOCUS-NAME . MONOMER-CONTEXT) pair MATCH can produce, in iteration order.

MATCH-ITERATOR hands back the focus name as its second value (context.lisp:107), so recovering
which name a context belongs to needs no parsing of the context symbol."
  (when match
    (let ((iterator (monomer-context:match-iterator match))
          (pairs nil))
      (loop (multiple-value-bind (context focus-name) (funcall iterator)
              (unless context (return))
              (pushnew (cons focus-name context) pairs :test #'equal)))
      (nreverse pairs))))

(defun group-rotamer-bound (context-pairs group rotamers-database)
  "Largest rotamer count over the contexts in CONTEXT-PAIRS whose focus name is in GROUP.

0 when the group reaches no context at all - which gives it a single slot, the reserved
refinement slot and nothing else."
  (let ((counts (loop for (focus-name . context) in context-pairs
                      when (member focus-name group)
                        collect (context-max-rotamer-count rotamers-database context))))
    (if counts (reduce #'max counts) 0)))

(defun compute-blueprint-layout (assembler &key (moveable-p-fn #'locus-moveable-p))
  "PASS 1.  Return (values LOCI LOCUS-MAP TOTAL-SLOTS MONOMER-COUNT) for the design ASSEMBLER.

LOCI is a vector of BLUEPRINT-LOCUS in mrkindex order, each holding its BLUEPRINT-MONOMERs with
ROTAMER-SLOT-COUNT, BASE and MONOMER-ORDINAL assigned.  Nothing is materialized - no residues,
no atresidues, no energy components, no ROTAMER-SCANs.  Those need the aggregate and the energy
function, which pass 2 builds.

The numbering established here is permanent.  MRKINDEX is (+ (base locus) (base monomer) slot),
and no cached energy stays addressable across a backbone change unless these values never move.

Loci come in two tiers.  A FIXED locus - non-moveable receptor backbone or an immobile sidechain
- is recorded so pass 2 materializes its residue, but it gets no monomers, no slots and no
mrkindexes; its atoms reach the energy through the backbone bucket.  Only MOVEABLE loci are
choices, and only choices belong in the pair table.  MOVEABLE-P-FN decides, and defaults to
LOCUS-MOVEABLE-P."
  (let ((loci (make-array 16 :adjustable t :fill-pointer 0))
        (locus-map (make-hash-table :test 'eq))
        (locus-number 0)
        (monomer-ordinal 0)
        (next-mrkindex 0))
    (loop for oligomer-shape in (oligomer-shapes assembler)
          for oligomer = (oligomer oligomer-shape)
          for space = (oligomer-space oligomer)
          ;; Derived per oligomer-shape, never passed in.  Each oligomer-shape reaches its own
          ;; foldamer's database (shape.lisp:251-254); handing one database to the whole
          ;; assembler silently gave every ligand locus a single slot, because the spiro
          ;; contexts are not in the amber-protein database.
          for shape-rotamers-database = (rotamers-database oligomer-shape)
          for monomer-shape-map = (monomer-shape-map oligomer-shape)
          for foldamer = (foldamer space)
          do (loop for locus-monomer across (monomers space)
                   for allowed-names = (monomers locus-monomer)
                   for monomer-shape = (gethash locus-monomer monomer-shape-map)
                   if (not (funcall moveable-p-fn monomer-shape))
                     ;; FIXED.  BASE is still recorded, and equals the next moveable locus's
                     ;; base since EXTENT is 0 - the mrkindex space simply skips it.
                     do (let ((blueprint-locus
                                (make-instance 'blueprint-locus
                                               :locus locus-number
                                               :kind :fixed
                                               :original-monomer locus-monomer
                                               :oligomer-shape oligomer-shape
                                               :residue-shape monomer-shape
                                               :monomers #()
                                               :base next-mrkindex
                                               :extent 0)))
                          (vector-push-extend blueprint-locus loci)
                          (setf (gethash locus-monomer locus-map) blueprint-locus)
                          (incf locus-number))
                   else
                     do (multiple-value-bind (locus-matcher locus-match)
                            (best-monomer-context-matcher locus-monomer space foldamer)
                       (let* ((locus-base next-mrkindex)
                             (blueprint-monomers (make-array (length allowed-names)
                                                             :fill-pointer 0))
                             (offset-within-locus 0)
                             ;; Every context reachable here, with its focus name.  Enumerated
                             ;; from the SPACE, so this is a function of the design space alone
                             ;; and reproducible across runs - see "The bound" above.
                             (context-pairs (reachable-monomer-contexts locus-match))
                             ;; Group by constitution ONLY for backbone.  Sidechain names never
                             ;; share a constitution, so grouping them is wasted lookups - and
                             ;; if two ever did, merging them would silently delete a design
                             ;; alternative that MCSTATE is supposed to choose between.  Only
                             ;; backbone stereoisomers are meant to collapse into one
                             ;; materialization.  The rotamer container's class is the
                             ;; discriminator and is already what MAX-ROTAMER-COUNT dispatches
                             ;; on (internals.lisp:452-462).
                             ;; Any reachable context answers the backbone-or-sidechain
                             ;; question - a locus does not mix the two kinds - so the first
                             ;; one's container is the discriminator.
                             (container (gethash (cdr (first context-pairs))
                                                 (context-to-rotamers
                                                  shape-rotamers-database)))
                             (groups (if (typep container 'backbone-rotamers-base)
                                         (group-names-by-constitution allowed-names)
                                         (mapcar #'list allowed-names))))
                        (loop for group in groups
                              ;; One BLUEPRINT-MONOMER per group.  A backbone group is the
                              ;; stereoisomers sharing a constitution - same atoms, same bonding
                              ;; - so they get ONE materialization and a STEREOISOMERS set the
                              ;; pose selects from: the four PRO4** at a locus become a single
                              ;; group of four.  A sidechain group is always one name, so each
                              ;; gets its own materialization and its own slot range.
                              ;; The bound: max over every reachable context whose focus is in
                              ;; this group, so a switch to a richer stereoisomer - or a
                              ;; neighbouring backbone change that alters this monomer's context
                              ;; - fits in slots that already exist and no mrkindex moves.
                              for rotamers = (group-rotamer-bound context-pairs group
                                                                  shape-rotamers-database)
                              ;; +1 for the reserved refinement slot - see ROTAMER-SLOT-COUNT.
                              ;; Maxed over the whole group so a SWITCH-STEREOISOMER to a richer
                              ;; member fits in slots that already exist and no mrkindex moves.
                              for slots = (1+ rotamers)
                              do (vector-push-extend
                                  (make-instance 'blueprint-monomer
                                                 ;; Representative of the group.  Which member is
                                                 ;; actually instantiated is the persona's
                                                 ;; ACTIVE-STEREOISOMER, not this.
                                                 :monomer-name (first group)
                                                 :stereoisomers (when (cdr group) group)
                                                 :topology (chem:find-topology (first group) t)
                                                 :matcher locus-matcher
                                                 :rotamer-slot-count slots
                                                 :base offset-within-locus
                                                 :monomer-ordinal monomer-ordinal)
                                  blueprint-monomers)
                                 (incf offset-within-locus slots)
                                 (incf monomer-ordinal))
                        ;; A :BACKBONE locus keeps its BLUEPRINT-MONOMERs - the persona needs
                        ;; them to track rotamers, contexts and stereoisomers - but claims NO
                        ;; mrkindex.  Only one backbone conformation exists at a time, so there
                        ;; is nothing for a range of slots to hold, and its energy belongs to the
                        ;; backbone bucket rather than the pair table.  EXTENT 0 is what stops
                        ;; BLUEPRINT-MRKINDEX from ever being called on it; see OWNS-SLOTS-P.
                        (let* ((backbone-p (typep container 'backbone-rotamers-base))
                               (extent (if backbone-p 0 offset-within-locus))
                               (blueprint-locus
                                 (make-instance 'blueprint-locus
                                                :locus locus-number
                                                :kind (if backbone-p :backbone :sidechain)
                                                :original-monomer locus-monomer
                                                :oligomer-shape oligomer-shape
                                                :monomers (copy-seq blueprint-monomers)
                                                :base locus-base
                                                :extent extent)))
                          (vector-push-extend blueprint-locus loci)
                          (setf (gethash locus-monomer locus-map) blueprint-locus)
                          (incf next-mrkindex extent))
                        (incf locus-number)))))
    ;; MONOMER-ORDINAL was incremented once per blueprint-monomer, so it IS the count.
    (values (copy-seq loci) locus-map next-mrkindex monomer-ordinal)))


;;; ------------------------------------------------------------------
;;; Constructing the blueprint
;;;
;;; What can be carried over from the source assembler and what cannot is settled by one
;;; question: is the slot keyed by the DESIGN SPACE or by the MATERIALIZATION?
;;;
;;;   MONOMER-CONTEXTS is keyed by MONOMER (assembler.lisp:363-365), an oligomer-space object
;;;   the blueprint shares with the assembler because both describe the same design space.  It
;;;   carries over.
;;;
;;;   ADJUSTMENTS is keyed by ATRESIDUE (assembler.lisp:793, 1038, 1074) and MONOMER-POSITIONS
;;;   locates monomers within an aggregate.  Pass 2 builds the blueprint its OWN atresidues and
;;;   aggregate, so the assembler's tables would miss on every lookup here - and miss silently,
;;;   since GETHASH answers NIL instead of signalling.  They are left UNBOUND, which turns that
;;;   silent wrong answer into an UNBOUND-SLOT error the moment pass 2 is skipped.
;;;
;;; AGGREGATE, ENERGY-FUNCTION, ATAGGREGATE and JOINT-TREE are unbound for the same reason.
;;; ------------------------------------------------------------------

(defun make-blueprint (assembler &key (moveable-p-fn #'locus-moveable-p))
  "Build the BLUEPRINT for the design ASSEMBLER.  PASS 1 ONLY.

The returned blueprint has its numbering final and permanent - LOCI, LOCUS-MAP and TOTAL-SLOTS -
and nothing else.  It is not yet usable for energy: no aggregate, no energy-function, no
ataggregate, no ROTAMER-SCANs.  Pass 2 materializes those, and until it runs the inherited
ASSEMBLER-BASE slots are unbound on purpose.

MOVEABLE-P-FN classifies each locus; see LOCUS-MOVEABLE-P for why the default reads the answer
off the oligomer-shape rather than deciding it here."
  (unless (typep assembler 'assembler-base)
    (error "MAKE-BLUEPRINT needs an assembler-base, got ~s" (type-of assembler)))
  (let ((oligomer-shapes (oligomer-shapes assembler)))
    (unless oligomer-shapes
      (error "~s has no oligomer-shapes - there is nothing to lay out." assembler))
    (multiple-value-bind (loci locus-map total-slots monomer-count)
        (compute-blueprint-layout assembler :moveable-p-fn moveable-p-fn)
      (make-instance 'blueprint
                     ;; The LIST is copied so pushing onto the assembler's cannot reach us; the
                     ;; SHAPES themselves are shared deliberately, since LOCUS-MAP keys off
                     ;; monomers living inside them and every fixed locus holds one of their
                     ;; RESIDUE-SHAPEs.
                     :oligomer-shapes (copy-list oligomer-shapes)
                     ;; Read only if present - ASSEMBLER-BASE does not define these, and a
                     ;; receptor-free assembler leaves both NIL anyway.
                     :ligand-oligomer-shape (when (typep assembler 'assembler)
                                              (%ligand-oligomer-shape assembler))
                     :receptor-oligomer-shape (when (typep assembler 'assembler)
                                                (%receptor-oligomer-shape assembler))
                     :monomer-contexts (monomer-contexts assembler)
                     :loci loci
                     :locus-map locus-map
                     :total-slots total-slots
                     :monomer-count monomer-count
                     :rotamer-scans (make-array total-slots :initial-element nil)))))


(defun make-blueprint-and-materialize (assembler &rest args &key (moveable-p-fn #'locus-moveable-p))
  (let ((bp (apply 'make-blueprint assembler args)))
    (materialize-blueprint bp)
    bp))

(defun moveable-loci (blueprint)
  "The loci that are design choices - the ones that own slots and mrkindexes."
  (remove-if-not #'moveable-p (loci blueprint)))

(defun fixed-loci (blueprint)
  "The loci fixed in space, whose RESIDUE-SHAPE already holds their coordinates.  Pass 2
materializes them; they never appear in MCSTATE or the pair table."
  (remove-if #'moveable-p (loci blueprint)))


;;; ------------------------------------------------------------------
;;; Persona - which conformations are loaded
;;;
;;; The third lifetime.  There were always three, and squeezing them into two objects is what
;;; made SLOT-TO-ROTAMER look like it belonged first on the blueprint and then on the pose:
;;;
;;;   BLUEPRINT   never changes            structure and addressing
;;;   PERSONA     backbone move, refinement, stereoisomer switch
;;;   POSE        every Monte Carlo step   coordinates and the selection producing them
;;;
;;; This is not a new layering.  It is the one the kinematics layer already makes, and the
;;; reason a persona can be a VALUE rather than an installed state.  A joint stores no
;;; coordinates: Joint_O holds _PositionIndexX3 and _EndPositionIndexX3 - indices - and
;;; UPDATE-ATAGGREGATE-JOINT-TREE-INTERNAL-COORDINATES (assembler.lisp:975) takes the internals
;;; and coordinates as ARGUMENTS.  So several personas coexist against one joint tree, each
;;; folding with its own internals, and nothing shared is written.  Map:
;;;
;;;   blueprint <-> joint tree      persona <-> internals vector      pose <-> coordinates vector
;;;
;;; (XyzJoint_O is the one joint that can hold a value inline, in _AtomPosition.  It only does
;;; so when _AtomPositionIndex3 is -1; at >= 0 it refreshes from a coords vector like everything
;;; else.  The fixed tier can use either, since those positions never change.)
;;;
;;; Two caveats for multiple live personas: ATRESIDUE's STEREOISOMER-NAME and ROTAMER-INDEX
;;; (fold.lisp:172,176) are per-atresidue and would reflect whichever persona wrote last.  The
;;; second is documented debugging-only; the first is real state that ACTIVE-STEREOISOMER here
;;; supersedes.
;;; ------------------------------------------------------------------

(defclass persona ()
  ((blueprint :initarg :blueprint :reader blueprint
              :documentation "The blueprint this persona describes conformations for.  Every
vector below is sized from it, so a persona is meaningless against any other.")
   (slot-to-rotamer :initarg :slot-to-rotamer :reader slot-to-rotamer
                    :documentation "Rotamers holding the internal coordinates for each slot,
indexed by MRKINDEX - flat, length TOTAL-SLOTS, because an mrkindex already addresses
(locus, monomer, slot) densely and there is no reason to nest what is already flat.

The rotamers live here rather than indexes into the database so placement never resolves a
lookup - the slot IS the conformation.  Entries are GC-managed pointers shared with the
database, not copies; a refinement-derived rotamer is owned here and sits at the reserved slot.

Rewritten on a backbone change or refinement.  Its counterpart in the BLUEPRINT - mrkindexes,
residues, atresidues, energy components - never changes, which is what keeps cached energies
ADDRESSABLE across such a change.  Addressable is not valid: a slot whose rotamer changed has a
stale energy at an index that is still correct.")
   (active-rotamer-count :initarg :active-rotamer-count :reader active-rotamer-count
                         :documentation "Per BLUEPRINT-MONOMER, indexed by MONOMER-ORDINAL: how
many rotamers are active under the current shape-key.  Active slots are [0, count); slots above
are materialized but unaddressed and must be kept out of the energy path.  They still hold
joint trees and residues, so they are ignored rather than cleared.")
   (active-stereoisomer :initarg :active-stereoisomer :reader active-stereoisomer
                        :documentation "Per BLUEPRINT-MONOMER: which of that monomer's
STEREOISOMERS is instantiated, NIL where it has no alternatives.  Backbone stereochemistry needs
a home here because backbone is not in the mrkindex space, so MCSTATE cannot carry it the way it
carries sidechain selection.  Written only by SWITCH-STEREOISOMER.")
   (refined-rotamer-index-or-nil :initarg :refined-rotamer-index-or-nil
                                 :reader refined-rotamer-index-or-nil
                                 :documentation "Per BLUEPRINT-MONOMER: NIL when no
refinement-derived rotamer exists, otherwise the slot holding it.  Moves with
ACTIVE-ROTAMER-COUNT because it is meaningless apart from it - a refinement changes both at the
same instant.")
   (backbone-rotamer-index :initarg :backbone-rotamer-index :reader backbone-rotamer-index
                           :documentation "Per BLUEPRINT-MONOMER of a :BACKBONE locus: the index
of the rotamer currently INSTALLED.  NIL for :SIDECHAIN monomers, which have no single installed
rotamer - all of theirs are materialized at once and SLOT-TO-ROTAMER holds them.

One integer rather than a range of slots because only one backbone conformation exists at a time.
A backbone rotamer moves everything downstream of it, so alternatives cannot coexist in one
structure the way sidechain rotamers can, and there is nothing for a slot range to hold.

This is why backbone owns no mrkindex.  ACTIVE-ROTAMER-COUNT still says how many are AVAILABLE
here - the index is valid in [0, that) - but the choice is made and evaluated by installing it,
not by indexing a pair table.")
   (refined-backbone-mask :initarg :refined-backbone-mask :reader refined-backbone-mask
                          :documentation "Per LOCUS bit vector: 1 means that this persona's
INTERNALS vector contains refinement-derived geometry for the backbone at that locus.  The bit
is separate from BACKBONE-ROTAMER-INDEX because NIL there can also mean that a catalogue index
could not be resolved.  No geometry is duplicated here; this is only the marker that prevents a
later UPDATE-BLUEPRINT-INTERNALS from overwriting refined values with the oligomer-shape.")
   (monomer-contexts :initarg :monomer-contexts :reader monomer-contexts
                     :documentation "Per BLUEPRINT-MONOMER, indexed by MONOMER-ORDINAL: the
monomer-context that applies right now.  With SHAPE-KEYS below it is the exact key pair that
produced SLOT-TO-ROTAMER -

    context-to-rotamers[monomer-context]     -> container
    container.shape-key-to-index[shape-key]  -> rotamer indices

- so a refresh can ask 'do these rotamers still apply?' without recomputing them.

Here rather than on BLUEPRINT-MONOMER because a context is not structural: it names the focus
monomer AND its neighbours, so it changes when a neighbouring backbone name or conformation
changes.  The blueprint holds the MATCHER, which does not change; the persona holds what the
matcher currently produces.")
   (shape-keys :initarg :shape-keys :reader shape-keys
               :documentation "Per LOCUS, the backbone shape-key that determined everything
above.  Stored because it is the INPUT the rest of this persona was derived from: without it,
'do these rotamers still apply?' cannot be answered without recomputing them.

Never interpreted here.  MAKE-SHAPE-KEY is (&rest pparts) -> pparts (internals.lisp:26-27), so a
shape-key is whatever the foldamer says: backbone dihedrals for amino acids, Cremer-Pople ring
pucker parameters for spiroligomers.")
   (internals :initarg :internals :accessor internals :initform nil
              :documentation "The internal-coordinates vector this persona folds with - the
ASSEMBLER-INTERNALS argument of UPDATE-ATAGGREGATE-JOINT-TREE-INTERNAL-COORDINATES.  This is what
makes a persona a value rather than an installed state: it is passed to the fold, not written
into the shared joint tree.  NIL until pass 2 builds the joint tree that sizes it.")
   (generation :initform 0 :accessor generation
               :documentation "Bumped whenever this persona's conformations are rewritten.  The
invalidation handle for anything cached against it - a consumer holding a generation can tell
that its cache is stale without comparing conformations.  HREMC and DEE do not change backbone
conformation and do not bump it."))
  (:documentation "Which conformations are currently loaded into a blueprint's slots.

Everything here answers a question the blueprint deliberately does not: the blueprint says slot
1234 EXISTS and always will, the persona says what is in it right now.

Cached energies stay a separate object - they have their own lifetime, expensive to compute and
cheap to reweight, which is why they are stored rep-weight-independently as interleaved
fa_rest/fa_rep and one cache serves every HREMC replica.  But they are SCOPED to a persona, and
must record which persona and which GENERATION they were computed against.  A backbone change
leaves every mrkindex addressable while making the values at those indices wrong, because the
slot now denotes a different rotamer, and nothing about that failure is loud: the index is still
correct.

Fixed indexing was never about keeping values valid.  It is about not reallocating arrays, not
moving the atom table, not rebuilding the energy function, and letting slots whose rotamer did
NOT change keep their values across a partial refresh."))

(defmethod print-object ((obj persona) stream)
  ;; A persona has no name - what identifies it is HOW LOADED it is and WHICH GENERATION, since
  ;; several can be live against one blueprint at once and the interesting question is always
  ;; "is this one current, and is it filled in?"
  ;;
  ;; SLOT-BOUNDP rather than the accessor: a print-object that errors on a half-built object is
  ;; worst exactly when you most need to look at it.
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (if (slot-boundp obj 'slot-to-rotamer)
            (let ((v (slot-to-rotamer obj)))
              (format stream "gen ~d ~d/~d loaded~@[ ~a~]"
                      (generation obj)
                      (count-if-not #'null v)
                      (length v)
                      (when (internals obj) "+internals")))
            (format stream "gen ~d unallocated" (generation obj))))))

(defun make-persona (blueprint &key (fill t) verbose)
  "An empty persona sized for BLUEPRINT - all vectors allocated, no conformations loaded.

Allocating here rather than on demand means every vector has its final length from the start, so
a later refresh overwrites in place and no index can move.  INTERNALS stays NIL until pass 2
builds the joint tree that determines its size."
  (let ((persona (make-instance 'persona
                                :blueprint blueprint
                                :slot-to-rotamer (make-array (total-slots blueprint) :initial-element nil)
                                :active-rotamer-count (make-array (monomer-count blueprint) :initial-element 0)
                                :active-stereoisomer (make-array (monomer-count blueprint) :initial-element nil)
                                :refined-rotamer-index-or-nil (make-array (monomer-count blueprint)
                                                                          :initial-element nil)
                                :monomer-contexts (make-array (monomer-count blueprint) :initial-element nil)
                                :backbone-rotamer-index (make-array (monomer-count blueprint)
                                                                    :initial-element nil)
                                :refined-backbone-mask (make-array (length (loci blueprint))
                                                                   :element-type 'bit
                                                                   :initial-element 0)
                                :shape-keys (make-array (length (loci blueprint)) :initial-element nil))))
    (if fill
        (multiple-value-bind (loaded overflows)
            (fill-persona persona :verbose verbose)
          (values persona loaded overflows))
        (values persona 0 nil))))

(defun copy-persona (persona)
  "Return an independent conformational-state copy of PERSONA.

The blueprint and the rotamer objects installed in its slots are immutable and remain shared.
Every mutable persona vector is copied, including INTERNALS when it has been built.  The copy has
the same GENERATION because copying does not change a conformation; a subsequent refinement or
backbone refresh is responsible for incrementing it.

Geometry movers use this operation before installing refinement-derived rotamers so that applying
a mover cannot change the input pose through its shared PERSONA."
  (let ((copy (make-instance 'persona
                             :blueprint (blueprint persona)
                             :slot-to-rotamer (copy-seq (slot-to-rotamer persona))
                             :active-rotamer-count (copy-seq (active-rotamer-count persona))
                             :active-stereoisomer (copy-seq (active-stereoisomer persona))
                             :refined-rotamer-index-or-nil
                             (copy-seq (refined-rotamer-index-or-nil persona))
                             :backbone-rotamer-index (copy-seq (backbone-rotamer-index persona))
                             :monomer-contexts (copy-seq (monomer-contexts persona))
                             :shape-keys (copy-seq (shape-keys persona))
                             :refined-backbone-mask (copy-seq (refined-backbone-mask persona))
                             :internals (when (internals persona)
                                          (copy-seq (internals persona))))))
    (setf (generation copy) (generation persona))
    copy))


;;; ------------------------------------------------------------------
;;; Filling a persona
;;;
;;; This needs NO materialization - no aggregate, no joints, no residues.  SLOT-TO-ROTAMER holds
;;; ROTAMER objects straight from the database, so a persona can be filled and checked while the
;;; blueprint's ASSEMBLER-BASE slots are still unbound.
;;;
;;; That is the point of doing it first.  Nothing has yet verified that ROTAMER-SLOT-COUNT is
;;; actually big enough - the bound says GLU gets 55 and PRO4SS 21, and those numbers have never
;;; met a real rotamer set.  Filling a persona checks them for the cost of some hash lookups,
;;; before 526 residues get built on slot counts that might be too small.
;;; ------------------------------------------------------------------

(defun persona-focus-name (blueprint-monomer blueprint-locus)
  "The monomer-name BLUEPRINT-MONOMER currently stands for - which is not always MONOMER-NAME.

MONOMER-NAME is (FIRST GROUP).  For an ungrouped monomer that is the whole story, and it holds
whether or not the oligomer has that name selected: METHYL and BNZ at one locus are separately
materialized and BOTH need their own rotamers loaded, so each uses its own name regardless of
the selection.

For a GROUPED monomer - the four PRO4** sharing a constitution - only one member is instantiated
at a time, so the answer is whichever the oligomer has selected.  Using (FIRST GROUP) there
computed every context as though PRO4SS were live: a wrong context that the bound cannot catch,
because a wrong context is still a REACHABLE one and its rotamer count stays under the maximum.

The MEMBER test is what makes both cases right, including the awkward one - a locus holding two
distinct groups, where the group that does not own the selected name must fall back to its own
representative rather than adopting the other group's."
  (let* ((oligomer (oligomer (oligomer-shape blueprint-locus)))
         (selected (current-stereoisomer-name (original-monomer blueprint-locus) oligomer)))
    (if (member selected (stereoisomers blueprint-monomer))
        selected
        (monomer-name blueprint-monomer))))

(defun persona-monomer-context (blueprint-monomer blueprint-locus)
  "The monomer-context for BLUEPRINT-MONOMER under the oligomer's CURRENT selections, or NIL.

Matched against the OLIGOMER, not the space: the bound needed every context this locus could
reach, this needs the one that applies now.

The match cache is primed with this blueprint-monomer's own name for the FOCUS monomer.  Without
it, CACHED-STEREOISOMER-NAME resolves the focus through CURRENT-STEREOISOMER-NAME and returns
whichever name the oligomer has SELECTED there - so METHYL and BNZ at one locus would both get
the selected name's context, which is exactly the bug the bound had.  Neighbours are deliberately
not primed: their names ARE the current selection, which is what the context should say.

Reaches two MONOMER-CONTEXT internals, since neither the cache struct nor *MATCH-CACHE* is
exported.  An exported priming entry point would be the honest fix.

Returns (values CONTEXT FOCUS-NAME)."
  (let ((oligomer (oligomer (oligomer-shape blueprint-locus)))
        (focus (original-monomer blueprint-locus))
        (focus-name (persona-focus-name blueprint-monomer blueprint-locus)))
    (monomer-context:with-match-cache
        (setf (gethash focus (monomer-context::match-cache-monomer-names
                              monomer-context::*match-cache*))
              focus-name)
      (let ((match (monomer-context:match (matcher blueprint-monomer) focus oligomer)))
        (values (when match
                  (monomer-context:match-as-symbol match))
                focus-name)))))

(defun active-rotamers-for-monomer (container monomer-context oligomer-shape blueprint-locus
                                    foldamer)
  "Return (values ROTAMERS SHAPE-KEY INDEXES): the rotamers active for this monomer right now.

The same split MAX-ROTAMER-COUNT dispatches on.  BACKBONE rotamers are not keyed by anything -
they ARE the backbone - so all of them are active and there is no shape-key.  SIDECHAIN rotamers
are backbone-dependent, so the shape-key selects which apply.

INDEXES is the SLOT -> DATABASE-INDEX table: slot k holds the rotamer at database index
(AREF INDEXES k).  NIL means the two numberings coincide - the backbone branch takes the whole
rotamer vector in order, so slot k IS database index k.

It is returned because those numberings are NOT interchangeable and reading one as the other is
silent: a rotamer index taken from a MONOMER-SHAPE (see BACKBONE-INSTALLED-ROTAMER-INDEX) is a
database index, and using it as a slot number selects a real, valid, WRONG rotamer whenever it
happens to fall below the slot count.  BLUEPRINT-SLOT-FOR-ROTAMER-INDEX inverts this."
  (etypecase container
    (null (values #() nil nil))
    (backbone-rotamers-base (values (rotamer-vector container) nil nil))
    (sidechain-rotamers
     (let* ((shape-key (shape-key-for-sidechain-monomer
                        foldamer monomer-context oligomer-shape
                        (original-monomer blueprint-locus)))
            (indexes (gethash shape-key (shape-key-to-index container)))
            (vector (rotamer-vector container)))
       (values (if indexes
                   (map 'vector (lambda (i) (aref vector i)) indexes)
                   #())
               shape-key
               indexes)))))

(defun blueprint-rotamer-container (blueprint-locus blueprint-monomer)
  "The ROTAMERS container for one (locus, monomer) - the object holding this monomer-context's
ROTAMER-VECTOR and, for a sidechain, its SHAPE-KEY-TO-INDEX map.  NIL when the context has none.

This is the vector every rotamer index in the system is an index INTO, so it is the one place that
makes 'which rotamer does this slot hold' answerable as a lookup rather than an inference."
  (let* ((shape (oligomer-shape blueprint-locus))
         (database (rotamers-database shape))
         (context (persona-monomer-context blueprint-monomer blueprint-locus)))
    (and context (gethash context (context-to-rotamers database)))))

(defun slot-rotamer-database-index (blueprint persona mrkindex)
  "Which rotamer, in its monomer-context's ROTAMER-VECTOR, the slot at MRKINDEX was built from.

Returns an INDEX into that vector, or:

  :UNUSED              nothing is loaded here - the reserved refinement slot, or the slack between
                       this monomer's rotamer count and the locus-wide bound.  FILL-PERSONA fills
                       slots 0..FITS-1 and leaves the rest NIL.
  :GEOMETRY-OPTIMIZED  a rotamer object that is NOT in the database vector, so it came from
                       somewhere else - refinement rather than the library.
  :NO-CONTAINER        the monomer-context resolves to no rotamers container at all.

Compares by IDENTITY: SLOT-TO-ROTAMER holds the very objects ACTIVE-ROTAMERS-FOR-MONOMER selected
out of ROTAMER-VECTOR (it MAPs over the vector, it does not copy), so EQ is meaningful and is the
only test that cannot confuse two rotamers with equal contents."
  (multiple-value-bind (bp-locus bp-monomer slot)
      (blueprint-locus-of-mrkindex blueprint mrkindex)
    (declare (ignore slot))
    (if (null bp-locus)
        :unused
        (let* ((rot (aref (slot-to-rotamer persona) mrkindex))
               (container (blueprint-rotamer-container bp-locus bp-monomer))
               (vec (and container (rotamer-vector container))))
          (cond ((null rot) :unused)
                ((null vec) :no-container)
                ((position rot vec))          ; 0 is a valid answer and is TRUE in CL
                (t :geometry-optimized))))))

(defun active-rotamers-for-locus-monomer (blueprint-locus blueprint-monomer)
  "(values ROTAMERS SHAPE-KEY INDEXES) for one (locus, monomer), resolving the container itself.

Exists so a diagnostic can see INDEXES - the SLOT -> DATABASE-INDEX table that
BLUEPRINT-SLOT-FOR-ROTAMER-INDEX inverts.  NIL there means the two numberings coincide, and a
diagnostic that cannot tell that apart from a successful inversion cannot tell a working conversion
from a no-op."
  (let* ((shape (oligomer-shape blueprint-locus))
         (foldamer (foldamer (oligomer-space (oligomer shape))))
         (database (rotamers-database shape))
         (context (persona-monomer-context blueprint-monomer blueprint-locus))
         (container (and context (gethash context (context-to-rotamers database)))))
    (active-rotamers-for-monomer container context shape blueprint-locus foldamer)))

(defun blueprint-slot-for-rotamer-index (blueprint-locus blueprint-monomer rotamer-index)
  "Which SLOT of BLUEPRINT-MONOMER holds the rotamer at DATABASE index ROTAMER-INDEX.

Returns (values SLOT NIL) on success, or (values NIL REASON) where REASON is

  :NOT-ACTIVE       the rotamer is not in this monomer's active set for the current shape-key -
                    the blueprint never materialized it, so no slot can describe it
  :BEYOND-CAPACITY  it is active but past the slots available, so FILL-PERSONA dropped it (the
                    same overflow FILL-PERSONA reports; the last slot is reserved for refinement)

WHY THIS EXISTS: a rotamer index read off a MONOMER-SHAPE indexes the rotamers DATABASE, while a
slot indexes this monomer's ACTIVE subset.  For a sidechain those differ, and using one as the
other picks a different rotamer without erroring - it only looks wrong once the energies disagree."
  (let* ((shape (oligomer-shape blueprint-locus))
         (foldamer (foldamer (oligomer-space (oligomer shape))))
         (database (rotamers-database shape))
         (context (persona-monomer-context blueprint-monomer blueprint-locus))
         (container (and context (gethash context (context-to-rotamers database)))))
    (multiple-value-bind (rotamers shape-key indexes)
        (active-rotamers-for-monomer container context shape blueprint-locus foldamer)
      (declare (ignore shape-key))
      (let* ((capacity (1- (rotamer-slot-count blueprint-monomer)))
             (fits (min (length rotamers) capacity))
             (slot (if indexes
                       (position rotamer-index indexes)
                       rotamer-index)))
        (cond ((null slot) (values nil :not-active))
              ((>= slot fits) (values nil :beyond-capacity))
              (t (values slot nil)))))))

(defun backbone-installed-rotamer-index (blueprint-locus)
  "The rotamer index currently installed at a :BACKBONE locus, or NIL.

Read from the oligomer-shape's MONOMER-SHAPE, which is where the selection actually lives - the
blueprint does not choose it and the persona only records it.  ROTAMER-INDEX is defined on
ROTAMER-SHAPE and, via CLOSEST-ROTAMER-INDEX, on RESIDUE-SHAPE (shape.lisp:139-140); anything
else, or an unfilled shape, gives NIL rather than an error."
  (ignore-errors
   (let* ((shape (oligomer-shape blueprint-locus))
          (monomer-shape (gethash (original-monomer blueprint-locus)
                                  (monomer-shape-map shape))))
     (when monomer-shape
       (rotamer-index monomer-shape)))))

(defun fill-persona (persona &key verbose)
  "Load PERSONA's slots with the rotamers active under the blueprint's current backbone state.

Fills SLOT-TO-ROTAMER, ACTIVE-ROTAMER-COUNT, MONOMER-CONTEXTS and SHAPE-KEYS, and bumps
GENERATION.  Fixed loci are skipped - they own no mrkindex and have no rotamers.

Returns (values LOADED OVERFLOWS).  LOADED is how many slots were filled.  OVERFLOWS is a list of
(MONOMER-NAME NEEDED CAPACITY) for any monomer whose active rotamers do not FIT its slots, and it
must be empty: an overflow means the bound was computed too small, so rotamers are being silently
dropped and ACTIVE-ROTAMER-COUNT is a lie.  It is returned rather than signalled so one run
reports every offender instead of dying on the first."
  (let* ((blueprint (blueprint persona))
         (slots (slot-to-rotamer persona))
         (contexts (monomer-contexts persona))
         (counts (active-rotamer-count persona))
         (stereo (active-stereoisomer persona))
         (installed (backbone-rotamer-index persona))
         (keys (shape-keys persona))
         (loaded 0)
         (overflows nil))
    (loop for bp-locus across (loci blueprint)
          when (moveable-p bp-locus)
            do (let* ((shape (oligomer-shape bp-locus))
                      (foldamer (foldamer (oligomer-space (oligomer shape))))
                      (database (rotamers-database shape)))
                 (loop for bp-monomer across (monomers bp-locus)
                       for ordinal = (monomer-ordinal bp-monomer)
                       do (multiple-value-bind (context focus-name)
                              (persona-monomer-context bp-monomer bp-locus)
                          (let ((container (and context
                                                (gethash context
                                                         (context-to-rotamers database)))))
                           (multiple-value-bind (rotamers shape-key)
                              (active-rotamers-for-monomer container context shape
                                                           bp-locus foldamer)
                            (setf (aref contexts ordinal) context)
                            ;; Only for a GROUPED monomer, so NIL keeps its documented meaning:
                            ;; this monomer has no alternatives to be active among.
                            (setf (aref stereo ordinal)
                                  (when (cdr (stereoisomers bp-monomer)) focus-name))
                            ;; Per LOCUS, not per monomer: the shape-key comes from the flanking
                            ;; BACKBONE, so every monomer at one locus derives the same one.
                            ;; Backbone monomers have none and leave it NIL.
                            (when shape-key
                              (setf (aref keys (locus bp-locus)) shape-key))
                            ;; The reserved refinement slot is not available to database
                            ;; rotamers - see ROTAMER-SLOT-COUNT - so capacity is one less.
                            (let* ((capacity (1- (rotamer-slot-count bp-monomer)))
                                   (found (length rotamers))
                                   (fits (min found capacity)))
                              (when (> found capacity)
                                (push (list (monomer-name bp-monomer) found capacity)
                                      overflows))
                              ;; Only a :SIDECHAIN locus owns mrkindexes.  BLUEPRINT-MRKINDEX on
                              ;; a :BACKBONE locus would return an index belonging to the next
                              ;; slot-owning locus, since its EXTENT is 0 and its BASE is that
                              ;; locus's BASE.
                              (if (owns-slots-p bp-locus)
                                  (loop for slot below fits
                                        do (setf (aref slots
                                                       (blueprint-mrkindex bp-locus
                                                                           bp-monomer slot))
                                                 (aref rotamers slot))
                                           (incf loaded))
                                  ;; :BACKBONE - record which one is installed, not all of them.
                                  (setf (aref installed ordinal)
                                        (backbone-installed-rotamer-index bp-locus)))
                              (setf (aref counts ordinal) fits))))))))
    (when verbose
      (format t "~&fill-persona: ~d slot~:p loaded of ~d~%" loaded (total-slots blueprint))
      (if overflows
          (progn
            (format t "~&  ~d OVERFLOW~:p - the bound is TOO SMALL:~%" (length overflows))
            (loop for (name found capacity) in overflows
                  do (format t "     ~a needs ~d, has ~d~%" name found capacity)))
          (format t "~&  no overflows - every monomer's rotamers fit its slots~%")))
    (incf (generation persona))
    (values loaded overflows)))


(defun reload-persona-sidechain-rotamers (persona &key verbose)
  "Reload PERSONA's catalogue sidechain slots using its preserved SHAPE-KEYS.

Unlike FILL-PERSONA, this function does not derive shape keys from the blueprint's
OLIGOMER-SHAPEs.  It is the refinement path: optimized backbone internals remain in PERSONA,
while the discrete shape-key inherited from the catalogue backbone selects the sidechain
rotamers that are still compatible with that state.

For each :SIDECHAIN locus, catalogue rotamers are written contiguously starting at slot zero.
Every remaining slot is cleared, including any previously installed refined rotamer, and
REFINED-ROTAMER-INDEX-OR-NIL is reset.  MONOMER-CONTEXTS, SHAPE-KEYS and
ACTIVE-STEREOISOMER are preserved.

Returns (values LOADED OVERFLOWS), with the same meanings as FILL-PERSONA, and increments
GENERATION because the slot conformations have been rewritten."
  (let* ((blueprint (blueprint persona))
         (slots (slot-to-rotamer persona))
         (contexts (monomer-contexts persona))
         (counts (active-rotamer-count persona))
         (refined (refined-rotamer-index-or-nil persona))
         (keys (shape-keys persona))
         (loaded 0)
         (overflows nil))
    (loop for bp-locus across (loci blueprint)
          when (owns-slots-p bp-locus)
            do (let* ((shape (oligomer-shape bp-locus))
                      (database (rotamers-database shape))
                      (shape-key (aref keys (locus bp-locus))))
                 (loop for bp-monomer across (monomers bp-locus)
                       for ordinal = (monomer-ordinal bp-monomer)
                       for context = (aref contexts ordinal)
                       for container = (and context
                                            (gethash context
                                                     (context-to-rotamers database)))
                       for indexes = (and container
                                          (gethash shape-key
                                                   (shape-key-to-index container)))
                       for rotamer-vector = (and container (rotamer-vector container))
                       for found = (length indexes)
                       for capacity = (1- (rotamer-slot-count bp-monomer))
                       for fits = (min found capacity)
                       do (when (> found capacity)
                            (push (list (monomer-name bp-monomer) found capacity)
                                  overflows))
                          ;; Clear the entire fixed-capacity range first.  Besides removing a
                          ;; previous refinement, this prevents stale catalogue rotamers from a
                          ;; larger active set remaining visible above the new active count.
                          (loop for slot below (rotamer-slot-count bp-monomer)
                                do (setf (aref slots
                                               (blueprint-mrkindex bp-locus
                                                                   bp-monomer slot))
                                         nil))
                          (loop for slot below fits
                                for database-index = (aref indexes slot)
                                do (setf (aref slots
                                               (blueprint-mrkindex bp-locus
                                                                   bp-monomer slot))
                                         (aref rotamer-vector database-index))
                                   (incf loaded))
                          (setf (aref counts ordinal) fits
                                (aref refined ordinal) nil))))
    (when verbose
      (format t "~&reload-persona-sidechain-rotamers: ~d catalogue slot~:p loaded~%"
              loaded)
      (if overflows
          (loop initially
                  (format t "~&  ~d OVERFLOW~:p - the bound is TOO SMALL:~%"
                          (length overflows))
                for (name found capacity) in overflows
                do (format t "     ~a needs ~d, has ~d~%" name found capacity))
          (format t "~&  no overflows - every monomer's rotamers fit its slots~%")))
    (incf (generation persona))
    (values loaded (nreverse overflows))))


(defun report-persona (persona &key (stream *standard-output*) zeros-only)
  "Print ACTIVE/CAPACITY per BLUEPRINT-MONOMER; return the list that loaded NOTHING.

A monomer with zero active rotamers is invisible in FILL-PERSONA's totals - the loaded count
simply comes out lower - but it means that monomer can never be selected, so the design has
silently lost an alternative.  That is the opposite failure from an overflow: an overflow drops
rotamers it knows about, a zero never finds any.

CAPACITY excludes the reserved refinement slot, so it is ROTAMER-SLOT-COUNT minus one - the same
figure FILL-PERSONA checks against.  ACTIVE below CAPACITY is normal and expected: the bound
covers every reachable context, the active count is for the current one.

ZEROS-ONLY prints just the offenders."
  (let ((blueprint (blueprint persona))
        (counts (active-rotamer-count persona))
        (contexts (monomer-contexts persona))
        (stereo (active-stereoisomer persona))
        (zeros nil))
    (loop for bp-locus across (loci blueprint)
          when (moveable-p bp-locus)
            do (loop for bp-monomer across (monomers bp-locus)
                     for ordinal = (monomer-ordinal bp-monomer)
                     for active = (aref counts ordinal)
                     for capacity = (1- (rotamer-slot-count bp-monomer))
                     do (when (zerop active)
                          (push bp-monomer zeros))
                        (unless (and zeros-only (plusp active))
                          (format stream
                                  "~&  locus ~3d ~10a ~14a ~4d/~4d~@[ ~a~]~@[ ~a~]~@[ ~a~]~%"
                                  (locus bp-locus)
                                  ;; :BACKBONE shows the installed index instead of a slot range
                                  ;; - it owns no mrkindexes, so there is no range to show.
                                  (if (owns-slots-p bp-locus)
                                      :sidechain
                                      (format nil "bb@~a"
                                              (aref (backbone-rotamer-index persona) ordinal)))
                                  (monomer-name bp-monomer)
                                  active
                                  capacity
                                  ;; Only when the group has more than one member, so its
                                  ;; presence is the signal.  MONOMER-NAME alone cannot answer
                                  ;; this - it is (FIRST GROUP), which reads the same whether
                                  ;; the group collapsed four stereoisomers or held one name.
                                  ;; That distinction decides whether idle slots are reserve
                                  ;; for SWITCH-STEREOISOMER or headroom for neighbour contexts.
                                  ;; The ACTIVE member is starred.  No star on a grouped monomer
                                  ;; means ACTIVE-STEREOISOMER was never written, or was written
                                  ;; with a name outside the group - both silent otherwise,
                                  ;; since neither changes any count.
                                  (let ((names (stereoisomers bp-monomer))
                                        (active (aref stereo ordinal)))
                                    (when (cdr names)
                                      (format nil "stereo[~{~a~^ ~}]"
                                              (mapcar (lambda (name)
                                                        (if (eq name active)
                                                            (format nil "*~a" name)
                                                            name))
                                                      names))))
                                  (when (zerop active) "*** ZERO ***")
                                  ;; A zero WITH a context means the context resolved but the
                                  ;; database had nothing for it at this shape-key.  A zero
                                  ;; WITHOUT one means the match itself failed - a different
                                  ;; bug, and worth telling apart on sight.
                                  (unless (aref contexts ordinal) "<no context>")))))
    (setf zeros (nreverse zeros))
    (format stream "~&  ~d of ~d monomer~:p loaded zero rotamers~%"
            (length zeros) (monomer-count blueprint))
    zeros))


;;; ------------------------------------------------------------------
;;; Testing the bound
;;;
;;; FILL-PERSONA proves the bound is big enough for ONE backbone state.  The bound claims to
;;; cover EVERY state, and that claim is what every mrkindex rests on - if a single reachable
;;; selection overflows, slots are too small, rotamers are dropped, and the numbering is wrong.
;;;
;;; So drive the selection through the combinations the space allows and re-fill at each one.
;;; This mutates the oligomer's MONOMER-INDEXES directly, which is exactly what a persona is
;;; meant to insulate callers from - acceptable in a test that restores what it changed, and not
;;; a pattern to copy into the search.
;;; ------------------------------------------------------------------

(defun locus-allowed-names (blueprint-locus)
  "The monomer-names the SPACE allows at this locus - all of them, not the one selected."
  (monomers (original-monomer blueprint-locus)))

(defun locus-name-index (blueprint-locus)
  "Index into LOCUS-ALLOWED-NAMES of the name currently selected in the oligomer."
  (let* ((oligomer (oligomer (oligomer-shape blueprint-locus)))
         (space (oligomer-space oligomer))
         (position (position (original-monomer blueprint-locus) (monomers space))))
    (aref (monomer-indexes oligomer) position)))

(defun (setf locus-name-index) (index blueprint-locus)
  "Select the INDEXth allowed name at this locus by writing the oligomer's MONOMER-INDEXES -
the same slot CURRENT-STEREOISOMER-NAME reads (topology-classes.lisp:326-330)."
  (let* ((oligomer (oligomer (oligomer-shape blueprint-locus)))
         (space (oligomer-space oligomer))
         (position (position (original-monomer blueprint-locus) (monomers space))))
    (setf (aref (monomer-indexes oligomer) position) index)))

(defun check-bound-over-selections (persona &key (limit 4096) (stream *standard-output*))
  "Re-fill PERSONA under every combination of monomer-names the space allows, and report any
overflow.  Returns the list of failing combinations - empty means the bound survived.

Varies every moveable locus that has more than one allowed name, which is more than just the
backbone stereoisomers: a sidechain alternative changes contexts too, and the bound claims to
cover those as well.

LIMIT caps the number of combinations tried, since the count is a product.  A truncated sweep is
reported as such - a partial pass must not read as a clean one.

Restores the original selection and re-fills before returning, so the persona is left as it was
found."
  (let* ((blueprint (blueprint persona))
         (varying (loop for bp-locus across (loci blueprint)
                        when (and (moveable-p bp-locus)
                                  (cdr (locus-allowed-names bp-locus)))
                          collect bp-locus))
         (sizes (mapcar (lambda (l) (length (locus-allowed-names l))) varying))
         (total (reduce #'* sizes :initial-value 1))
         (saved (mapcar #'locus-name-index varying))
         (budget (min total limit))
         (counters (make-list (length varying) :initial-element 0))
         (failures nil)
         (tested 0))
    (format stream "~&check-bound: ~d varying loci ~s, ~:d combination~:p~@[ (testing ~:d)~]~%"
            (length varying)
            (mapcar #'locus varying)
            total
            (when (> total limit) budget))
    (when (null varying)
      (format stream "~&  nothing varies - the space allows one name at every moveable locus~%")
      (return-from check-bound-over-selections nil))
    (loop while (< tested budget)
          do (loop for bp-locus in varying
                   for index in counters
                   do (setf (locus-name-index bp-locus) index))
             (multiple-value-bind (loaded overflows) (fill-persona persona)
               (declare (ignore loaded))
               (incf tested)
               (when overflows
                 (push (cons (mapcar (lambda (l i) (elt (locus-allowed-names l) i))
                                     varying counters)
                             overflows)
                       failures)))
             ;; Odometer over the allowed-name lists.
             (loop for tail on counters
                   for size in sizes
                   do (incf (car tail))
                      (when (< (car tail) size) (return))
                      (setf (car tail) 0)))
    ;; Put the oligomer back the way it was, and leave the persona consistent with it.
    (loop for bp-locus in varying
          for index in saved
          do (setf (locus-name-index bp-locus) index))
    (fill-persona persona)
    (setf failures (nreverse failures))
    (cond (failures
           (format stream "~&  *** THE BOUND FAILS *** ~d of ~:d combination~:p overflowed~%"
                   (length failures) tested)
           (loop for (names . overflows) in failures
                 repeat 10
                 do (format stream "~&    ~{~a~^ ~}~%" names)
                    (loop for (name found capacity) in overflows
                          do (format stream "        ~a needs ~d, has ~d~%"
                                     name found capacity))))
          ((< tested total)
           (format stream "~&  no overflows in ~:d of ~:d combinations - PARTIAL, not a pass~%"
                   tested total))
          (t
           (format stream "~&  no overflows in all ~:d combinations - the bound holds~%" tested)))
    failures))


;;; ------------------------------------------------------------------
;;; Switching a backbone stereoisomer
;;;
;;; The four stereoisomers at a backbone locus (PRO4SS / PRO4RR / PRO4SR / PRO4RS) share ONE
;;; constitution - same atoms, same bonding - so a switch preserves everything structural: the
;;; residue, the atresidue, the atom table, and every bonded term.  Atom TYPES are identical
;;; too, so the nonbond, LK and bonded parameters are untouched.  That makes a switch a REFRESH
;;; rather than a rebuild, and it is why ROTAMER-SLOT-COUNT must be sized over the whole
;;; switchable set rather than over the stereoisomer the blueprint was handed - a switch to a
;;; richer stereoisomer must fit in slots that already exist, or every mrkindex moves.
;;;
;;; Three things DO change, and none of them errors if skipped - each just yields silently wrong
;;; energies.  That is the entire reason this is one entry point instead of three calls:
;;;
;;;   1. the rotamer set          -> rewrite SLOT-TO-ROTAMER / ACTIVE-ROTAMER-COUNT
;;;   2. the chiral restraints    -> the atoms' CONFIGURATION changes, and
;;;                                  energyFunction.cc:1832-1844 reads exactly that to pick
;;;                                  `side`.  Terms built for the old stereoisomer will fight
;;;                                  the new geometry at 100000 force constant per term.
;;;   3. the partial charges      -> stereoisomers may differ in charge (not in atom type).
;;;                                  EnergyRosettaElec caches them per atom-table index in
;;;                                  _CachedCharge; since the blueprint's atom table never
;;;                                  moves, a stale cache stays the right LENGTH with the right
;;;                                  INDICES and the wrong NUMBERS.

(defun switch-stereoisomer (persona blueprint-monomer new-monomer-name)
  "Switch BLUEPRINT-MONOMER to NEW-MONOMER-NAME, updating everything that depends on it.

Operates on a PERSONA, not a pose: a stereoisomer switch changes which conformations are loaded,
which is the persona's rate.  The blueprint comes from the persona rather than as its own
argument, so a mismatched pair cannot be passed.

Takes the BLUEPRINT-MONOMER rather than a locus: it owns STEREOISOMERS, so the call validates
itself, and it carries MONOMER-ORDINAL, which is exactly the index into the persona's
ACTIVE-STEREOISOMER / ACTIVE-ROTAMER-COUNT / REFINED-ROTAMER-INDEX-OR-NIL vectors.  A locus
would be ambiguous - sidechain loci have several monomers.

A very rare search step - not something HREMC or DEE do.  See the comment above for why all
three updates must happen together."
  (declare (ignorable persona))
  (unless (member new-monomer-name (stereoisomers blueprint-monomer))
    (error "~a is not among ~a's stereoisomers ~a - a stereoisomer switch may only move within~@
            names sharing one constitution."
           new-monomer-name (monomer-name blueprint-monomer)
           (stereoisomers blueprint-monomer)))
  ;; 1. rotamer set
  (error "switch-stereoisomer: not implemented.

Wiring needed before this can work:

  (1) ROTAMER SET.  Needs the blueprint-monomer to hold the switchable stereoisomer SET and the
      pose to hold which one is active - MONOMER-NAME is currently a :READER on
      BLUEPRINT-MONOMER, i.e. structural, which is wrong if it can change.  Deciding that is a
      prerequisite, not a detail: it determines this function's signature.

  (2) CHIRAL RESTRAINTS.  EnergyChiralRestraint exposes ADD-CHIRAL-RESTRAINT-TERM but nothing
      to clear or rebuild one residue's terms.  Either add that, or make `side` a value read at
      evaluation rather than baked in when the term is built.

  (3) CHARGES.  EnergyRosettaElec::invalidateParameterCache (energyRosettaElec.h:166) has no
      CL_DEFMETHOD, so it cannot be called from Lisp.  One line in the header fixes that; then
      this becomes CHEM:ATOM/SET-CHARGE over the residue's atoms followed by the invalidate,
      and the next pairlist build re-reads _CachedCharge.

Deliberately an ERROR rather than a partial implementation: doing 1 of 3 produces wrong
energies with no error, which is the failure mode this entry point exists to prevent."))


 (defun locus-backbone-attachment (blueprint-locus)
    "Return (values BACKBONE-MONOMER COUPLING IN-PLUG) - the backbone monomer this sidechain
  locus hangs off, and how.
      
  *** THIS IS THE LOOKUP TO VERIFY. ***  The assumption is that a sidechain monomer's coupling is
  keyed in its COUPLINGS table by its IN-PLUG's name.  topology-classes.lisp:472 registers a
  coupling under (source-plug-name coupling) on the source monomer and, presumably on the next
  line, under the target plug name on the target monomer - so this holds only if the sidechain is
  the TARGET of the coupling and the target-plug-name is the in-plug's name.  If either is false
  this silently returns the wrong monomer and 470 sidechains bond to the wrong residue.
             
  Every blueprint-monomer at a locus shares one attachment point, so the first one's topology
  answers for all of them."
    (let* ((monomer (original-monomer blueprint-locus))
           (topology (topology (aref (monomers blueprint-locus) 0)))
           (in-plug (find-in-plug topology))
           (coupling (and in-plug (gethash (name in-plug) (couplings monomer)))))
      (unless in-plug
        (error "Locus ~d monomer ~a has no in-plug - it cannot attach to a backbone."
               (locus blueprint-locus) (monomer-name (aref (monomers blueprint-locus) 0))))
      (unless coupling
        (error "Locus ~d has no coupling under its in-plug name ~s.  Couplings present: ~s"
               (locus blueprint-locus) (name in-plug)
               (alexandria:hash-table-keys (couplings monomer))))
      (values (other-monomer coupling monomer) coupling in-plug)))


(defun slot-attachment-atom-names (blueprint-monomer blueprint-locus)
    "Return (values SIDECHAIN-ATOM-NAME BACKBONE-ATOM-NAME) for a slot's attachment bond.
 
  Both come from PLUG-BONDS, the data joint-tree.lisp:100-114 reads to find a plug's joint: an
  in-plug's first plug-bond names the atom on the sidechain side, an out-plug's the backbone side."
    (multiple-value-bind (backbone-monomer coupling in-plug)
        (locus-backbone-attachment blueprint-locus)
      (let* ((sidechain-atom (atom-name (elt (plug-bonds in-plug) 0)))
             (backbone-topology (monomer-topology backbone-monomer
                                                  (oligomer (oligomer-shape blueprint-locus))))
             (out-plug (plug-named backbone-topology (source-plug-name coupling)))
             (backbone-atom (atom-name (elt (plug-bonds out-plug) 0))))
        (declare (ignore blueprint-monomer))
        (values sidechain-atom backbone-atom))))


(defun build-slot-residue (blueprint-monomer blueprint-locus molecule backbone-residue
                           foldamer constitution-context)
  "Build one slot's residue, add it to MOLECULE, and bond it to BACKBONE-RESIDUE.
 
  The bond is marked BONDMASKBLUEPRINT so spanning trees can decline to cross it - see
  NEXTNOTCROSSING, which is the only thing keeping a tree walk out of every rotamer at this locus."
  ;; One in-plug, checked here rather than counting masked bonds afterwards.  A two-point
  ;; attachment - proline's :-side.0 / :-side.1 - would break the count rule that keeps
  ;; cross-slot energy terms out, and it would break it silently, as missing terms.
  (let ((in-plugs 0))
    (maphash (lambda (name plug)
               (declare (ignore name))
               (when (typep plug 'in-plug) (incf in-plugs)))
             (plugs (topology blueprint-monomer)))
    (unless (= in-plugs 1)
      (error "Slot monomer ~a has ~d in-plugs; the fan-out assumes exactly one attachment point."
             (monomer-name blueprint-monomer) in-plugs)))
  (multiple-value-bind (sidechain-atom-name backbone-atom-name)
      (slot-attachment-atom-names blueprint-monomer blueprint-locus)
    (let* ((residue (build-residue-for-monomer-name (topology blueprint-monomer)
                                                    (monomer-name blueprint-monomer)))
           (sidechain-atom (chem:atom-with-name residue sidechain-atom-name))
           (backbone-atom (chem:atom-with-name backbone-residue backbone-atom-name)))
      (chem:add-matter molecule residue)
      ;; BOND-TO returns the bond (atom.cc:564-566), so no second lookup.
      (let ((bond (chem:bond-to sidechain-atom backbone-atom :single-bond :error-if-exceed-valence nil)))
        (chem:bond/setf-mask bond (chem:bond-mask :blueprint)))
      ;; Stamp :given-atom-type = (atom-name . constitution-context).  MAKE-ASSEMBLER does this
      ;; for the base residues (assembler.lisp:743) but cannot reach a slot residue: it walks
      ;; MONOMERS-TO-RESIDUES, and a slot residue has no monomer of its own.
      ;;
      ;; This is what lets the cached SMIRNOFF force field type these atoms by hash lookup
      ;; (smirnoff-cache.lisp:318) instead of running SMIRKS over the fanned-out graph - where a
      ;; CA carrying 58-79 substituents is not a real chemical environment and the types real
      ;; SMIRNOFF would infer are not the ones any actual molecule has.
      (assign-residue-given-atom-types foldamer residue constitution-context)
      residue)))


(defun build-blueprint-fan-out (blueprint oligomer-shapes-molecules monomers-to-residues)
  "Build every sidechain slot's residue and bond it to its backbone.  Returns how many.

  Called from MATERIALIZE-BLUEPRINT-BASE's AFTER-MOLECULES-FN - after the molecules exist and
  before MAKE-ENERGY-FUNCTION, which is the last moment a residue can be added and still appear in
  the atom table.  Joints are NOT built here: they need the atom table, which does not exist yet."
  (let ((built 0))
    (loop for bp-locus across (loci blueprint)
          when (owns-slots-p bp-locus)
            do (let* ((molecule (cdr (assoc (oligomer-shape bp-locus) oligomer-shapes-molecules)))
                      (backbone-monomer (locus-backbone-attachment bp-locus))
                      (backbone-residue (gethash backbone-monomer monomers-to-residues)))
                 (unless molecule
                   (error "No molecule for locus ~d's oligomer-shape" (locus bp-locus)))
                 (unless backbone-residue
                   (error "Locus ~d's backbone monomer has no residue - it may have been ~
                             declined by BUILD-MONOMER-P, which would be a bug: only :SIDECHAIN ~
                             loci are declined." (locus bp-locus)))
                 (let ((foldamer (foldamer (oligomer-space
                                            (oligomer (oligomer-shape bp-locus))))))
                   (loop for bp-monomer across (monomers bp-locus)
                         for cc = (blueprint-monomer-constitution-context bp-monomer bp-locus)
                         do (format t "~&  locus ~3d ~12a foldamer ~a  cc ~s~%"
                                    (locus bp-locus) (monomer-name bp-monomer) (type-of foldamer) cc)
                            ;; A NIL context has to be loud.  ASSIGN-RESIDUE-GIVEN-ATOM-TYPES is a
                            ;; no-op on NIL, so this slot's residues would go unstamped, coverage
                            ;; would fail for the whole molecule, and the cached force field would
                            ;; fall back to real SMIRNOFF over the fanned-out graph - typing a CA
                            ;; that carries 58-79 substituents as though it were a real chemical
                            ;; environment.  That yields a wrong energy function, not a slow one,
                            ;; and nothing downstream would flag it.
                            (unless cc
                              (error "No constitution-context for locus ~d monomer ~a (matcher ~a).~%~
                                      Its slot residues would go unstamped and the cached SMIRNOFF ~
                                      force field would silently fall back to real SMIRNOFF on the ~
                                      fanned-out graph."
                                     (locus bp-locus) (monomer-name bp-monomer)
                                     (matcher bp-monomer)))
                         do (loop for slot below (rotamer-slot-count bp-monomer)
                                  for residue = (build-slot-residue bp-monomer bp-locus
                                                                    molecule backbone-residue foldamer cc)
                                  do (setf (aref (rotamer-scans blueprint)
                                                 (blueprint-mrkindex bp-locus bp-monomer slot))
                                           residue)
                                     (incf built)))
                   )))
    built))

(defun blueprint-slot-parent-joint (blueprint blueprint-locus)
  "The backbone joint a slot's root joint hangs from.

  Same derivation joint-tree.lisp:100-114 uses for an out-plug: find the attachment atom's index in
  the backbone's constitution-atoms, then take that joint from the backbone atresidue."
  (multiple-value-bind (sidechain-atom-name backbone-atom-name)
      (slot-attachment-atom-names (aref (monomers blueprint-locus) 0) blueprint-locus)
    (declare (ignore sidechain-atom-name))
    (let* ((backbone-monomer (locus-backbone-attachment blueprint-locus))
           (position (gethash backbone-monomer (monomer-positions blueprint)))
           (backbone-atresidue (at-position (ataggregate blueprint) position))
           (constitution-atoms (constitution-atoms
                                (constitution (topology backbone-atresidue))))
           (index (position backbone-atom-name constitution-atoms :key #'atom-name)))
      (unless index
        (error "Backbone atom ~a not found in ~a's constitution atoms"
               backbone-atom-name (topology backbone-atresidue)))
      (values (aref (joints backbone-atresidue) index) index))))

(defun build-blueprint-fan-out-joints (blueprint)
  "PASS 2, STAGE 2b.  Build an atresidue and joints for every slot, and replace each staged
  residue in ROTAMER-SCANS with a real ROTAMER-SCAN.
                                        
  Runs AFTER MAKE-ASSEMBLER, because WRITE-INTO-JOINT-TREE takes atomids from the atom table, which
  the energy function owns and which does not exist until then.
                    
  MONOMER and MONOMER-SUBSET are NIL deliberately: WRITE-INTO-JOINT-TREE's bonded and
  in-plug-bonded methods ignore MONOMER on this path, and a NIL subset keeps it out of
  DEAL-WITH-PAIR-SCAN-SIDECHAIN and CREATE-ANCHORED-JOINT - the piece-assembler machinery the
  blueprint exists to replace."         
  (let ((atom-table (chem:atom-table (energy-function blueprint)))
        (adjustments (adjustments blueprint))
        (scans (rotamer-scans blueprint))
        (built 0))                         
    (loop for bp-locus across (loci blueprint)            
          when (owns-slots-p bp-locus)                    
            do (let* ((molecule-index (position (oligomer-shape bp-locus)
                                                (oligomer-shapes blueprint)))
                      (atmolecule (aref (atmolecules (ataggregate blueprint)) molecule-index))
                      (parent-joint (blueprint-slot-parent-joint blueprint bp-locus)))
                 (loop for bp-monomer across (monomers bp-locus) 
                       for topology = (topology bp-monomer)     
                       for constitution-atoms = (constitution-atoms (constitution topology))
                       do (loop for slot below (rotamer-slot-count bp-monomer)
                                for mrkindex = (blueprint-mrkindex bp-locus bp-monomer slot)
                                for residue = (aref scans mrkindex)
                                for atresidue-index = (length (atresidues atmolecule))
                                do (let ((atresidue (make-instance 'atresidue
                                                                   :name (chem:get-name residue)
                                                                   :residue residue
                                                                   :topology topology)))
                                     (resize-atatoms atresidue (length constitution-atoms))
                                     (recursively-write-into-atresidue
                                      nil ; monomer
                                      nil ; nil-or-monomer-shape -> BondedJoints
                                      nil ; monomer-subset
                                      (joint-template topology)
                                      parent-joint
                                      atresidue
                                      molecule-index
                                      atresidue-index
                                      atom-table
                                      adjustments)
                                     (vector-push-extend atresidue (atresidues atmolecule))
                                     (setf (aref scans mrkindex)
                                           (make-instance 'rotamer-scan
                                                          :mrkindex mrkindex
                                                          :residue residue
                                                          :atresidue atresidue))
                                     (incf built))))))
    ;; Initialize the adjustments, mirroring assembler.lisp:816-826.
    ;;
    ;; RECURSIVELY-WRITE-INTO-ATRESIDUE above REGISTERED an adjustment for every slot joint whose
    ;; template carries one (joint-templates.lisp:541,585) - but registration only pushes it onto
    ;; the adjustments object.  MAKE-ASSEMBLER initializes its own atresidues' adjustments as its
    ;; last act, and these atresidues are created after it has returned, so without this loop the
    ;; slot adjustments exist and are applied by ADJUST-ATOM-TREE-EXTERNAL-COORDINATES while
    ;; uninitialized - which MISPLACES the fan rather than leaving it alone.
    ;;
    ;; The whole ataggregate is walked, not just the slots, so the base atresidues' adjustments are
    ;; (re)initialized against the BLUEPRINT.  MAKE-ASSEMBLER initialized them against the
    ;; throwaway base assembler, which is not the object anything folds with afterwards.
    (let ((adjustments (adjustments blueprint)))
      (loop for atmol across (atmolecules (ataggregate blueprint))
            do (loop for atres across (atresidues atmol)
                     do (loop for adjust in (gethash atres (internal-adjustments adjustments))
                              do (initialize-adjustment adjust blueprint))
                        (loop for adjust in (gethash atres (external-adjustments adjustments))
                              do (initialize-adjustment adjust blueprint)))))
    built))
                              
;;; ------------------------------------------------------------------
  ;;; Pass 2, stage 1 - the structure the slots hang off
  ;;;            
  ;;; MAKE-ASSEMBLER does more of this than expected.  It builds one residue per locus for every
  ;;; oligomer-shape, and the FIXED tier comes out right for free: BUILD-ATMOLECULE-USING-OLIGOMER
  ;;; consults the monomer-shape-map and overrides the joint-templates when the shape is a
  ;;; RESIDUE-SHAPE (fold.lisp:135), so those 62 loci get XyzJoints with no blueprint code.
  ;;;       
  ;;; What it does NOT build is the fan-out.  It makes ONE residue at each sidechain locus where
  ;;; the blueprint needs one per slot, so those loci are declined here and materialized in stage
  ;;; 2.  After this runs, the aggregate has a backbone, the fixed tier, and nothing at the 13
  ;;; sidechain loci.
  ;;; ------------------------------------------------------------------

(defparameter *blueprint-rep-weight* 1.0d0
  "The rep-weight baked into the blueprint's rosetta nonbond setup.

Stored energies are rep-weight-INDEPENDENT - they are kept as interleaved (fa_rest, fa_rep) per
mrkindex precisely so one scan serves every weight - so this value does not affect what
FILL-BLUEPRINT-SINGLE-ENERGIES writes.  It matters only to anything reading a COMBINED energy off
this function.")

(defun blueprint-bare-energy-function (aggregate &key (rep-weight *blueprint-rep-weight*))
  "An energy function with an atom table and nothing else.

  MAKE-ASSEMBLER needs one - the joints take their atomids from (CHEM:ATOM-TABLE ...) - but the
  blueprint does not want the terms it would generate: they would pool every slot's contributions
  into one component set, including the cross-slot junk.  The real terms are built afterwards, per
  slot, into EnergyComponentGroups.

  :SETUP IS NOT DECORATION HERE even though nothing is generated.  DEFINE-FOR-AGGREGATE stores the
  full setup on the energy function BEFORE GENERATE-NONBOND-ENERGY-FUNCTION-TABLES takes its
  keep-interaction-factory-is-NIL early return, and GENERATE-INTO-GROUP later REPLAYS that stored
  setup for every slot group.  Omitting it defaults the setup to (:AMBER), so all 470 slot groups
  would come back with amber nonbond terms: different units, no fa_rep channel at all, and
  silently incomparable with the piece-assembler numbers they have to be checked against.  The
  list matches PROTEIN-BINDER.LISP:265, which is what the design movers score with."
  (chem:make-energy-function :matter aggregate
                             :assign-types nil
                             :keep-interaction-factory nil
                             :setup (list :rosetta
                                          (list 'chem:energy-rosetta-nonbond
                                                :rep-weight rep-weight))))

(defun materialize-blueprint (blueprint &key (energy-function-factory #'blueprint-bare-energy-function)
                                                  (ensure-trained t) (verbose nil))
  "PASS 2, STAGE 1.  Build the base structure and bind the inherited ASSEMBLER-BASE slots.

  Returns BLUEPRINT.  AGGREGATE, ENERGY-FUNCTION, ATAGGREGATE, JOINT-TREE, ADJUSTMENTS and
  MONOMER-POSITIONS are bound afterwards - they are deliberately unbound before, so reaching for
  them early signals rather than returning a stale answer.

  ENSURE-TRAINED T (the default) trains the SMIRNOFF parameter cache on this blueprint's trainers
  FIRST.  It has to happen here rather than being left to the caller: MAKE-ASSEMBLER builds the
  energy function, and an untrained cache would be taught by the fanned-out aggregate - a CA with
  58-79 substituents - poisoning every later build silently.  Pass NIL only when the cache is
  known to be trained already, and expect a slow build with no coverage if you are wrong."
  (unless (cdr (oligomer-shapes blueprint))
    (error "MAKE-ASSEMBLER wants :RECEPTOR-ONLY when there is a single oligomer-shape - decide ~
              what that means for a blueprint before calling this."))
  (when ensure-trained
    (ensure-blueprint-trained blueprint :verbose verbose))
  (let ((base (make-assembler
               (oligomer-shapes blueprint) 
               ;; Decline exactly the loci that own mrkindexes.  :FIXED and :BACKBONE build
               ;; normally - one residue each, which is all either needs.
               :build-monomer-p (lambda (monomer oligomer)
                                  (declare (ignore oligomer))
                                  (let ((bp-locus (gethash monomer (locus-map blueprint))))
                                    (not (and bp-locus (owns-slots-p bp-locus)))))
               :after-molecules-fn
               (lambda (aggregate oligomer-shapes-molecules monomers-to-residues monomer-positions)
                 (declare (ignore aggregate monomer-positions))
                 (let ((n (build-blueprint-fan-out blueprint oligomer-shapes-molecules
                                                   monomers-to-residues)))
                   (format t "~&fan-out: built ~d slot residues~%" n)))
               :monomer-contexts (monomer-contexts blueprint)
               :ligand-oligomer-shape (%ligand-oligomer-shape blueprint)
               :receptor-oligomer-shape (%receptor-oligomer-shape blueprint)
               :energy-function-factory energy-function-factory)))
    (setf (aggregate blueprint)         (aggregate base)
          (energy-function blueprint)   (energy-function base)
          (ataggregate blueprint)       (ataggregate base)
          (joint-tree blueprint)        (joint-tree base)
          (adjustments blueprint)       (adjustments base)
          (monomer-positions blueprint) (monomer-positions base))
    (let ((n (build-blueprint-fan-out-joints blueprint)))
      (format t "~&fan-out joints: ~d atresidues~%" n))
    blueprint))


(defun blueprint-materialized-counts (blueprint)
  "Return (values RESIDUES ATRESIDUES) actually built.  They MUST be equal -
  RECURSIVELY-BUILD-ATMOLECULE checks (/= next-atresidue-index next-residue-index) as it goes, so a
  disagreement between the molecule half and the atmolecule half errors there rather than here;
  this is the after-the-fact confirmation."
  (let ((residues 0)
        (atresidues 0))
    (chem:do-molecules (m (aggregate blueprint))
      (incf residues (chem:content-size m)))
    (loop for atmol across (atmolecules (ataggregate blueprint))
          do (incf atresidues (length (atresidues atmol))))
    (values residues atresidues)))


(defun blueprint-monomer-constitution-context (blueprint-monomer blueprint-locus)
  "The constitution-context for this monomer's slot residues.
  
  Needed because the slot residues are built outside the normal path, so MAKE-ASSEMBLER's
  ASSIGN-GIVEN-ATOM-TYPES never sees them - it walks MONOMERS-TO-RESIDUES, and a slot residue has
  no monomer of its own.  Same problem the piece assembler solved with
  ASSIGN-PIECE-GIVEN-ATOM-TYPES (piece-assembler.lisp:406), which notes that
  TOPOLOGY:ASSIGN-GIVEN-ATOM-TYPES cannot be used when the residue's oligomer has no backbone to
  match against.

  MATCH-AS-CONSTITUTION-SYMBOL rather than MATCH-AS-SYMBOL - otherwise identical to
  PERSONA-MONOMER-CONTEXT, including the cache priming that makes the focus resolve to THIS
  monomer's name rather than whichever the oligomer has selected."
  (let ((oligomer (oligomer (oligomer-shape blueprint-locus)))
        (focus (original-monomer blueprint-locus))
        (focus-name (persona-focus-name blueprint-monomer blueprint-locus)))
    (monomer-context:with-match-cache
        (setf (gethash focus (monomer-context::match-cache-monomer-names
                              monomer-context::*match-cache*))
              focus-name)
      (let ((match (monomer-context:match (matcher blueprint-monomer) focus oligomer)))
        (when match
          (monomer-context:match-as-constitution-symbol match))))))


;;; ------------------------------------------------------------------
;;; Internals and coordinates — the blueprint analogue of
;;; MAKE-INTERNALS-FOR-ASSEMBLER / MAKE-COORDINATES-FOR-ASSEMBLER and the
;;; UPDATE-INTERNALS / UPDATE-EXTERNALS pair (assembler.lisp:568-576, 1212, 1468).
;;;
;;; REBUILDS EVERYTHING.  One backbone change invalidates only the loci downstream of it and one
;;; rotamer change invalidates only its own slot, so an incremental version is the real goal - but
;;; correctness first, and a full rebuild is the reference an incremental one gets checked against.
;;;
;;; Two things stop the assembler's own functions being reused verbatim:
;;;
;;;   1. FILL-INTERNALS-FROM-OLIGOMER-SHAPE walks ORDERED-MONOMERS and immediately does
;;;      (molecule-index (gethash monomer (monomer-positions ...))).  The blueprint DECLINED its
;;;      sidechain monomers via BUILD-MONOMER-P, so they have no position and that is a NIL
;;;      dereference.  The loop below skips them.
;;;   2. A slot residue has no MONOMER-SHAPE, so APPLY-MONOMER-SHAPE-TO-ATRESIDUE-INTERNALS has
;;;      nothing to dispatch on.  Slots take the same path by hand: context + rotamer-index ->
;;;      rotamer -> APPLY-FRAGMENT-INTERNALS-TO-ATRESIDUE, exactly as shape.lisp:98-111 does for a
;;;      rotamer-shape and piece-assembler.lisp:432-439 for a piece.
;;; ------------------------------------------------------------------

(defun make-internals-for-blueprint (blueprint)
  "A vector sized to hold every joint's internal coordinates for BLUEPRINT."
  (make-coordinates-for-number-of-atoms (chem:number-of-atoms (aggregate blueprint))))

(defun make-coordinates-for-blueprint (blueprint)
  "A vector sized to hold every atom's external coordinates for BLUEPRINT."
  (make-coordinates-for-number-of-atoms (chem:number-of-atoms (aggregate blueprint))))

(defun blueprint-locus-atresidue (blueprint blueprint-locus)
  "Return the assembler ATRESIDUE represented by a fixed or backbone BLUEPRINT-LOCUS.

Sidechain loci have fan-out atresidues per mrkindex rather than one base atresidue, so accepting
one here would hide an addressing error."
  (unless (member (kind blueprint-locus) '(:fixed :backbone))
    (error "Locus ~d is ~s, not a fixed or backbone locus"
           (locus blueprint-locus) (kind blueprint-locus)))
  (unless (and (< (locus blueprint-locus) (length (loci blueprint)))
               (eq blueprint-locus (aref (loci blueprint) (locus blueprint-locus))))
    (error "~s does not belong to ~s" blueprint-locus blueprint))
  (let ((position (gethash (original-monomer blueprint-locus)
                           (monomer-positions blueprint))))
    (unless position
      (error "No assembler position for locus ~d monomer ~s"
             (locus blueprint-locus) (original-monomer blueprint-locus)))
    (at-position (ataggregate blueprint) position)))

(defun install-refined-backbone-internals
    (persona blueprint-locus refined-internals defined-mask &key verbose)
  "Install one optimized backbone's residue-local internals into PERSONA.

REFINED-INTERNALS has three values per joint of BLUEPRINT-LOCUS's ATRESIDUE.  DEFINED-MASK has
one bit per joint and prevents undefined values from replacing valid geometry.  The full geometry
is stored only in PERSONA's blueprint-sized INTERNALS vector; REFINED-BACKBONE-MASK merely tells a
later rebuild to preserve it.

The caller must subsequently reload sidechain catalogue rotamers using the preserved persona
shape keys, then append the optimized sidechain rotamers.  It must not replace the discrete keys
with classifications derived from the optimized geometry; shape-key drift is diagnostic only.
This function increments GENERATION so energies cached against the old geometry are detectably
stale."
  (let* ((blueprint (blueprint persona))
         (atresidue (blueprint-locus-atresidue blueprint blueprint-locus))
         (joint-count (length (joints atresidue))))
    (unless (eq (kind blueprint-locus) :backbone)
      (error "Cannot install refined backbone internals at ~s" blueprint-locus))
    (unless (internals persona)
      (update-blueprint-internals blueprint persona :verbose verbose))
    (unless (= (length refined-internals) (* 3 joint-count))
      (error "Refined internals have length ~d; locus ~d needs ~d"
             (length refined-internals) (locus blueprint-locus) (* 3 joint-count)))
    (unless (or (null defined-mask) (= (length defined-mask) joint-count))
      (error "Defined mask has length ~d; locus ~d needs ~d"
             (length defined-mask) (locus blueprint-locus) joint-count))
    (write-internals blueprint (internals persona) atresidue
                     refined-internals defined-mask :verbose verbose)
    (setf (aref (refined-backbone-mask persona) (locus blueprint-locus)) 1)
    (loop for blueprint-monomer across (monomers blueprint-locus)
          for ordinal = (monomer-ordinal blueprint-monomer)
          do (setf (aref (backbone-rotamer-index persona) ordinal) nil
                   (aref (active-rotamer-count persona) ordinal) 1))
    (incf (generation persona))
    persona))

(defun install-refined-sidechain-rotamer
    (persona blueprint-locus blueprint-monomer refined-rotamer)
  "Append REFINED-ROTAMER to one sidechain monomer's contiguous active slot range.

The first installation uses the slot at ACTIVE-ROTAMER-COUNT, which is the extra capacity
reserved when the blueprint layout was computed, then increments that count.  If this monomer
already has a refined rotamer, it must be the final active slot and is replaced in place; repeated
refinement therefore does not grow the range.

Returns (values PERSONA SLOT MRKINDEX).  The caller must rebuild the persona internals before
folding coordinates and must refresh energies cached against the previous generation."
  (check-type refined-rotamer sidechain-rotamer)
  (let ((blueprint (blueprint persona)))
    (unless (and (< (locus blueprint-locus) (length (loci blueprint)))
                 (eq blueprint-locus
                     (aref (loci blueprint) (locus blueprint-locus))))
      (error "~s does not belong to ~s" blueprint-locus blueprint))
    (unless (owns-slots-p blueprint-locus)
      (error "Cannot install a refined sidechain rotamer at ~a locus ~d"
             (kind blueprint-locus) (locus blueprint-locus)))
    (unless (find blueprint-monomer (monomers blueprint-locus) :test #'eq)
      (error "~s does not belong to blueprint locus ~d"
             blueprint-monomer (locus blueprint-locus)))
    (let* ((ordinal (monomer-ordinal blueprint-monomer))
           (counts (active-rotamer-count persona))
           (refined (refined-rotamer-index-or-nil persona))
           (count (aref counts ordinal))
           (old-refined-slot (aref refined ordinal))
           (slot-count (rotamer-slot-count blueprint-monomer))
           (slot (or old-refined-slot count)))
      (when old-refined-slot
        (unless (and (plusp count)
                     (= old-refined-slot (1- count))
                     (< old-refined-slot slot-count))
          (error "Locus ~d monomer ~a has inconsistent refined slot ~s and active count ~d"
                 (locus blueprint-locus) (monomer-name blueprint-monomer)
                 old-refined-slot count)))
      (unless (< slot slot-count)
        (error "Locus ~d monomer ~a has no capacity after ~d active rotamers"
               (locus blueprint-locus) (monomer-name blueprint-monomer) count))
      (let ((mrkindex (blueprint-mrkindex blueprint-locus blueprint-monomer slot)))
        (setf (aref (slot-to-rotamer persona) mrkindex) refined-rotamer
              (aref refined ordinal) slot)
        (unless old-refined-slot
          (incf (aref counts ordinal)))
        (incf (generation persona))
        (values persona slot mrkindex)))))

(defun build-refined-persona
    (source-persona backbone-refinements sidechain-refinements &key verbose)
  "Build a complete refined persona without modifying SOURCE-PERSONA.

BACKBONE-REFINEMENTS is a list of

  (BLUEPRINT-LOCUS INTERNALS DEFINED-MASK)

and SIDECHAIN-REFINEMENTS is a list of

  (BLUEPRINT-LOCUS BLUEPRINT-MONOMER SIDECHAIN-ROTAMER).

All records are validated before a private persona copy is made.  Refined backbones are installed
first; catalogue sidechains are then reloaded from the copy's preserved shape keys; refined
sidechains are appended last.  The full blueprint internals vector is rebuilt only after every
conformation is installed.  SOURCE-PERSONA remains unchanged if validation or construction
signals an error.

Returns (values REFINED-PERSONA MAPPINGS).  Each entry of MAPPINGS is

  (:LOCUS BLUEPRINT-LOCUS :MONOMER BLUEPRINT-MONOMER :SLOT SLOT :MRKINDEX MRKINDEX)

and provides the design layer with the new selection needed to remap an MC state.  Energies
cached against SOURCE-PERSONA are not refreshed here and must not be used with the result."
  (check-type source-persona persona)
  (let* ((blueprint (blueprint source-persona))
         (backbone-loci (make-hash-table :test #'eq))
         (sidechain-loci (make-hash-table :test #'eq)))
    (labels ((validate-locus (bp-locus expected-kind)
               (check-type bp-locus blueprint-locus)
               (unless (and (< (locus bp-locus) (length (loci blueprint)))
                            (eq bp-locus (aref (loci blueprint) (locus bp-locus))))
                 (error "~s does not belong to ~s" bp-locus blueprint))
               (unless (eq (kind bp-locus) expected-kind)
                 (error "Locus ~d is ~s, not ~s"
                        (locus bp-locus) (kind bp-locus) expected-kind))))
      ;; Validate the complete request before allocating or changing the result persona.  One
      ;; concrete optimized structure can contribute at most one conformation per locus.
      (dolist (refinement backbone-refinements)
        (destructuring-bind (bp-locus refined-internals defined-mask) refinement
          (validate-locus bp-locus :backbone)
          (when (gethash bp-locus backbone-loci)
            (error "Backbone locus ~d occurs more than once in the refinement request"
                   (locus bp-locus)))
          (setf (gethash bp-locus backbone-loci) t)
          (let ((joint-count
                  (length (joints (blueprint-locus-atresidue blueprint bp-locus)))))
            (unless (= (length refined-internals) (* 3 joint-count))
              (error "Refined internals have length ~d; backbone locus ~d needs ~d"
                     (length refined-internals) (locus bp-locus) (* 3 joint-count)))
            (unless (or (null defined-mask) (= (length defined-mask) joint-count))
              (error "Defined mask has length ~d; backbone locus ~d needs ~d"
                     (length defined-mask) (locus bp-locus) joint-count)))))
      (dolist (refinement sidechain-refinements)
        (destructuring-bind (bp-locus bp-monomer refined-rotamer) refinement
          (validate-locus bp-locus :sidechain)
          (check-type bp-monomer blueprint-monomer)
          (check-type refined-rotamer sidechain-rotamer)
          (unless (find bp-monomer (monomers bp-locus) :test #'eq)
            (error "~s does not belong to blueprint locus ~d"
                   bp-monomer (locus bp-locus)))
          (when (gethash bp-locus sidechain-loci)
            (error "Sidechain locus ~d occurs more than once in the refinement request"
                   (locus bp-locus)))
          (setf (gethash bp-locus sidechain-loci) t)
          (let* ((mrkindex (blueprint-mrkindex bp-locus bp-monomer 0))
                 (atresidue (atresidue (aref (rotamer-scans blueprint) mrkindex)))
                 (needed (* 3 (length (joints atresidue))))
                 (provided (length (internals-values refined-rotamer))))
            (unless (= provided needed)
              (error "Refined rotamer internals have length ~d; locus ~d monomer ~a needs ~d"
                     provided (locus bp-locus) (monomer-name bp-monomer) needed)))))
      (let ((result (copy-persona source-persona))
            (mappings nil))
        (dolist (refinement backbone-refinements)
          (destructuring-bind (bp-locus refined-internals defined-mask) refinement
            (install-refined-backbone-internals
             result bp-locus refined-internals defined-mask :verbose verbose)))
        (multiple-value-bind (loaded overflows)
            (reload-persona-sidechain-rotamers result :verbose verbose)
          (declare (ignore loaded))
          (when overflows
            (error "Refined persona sidechain reload overflowed: ~s" overflows)))
        (dolist (refinement sidechain-refinements)
          (destructuring-bind (bp-locus bp-monomer refined-rotamer) refinement
            (multiple-value-bind (persona slot mrkindex)
                (install-refined-sidechain-rotamer
                 result bp-locus bp-monomer refined-rotamer)
              (declare (ignore persona))
              (push (list :locus bp-locus :monomer bp-monomer
                          :slot slot :mrkindex mrkindex)
                    mappings))))
        (update-blueprint-internals blueprint result :verbose verbose)
        (values result (nreverse mappings))))))

(defun update-blueprint-base-internals (blueprint persona internals &key verbose)
  "Fill internals for the :FIXED and :BACKBONE residues - everything MAKE-ASSEMBLER built.

  Mirrors FILL-INTERNALS-FROM-OLIGOMER-SHAPE (assembler.lisp:1142-1155), including its ordering
  requirement: ORDERED-MONOMERS runs root-outwards so a monomer's internals are installed after
  those of everything it hangs from.  The one difference is the MONOMER-POSITION test - a declined
  sidechain monomer has no position, and is skipped rather than dereferenced.

  A backbone marked in REFINED-BACKBONE-MASK is already authoritative in PERSONA's INTERNALS and is
  preserved instead of being overwritten from its old OLIGOMER-SHAPE."
  (let ((atagg (ataggregate blueprint))
        (filled 0))
    (loop for oligomer-shape in (oligomer-shapes blueprint)
          for oligomer = (oligomer oligomer-shape)
          do (loop for monomer in (ordered-monomers oligomer)
                   for monomer-position = (gethash monomer (monomer-positions blueprint))
                   when monomer-position
                     do (let* ((atmol (elt (atmolecules atagg) (molecule-index monomer-position)))
                               (atres (elt (atresidues atmol) (residue-index monomer-position)))
                               (monomer-context (gethash monomer (monomer-contexts blueprint)))
                               (monomer-shape (gethash monomer (monomer-shape-map oligomer-shape)))
                               (blueprint-locus (gethash monomer (locus-map blueprint)))
                               (refined-p (and blueprint-locus
                                               (eq (kind blueprint-locus) :backbone)
                                               (= 1 (aref (refined-backbone-mask persona)
                                                          (locus blueprint-locus))))))
                          (when (and monomer-shape (not refined-p))
                            (apply-monomer-shape-to-atresidue-internals
                             blueprint internals oligomer-shape monomer-shape
                             monomer-context atres :verbose verbose))
                          (when (or monomer-shape refined-p)
                            (incf filled)))))
    filled))

(defun slot-rotamer-internals (persona mrkindex)
  "The ROTAMER object PERSONA has installed at MRKINDEX, or NIL when that slot is unloaded.

  SLOT-TO-ROTAMER holds rotamer OBJECTS, not indices - FILL-PERSONA puts them there straight from
  the database - so there is no monomer-context lookup and no ROTAMER-VECTOR indexing to do here.
  That also means a slot's geometry cannot go wrong through a mis-derived context: the persona
  already resolved it."
  (aref (slot-to-rotamer persona) mrkindex))

(defun update-blueprint-slot-internals (blueprint persona internals &key verbose)
  "Fill internals for every fan-out slot from the rotamer PERSONA selected for it.

  Runs AFTER the base fill: a slot's root joint hangs off a backbone joint, so its internals are
  only meaningful once that backbone joint's are in place."
  (let ((scans (rotamer-scans blueprint))
        (filled 0)
        (fallbacks 0))
    (loop for bp-locus across (loci blueprint)
          when (owns-slots-p bp-locus)
            do (loop for bp-monomer across (monomers bp-locus)
                     ;; EVERY slot must get internals, including the unloaded ones.  A slot's
                     ;; residue and joints exist whether or not a rotamer was loaded into it - that
                     ;; is the fixed-addressing invariant - and the fold walks the whole tree, so a
                     ;; slot left at zero internals puts its first atom exactly on its parent and
                     ;; KIN:UPDATE-XYZ-COORDS dies normalizing a zero-length b-c.
                     ;;
                     ;; Unloaded slots are spare capacity: the reserved refinement slot, plus the
                     ;; slack between this monomer's rotamer count and the locus-wide bound.  They
                     ;; are never selected, so their geometry is arbitrary - any valid rotamer of
                     ;; the same monomer makes the tree foldable and costs nothing.
                     for fallback = (loop for slot below (rotamer-slot-count bp-monomer)
                                          for r = (slot-rotamer-internals
                                                   persona
                                                   (blueprint-mrkindex bp-locus bp-monomer slot))
                                          when r return r)
                     do (unless fallback
                          (error "Locus ~d monomer ~a has no loaded rotamer in any of its ~d ~
                                  slots, so there is no geometry to fold its residues with.  ~
                                  Either FILL-PERSONA found no rotamers for its context, or the ~
                                  context resolved wrongly."
                                 (locus bp-locus) (monomer-name bp-monomer)
                                 (rotamer-slot-count bp-monomer)))
                        ;; ONE catalogue lookup table per monomer, not per slot.  ROTAMER-VECTOR is
                        ;; thousands of entries for an amino acid and there are up to 76 slots, so a
                        ;; POSITION per slot would be millions of comparisons for a value that is
                        ;; pure bookkeeping.  EQ because the persona holds the very objects
                        ;; ACTIVE-ROTAMERS-FOR-MONOMER selected out of the vector.
                        (let* ((container (blueprint-rotamer-container bp-locus bp-monomer))
                               (catalogue (when container
                                            (let ((h (make-hash-table :test 'eq))
                                                  (v (rotamer-vector container)))
                                              (dotimes (i (length v) h)
                                                (setf (gethash (aref v i) h) i))))))
                          (loop for slot below (rotamer-slot-count bp-monomer)
                                for mrkindex = (blueprint-mrkindex bp-locus bp-monomer slot)
                                for scan = (aref scans mrkindex)
                                for rotamer = (slot-rotamer-internals persona mrkindex)
                                do (unless rotamer (incf fallbacks))
                                   ;; SLOT is the bookkeeping rotamer-index: which of this monomer's
                                   ;; rotamers sits here.  APPLY-FRAGMENT-INTERNALS-TO-ATRESIDUE only
                                   ;; stamps it on the atresidue - the geometry comes from the
                                   ;; rotamer's INTERNALS-VALUES.
                                   (apply-fragment-internals-to-atresidue
                                    blueprint internals (or rotamer fallback) slot (atresidue scan))
                                   ;; Recorded HERE because this is where the fold commits to a
                                   ;; rotamer, and it records what was ACTUALLY applied - so an
                                   ;; unloaded slot reads :UNUSED rather than reporting the stand-in
                                   ;; it was folded with as though it had been chosen.
                                   (setf (rotamer-catalogue-index scan)
                                         (cond ((null rotamer) :unused)
                                               ((null catalogue) :no-container)
                                               ((gethash rotamer catalogue))
                                               (t :geometry-optimized)))
                                   (incf filled)))))
    (when (and verbose (plusp fallbacks))
      (format t "~&  ~d unloaded slots folded with a stand-in rotamer~%" fallbacks))
    (values filled fallbacks)))

(defun update-blueprint-internals (blueprint persona &key internals verbose)
  "Fill and retain PERSONA's internals vector for the WHOLE blueprint.

  PERSONA supplies the rotamer selection.  When INTERNALS is supplied, it becomes this persona's
  vector; otherwise the persona's existing vector is reused or a new one is allocated.  Nothing is
  written into the joints - internals are a value passed to the fold, which is what lets several
  personas coexist over one blueprint.

  The ADJUST-ALL-INTERNALS call is not optional and it is the reason the fan used to land in the
  wrong place.  The normal path is UPDATE-INTERNALS = fill THEN adjust (assembler.lisp:1212-1223);
  filling alone leaves every adjusted dihedral at whatever the raw rotamer internals happened to
  hold.  For an amino-acid sidechain that is CB, whose real dihedral is not stored in the rotamer
  at all - AMBER-PROTEIN's BETA-CARBON adjustment derives it from the sibling C joint as
  c-dih - 120 degrees (aminoacids.lisp:84-94).  Unadjusted, CB lands on top of the carbonyl carbon.

  ADJUST-INTERNALS walks every atresidue of each oligomer-shape's atmolecule (assembler.lisp:1100),
  and the 470 slot atresidues are appended to those same atmolecules, so the fan is reached without
  a separate traversal.  Their adjustments were registered at fan-out and initialized by the
  INITIALIZE-ADJUSTMENT loop in BUILD-BLUEPRINT-FAN-OUT-JOINTS.

  Stores and returns the internals vector."
  (let ((internals (or internals
                       (internals persona)
                       (make-internals-for-blueprint blueprint))))
    (let ((base (update-blueprint-base-internals blueprint persona internals :verbose verbose)))
      (multiple-value-bind (slots fallbacks)
          (update-blueprint-slot-internals blueprint persona internals :verbose verbose)
        (adjust-all-internals blueprint internals)
        (format t "~&blueprint internals: ~d base residues, ~d slots filled~
                   ~@[ (~d with a stand-in rotamer)~], adjusted~%"
                base slots (and (plusp fallbacks) fallbacks))))
    (setf (internals persona) internals)
    internals))

(defun build-blueprint-coordinates (blueprint persona &key verbose internals coords
                                                          (ligand-orientation nil ligand-orientation-p))
  "Rebuild EVERYTHING: internals from PERSONA, then coordinates, then write them into the atoms.

  INTERNALS and COORDS are reused when supplied, and freshly allocated when not.  They must be
  sized for THIS blueprint; nothing checks.  Reuse also means a bug that skips atoms leaves the
  PREVIOUS fold's positions in place rather than zeros, which looks plausible.

  The fold is UPDATE-EXTERNALS, inherited from ASSEMBLER now that BLUEPRINT subclasses it.  Do not
  hand-roll it: for a two-oligomer-shape system it dispatches ligand and receptor separately
  (assembler.lisp:1500-1522), and only the ligand branch applies TRANSFORM-EXTERNALS-TO-GLOBAL-
  FRAME.  Without that the ligand folds in its own local frame and lands nowhere near the receptor.
  Both branches also run ADJUST-ATOM-TREE-EXTERNAL-COORDINATES, which a bare root-joint walk skips.

  The fan comes along free either way: BUILD-ATOM-TREE-EXTERNAL-COORDINATES* recurses from root
  joints (assembler.lisp:1040-1046) and a slot's root joint is a CHILD of its backbone joint, so
  every rotamer is reached without a separate traversal.  ADJUST-ATOM-TREE-EXTERNAL-COORDINATES
  walks all atresidues including the 470 slot ones, which simply have no adjustments registered.

  LIGAND-ORIENTATION IS REQUIRED and deliberately has no default.  It used to default to a fresh
  identity (MAKE-ORIENTATION), which is not a harmless placeholder: TRANSFORM-EXTERNALS-TO-GLOBAL-
  FRAME still ran, still applied a transform, and still reported success - it just applied the
  wrong one.  The ligand folded in its own local frame and landed at the origin, ~150 A from a
  receptor sitting near (-11 28 155), which reads as \"the ligand is missing\" in a viewer even
  though the molecule is present in the aggregate and written to the mol2.  UPDATE-EXTERNALS makes
  the same call and errors rather than defaulting (assembler.lisp:1513-1514); match it.

  The orientation is POSE-rate, not blueprint-rate - a rigid-body move of the ligand changes it
  every MC step - so it cannot be captured on the blueprint at construction.  It comes from the
  pose: (DESIGN:ORIENTATION pose), the same value POSE.LISP:158 threads into UPDATE-EXTERNALS, and
  the sibling slot of the ASSEMBLER the blueprint was built from (POSE.LISP:34 and :52), so the two
  are a matched pair by construction.

  CAVEAT for a pose mid-trajectory: (DESIGN:ORIENTATION pose) is the REFERENCE placement, not
  necessarily the current one.  MONTECARLO.LISP:215-218 folds the ligand through a FRESH
  orientation built from that one's LOCAL-FRAME-SPECS and GLOBAL-POSITIONING-TRANSFORM plus the
  rigid state's ADJUSTMENT-TRANSFORM, which lives in the MCSTATE and not in the orientation slot.
  A pose that has taken rigid-body moves needs that same combination reconstructed here; passing
  the bare slot rebuilds the ligand at its pre-move placement.  Poses straight out of
  COMPUTE-ENERGIES-MOVER have an identity adjustment, which is why the bench is unaffected.

  Returns (values COORDS INTERNALS)."
  (unless ligand-orientation-p
    (error "BUILD-BLUEPRINT-COORDINATES requires :LIGAND-ORIENTATION.  Pass the pose's own ~
            (DESIGN:ORIENTATION pose); a fresh (MAKE-ORIENTATION) silently builds the ligand at ~
            the origin instead of on the receptor."))
  (let* ((internals (update-blueprint-internals blueprint persona
                                                :internals internals
                                                :verbose verbose))
         (coords (update-externals blueprint internals
                                   :coords (or coords (make-coordinates-for-assembler blueprint))
                                   :ligand-orientation ligand-orientation)))
    (copy-all-joint-positions-into-atoms blueprint coords)
    (values coords internals)))


;;; ------------------------------------------------------------------
;;; Atom maps for the per-slot generation passes
;;;
;;; A slot's own energy is the terms INSIDE its fan-out residue plus the terms between that
;;; residue and its backbone.  Two EQ atom maps express exactly that, and the partition falls
;;; out of them without any "not in some other slot" test:
;;;
;;;   at least one atom in the SIDECHAIN map  - otherwise it is a backbone-only term and belongs
;;;                                             to the backbone bucket, not to any slot
;;;   every other atom in the BACKBONE map    - a CB from a DIFFERENT rotamer is in neither map,
;;;                                             so the cross-slot junk CB(i)-CA-CB(j) is rejected
;;;                                             by absence rather than by a special rule
;;;
;;; The maps are keyed by ATOM, not by atom-table index, because the keep-interaction-factory's
;;; predicates are funcalled with Atom_sp objects (energyFunction.cc:1040) - a bit vector would
;;; need an atom->index lookup on every test.
;;;
;;; Per SLOT, not per locus: two rotamers of the same monomer at the same locus are separate
;;; residues holding separate atoms that merely share names.
;;; ------------------------------------------------------------------

(defun slot-residue-of (scan-or-residue)
  "ROTAMER-SCANS holds a bare residue between stage 2a and 2b and a ROTAMER-SCAN after.
Accept either, so map construction is not pinned to one side of that swap."
  (if (typep scan-or-residue 'rotamer-scan)
      (residue scan-or-residue)
      scan-or-residue))

(defun blueprint-slot-atom-map (blueprint mrkindex)
  "An EQ set of the atoms of the fan-out residue at MRKINDEX."
  (let ((map (make-hash-table :test 'eq))
        (residue (slot-residue-of (aref (rotamer-scans blueprint) mrkindex))))
    (unless residue
      (error "No slot residue staged at mrkindex ~d - MATERIALIZE-BLUEPRINT-BASE has not run."
             mrkindex))
    (chem:do-atoms (atom residue)
      (setf (gethash atom map) t))
    map))

(defun blueprint-backbone-atom-maps (blueprint)
  "A vector indexed by molecule-index: the atoms of that molecule that are NOT in any fan-out
residue.

  One map per MOLECULE, which is what keeps a ligand fan-out from generating terms against the
  receptor and vice versa - the intermolecular half is `setMatters`' job, not this pass's.

  'Backbone' here means everything in the molecule that was not fanned out, so it includes
  immobile sidechains.  Those are real intramolecular partners at fixed positions, so their terms
  belong in the slot's energy; narrow this if you want strictly main-chain atoms.

  Call after MATERIALIZE-BLUEPRINT-BASE - it reads the aggregate and the staged slot residues."
  (let* ((agg (aggregate blueprint))
         (slot-atoms (make-hash-table :test 'eq))
         (molecule-count (length (oligomer-shapes blueprint)))
         (maps (make-array molecule-count)))
    (loop for scan across (rotamer-scans blueprint)
          for residue = (slot-residue-of scan)
          when residue
            do (chem:do-atoms (atom residue)
                 (setf (gethash atom slot-atoms) t)))
    (dotimes (molecule-index molecule-count)
      (let ((map (make-hash-table :test 'eq)))
        (chem:map-atoms nil
                        (lambda (atom)
                          (unless (gethash atom slot-atoms)
                            (setf (gethash atom map) t)))
                        (chem:content-at agg molecule-index))
        (setf (aref maps molecule-index) map)))
    maps))

(defun blueprint-locus-molecule-index (blueprint blueprint-locus)
  "Which molecule of the aggregate this locus lives in.  Same derivation
BUILD-BLUEPRINT-FAN-OUT-JOINTS uses."
  (position (oligomer-shape blueprint-locus) (oligomer-shapes blueprint)))

;;; ------------------------------------------------------------------
;;; Scope aggregates - bounding what a parameterization pass walks
;;;
;;; A blueprint's aggregate holds EVERY rotamer of EVERY slot at once.  Generating one slot's bonded
;;; terms walks all of them: the Loop hierarchy behind map-bonds/angles/dihedrals/impropers, ring
;;; perception and SMIRKS matching whenever the SMIRNOFF cache misses, and construct14InteractionTerms.
;;; Every one of those scales with the fan-out while the answer scales with the slot.
;;;
;;; A SCOPE AGGREGATE is a throwaway aggregate holding one throwaway molecule holding BORROWED
;;; residues - the molecule's non-slot residues plus whichever rotamer residue is being parameterized.
;;; Handed to GENERATE-INTO-GROUP as :SCOPE-AGGREGATE it replaces the aggregate for that pass alone.
;;;
;;; WHY BORROWING IS SAFE: ADD-MATTER-DONT-UPDATE-ID appends without touching the child, and
;;; Matter_O::addMatter sets no parent pointer (the setContainedBy call is commented out), so a
;;; residue can sit in two containers at once and nothing points back up to disagree.
;;;
;;; WHY THE BACKBONE MUST BE IN IT, not just the slot residue: Loop emits an angle or dihedral only
;;; from the end with lower unique-atom-order (loop.cc:346).  A term spanning the container boundary
;;; is therefore DROPPED, not duplicated - the CA-CB attachment bond, its angles and the
;;; backbone-CA-CB-CG dihedrals would vanish with no error at all.
;;;
;;; AND THE OTHER FAN-OUT BONDS MUST STILL BE BROKEN.  The container bounds only where a loop STARTS;
;;; every step outward is `bondedNeighbor`, the bond graph (loop.cc:338).  A backbone atom with a live
;;; bond to a detached rotamer reaches into it regardless of what the container holds.  The detach
;;; scan already breaks them - the two mechanisms are complementary, not alternatives.
;;;
;;; NEVER add a scope aggregate's molecule to the REAL aggregate: the residues are already there via
;;; the real molecule, so they would appear twice and the nonbond pass would double-count them.
;;; ------------------------------------------------------------------

(defun blueprint-slot-residue-set (blueprint)
  "EQ hash set of every residue holding fanned-out rotamer atoms.

Includes the 85-odd slots folded with a stand-in rotamer - they are slot residues like any other and
belong out of the backbone just the same."
  (let ((set (make-hash-table :test 'eq)))
    (loop for scan across (rotamer-scans blueprint)
          for residue = (and scan (slot-residue-of scan))
          when residue do (setf (gethash residue set) t))
    set))

(defun make-blueprint-scope-aggregate (blueprint molecule-index slot-residue-set)
  "A scope aggregate for molecule MOLECULE-INDEX - see the section comment.

Returns (values AGGREGATE MOLECULE SLOT-POSITION), where SLOT-POSITION is the index in MOLECULE's
contents reserved for the rotamer residue being parameterized.  Swap it with

    (chem:set-content-at MOLECULE SLOT-POSITION residue)

which costs nothing and, unlike ADD-MATTER, does not renumber the residue.

NIL SLOT-POSITION means this molecule owns no slots, so nothing here will ever be parameterized."
  (let* ((agg (aggregate blueprint))
         (real (chem:content-at agg molecule-index))
         (scope-agg (chem:make-aggregate))
         (scope-mol (chem:make-molecule (chem:get-name real)))
         (first-slot nil))
    ;; The force-field name is not cosmetic: resolveMoleculeForceFields runs over whatever aggregate
    ;; generate-into-group is handed, and a molecule without one cannot be parameterized at all.
    (chem:setf-force-field-name scope-mol (chem:force-field-name real))
    (loop for i below (chem:content-size real)
          for res = (chem:content-at real i)
          if (gethash res slot-residue-set)
            do (unless first-slot (setf first-slot res))
          else
            do (chem:add-matter-dont-update-id scope-mol res))
    (chem:add-matter scope-agg scope-mol)
    (let ((slot-position (when first-slot
                           ;; Reserve the trailing position with any slot residue; every later pass
                           ;; overwrites it before use.
                           (let ((pos (chem:content-size scope-mol)))
                             (chem:add-matter-dont-update-id scope-mol first-slot)
                             pos))))
      (values scope-agg scope-mol slot-position))))

(defun blueprint-locus-of-mrkindex (blueprint mrkindex)
  "The (values LOCUS MONOMER SLOT) that MRKINDEX names, or NIL when nothing owns it.

A linear walk of the loci.  MRKINDEX is (+ (base locus) (base monomer) slot) and is frozen at
construction, so this could be a stored reverse index - it is not, because the only callers are
diagnostics that run over a handful of slots."
  (loop for bp-locus across (loci blueprint)
        when (owns-slots-p bp-locus)
          do (loop for bp-monomer across (monomers bp-locus)
                   do (loop for slot below (rotamer-slot-count bp-monomer)
                            when (= mrkindex (blueprint-mrkindex bp-locus bp-monomer slot))
                              do (return-from blueprint-locus-of-mrkindex
                                   (values bp-locus bp-monomer slot))))))

(defun mrkindex-info (blueprint mrkindex)
  "Everything MRKINDEX names, as a plist.  NIL when no locus owns it.

  THREE INDEX SPACES meet here and conflating them is the recurring bug in this area:

    :BLUEPRINT-LOCUS   counts loci across ALL oligomer-shapes, fixed and moveable alike.
                       COMPUTE-BLUEPRINT-LAYOUT never resets its counter per shape.
    :LOCUS-IN-SHAPE    position within one shape's (MONOMERS SPACE) - what an oligomer-shape
                       itself would call a locus.
    :MCSTATE-INDEX     position among SLOT-OWNING loci only.  This is the mcstate's numbering,
                       compact because sampler.lisp derives ranges by subtraction, and it skips
                       every fixed locus.

  :ROTAMER-SLOT is likewise NOT a rotamers-database index - it is the position within this
  monomer's ACTIVE subset, which is what MRKINDEX is built from.  The database index is what
  BACKBONE-INSTALLED-ROTAMER-INDEX returns; BLUEPRINT-SLOT-FOR-ROTAMER-INDEX converts."
  (multiple-value-bind (bp-locus bp-monomer slot)
      (blueprint-locus-of-mrkindex blueprint mrkindex)
    (when bp-locus
      (let* ((shape (oligomer-shape bp-locus))
             ;; (MONOMERS SPACE) is the space's monomer vector.  MONOMERS is overloaded - on a
             ;; locus-monomer it returns the ALLOWED NAMES instead - so reaching for the wrong
             ;; object here yields a plausible position into an unrelated vector.
             (space (oligomer-space (oligomer shape))))
        (list :mrkindex        mrkindex
              :oligomer-shape  shape
              :shape-index     (position shape (oligomer-shapes blueprint))
              :blueprint-locus (locus bp-locus)
              :locus-in-shape  (position (original-monomer bp-locus) (monomers space))
              ;; Computed here rather than through DESIGN::BLUEPRINT-SLOT-OWNING-LOCI, which is in
              ;; another package and loaded later.  Same enumeration - slot-owning loci in LOCI
              ;; order - and the two must not drift apart.
              :mcstate-index   (position bp-locus (remove-if-not #'owns-slots-p (loci blueprint)))
              :monomer-name    (monomer-name bp-monomer)
              :rotamer-slot    slot
              :slot-count      (rotamer-slot-count bp-monomer)
              ;; The number the pose actually installed, in DATABASE numbering - printed beside
              ;; ROTAMER-SLOT so the two spaces can be compared at a glance.
              :database-index  (backbone-installed-rotamer-index bp-locus))))))

(defparameter +remappable-component-types+
  '(or chem:energy-stretch chem:energy-angle chem:energy-linear-angle
    chem:energy-dihedral chem:energy-nonbond14)
  "The component classes COPY-REMAPPED-BONDED-TERMS reproduces.

ENERGY-NONBOND14 is in here, which looks wrong until you see where its terms come from:
construct14InteractionTerms loops PROPERS over the matter and takes each proper dihedral's terminal
atoms.  Same enumeration and same source as the dihedral terms, so it renames across a monomer's
rotamers for the same reason they do - only the ENERGY is nonbond-shaped.

This list and the C++ dispatch in COPY-REMAPPED-BONDED-TERMS must agree.  If a class is remapped but
missing here, the count assertion compares different things and fails on every slot; if a class is
here but not remapped, the assertion fails the same way.  Either direction is loud, which is the
point - the dangerous version is a class in NEITHER, which is silently dropped from every
copy-remapped rotamer.")

(defun group-bonded-term-count (group)
  "Terms in GROUP's remappable components - see +REMAPPABLE-COMPONENT-TYPES+.

The rosetta nonbond/elec/LK components are excluded deliberately: they are rebuilt from the query set
at evaluation, not copied, so counting them would make the detach scan's equality assertion compare
two different things and always fail."
  (loop for c in (chem:energy-component-group/components group)
        when (typep c +remappable-component-types+)
          sum (chem:number-of-terms c)))

(defun blueprint-locus-molecule (blueprint blueprint-locus)
  "The MOLECULE this locus lives in - what :ONLY-MOLECULE wants.

A slot's keep-interaction factory admits only its own molecule's atoms, and no bonded term spans
molecules, so restricting the generation pass to this one molecule discards nothing.  Without it
every other molecule is re-parameterized once per slot purely to have its terms thrown away."
  (chem:content-at (aggregate blueprint)
                   (blueprint-locus-molecule-index blueprint blueprint-locus)))

(defun make-slot-keep-interaction-factory (sidechain-map backbone-map)
  "The KEEP-INTERACTION-FACTORY for one slot's generation pass.

  Returns a function of a component CLASS returning T, NIL, or a per-tuple predicate - the
  protocol `specializeKeepInteractionFactory` expects (energyFunction.cc:198).

  Fixed-arity predicates for the three bonded classes rather than one &REST version: this is
  called once per candidate tuple across every bond, angle, proper and improper in the molecule,
  470 times over, and consing a fresh list per call there is not free.

  Chiral restraints are refused outright.  The blueprint starts from library geometries so it does
  not need them, and a chiral restraint requires EXACTLY four neighbours while a fanned CA has
  58-79."
  (flet ((side (atom) (gethash atom sidechain-map))
         (scope (atom) (or (gethash atom sidechain-map)
                           (gethash atom backbone-map))))
    (let ((two (lambda (a1 a2)
                 (and (or (side a1) (side a2))
                      (scope a1) (scope a2)
                      t)))
          (three (lambda (a1 a2 a3)
                   (and (or (side a1) (side a2) (side a3))
                        (scope a1) (scope a2) (scope a3)
                        t)))
          (four (lambda (a1 a2 a3 a4)
                  (and (or (side a1) (side a2) (side a3) (side a4))
                       (scope a1) (scope a2) (scope a3) (scope a4)
                       t)))
          (any (lambda (&rest atoms)
                 (and (some #'side atoms) (every #'scope atoms) t))))
      (lambda (component-class)
        (cond ((eq component-class (load-time-value
                                    (find-class 'chem:energy-chiral-restraint)))
               nil)
              ((eq component-class (load-time-value (find-class 'chem:energy-stretch)))  two)
              ((eq component-class (load-time-value (find-class 'chem:energy-angle)))    three)
              ((eq component-class (load-time-value (find-class 'chem:energy-dihedral))) four)
              (t any))))))

(defun blueprint-slot-keep-interaction-factory (blueprint blueprint-locus mrkindex
                                                &optional backbone-maps)
  "Convenience: build the factory for one slot.  Pass BACKBONE-MAPS from
BLUEPRINT-BACKBONE-ATOM-MAPS when doing many slots - it walks the whole aggregate, so rebuilding
it 470 times would dominate."
  (let ((maps (or backbone-maps (blueprint-backbone-atom-maps blueprint))))
    (make-slot-keep-interaction-factory
     (blueprint-slot-atom-map blueprint mrkindex)
     (aref maps (blueprint-locus-molecule-index blueprint blueprint-locus)))))


(defun build-blueprint-slot-energy-groups (blueprint &key verbose)
  "PASS 2, STAGE 3.  Generate every slot's own bonded terms into its own ENERGY-COMPONENT-GROUP
and store it on that slot's ROTAMER-SCAN.

  Runs AFTER BUILD-BLUEPRINT-FAN-OUT-JOINTS: it needs real ROTAMER-SCANs to store into, and it
  reads the atom table, which the energy function owns.

  BONDED ONLY for now.  GENERATE-INTO-GROUP runs pass 2 - stretch/angle/dihedral/improper -
  and `generateNonbondEnergyFunctionTables` does not yet take a group, so a slot's group holds no
  fa_rep.  Since every rotamer of a monomer has nearly identical bonded energy, these groups are
  good for verifying the plumbing and NOT yet usable for choosing rotamers.

  Groups are left ENABLED.  A disabled component returns 0 from `evaluateAllComponent` even when
  evaluated directly, so disabling them would make the per-slot energies this exists to compute
  come back as zero.  The cost is that evaluating the whole energy function is now meaningless -
  it would sum 470 mutually exclusive rotamers into a plausible-looking number.  Nothing should
  do that; an MC step reads cached per-slot energies selected by MCSTATE.

  Returns (values GROUPS-BUILT EMPTY-GROUPS)."
  (let* ((ef (energy-function blueprint))
         (backbone-maps (blueprint-backbone-atom-maps blueprint))
         (scans (rotamer-scans blueprint))
         ;; Same scoping as the detach scan - see MAKE-BLUEPRINT-SCOPE-AGGREGATE.  Built once per
         ;; molecule; only the reserved trailing residue is swapped per slot.
         (slot-residue-set (blueprint-slot-residue-set blueprint))
         (scope-aggregates (make-array (length backbone-maps)))
         (scope-molecules (make-array (length backbone-maps)))
         (scope-positions (make-array (length backbone-maps)))
         (start (get-internal-real-time))
         (built 0)
         (empty 0))
    (dotimes (m (length backbone-maps))
      (multiple-value-bind (agg mol pos)
          (make-blueprint-scope-aggregate blueprint m slot-residue-set)
        (setf (aref scope-aggregates m) agg
              (aref scope-molecules m) mol
              (aref scope-positions m) pos)))
    (loop for bp-locus across (loci blueprint)
          when (owns-slots-p bp-locus)
            do (let* ((locus-mol-index (blueprint-locus-molecule-index blueprint bp-locus))
                      (backbone-map (aref backbone-maps locus-mol-index))
                      (scope (aref scope-aggregates locus-mol-index))
                      (scope-mol (aref scope-molecules locus-mol-index))
                      (scope-pos (aref scope-positions locus-mol-index))
                      (locus-terms 0))
                 (loop for bp-monomer across (monomers bp-locus)
                       do (loop for slot below (rotamer-slot-count bp-monomer)
                                for mrkindex = (blueprint-mrkindex bp-locus bp-monomer slot)
                                for scan = (aref scans mrkindex)
                                do (unless (typep scan 'rotamer-scan)
                                     (error "mrkindex ~d holds ~s, not a ROTAMER-SCAN - ~
                                             BUILD-BLUEPRINT-FAN-OUT-JOINTS has not run."
                                            mrkindex scan))
                                   (let ((factory (make-slot-keep-interaction-factory
                                                   (blueprint-slot-atom-map blueprint mrkindex)
                                                   backbone-map))
                                         (group (chem:make-energy-component-group)))
                                     (chem:setf-energy-component-group-name
                                      group
                                      (list :source :blueprint
                                            :scope :slot
                                            :channel :intramolecular
                                            :terms :bonded
                                            :mrkindex mrkindex
                                            :blueprint-locus (locus bp-locus)
                                            :monomer (monomer-name bp-monomer)
                                            :slot slot))
                                     (chem:set-content-at scope-mol scope-pos
                                                          (slot-residue-of scan))
                                     (chem:energy-function/generate-into-group
                                      ef factory group :scope-aggregate scope)
                                     (setf (intramolecular-energy-components scan) group)
                                     (incf built)
                                     (let ((terms (chem:number-of-terms group)))
                                       (incf locus-terms terms)
                                       (when (zerop terms) (incf empty))))))
                 (when verbose
                   (format t "~&  locus ~3d  ~4d terms over its slots~%"
                           (locus bp-locus) locus-terms)
                   (finish-output))))
    (format t "~&slot energy groups: ~d built in ~,1f s~%"
            built (/ (- (get-internal-real-time) start)
                     internal-time-units-per-second))
    ;; An empty group cannot be right - a rotamer always has bonds among its own atoms - so it
    ;; means the keep-interaction-factory rejected everything.  Reported rather than signalled so
    ;; one run shows how many and at which loci, instead of dying on the first.
    (when (plusp empty)
      (format t "~&*** ~d OF ~d SLOT GROUPS ARE EMPTY - the keep-interaction-factory matched ~
                 nothing.  Check that the slot atom maps hold the fan-out residue's atoms and ~
                 that the backbone map is the right molecule's. ***~%"
              empty built))
    (values built empty)))


;;; ------------------------------------------------------------------
;;; Detach-and-reattach generation
;;;
;;; Removes the CAUSE of the cost rather than filtering its consequences.  With all 470 fan-out
;;; bonds broken and exactly one restored at a time, the CA carrying 58-79 substituents is back to
;;; ordinary valence, so every generation pass enumerates a chemically ordinary graph.  The ~3,000
;;; angles and tens of thousands of dihedrals through that hub - the terms the two-map filter spent
;;; all its time rejecting, 470 times over - simply never exist.
;;;
;;; A keep-interaction-factory is still needed, but a far simpler one.  Breaking a bond DISCONNECTS
;;; a slot residue without removing it from the molecule, so its internal bonds/angles/dihedrals
;;; are still enumerated by LOOPTOPGOAL and would otherwise land in whichever pass is running.
;;; What disappears is the cross-slot case: with one sidechain attached there is no CB(j) to
;;; reject, so the "every other atom in the backbone map" condition is gone and each test is a
;;; single lookup.
;;; ------------------------------------------------------------------

(defun residue-atom-set (residue)
  "An EQ set of RESIDUE's atoms."
  (let ((set (make-hash-table :test 'eq)))
    (chem:do-atoms (atom residue) (setf (gethash atom set) t))
    set))

(defun blueprint-fan-out-bonds (blueprint)
  "Record every fan-out bond as (SIDECHAIN-ATOM BACKBONE-ATOM ORDER-INT MASK), one per slot.

  Found structurally rather than by re-deriving attachment names: a slot residue's fan-out bond is
  its ONLY bond leaving the residue, which BUILD-SLOT-RESIDUE guarantees by refusing any monomer
  without exactly one in-plug.

  ORDER and MASK are both recorded because restoring has to reproduce the bond exactly - the mask
  carries BONDMASKBLUEPRINT, which is what stops a spanning tree walking into every rotamer
  (NEXTNOTCROSSING).  Losing it would not show up until something walked a tree."
  (let ((records nil))
    (loop for scan across (rotamer-scans blueprint)
          for residue = (slot-residue-of scan)
          when residue
            do (let ((mine (residue-atom-set residue))
                     (found nil))
                 (chem:do-atoms (atom residue)
                   (dotimes (i (chem:number-of-bonds atom))
                     (let ((neighbor (chem:bonded-neighbor atom i)))
                       (unless (gethash neighbor mine)
                         (when found
                           (error "Slot residue ~a has more than one bond leaving it - the ~
                                   detach/reattach scan assumes exactly one." residue))
                         (let ((bond (chem:get-bond-to atom neighbor)))
                           (setf found (list atom neighbor
                                             (chem:get-order-from-atom-as-int bond atom)
                                             (chem:bond/mask bond))))))))
                 (unless found
                   (error "Slot residue ~a has no bond leaving it - it was never attached."
                          residue))
                 (push found records)))
    (nreverse records)))

(defun break-fan-out-bond (record)
  (chem:remove-bond-to (first record) (second record)))

(defun restore-fan-out-bond (record)
  "Re-create the bond RECORD describes, mask included.

  :ERROR-IF-EXCEED-VALENCE NIL because the final restore puts all 79 back on one CA.  The recorded
  order is asserted rather than translated: every fan-out bond is created single by
  BUILD-SLOT-RESIDUE, and silently downgrading a non-single one would be worse than stopping."
  (destructuring-bind (sidechain-atom backbone-atom order mask) record
    (unless (= order 1)
      (error "Fan-out bond ~a-~a has order ~d; only single bonds can be restored."
             sidechain-atom backbone-atom order))
    (let ((bond (chem:bond-to sidechain-atom backbone-atom :single-bond
                              :error-if-exceed-valence nil)))
      (chem:bond/setf-mask bond mask))))

;;; ------------------------------------------------------------------
;;; I3-indexed scope tests for the nonbond predicate
;;;
;;; The nonbond predicate is called once per CANDIDATE ATOM PAIR from
;;; rebuildPairListImpl - roughly two million times per component, 470 components over - so what it
;;; costs per call is the whole story.  Hash lookups on the atom objects were 3-6 probes per pair.
;;;
;;; Both tests become integer arithmetic on the i3 values the predicate is already handed:
;;;
;;;   SIDECHAIN  a range test.  CONSTRUCT-FROM-MOLECULE walks residues in order and atoms within a
;;;              residue in order, incrementing idx and coordinateIndex together
;;;              (energyAtomTable.cc:595-597), so one residue's atoms are a CONTIGUOUS block and a
;;;              slot is [lo3, hi3).  Two comparisons, zero heap.
;;;
;;;   BACKBONE   a bit-vector indexed by i3 DIRECTLY, not by i3/3.  That wastes two thirds of the
;;;              bits - only every third one can ever be set - and costs 2.5 KB per molecule rather
;;;              than 846 bytes, to remove a division by 3 from the hot path.  Dividing by 3 is not
;;;              a shift; it compiles to a multiply-high and shift.  The memory is irrelevant here
;;;              and the division is not.
;;;
;;; Together they replace ~700 KB of hash tables with ~5 KB of bit-vectors and a few fixnums.
;;; ------------------------------------------------------------------

(defun blueprint-atom-i3 (atom-table atom)
  "The coordinate index times 3 of ATOM - its index into the coordinate vector, and the value the
nonbond predicate is handed."
  (chem:get-coordinate-index-times3 atom-table atom))

(defun blueprint-slot-i3-range (atom-table sidechain-map)
  "Return (values LO3 HI3), the half-open i3 range covering SIDECHAIN-MAP's atoms.

Correct only because a residue's atoms are contiguous in the atom table - see the section comment.
Computed by min/max over the slot's ~15 atoms rather than read from AtomTable's _ResiduePointers,
which carries an extra terminator entry per molecule (energyAtomTable.cc:611) and is awkward to
index across a two-molecule aggregate.

Signals when the range is not densely occupied by this slot: that would mean the contiguity
assumption has been broken - a residue split, virtual atoms interleaved, a map spanning two
residues - and the range test would then silently accept atoms belonging to something else."
  (let ((lo nil) (hi nil) (n 0))
    (maphash (lambda (atom present)
               (declare (ignore present))
               (let ((i3 (blueprint-atom-i3 atom-table atom)))
                 (incf n)
                 (when (or (null lo) (< i3 lo)) (setf lo i3))
                 (when (or (null hi) (> i3 hi)) (setf hi i3))))
             sidechain-map)
    (when (null lo)
      (error "Empty sidechain map - a slot with no atoms cannot have an i3 range"))
    (let ((span (+ (- hi lo) 3)))
      (unless (= span (* 3 n))
        (error "Slot atoms are not contiguous in the atom table: ~d atoms spanning i3 ~d..~d ~
                (~d slots' worth).  BLUEPRINT-SLOT-I3-RANGE assumes one residue's atoms occupy a ~
                consecutive block, which CONSTRUCT-FROM-MOLECULE guarantees only while the slot is ~
                exactly one residue."
               n lo hi (/ span 3)))
      (values lo (+ hi 3)))))

(defun blueprint-backbone-i3-bits (atom-table backbone-map)
  "A SIMPLE-BIT-VECTOR with bit i3 set for every atom of BACKBONE-MAP.

Indexed by i3 rather than by atom index, so the hot path needs no division - see the section
comment.  Length is 3x the atom count; the two intermediate bits of each triple are never set and
never read."
  (let* ((n3 (* 3 (chem:get-number-of-atoms atom-table)))
         (bits (make-array n3 :element-type 'bit :initial-element 0)))
    (maphash (lambda (atom present)
               (declare (ignore present))
               (let ((i3 (blueprint-atom-i3 atom-table atom)))
                 (unless (< i3 n3)
                   (error "i3 ~d is outside the atom table's ~d coordinate slots" i3 n3))
                 (setf (sbit bits i3) 1)))
             backbone-map)
    bits))

(defun make-detached-slot-factory (sidechain-map backbone-map
                                   &key slot-lo3 slot-hi3 backbone-bits)
  "Keep-interaction-factory for a slot whose bond is the ONLY one attached.

  Two conditions, the same pair MAKE-SLOT-KEEP-INTERACTION-FACTORY uses: at least one atom in this
  slot's sidechain map, and EVERY atom within slot-or-backbone scope.

  The SCOPE test used to be omitted here, on the argument that detaching makes it unnecessary - no
  other rotamer is bonded, so no cross-slot tuple can be enumerated.  That argument is sound for
  BONDED terms and false for NONBOND ones, which is the whole difference: nonbond enumeration does
  not walk bonds.  Without SCOPE every slot's nonbond pairs its atoms against all 469 OTHER slots'
  atoms, which is wrong twice over - rotamers at one locus are mutually exclusive alternatives, and
  two slots at DIFFERENT loci are a PAIR term belonging in the pair matrix rather than in this
  slot's single-scan energy - and it makes each of the 470 passes enumerate the whole atom table
  unbounded.

  Detaching is still what makes the bonded side cheap; SCOPE is what makes the nonbond side
  correct.  The two are not alternatives."
  (flet ((side (atom) (gethash atom sidechain-map))
         (scope (atom) (or (gethash atom sidechain-map)
                           (gethash atom backbone-map))))
    (let ((two   (lambda (a1 a2) (and (or (side a1) (side a2))
                                      (scope a1) (scope a2) t)))
          (three (lambda (a1 a2 a3) (and (or (side a1) (side a2) (side a3))
                                         (scope a1) (scope a2) (scope a3) t)))
          (four  (lambda (a1 a2 a3 a4)
                   (and (or (side a1) (side a2) (side a3) (side a4))
                        (scope a1) (scope a2) (scope a3) (scope a4) t)))
          ;; The NONBOND arity is (atomA atomB ti3A ti3B) - four values, only the first two of
          ;; which are atoms (skipInteraction_EnergyNonbond, energyFunction.cc:1035).
          ;;
          ;; This one is called from rebuildPairListImpl for EVERY candidate atom pair - millions
          ;; per component, 470 components over - so it uses the I3 values and never touches the
          ;; atoms at all.  No hashing, no GC object access, no allocation: a range test for the
          ;; slot and a bit test for the backbone, both on integers already in hand.  Earlier
          ;; versions cost 3-6 hash probes per pair, and the &REST one before that showed up as
          ;; SOME/1 under rebuildPairListImpl in the profile.
          ;;
          ;; Falls back to the hash predicate when the i3 forms were not supplied, so a caller that
          ;; has not built them still gets correct answers rather than a wrong-arity error.
          (nonbond (if (and slot-lo3 slot-hi3 backbone-bits)
                       (let ((lo3 slot-lo3) (hi3 slot-hi3)
                             (bits backbone-bits))
                         ;; CL:TYPE spelled out.  TOPOLOGY used to shadow #:TYPE, which turned a
                         ;; bare (declare (type fixnum ...)) into the type-declaration ABBREVIATION
                         ;; "declare FIXNUM and LO3 to be of type TOPOLOGY::TYPE" - silent at
                         ;; compile time, failing at runtime inside TYPEP against a symbol that
                         ;; appears nowhere in the source.  The shadow is gone now
                         ;; (topology/packages.lisp), so plain TYPE would work, but the explicit
                         ;; form is correct under either package definition and costs nothing -
                         ;; which matters while a compiled fasl may predate the package change.
                         (declare (cl:type fixnum lo3 hi3)
                                  (cl:type simple-bit-vector bits))
                         (lambda (a1 a2 i3a i3b)
                           (declare (ignore a1 a2)
                                    (cl:type fixnum i3a i3b)
                                    (optimize (speed 3) (safety 0)))
                           (let ((side-a (and (>= i3a lo3) (< i3a hi3)))
                                 (side-b (and (>= i3b lo3) (< i3b hi3))))
                             (and (or side-a side-b)
                                  (or side-a (= 1 (sbit bits i3a)))
                                  (or side-b (= 1 (sbit bits i3b)))
                                  t))))
                       (lambda (a1 a2 i3a i3b)
                         (declare (ignore i3a i3b))
                         (and (or (side a1) (side a2))
                              (scope a1) (scope a2) t))))
          (any   (lambda (&rest atoms)
                   (and (some #'side atoms) (every #'scope atoms) t))))
      (lambda (component-class)
        (cond ((eq component-class (load-time-value
                                    (find-class 'chem:energy-chiral-restraint)))
               nil)
              ((eq component-class (load-time-value (find-class 'chem:energy-stretch)))  two)
              ((eq component-class (load-time-value (find-class 'chem:energy-angle)))    three)
              ((eq component-class (load-time-value (find-class 'chem:energy-dihedral))) four)
              ;; Every class routed through skipInteraction_EnergyNonbond gets the fixed-arity
              ;; NONBOND predicate.  ANY still catches the restraint classes, whose arities differ.
              ((or (eq component-class (load-time-value (find-class 'chem:energy-rosetta-nonbond)))
                   (eq component-class (load-time-value (find-class 'chem:energy-rosetta-elec)))
                   (eq component-class (load-time-value (find-class 'chem:energy-rosetta-lksolvation)))
                   (eq component-class (load-time-value (find-class 'chem:energy-nonbond)))
                   (eq component-class (load-time-value (find-class 'chem:energy-nonbond14))))
               nonbond)
              (t any))))))

;;; ------------------------------------------------------------------
;;; The backbone in two channels
;;;
;;; One backbone group cannot serve both: DESIGN::ENERGIES keeps INTRAMOLECULAR-BACKBONE-ENERGY and
;;; INTERMOLECULAR-BACKBONE-ENERGY as two separate interleaved 2-vectors (protein-binder.lisp:
;;; 1703-1708), and the piece path fills them by evaluating two DIFFERENT energy functions
;;; (protein-binder.lisp:2704-2730).  A single group over "no atom in any slot" mixes each
;;; molecule's internal terms with the ligand/receptor cross terms and can fill neither slot
;;; correctly.
;;;
;;; INTRA is "both atoms in the SAME molecule's backbone" - summed over both molecules, because
;;; there is one intramolecular scalar, not one per molecule.
;;; INTER is "one atom in each molecule's backbone", nonbond only: nothing bonded spans molecules.
;;; ------------------------------------------------------------------

(defun make-backbone-intramolecular-factory (backbone-maps &key backbone-bits)
  "Keep-interaction-factory for the INTRAmolecular backbone: every atom in ONE molecule's backbone.

  Run with every fan-out bond broken, so the backbone is a clean molecule.  The slot residues are
  still present but disconnected, and the backbone maps - which exclude every fan-out atom by
  construction - are what keep their terms out.

  Requiring all atoms in the SAME map is what excludes the cross-molecule terms; testing merely
  'not in any slot' would sweep them in here as well as into the intermolecular group, and they
  would then be counted twice."
  (flet ((mol-of (atom)
           (position-if (lambda (m) (gethash atom m)) backbone-maps)))
    (macrolet ((same-mol (&rest atoms)
                 `(let ((m (mol-of ,(first atoms))))
                    (and m ,@(mapcar (lambda (a) `(eql m (mol-of ,a))) (rest atoms)) t))))
      (let ((two   (lambda (a1 a2) (same-mol a1 a2)))
            (three (lambda (a1 a2 a3) (same-mol a1 a2 a3)))
            (four  (lambda (a1 a2 a3 a4) (same-mol a1 a2 a3 a4)))
            ;; i3-indexed bit tests when available - see MAKE-DETACHED-SLOT-FACTORY.  This pass
            ;; runs the full-table rebuildPairListImpl walk (SET-MATTERS cannot express "same
            ;; molecule": it emits cross-set pairs only), so the per-pair cost is all there is to
            ;; optimise here.  It is ONE pass, not 470, so the hash fallback is tolerable.
            (nonbond (if backbone-bits
                         (let ((b0 (aref backbone-bits 0))
                               (b1 (aref backbone-bits 1)))
                           (declare (cl:type simple-bit-vector b0 b1))
                           (lambda (a1 a2 i3a i3b)
                             (declare (ignore a1 a2)
                                      (cl:type fixnum i3a i3b)
                                      (optimize (speed 3) (safety 0)))
                             (or (and (= 1 (sbit b0 i3a)) (= 1 (sbit b0 i3b)))
                                 (and (= 1 (sbit b1 i3a)) (= 1 (sbit b1 i3b))))))
                         (lambda (a1 a2 i3a i3b)
                           (declare (ignore i3a i3b))
                           (same-mol a1 a2))))
            (any   (lambda (&rest atoms)
                     (let ((m (mol-of (first atoms))))
                       (and m (every (lambda (a) (eql m (mol-of a))) (rest atoms)) t)))))
        (lambda (component-class)
          (cond ((eq component-class (load-time-value
                                      (find-class 'chem:energy-chiral-restraint)))
                 nil)
                ((eq component-class (load-time-value (find-class 'chem:energy-stretch)))  two)
                ((eq component-class (load-time-value (find-class 'chem:energy-angle)))    three)
                ((eq component-class (load-time-value (find-class 'chem:energy-dihedral))) four)
                ((or (eq component-class (load-time-value (find-class 'chem:energy-rosetta-nonbond)))
                     (eq component-class (load-time-value (find-class 'chem:energy-rosetta-elec)))
                     (eq component-class (load-time-value (find-class 'chem:energy-rosetta-lksolvation)))
                     (eq component-class (load-time-value (find-class 'chem:energy-nonbond)))
                     (eq component-class (load-time-value (find-class 'chem:energy-nonbond14))))
                 nonbond)
                (t any)))))))

(defun make-backbone-intermolecular-factory (backbone-maps &key backbone-bits)
  "Keep-interaction-factory for the INTERmolecular backbone: one atom in EACH molecule's backbone.

  Bonded classes are refused outright - nothing bonded spans two molecules - so this group holds
  only nonbond components, exactly as the piece path's intermolecular energy function does
  (protein-binder.lisp:1213 filters to the rosetta nonbonded components).

  The predicate is a backstop; SET-GROUP-NONBOND-MATTERS expresses the same restriction
  structurally and is what makes this affordable."
  (flet ((in0 (atom) (gethash atom (aref backbone-maps 0)))
         (in1 (atom) (gethash atom (aref backbone-maps 1))))
    (let ((nonbond (if backbone-bits
                       (let ((b0 (aref backbone-bits 0))
                             (b1 (aref backbone-bits 1)))
                         (declare (cl:type simple-bit-vector b0 b1))
                         (lambda (a1 a2 i3a i3b)
                           (declare (ignore a1 a2)
                                    (cl:type fixnum i3a i3b)
                                    (optimize (speed 3) (safety 0)))
                           (or (and (= 1 (sbit b0 i3a)) (= 1 (sbit b1 i3b)))
                               (and (= 1 (sbit b1 i3a)) (= 1 (sbit b0 i3b))))))
                       (lambda (a1 a2 i3a i3b)
                         (declare (ignore i3a i3b))
                         (and (or (and (in0 a1) (in1 a2))
                                  (and (in1 a1) (in0 a2)))
                              t)))))
      (lambda (component-class)
        (cond ((or (eq component-class (load-time-value (find-class 'chem:energy-rosetta-nonbond)))
                   (eq component-class (load-time-value (find-class 'chem:energy-rosetta-elec)))
                   (eq component-class (load-time-value (find-class 'chem:energy-rosetta-lksolvation)))
                   (eq component-class (load-time-value (find-class 'chem:energy-nonbond)))
                   (eq component-class (load-time-value (find-class 'chem:energy-nonbond14))))
               nonbond)
              (t nil))))))

(defun atom-set-to-vector (atom-set)
  "A vector of the atoms in ATOM-SET, the form COLLECT-ATOMS-FROM-MATTER accepts directly
(pairList.h:156-161 takes a SimpleVector of atoms as an alternative to a Matter)."
  (let ((v (make-array (hash-table-count atom-set)))
        (i 0))
    (maphash (lambda (atom present)
               (declare (ignore present))
               (setf (aref v i) atom)
               (incf i))
             atom-set)
    v))

(defun group-pairlist-components (group)
  "The components of GROUP that build a pair list.

  Tested with TYPEP against CHEM:ENERGY-PAIRLIST-COMPONENT rather than by enumerating class names.
  That base class exists precisely so this question has one answer: EnergyNonbond, the three rosetta
  components and anything added later inherit it, while EnergyNonbond14 - which builds its terms
  from the bond graph and has no pair list - does not."
  (remove-if-not (lambda (c) (typep c 'chem:energy-pairlist-component))
                 (chem:energy-component-group/components group)))

(defun set-group-query-cutoff (comps)
  "Set every component's query cutoff to the MAX rpairlist across COMPS.

  So the shared neighbour list is built once at the widest cutoff regardless of which component
  evaluates first.  The three rosetta components differ - nonbond 9.0, LK 8.0, elec 7.5 - and
  validity is `built >= requested`, so the widest answers all three while each still applies its
  own exact distance test per pair.

  Without this the first component to evaluate sets the width.  That works today only because
  components are added nonbond-first; overriding :rpairlist to make nonbond narrower than LK would
  make the first build too small and force a second full pass.  PROTEIN-BINDER.LISP:2328 takes the
  same max for the same reason."
  (when comps
    (let ((max-cutoff (reduce #'max comps :key #'chem:rpairlist)))
      (dolist (c comps) (chem:set-query-cutoff c max-cutoff))
      max-cutoff)))

(defun set-group-query-range (group lo3 hi3 target-bits)
  "Enumerate GROUP's pair lists from the i3 range [LO3,HI3), finding only atoms set in TARGET-BITS.

  THE RANGE FORM ALSO GENERATES THE INTRA-QUERY PAIRS.  REBUILD-PAIR-LIST-FROM-NEIGHBORS-IMPL runs
  a second, brute-force loop over the range for exactly this reason: a slot's own atoms are not in
  the target set, so they never appear in anyone's neighbour list, and the 1-5-and-beyond contacts
  inside a rotamer are real energy that would otherwise vanish.  There is nothing else to set up -
  choosing a RANGE is what asks for them; a MASK does not.

  For an INTERmolecular group those intra-slot pairs are generated and then rejected by the
  keep-interaction predicate, which demands one atom in the partner molecule.  ~105 wasted
  candidates per slot, and the alternative - a flag to suppress the loop - would be one more thing
  to get wrong for no measurable gain."
  (let ((comps (group-pairlist-components group)))
    (set-group-query-cutoff comps)
    (dolist (comp comps)
      (chem:set-query-range comp lo3 hi3)
      (chem:set-query-targets comp target-bits))))

(defun set-group-query-mask (group query-bits target-bits)
  "Enumerate GROUP's pair lists from the atoms set in QUERY-BITS, finding only those in TARGET-BITS.

  For the backbone groups, whose query atoms are not contiguous.  No intra-query loop runs here and
  none is wanted: a mask-driven caller's own atoms are in the target set, so those pairs already
  come through the neighbour lists - with the i3 < j3 guard in the impl keeping each to one copy."
  (let ((comps (group-pairlist-components group)))
    (set-group-query-cutoff comps)
    (dolist (comp comps)
      (chem:set-query-mask comp query-bits)
      (chem:set-query-targets comp target-bits))))

(defun blueprint-union-bits (a b)
  "Bitwise OR of two i3-indexed bit-vectors - the backbone of both molecules as one target set."
  (bit-ior a b))

(defun set-group-nonbond-matters (group slot-atom-vector partner-atom-vector)
  "Restrict GROUP's nonbond components to SLOT-ATOM-VECTOR x PARTNER-ATOM-VECTOR.

  FOR THE INTERMOLECULAR CHANNEL ONLY.  The between-matters loop discards every pair with both
  atoms in the same set (pairList.h:291), so a slot's own internal nonbond is not emitted.  Across
  two molecules that is exactly right - a rotamer's internal sterics are not intermolecular energy,
  and there is no bonded path between the sets.  Applied to the INTRAmolecular channel the same
  behaviour would silently drop each rotamer's internal nonbond, which is a real part of the
  quantity rotamers are ranked by.

  Without this a nonbond component rebuilds its pair list over the WHOLE atom table:
  REBUILD-PAIR-LIST-IMPL grids all 6762 atoms and enumerates every neighbour pair within rpairlist
  - about two million candidates - before the keep-interaction predicate sees any of them
  (pairList.h:390-440).  Measured at 11 s per slot, dead flat across slots, which is the signature
  of constant full-table work rather than of anything being amortised.

  SET-MATTERS switches the component to REBUILD-PAIR-LIST-BETWEEN-MATTERS-IMPL, which enumerates
  only matter1 x matter2 - a slot's ~15 atoms against the partner molecule.  Same mechanism the
  piece-assembler uses per scan (protein-binder.lisp:1223).

  THE TWO SETS MUST BE DISJOINT.  An atom in both lands at two indices in the concatenated vector,
  passes the cross-set test against ITSELF, and yields a self-pair at distance zero.  Two different
  molecules are disjoint by construction; do not reuse this with sets that can overlap.

  EXCLUSIONS-POSSIBLE is passed NIL: no bonded path spans two molecules, so the 1-2/1-3/1-4 probes
  would always miss and are pure cost per pair.  This mirrors protein-binder.lisp:1223-1226."
  (loop for comp in (chem:energy-component-group/components group)
        when (member (class-name (class-of comp))
                     '(chem:energy-rosetta-nonbond
                       chem:energy-rosetta-elec
                       chem:energy-rosetta-lksolvation))
          do (chem:set-matters comp slot-atom-vector partner-atom-vector nil)))

(defun blueprint-partner-molecule-index (blueprint molecule-index)
  "The index of the OTHER molecule - the intermolecular partner of MOLECULE-INDEX.

  Only meaningful for the two-oligomer-shape ligand/receptor system; anything else has no single
  partner and the caller has to say what it means."
  (let ((n (length (oligomer-shapes blueprint))))
    (unless (= n 2)
      (error "BLUEPRINT-PARTNER-MOLECULE-INDEX needs exactly 2 oligomer-shapes, found ~d" n))
    (- 1 molecule-index)))

(defun make-intermolecular-slot-factory (sidechain-map partner-map)
  "Keep-interaction-factory for one slot's INTERMOLECULAR pass.

  Bonded classes are refused outright: no bond, angle or dihedral spans two molecules, so anything
  those would enumerate here is junk by definition.  What remains is the nonbond set, restricted to
  one atom in this slot and the other in the PARTNER molecule.

  PARTNER-MAP is the partner's BACKBONE map, not all of its atoms.  A pair of atoms from two
  different fan-out slots is a PAIR term belonging in the pair matrix - it is a function of two
  rotamer choices, not one - so it must not land in either slot's single-scan energy.

  This predicate is a backstop rather than the mechanism: SET-GROUP-NONBOND-MATTERS expresses the
  same restriction structurally and is what makes it cheap.  It is kept because it is the only
  thing enforcing the rule if matters are ever not set."
  (flet ((side (atom) (gethash atom sidechain-map))
         (partner (atom) (gethash atom partner-map)))
    (let ((nonbond (lambda (a1 a2 i3a i3b)
                     (declare (ignore i3a i3b))
                     (and (or (and (side a1) (partner a2))
                              (and (side a2) (partner a1)))
                          t))))
      (lambda (component-class)
        (cond ((eq component-class (load-time-value (find-class 'chem:energy-stretch)))  nil)
              ((eq component-class (load-time-value (find-class 'chem:energy-angle)))    nil)
              ((eq component-class (load-time-value (find-class 'chem:energy-dihedral))) nil)
              ((eq component-class (load-time-value
                                    (find-class 'chem:energy-chiral-restraint)))         nil)
              ((or (eq component-class (load-time-value (find-class 'chem:energy-rosetta-nonbond)))
                   (eq component-class (load-time-value (find-class 'chem:energy-rosetta-elec)))
                   (eq component-class (load-time-value (find-class 'chem:energy-rosetta-lksolvation)))
                   (eq component-class (load-time-value (find-class 'chem:energy-nonbond)))
                   (eq component-class (load-time-value (find-class 'chem:energy-nonbond14))))
               nonbond)
              (t nil))))))

;;; ------------------------------------------------------------------
;;; Pair-scan - ONE group, matters re-set per pair
;;;
;;; The single-scan groups are per slot because each holds different TERMS.  A pair group does not:
;;; the terms for slot i against slot j are entirely determined by the two atom sets, which
;;; SET-MATTERS supplies at evaluation time.  So there is ONE group, re-pointed 110,685 times,
;;; not 110,685 groups.  Same pattern the piece-assembler uses - the comment at pairList.h:150-153
;;; notes that fill-piece-pair-energies "sets new matters before every single evaluation".
;;;
;;; NOT the query set.  A query set is baked in at setup and the neighbour list is masked to the
;;; BACKBONE, which contains no slot atoms at all - it cannot answer "which of slot j's atoms are
;;; near slot i".  Matters are the right mechanism here and they are also mutually exclusive with a
;;; query set (rebuildPairListImpl errors if both are set).
;;;
;;; ONE energy matrix, not two.  Unlike backbone and single-scan, ENERGIES keeps a single
;;; PAIR-SCAN-ENERGY-LOWER-TRIANGULAR-MATRIX plus a per-pair INTERMOLECULAR-P flag
;;; (protein-binder.lisp:1719-1727).  So the intra/inter distinction is recorded per pair by the
;;; caller rather than by having two groups.
;;;
;;; Only CROSS-LOCUS pairs are meaningful: two rotamers at the SAME locus are mutually exclusive
;;; alternatives that never coexist, so their interaction energy is not a thing that can be
;;; observed.  Most cross-locus pairs are also too far apart to interact at any rotamer choice -
;;; that is what the locus bounding-sphere filter is for (see PROTEIN-BINDER.LISP:2326-2331, which
;;; already does exactly this for pieces).
;;; ------------------------------------------------------------------

(defun make-pair-scan-factory ()
  "Keep-interaction-factory for the pair-scan group: nonbond only, unrestricted.

  Bonded classes are refused - two rotamers at different loci are not bonded to each other, so
  anything a bond, angle or dihedral would enumerate between them is junk by definition.  Refusing
  ENERGY-DIHEDRAL also suppresses ENERGY-NONBOND14, whose creation is guarded on that class
  (energyFunction.cc:1903) rather than on its own; 1-4 terms are bond-graph derived and cannot span
  two unbonded rotamers.

  The nonbond predicate is T rather than an atom test.  It cannot be anything else: the predicate
  is fixed when the group is generated, while WHICH pair of slots is being scored changes on every
  evaluation.  SET-MATTERS is what bounds each evaluation, and the between-matters loop emits only
  cross-set pairs (pairList.h:291), so the restriction is structural and exact."
  (lambda (component-class)
    (cond ((or (eq component-class (load-time-value (find-class 'chem:energy-rosetta-nonbond)))
               (eq component-class (load-time-value (find-class 'chem:energy-rosetta-elec)))
               (eq component-class (load-time-value (find-class 'chem:energy-rosetta-lksolvation)))
               (eq component-class (load-time-value (find-class 'chem:energy-nonbond))))
           t)
          (t nil))))

(defun build-blueprint-pair-scan-group (blueprint &key verbose)
  "Generate the ONE pair-scan energy component group.  Returns it.

  Generated with :BONDED NIL, so no bonded terms are produced AND no per-molecule parameterization
  is run.  That matters beyond the wasted work: pass 2 parameterizes before any keep-predicate is
  consulted, so with the fan-out attached a spiro NG carries 16 bonds and SMIRNOFF's bond-order
  check kills the run.  With :BONDED NIL this is insensitive to bond state and can be generated
  anywhere.

  Exclusions are unaffected either way - they are computed once at ATOM TABLE construction
  (energyAtomTable.cc:576-593), not here.

  The group is left with NO matters set.  SET-PAIR-SCAN-MATTERS points it at a specific pair
  immediately before each evaluation; until then it would enumerate the whole atom table, so it
  must not be evaluated as-is."
  (let ((ef (energy-function blueprint))
        (group (chem:make-energy-component-group)))
    (chem:setf-energy-component-group-name
     group
     '(:blueprint :pair-scan))
    ;; :BONDED NIL - generate ONLY the nonbond components.  Refusing the bonded classes in the
    ;; factory would discard their terms but still re-run the whole per-molecule force-field
    ;; parameterization first, which is both wasted work and what made this sensitive to bond state.
    (chem:energy-function/generate-into-group ef (make-pair-scan-factory) group :bonded nil)
    (when verbose
      (format t "~&pair-scan group: ~d components~%"
              (length (chem:energy-component-group/components group))))
    group))

(defmacro with-pinned-coordinates ((atom-table) &body body)
  "Run BODY with the per-evaluation drift walk switched off.

  BODY MUST NOT change the coordinates, and must invalidate every component it evaluates -
  SET-PAIR-SCAN-MATTERS does the second by calling SET-MATTERS, which invalidates.  Both hold for a
  pair scan: one coordinate vector for the whole scan, and every component re-pointed each
  iteration.

  Worth having because the walk is O(natoms) and EVALUATE-ALL runs it once per evaluation
  (energyFunction.cc:625).  A pair scan evaluates n^2/2 times over coordinates that never move, so
  the walk recomputes zero for a quarter of the scan's runtime and every component then discards
  the answer, having just been invalidated.

  The snapshot is invalidated on the way out rather than left in place.  It should still be valid -
  the coordinates did not move - but that rests entirely on BODY having kept its promise, and one
  forced rebuild is a small price beside a silently stale pair list.

  RESTORES the previous pin state rather than clearing it, so nesting is safe.  Hard-setting NIL
  would leave an inner region silently unpinning its enclosing one, and the symptom would be a
  performance regression that looks like nothing at all."
  (let ((at (gensym "ATOM-TABLE"))
        (was (gensym "WAS-PINNED")))
    `(let* ((,at ,atom-table)
            (,was (chem:coordinates-pinned-p ,at)))
       (chem:set-coordinates-pinned ,at t)
       (unwind-protect (progn ,@body)
         (chem:set-coordinates-pinned ,at ,was)
         ;; Only the OUTERMOST region invalidates.  An inner one doing it would force every
         ;; component to rebuild in the middle of the outer scan - correct, but it would give back
         ;; exactly the cost this is here to remove.
         (unless ,was (chem:invalidate-drift-snapshot ,at))))))

(defun set-pair-scan-matters (group atoms-i atoms-j &key same-molecule-p components)
  "Point GROUP's components at one pair of slots.

  COMPONENTS is GROUP-PAIRLIST-COMPONENTS' answer, hoisted by the caller.  It does not change while
  the group is being re-pointed, and computing it here means a REMOVE-IF-NOT - a fresh list, from a
  fresh traversal - on every one of a pair scan's n^2/2 calls.  That showed up as 1.4% of a scan
  doing nothing but rebuilding the same three-element list.

  ATOMS-I and ATOMS-J are atom VECTORS - ATOM-SET-TO-VECTOR of each slot's atom map.  They must be
  DISJOINT, which holds for two different slots by construction: an atom in both would land at two
  indices in the concatenated vector, pass the cross-set test against itself, and yield a self-pair
  at distance zero.

  SAME-MOLECULE-P leaves exclusions ON.  Two sidechains at adjacent residues are separated by five
  bonds at the closest (CB-CA-C-N-CA-CB), so no 1-2/1-3/1-4 pair can actually span them and the
  probes will always miss - but the cost of being wrong is silently keeping a pair that should have
  been excluded, and the cost of the probes is three set lookups.  Across molecules there is no
  bonded path at all, so they are switched off."
  (dolist (comp (or components (group-pairlist-components group)))
    (chem:set-matters comp atoms-i atoms-j (if same-molecule-p t nil))))

(defun blueprint-all-slot-atoms (blueprint)
  "An EQ set of every atom in every fan-out residue."
  (let ((set (make-hash-table :test 'eq)))
    (loop for scan across (rotamer-scans blueprint)
          for residue = (slot-residue-of scan)
          when residue
            do (chem:do-atoms (atom residue) (setf (gethash atom set) t)))
    set))

#+(or)
(defmacro fformat (stream-sel fmt &rest args)
  `(progn
     (format ,stream-sel ,fmt ,@args)
     (finish-output ,stream-sel)))

(defmacro fformat (stream-sel fmt &rest args)
  `(progn))

(defun build-blueprint-energy-groups-by-detach (blueprint &key verbose)
  "Break every fan-out bond, generate the backbone group, then reattach one slot at a time to
generate that slot's group, and finally restore every bond.

  UNWIND-PROTECT around the whole scan is not optional: an error partway through would otherwise
  leave the aggregate with hundreds of bonds missing and every later result quietly wrong.

  Returns (values SLOT-GROUPS BACKBONE-INTRA-GROUP BACKBONE-INTER-GROUP PAIR-GROUP)."
  (let* ((ef (energy-function blueprint))
         (records (blueprint-fan-out-bonds blueprint))
         (all-slot-atoms (blueprint-all-slot-atoms blueprint))
         (scans (rotamer-scans blueprint))
         ;; Hoisted: BLUEPRINT-BACKBONE-ATOM-MAPS walks the whole aggregate, so building it per
         ;; slot would dominate the scan.  One map per molecule, which is also what keeps a ligand
         ;; slot from generating terms against the receptor.
         (backbone-maps (blueprint-backbone-atom-maps blueprint))
         ;; One atom VECTOR per molecule, built once - SET-GROUP-NONBOND-MATTERS wants a vector
         ;; rather than the hash set, and rebuilding a ~6200-element vector per slot would put back
         ;; a good part of the cost this exists to remove.
         (backbone-vectors (map 'vector #'atom-set-to-vector backbone-maps))
         (atom-table (chem:atom-table ef))
         ;; One i3-indexed bit-vector per molecule, built once.  ~2.5 KB each, replacing a ~6200
         ;; entry hash table that was probed 3-6 times per candidate atom pair.
         (backbone-bits (map 'vector
                             (lambda (m) (blueprint-backbone-i3-bits atom-table m))
                             backbone-maps))
         ;; ONE target mask for EVERY group.  Not per-molecule: ATOM-TABLE's neighbour list is
         ;; cached against the target mask by IDENTITY (neighborListValidFor), so two groups asking
         ;; with different masks invalidate each other.  Per-molecule masks would have made a
         ;; slot's intra and inter groups rebuild the list in turn - ~940 full rebuilds across the
         ;; fill, worse than no cache at all.
         ;;
         ;; Safe because the list is purely GEOMETRIC: which molecule an atom belongs to is decided
         ;; by the keep-interaction predicate's own per-molecule bit tests, not by what the list
         ;; contains.  Widening the target set can only add candidates the predicate then rejects.
         (all-backbone-bits (blueprint-union-bits (aref backbone-bits 0) (aref backbone-bits 1)))
         ;; ONE scope aggregate per molecule, built once.  Its backbone half never changes across the
         ;; scan; only the reserved trailing residue is swapped, per pass, by SET-CONTENT-AT.  See
         ;; the section comment on MAKE-BLUEPRINT-SCOPE-AGGREGATE for why the backbone has to be in
         ;; it and why the other fan-out bonds must still be broken.
         (slot-residue-set (blueprint-slot-residue-set blueprint))
         (scope-aggregates (make-array (length backbone-maps)))
         (scope-molecules (make-array (length backbone-maps)))
         (scope-positions (make-array (length backbone-maps)))
         (start (get-internal-real-time))
         (backbone-group nil)
         (backbone-inter-group nil)
         (pair-group nil)
         (built 0)
         (empty 0)
         ;; REPARAMETERIZED is one per (locus, monomer); REMAPPED is every other rotamer.  Their
         ;; ratio is the whole point of the copy path - if REMAPPED is near zero the monomers each
         ;; have one rotamer and this bought nothing.
         (reparameterized 0)
         (remapped 0))
    (declare (ignorable all-slot-atoms))
    (format t "~&detach scan: ~d fan-out bonds recorded~%" (length records))
    (dotimes (m (length backbone-maps))
      (multiple-value-bind (agg mol pos)
          (make-blueprint-scope-aggregate blueprint m slot-residue-set)
        (setf (aref scope-aggregates m) agg
              (aref scope-molecules m) mol
              (aref scope-positions m) pos)))
    (unwind-protect
         (progn
           (mapc #'break-fan-out-bond records)
           ;; ---- backbone, with nothing attached, in BOTH channels ----
           ;;
           ;; Two groups, not one: ENERGIES keeps INTRAMOLECULAR-BACKBONE-ENERGY and
           ;; INTERMOLECULAR-BACKBONE-ENERGY separately, and a single "no atom in any slot" group
           ;; mixes each molecule's internal terms with the ligand/receptor cross terms.
           (setf backbone-group (chem:make-energy-component-group))
           (chem:setf-energy-component-group-name backbone-group '(:blueprint :backbone :intramolecular))
           (chem:energy-function/generate-into-group
            ef (make-backbone-intramolecular-factory backbone-maps
                                                     :backbone-bits backbone-bits)
            backbone-group)
           ;; Query FROM every backbone atom of both molecules, find only backbone atoms.  The
           ;; same-molecule test stays in the predicate; the query set only bounds enumeration.
           (set-group-query-mask backbone-group all-backbone-bits all-backbone-bits)
           (format t "~&  backbone intra group: ~d terms~%"
                   (chem:number-of-terms backbone-group))
           (setf backbone-inter-group (chem:make-energy-component-group))
           (chem:setf-energy-component-group-name
            backbone-inter-group
            '(:blueprint :backbone :intermolecular))
           ;; :BONDED NIL - this factory refuses every bonded class, so the per-molecule pass
           ;; would re-parameterize each molecule only to have every term it produced discarded.
           (chem:energy-function/generate-into-group
            ef (make-backbone-intermolecular-factory backbone-maps
                                                     :backbone-bits backbone-bits)
            backbone-inter-group :bonded nil)
           ;; Query FROM all backbone, targets all backbone - the SAME mask as every other group,
           ;; so one neighbour list serves them all.  The predicate picks the cross-molecule pairs;
           ;; the query atoms are in the target set, so the i3 < j3 guard dedups.
           ;;
           ;; Querying only molecule 0 would look tighter and would be a TRAP: with a shared target
           ;; mask the query atoms are targets either way, so the i3 < j3 guard applies, and it
           ;; would then drop any pair whose molecule-1 partner sorted lower.  That happens to be
           ;; none of them today only because molecule 0 occupies the low end of the atom table.
           (set-group-query-mask backbone-inter-group all-backbone-bits all-backbone-bits)
           ;; Pair-scan group: generated HERE, while the bonds are broken, for the reason in
           ;; BUILD-BLUEPRINT-PAIR-SCAN-GROUP's docstring - pass 2 parameterizes every molecule
           ;; before any keep-predicate runs, and a fully attached spiro NG has 16 bonds.
           (setf pair-group (build-blueprint-pair-scan-group blueprint :verbose verbose))
           (format t "~&  backbone inter group: ~d terms~%"
                   (chem:number-of-terms backbone-inter-group))
           ;; ---- one slot at a time ----
           (loop for bp-locus across (loci blueprint)
                 when (owns-slots-p bp-locus)
                   do (let* ((locus-terms 0)
                             (locus-mol-index (blueprint-locus-molecule-index blueprint bp-locus))
                             ;; Built once per molecule and reused, since the backbone half never
                             ;; changes - only the reserved trailing residue does.
                             (scope (aref scope-aggregates locus-mol-index))
                             (scope-mol (aref scope-molecules locus-mol-index))
                             (scope-pos (aref scope-positions locus-mol-index)))
                        (fformat t "--- locus ~s~%" bp-locus)
                        (loop for bp-monomer across (monomers bp-locus)
                              do (fformat t "----  bp-monomer: ~s~%" bp-monomer)
                              ;; ONE force-field parameterization PER MONOMER, not per rotamer.
                              ;;
                              ;; Every rotamer of a monomer is the same molecule in a different
                              ;; conformation: same atoms, same bonds, same types.  Bonded parameters
                              ;; depend on all three and on nothing else - coordinates are not in the
                              ;; terms, only INDEXES INTO the coordinate vector - so the first
                              ;; rotamer's stretch/angle/dihedral terms become any other's by
                              ;; renaming.  These three carry that reference across the slot loop.
                              do (loop with ref-group = nil
                                       with ref-lo3 = 0
                                       with ref-hi3 = 0
                                       with ref-bonded = 0
                                       for slot below (rotamer-slot-count bp-monomer)
                                       for mrkindex = (blueprint-mrkindex bp-locus bp-monomer slot)
                                       for scan = (aref scans mrkindex)
                                       for record = (nth mrkindex records)
                                       ;; Still restored for EVERY slot, including the copied ones:
                                       ;; the nonbond half's 1-2/1-3/1-4 exclusions are read off the
                                       ;; bond graph, so the fan bond has to be on even when no
                                       ;; bonded term is being generated.
                                       do (fformat t "----- slot = ~s~%" slot)
                                       do (restore-fan-out-bond record)
                                          (unwind-protect
                                               (let* ((group (chem:make-energy-component-group))
                                                      (inter-group (chem:make-energy-component-group))
                                                      (map (blueprint-slot-atom-map blueprint mrkindex))
                                                      (mol (blueprint-locus-molecule-index blueprint bp-locus))
                                                      (partner (blueprint-partner-molecule-index blueprint mol))
                                                      (bb (aref backbone-maps mol))
                                                      ;; Set by the intramolecular block below and
                                                      ;; reused for the intermolecular group - the
                                                      ;; slot's range is the same for both, only
                                                      ;; the target set differs.
                                                      (slot-lo3 0) (slot-hi3 0))
                                                  (chem:setf-energy-component-group-name
                                                   group
                                                   (list :blueprint :slot :intramolecular
                                                                    :mrkindex mrkindex
                                                                    :locus (locus bp-locus)
                                                                    :monomer (monomer-name bp-monomer)
                                                                    :slot slot))
                                                  (chem:setf-energy-component-group-name
                                                   inter-group
                                                   (list :blueprint :slot :intermolecular
                                                                    :mrkindex mrkindex
                                                                    :locus (locus bp-locus)
                                                                    :monomer (monomer-name bp-monomer)
                                                                    :slot slot))
                                                 (multiple-value-bind (lo3 hi3)
                                                     (blueprint-slot-i3-range atom-table map)
                                                   (let ((factory (make-detached-slot-factory
                                                                   map bb
                                                                   :slot-lo3 lo3 :slot-hi3 hi3
                                                                   :backbone-bits
                                                                   (aref backbone-bits mol))))
                                                     (if (null ref-group)
                                                         ;; ---- reference rotamer: parameterize ----
                                                         ;;
                                                         ;; :SCOPE-AGGREGATE - swap THIS slot's
                                                         ;; residue into the reserved position and
                                                         ;; parameterize against backbone-plus-it,
                                                         ;; instead of against an aggregate holding
                                                         ;; every rotamer of every slot.  Bounds the
                                                         ;; Loop walk, ring perception, SMIRKS
                                                         ;; matching and the 1-4 walk at once.
                                                         (progn
                                                           (fformat t "------   generate-into-group ref-group NIL~%")
                                                           (chem:set-content-at
                                                            scope-mol scope-pos
                                                            (slot-residue-of scan))
                                                           (chem:energy-function/generate-into-group
                                                            ef factory group
                                                            :scope-aggregate scope)
                                                           (setf ref-group group
                                                                 ref-lo3 lo3
                                                                 ref-hi3 hi3
                                                                 ref-bonded
                                                                 (group-bonded-term-count group))
                                                           (incf reparameterized))
                                                         ;; ---- every later rotamer: rename ----
                                                         (progn
                                                           ;; The shift is only well defined if the
                                                           ;; two slots hold the same atoms in the
                                                           ;; same order.  Same monomer means they
                                                           ;; do - but "means they do" is how the fan
                                                           ;; misplacement survived, so check.
                                                           (unless (= (- hi3 lo3)
                                                                      (- ref-hi3 ref-lo3))
                                                             (error "mrkindex ~d spans ~d i3 slots ~
                                                                     but its monomer's reference ~
                                                                     spans ~d - rotamers of one ~
                                                                     monomer must have identical ~
                                                                     atoms"
                                                                    mrkindex (- hi3 lo3)
                                                                    (- ref-hi3 ref-lo3)))
                                                           ;; :BONDED NIL - and this is also what
                                                           ;; pushes GROUP onto the energy function.
                                                           ;; Remapping alone never would, and a
                                                           ;; group absent from _EnergyComponents is
                                                           ;; never visited by evaluateAll: it would
                                                           ;; score 0.0 rather than error.
                                                           ;; Scoped too: :BONDED NIL still runs
                                                           ;; construct14InteractionTerms, which
                                                           ;; walks the matter for 1-4 pairs.
                                                           (fformat t "------  !! generate-into-group  when ref-group = ~s  bonded: NIL~%" ref-group)
                                                           (chem:set-content-at
                                                            scope-mol scope-pos
                                                            (slot-residue-of scan))
                                                           (chem:energy-function/generate-into-group
                                                            ef factory group
                                                            :bonded nil
                                                            :scope-aggregate scope)
                                                           (fformat t "------  doing copy-remapped-bonded-terms~%")
                                                           (let ((copied
                                                                   (chem:copy-remapped-bonded-terms
                                                                    ef group ref-group
                                                                    ref-lo3 ref-hi3
                                                                    (- lo3 ref-lo3))))
                                                             (unless (= copied ref-bonded)
                                                               (error "mrkindex ~d copied ~d bonded ~
                                                                       terms from a reference ~
                                                                       holding ~d"
                                                                      mrkindex copied ref-bonded)))
                                                           (incf remapped))))
                                                   ;; Query FROM this slot.  Targets are ALL
                                                   ;; backbone, not just this molecule's - see
                                                   ;; ALL-BACKBONE-BITS.  The range form also
                                                   ;; brings the intra-slot pairs.
                                                   (set-group-query-range
                                                    group lo3 hi3 all-backbone-bits)
                                                   (setf slot-lo3 lo3 slot-hi3 hi3))
                                                 (setf (intramolecular-energy-components scan)
                                                       group)
                                                 ;; ---- intermolecular: this slot vs the OTHER
                                                 ;; molecule's backbone.  SET-MATTERS is what makes
                                                 ;; it affordable - without it this component would
                                                 ;; grid the whole atom table like the
                                                 ;; intramolecular one still does.
                                                 ;; :BONDED NIL - nothing bonded spans two
                                                 ;; molecules, so the per-molecule pass here was
                                                 ;; re-parameterizing every molecule 470 times to
                                                 ;; produce zero bonded terms.
                                                 (chem:energy-function/generate-into-group
                                                  ef
                                                  (make-intermolecular-slot-factory
                                                   map (aref backbone-maps partner))
                                                  inter-group :bonded nil)
                                                 ;; Query set, NOT set-matters - they are mutually
                                                 ;; exclusive (rebuildPairListImpl errors if both
                                                 ;; are set) and the query set is the faster path.
                                                 ;; Same slot range, but targets are the PARTNER
                                                 ;; molecule's backbone.
                                                 (set-group-query-range
                                                  inter-group slot-lo3 slot-hi3
                                                  all-backbone-bits)
                                                 (setf (intermolecular-energy-components scan)
                                                       inter-group)
                                                 (incf built)
                                                 (let ((terms (chem:number-of-terms group)))
                                                   (incf locus-terms terms)
                                                   (when (zerop terms) (incf empty))))
                                            ;; detach again even if generation threw, or the next
                                            ;; slot would be generated with two sidechains on.
                                            (break-fan-out-bond record))))
                        (when verbose
                          (format t "~&  locus ~3d  ~5d terms over its slots~%"
                                  (locus bp-locus) locus-terms)
                          (finish-output)))))
      (mapc #'restore-fan-out-bond records))
    (format t "~&detach scan: ~d slot groups in ~,1f s~@[  *** ~d EMPTY ***~]~%~
               ~&  bonded: ~d parameterized, ~d copy-remapped (~,1fx fewer passes)~%"
            built (/ (- (get-internal-real-time) start) internal-time-units-per-second)
            (when (plusp empty) empty)
            reparameterized remapped
            (if (plusp reparameterized)
                (/ (+ reparameterized remapped) (float reparameterized))
                1.0))
    (when (or (null backbone-group) (null backbone-inter-group) (null pair-group))
      (error "There are missing groups"))
    (values built backbone-group backbone-inter-group pair-group)))


;;; ------------------------------------------------------------------
;;; Is the copy-remap equivalent to a real parameterization?
;;;
;;; COPY-REMAPPED-BONDED-TERMS renames one rotamer's bonded terms into another's atoms instead of
;;; re-running the force field.  That is exact ONLY IF rotamers of a monomer really are the same
;;; molecule in different conformations - same atoms, same order, same types.  They are, by
;;; construction, but a wrong remap does not error: it produces a full set of plausible terms with
;;; some indexes pointing at the wrong atoms, and the energy that comes out is merely wrong.
;;;
;;; So: rebuild a sample of slots the expensive way and diff.
;;; ------------------------------------------------------------------

(defun group-bonded-term-tuples (group)
  "GROUP's bonded terms as a SORTED list of numeric tuples - parameters and i3 indexes, no atoms.

SORTED, not in term order.  Improper ordering comes out of FFItorDb_O::improperAtomSort and there is
no guarantee two builds enumerate them in the same sequence; the SET is what has to match, and a
positional diff would report differences that do not exist.

Numeric only.  The atoms are recoverable from the indexes, so including them would only add a second
way to say the same thing - and make a mismatch harder to read."
  (let ((tuples '()))
    (dolist (c (chem:energy-component-group/components group))
      (typecase c
        (chem:energy-stretch
         (dotimes (i (chem:number-of-terms c))
           (multiple-value-bind (kb r0 i1 i2) (chem:safe-amber-energy-stretch-term c i)
             (push (list :stretch i1 i2 kb r0) tuples))))
        (chem:energy-angle
         (dotimes (i (chem:number-of-terms c))
           (multiple-value-bind (kt t0 i1 i2 i3) (chem:safe-amber-energy-angle-term c i)
             (push (list :angle i1 i2 i3 kt t0) tuples))))
        ;; A linear angle carries only KT - there is no equilibrium angle to store, which is what
        ;; makes it linear - so its tuple is one shorter than an ordinary angle's.
        (chem:energy-linear-angle
         (dotimes (i (chem:number-of-terms c))
           (multiple-value-bind (kt i1 i2 i3) (chem:safe-amber-energy-linear-angle-term c i)
             (push (list :linear-angle i1 i2 i3 kt) tuples))))
        (chem:energy-dihedral
         (dotimes (i (chem:number-of-terms c))
           (multiple-value-bind (cosphase sinphase v n i1 i2 i3 i4)
               (chem:safe-amber-energy-dihedral-term c i)
             (push (list :dihedral i1 i2 i3 i4 v n cosphase sinphase) tuples))))
        ;; INDEXES ONLY, no parameters - EnergyNonbond14 has no safe-amber-...-term accessor and
        ;; CALL-FOR-EACH-TERM reports the atoms and i3 pair and nothing else.  That is the weaker
        ;; half of the check, and it is the half that matters here: the remap copies parameters
        ;; verbatim by struct assignment, so the failure mode it can actually have is a wrong INDEX.
        (chem:energy-nonbond14
         (chem:call-for-each-term
          c (lambda (a1 a2 i1 i2)
              (declare (ignore a1 a2))
              (push (list :nonbond14 i1 i2) tuples))))))
    (sort tuples (lambda (a b)
                   (loop for x in a for y in b
                         do (cond ((and (symbolp x) (symbolp y))
                                   (unless (eq x y)
                                     (return (string< (symbol-name x) (symbol-name y)))))
                                  ((< x y) (return t))
                                  ((> x y) (return nil)))
                         finally (return nil))))))

(defun check-remapped-bonded-terms (blueprint mrkindexes &key (stream *standard-output*))
  "Rebuild each MRKINDEX's bonded terms by real parameterization and diff against what the detach
scan left on the slot.  Returns T when every sampled slot agrees.

Pick indexes that were COPY-REMAPPED, not the first rotamer of a monomer - the first is the
reference and was parameterized for real, so it agrees with itself and proves nothing.

DESTRUCTIVE ONLY TRANSIENTLY: all fan-out bonds are broken for the duration and restored on the way
out, including on a non-local exit."
  (let* ((ef (energy-function blueprint))
         (records (blueprint-fan-out-bonds blueprint))
         (scans (rotamer-scans blueprint))
         (atom-table (chem:atom-table ef))
         (backbone-maps (blueprint-backbone-atom-maps blueprint))
         (backbone-bits (map 'vector
                             (lambda (m) (blueprint-backbone-i3-bits atom-table m))
                             backbone-maps))
         ;; The check must parameterize the SAME WAY the scan did, scope aggregate included - a
         ;; difference in what was walked would show up as a term difference and be blamed on the
         ;; remap.
         (slot-residue-set (blueprint-slot-residue-set blueprint))
         (scope-aggregates (make-array (length backbone-maps)))
         (scope-molecules (make-array (length backbone-maps)))
         (scope-positions (make-array (length backbone-maps)))
         (all-ok t))
    (dotimes (m (length backbone-maps))
      (multiple-value-bind (agg mol pos)
          (make-blueprint-scope-aggregate blueprint m slot-residue-set)
        (setf (aref scope-aggregates m) agg
              (aref scope-molecules m) mol
              (aref scope-positions m) pos)))
    (unwind-protect
         (progn
           (mapc #'break-fan-out-bond records)
           (format stream "~&remap check - real parameterization vs copy-remap~%")
           (dolist (mrkindex mrkindexes)
             (let* ((scan (aref scans mrkindex))
                    (have (and scan (intramolecular-energy-components scan))))
               (if (null have)
                   (format stream "~&  mrkindex ~5d  no group - SKIPPED~%" mrkindex)
                   (let* ((bp-locus (blueprint-locus-of-mrkindex blueprint mrkindex))
                          (mol (blueprint-locus-molecule-index blueprint bp-locus))
                          (map (blueprint-slot-atom-map blueprint mrkindex))
                          (record (nth mrkindex records)))
                     (restore-fan-out-bond record)
                     (unwind-protect
                          (multiple-value-bind (lo3 hi3)
                              (blueprint-slot-i3-range atom-table map)
                            (let ((fresh (chem:make-energy-component-group)))
                              (chem:set-content-at (aref scope-molecules mol)
                                                   (aref scope-positions mol)
                                                   (slot-residue-of scan))
                              (chem:energy-function/generate-into-group
                               ef
                               (make-detached-slot-factory
                                map (aref backbone-maps mol)
                                :slot-lo3 lo3 :slot-hi3 hi3
                                :backbone-bits (aref backbone-bits mol))
                               fresh
                               :scope-aggregate (aref scope-aggregates mol))
                              (let* ((a (group-bonded-term-tuples fresh))
                                     (b (group-bonded-term-tuples have))
                                     (ok (equal a b)))
                                (unless ok (setf all-ok nil))
                                (format stream "~&  mrkindex ~5d  real ~5d  remapped ~5d  ~a~%"
                                        mrkindex (length a) (length b)
                                        (if ok
                                            "AGREE"
                                            (format nil "*** DIFFER - ~d only in real, ~d only in remapped ***"
                                                    (length (set-difference a b :test #'equal))
                                                    (length (set-difference b a :test #'equal)))))
                                ;; The counts alone cannot tell a DROPPED term from a term that came
                                ;; out naming different atoms - and those have completely different
                                ;; causes.  Print the tuples, and print the offending indexes
                                ;; relative to the slot range, since a remap bug shows up as an
                                ;; index that should have shifted and did not (or the reverse).
                                (unless ok
                                  (let ((only-real (set-difference a b :test #'equal))
                                        (only-remapped (set-difference b a :test #'equal)))
                                    (format stream "~&      slot i3 range [~d,~d)~%" lo3 hi3)
                                    (flet ((show (label tuples)
                                             (dolist (tu (subseq tuples 0 (min 6 (length tuples))))
                                               (format stream "~&      ~a ~a~{ ~a~}~%"
                                                       label (first tu)
                                                       (mapcar (lambda (x)
                                                                 (if (and (integerp x)
                                                                          (>= x 0)
                                                                          (< x 1000000))
                                                                     (format nil "~d~a" x
                                                                             (cond ((and (>= x lo3)
                                                                                         (< x hi3))
                                                                                    "[slot]")
                                                                                   (t "[bb]")))
                                                                     (format nil "~a" x)))
                                                               (rest tu))))))
                                      (show "only-real    " only-real)
                                      (show "only-remapped" only-remapped)))))
                              ;; FRESH was pushed onto the energy function by
                              ;; GENERATE-INTO-GROUP.  Leaving it there would add a duplicate of
                              ;; this slot to every later whole-function walk.
                              (chem:remove-energy-component ef fresh)))
                       (break-fan-out-bond record))))))
           (format stream "~&remap check: ~a~%" (if all-ok "PASS" "*** FAIL ***")))
      (mapc #'restore-fan-out-bond records))
    all-ok))


;;; ------------------------------------------------------------------
;;; Training the SMIRNOFF parameter cache for THIS blueprint
;;;
;;; The blueprint must never teach the cache.  SMIRNOFF types and parameterizes by SMIRKS
;;; matching against the local graph, and the fan-out leaves a CA carrying 58-79 substituents -
;;; not a chemical environment, so what real SMIRNOFF infers there is not what any actual
;;; molecule has.  Harvesting it would be inherited by every later build SILENTLY, because a
;;; covered molecule never re-runs SMIRNOFF to notice the disagreement.
;;;
;;; Instead: train on small chemically-real trainers first, then the blueprint hits the cache
;;; and SMIRKS never sees the fanned graph at all.  The cache key is
;;; (atom-name . constitution-context) - nothing in it mentions molecule, oligomer or pose - so
;;; a trainer's parameters transfer verbatim.
;;;
;;; Enumerating trainers is cheap (MAKE-OLIGOMER + FOLDAMER-MONOMER-CONTEXT); only building the
;;; assembler is expensive.  So the whole index is computed up front and only the missing
;;; contexts are actually built.
;;; ------------------------------------------------------------------

(defun blueprint-required-contexts-by-foldamer (blueprint)
  "FOLDAMER -> the constitution-contexts this BLUEPRINT needs from it, deduplicated.

  Grouped by foldamer rather than pooled, because only SOME foldamers have a cache to train.  The
  ligand is a spiroligomer on cached SMIRNOFF; the receptor's mobile sidechains are amino acids on
  AMBER, which types from residue templates and needs no parameter cache at all.  Pooling them
  reported five AA contexts (:LYN_-SIDE_AA, :ARG_-SIDE_AA, ...) as untrainable when in truth they
  were never anything a SMIRNOFF trainer would produce."
  (let ((table (make-hash-table :test 'eq)))
    (loop for bp-locus across (loci blueprint)
          when (moveable-p bp-locus)
            do (let* ((olig (oligomer (oligomer-shape bp-locus)))
                      (fold (foldamer (oligomer-space olig))))
                 (flet ((need (cc)
                          (when cc
                            (pushnew cc (gethash fold table) :test #'equal))))
                   (if (owns-slots-p bp-locus)
                       (loop for bp-monomer across (monomers bp-locus)
                             do (need (blueprint-monomer-constitution-context
                                       bp-monomer bp-locus)))
                       (multiple-value-bind (monomer-context cc)
                           (foldamer-monomer-context (original-monomer bp-locus) olig fold)
                         (declare (ignore monomer-context))
                         (need cc))))))
    table))


(defgeneric train-foldamer-contexts (foldamer contexts &key verbose)
  (:documentation
   "Train FOLDAMER's parameter cache on CONTEXTS.  Returns (values TRAINED MISSING), or
(values 0 NIL) when this foldamer has no cache to train.

  A foldamer whose force field is not a CACHED-SMIRNOFF-FORCE-FIELD is skipped, not an error:
  AMBER types from residue templates, so there is nothing a trainer could contribute.

  The name need not be a registered force field at all - AMBER-PROTEIN's is the cons
  (:DEFAULT . :USE-GIVEN-TYPES), and CHEM:FIND-FORCE-FIELD SIGNALS on an unknown name rather than
  returning NIL.  An unregistered name is just a louder way of saying there is no cache here."))
  

(defun blueprint-required-constitution-contexts (blueprint)
  "Every constitution-context this BLUEPRINT will stamp as :GIVEN-ATOM-TYPE, pooled across
foldamers.  Diagnostic only - training goes foldamer by foldamer, see
BLUEPRINT-REQUIRED-CONTEXTS-BY-FOLDAMER.

  :SIDECHAIN loci contribute one per blueprint-monomer - every rotamer alternative is materialized,
  so all of them are needed.  :BACKBONE loci contribute the context of their CURRENT selection,
  which is what MAKE-ASSEMBLER's ASSIGN-GIVEN-ATOM-TYPES will stamp (assembler.lisp:743).

  :FIXED loci are skipped: they are the receptor, whose foldamer's ASSIGN-GIVEN-ATOM-TYPES is a
  no-op, so they contribute no cache keys.  If that ever stops being true the coverage check
  catches it - an unstamped atom fails BONDED-CACHE-COVERS-MOLECULE-P, which is loud."
  (let ((contexts nil))
    (loop for bp-locus across (loci blueprint)
          when (moveable-p bp-locus)
            do (let* ((olig (oligomer (oligomer-shape bp-locus)))
                      (fold (foldamer (oligomer-space olig))))
                 (if (owns-slots-p bp-locus)
                     (loop for bp-monomer across (monomers bp-locus)
                           for cc = (blueprint-monomer-constitution-context bp-monomer bp-locus)
                           when cc do (pushnew cc contexts :test #'equal))
                     (multiple-value-bind (monomer-context cc)
                         (foldamer-monomer-context (original-monomer bp-locus) olig fold)
                       (declare (ignore monomer-context))
                       (when cc (pushnew cc contexts :test #'equal))))))
    (nreverse contexts)))


(defgeneric foldamer-trainer-index (foldamer)
  (:documentation "CONSTITUTION-CONTEXT -> (OLIGOMER . FOCUS-MONOMER) for every trainer FOLDAMER can build.

  Cheap by design: MAKE-OLIGOMER and FOLDAMER-MONOMER-CONTEXT only, no assembler and no SMIRNOFF.
  That is what makes demand-driven training practical - the enumeration costs nothing, so only the
  parameterizations you actually need get paid for."))


(defun ensure-blueprint-trained (blueprint &key verbose)
  "Train the SMIRNOFF parameter cache on exactly the trainers BLUEPRINT needs, and no others.

  CALL THIS BEFORE MATERIALIZE-BLUEPRINT-BASE.  An untrained cache learns from whatever is built
  first, and for a blueprint that is the fanned-out aggregate - see the commentary above.

  Errors listing any required context that no trainer produces.  That is the useful failure: it
  means the foldamer's training oligomer spaces do not cover this blueprint's design space, and no
  amount of building will fix it.

  Returns (values TRAINED-COUNT REQUIRED-COUNT)."
  (let ((by-foldamer (blueprint-required-contexts-by-foldamer blueprint))
        (trained 0)
        (required 0)
        (missing nil))
    (maphash (lambda (foldamer contexts)
               (incf required (length contexts))
               (multiple-value-bind (n gaps)
                   (train-foldamer-contexts foldamer contexts :verbose verbose)
                 (incf trained n)
                 (when gaps (push (cons foldamer gaps) missing))))
             by-foldamer)
    (when missing
      (error "No trainer produces these constitution-contexts the blueprint requires:~%~
              ~{~{  ~a:~%~{    ~s~%~}~}~}~
              That foldamer's TRAINING-OLIGOMER-SPACES do not cover this blueprint's design ~
              space.  Without them the cached force field cannot cover the blueprint aggregate ~
              and would fall back to real SMIRNOFF over the fanned-out graph."
             (mapcar (lambda (entry)
                       (list (type-of (car entry)) (cdr entry)))
                     missing)))
    (format t "~&ensure-blueprint-trained: ~d contexts required across ~d foldamer~:p, ~
               ~d newly trained~%"
            required (hash-table-count by-foldamer) trained)
    (values trained required)))
