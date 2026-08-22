(in-package :design)

;;; ==================================================================
;;; STAGE 1 - do the blueprint and the pose hold the SAME BONDED TERMS?
;;;
;;; The blueprint scores a structure by summing buckets:
;;;
;;;     backbone  +  the selected slots' singles  +  the selected pairs
;;;
;;; and that is supposed to partition exactly the interactions the pose's own energy function holds.
;;; CHECK-TOTAL-ENERGY compares the SUMS; when they disagree it cannot say whether the blueprint is
;;; missing terms, holding extra ones, or holding the right ones at the wrong coordinates.  This
;;; compares the term SETS and answers the first two.
;;;
;;; BONDED ONLY, AND THAT IS THE POINT.  A stretch/angle/dihedral/1-4 term exists because of the
;;; BOND GRAPH; where the atoms sit is irrelevant.  So the two sets must match EXACTLY even if every
;;; coordinate differs, which makes this the one comparison that isolates topology from geometry.
;;; The nonbond components cannot be compared this way at all: a rosetta-nonbond term exists only if
;;; the pair falls within RPAIRLIST, so its term set is a function of the coordinates, and two
;;; different folds legitimately produce different sets.  Comparing those before the coordinates are
;;; known to agree reports differences that are not bugs.
;;;
;;; Hence the order: bonded terms (here) -> coordinates -> nonbond terms.  If the coordinates turn
;;; out to match exactly, the nonbond sets must match too and the third stage is a formality.
;;;
;;; THE PAIR BUCKET CONTRIBUTES NOTHING HERE.  Pair-scan groups are generated :BONDED NIL - two
;;; rotamers at different loci are not bonded to each other - so stage 1 is backbone plus the
;;; selected singles, against the pose's bonded components.
;;; ==================================================================

;;; ------------------------------------------------------------------
;;; The blueprint -> pose atom-index map
;;;
;;; There is no shared atom identity to exploit.  MATERIALIZE-BLUEPRINT-BASE calls MAKE-ASSEMBLER
;;; afresh (blueprint.lisp:1666), so the blueprint's aggregate, residues, atoms and atom table are
;;; all new objects; only the OLIGOMER-SHAPES and the monomers inside them are shared with the pose's
;;; assembler.  The correspondence therefore has to be constructed, by name, within paired residues.
;;;
;;; It maps INDEX SPACES rather than objects - blueprint i3 -> pose i3 - because everything
;;; downstream is already expressed in i3: DUMP-COMPONENT-TERMS keys on sorted i3 lists, and a
;;; coordinate check is COORDS[bp-i3] against POSE-COORDS[pose-i3].  Translating the key is then the
;;; whole of the cross-structure machinery.
;;;
;;; BUILT AT EXTRACTION TIME, NEVER CACHED ON THE BLUEPRINT.  It is a function of the blueprint AND
;;; the pose AND the mcstate, which have three different lifetimes; a copy held on the blueprint
;;; would go stale exactly the way BLUEPRINT-MCSTATE-FOR-POSE's aliasing did, and just as silently.
;;; ------------------------------------------------------------------

(defun pose-monomer-index-to-residue (assembler)
  "Hash (SHAPE-INDEX . MONOMER-INDEX) -> RESIDUE over ASSEMBLER's ataggregate.

MONOMER-INDEX is documented on ATRESIDUE as \"only for debugging\" (fold.lisp:181), so this verifies
rather than trusts: a duplicate key means two atresidues claim the same monomer, and the map would
silently keep one - which reads downstream as a whole residue's worth of missing terms."
  (let ((map (make-hash-table :test 'equal))
        (atmolecules (topology::atmolecules (topology::ataggregate assembler))))
    (dotimes (shape-index (length atmolecules))
      (loop for atres across (topology::atresidues (aref atmolecules shape-index))
            when atres
              do (let ((key (cons shape-index (topology::monomer-index atres))))
                   (when (gethash key map)
                     (error "Two atresidues in shape ~d claim monomer-index ~d - MONOMER-INDEX ~
                             cannot identify residues in this assembler."
                            shape-index (topology::monomer-index atres)))
                   (setf (gethash key map) (topology::residue atres)))))
    map))

(defun pair-residue-atoms (bp-residue pose-residue bp-table pose-table map names residue-key what)
  "Fill MAP: for each atom of BP-RESIDUE, MAP[bp-i3] = the pose i3 of the same-named atom.
Also fill NAMES: NAMES[pose-i3] = a printable (RESIDUE-KEY . ATOM-NAME).

Atom NAMES are the correspondence.  They are unique within a residue, which is what lets one atom in
the blueprint be identified with one in the pose when the two structures share no objects at all.

NAMES exists so a reported difference can be READ.  A term key is a list of i3 integers, and a diff
that prints (ENERGY-ANGLE 1023 1026 1029) localizes nothing - acting on it would mean translating
three indices through the atom table by hand.  Built here because this is the one place that has an
atom and the residue it belongs to at the same moment; atoms carry no back-pointer to their residue.

Signals on a name that is missing or duplicated on either side rather than skipping it: an unpaired
atom means the two residues are not the same chemistry, and every term touching it would then be
reported as a difference when the real fault is here."
  (let ((pose-by-name (make-hash-table :test 'equal)))
    (chem:do-atoms (a pose-residue)
      (let ((n (string (chem:get-name a))))
        (when (gethash n pose-by-name)
          (error "~a: pose residue ~a has two atoms named ~a" what pose-residue n))
        (setf (gethash n pose-by-name) a)))
    (let ((seen (make-hash-table :test 'equal)))
      (chem:do-atoms (a bp-residue)
        (let* ((n (string (chem:get-name a)))
               (pose-atom (gethash n pose-by-name)))
          (when (gethash n seen)
            (error "~a: blueprint residue ~a has two atoms named ~a" what bp-residue n))
          (setf (gethash n seen) t)
          (unless pose-atom
            (error "~a: blueprint atom ~a has no counterpart in the pose residue ~a.  The two ~
                    residues are not the same chemistry, so no term comparison between them means ~
                    anything."
                   what n pose-residue))
          (let ((pose-i3 (chem:get-coordinate-index-times3 pose-table pose-atom)))
            (setf (aref map (topology::blueprint-atom-i3 bp-table a)) pose-i3)
            (setf (gethash pose-i3 names) (cons residue-key n)))))
      ;; Both directions: a pose atom with no blueprint counterpart is equally fatal, and only this
      ;; check catches it - the loop above would simply never visit it.
      (maphash (lambda (n a)
                 (declare (ignore a))
                 (unless (gethash n seen)
                   (error "~a: pose atom ~a has no counterpart in the blueprint residue ~a"
                          what n bp-residue)))
               pose-by-name))))

(defun blueprint-to-pose-atom-index-map (blueprint pose mcstate)
  "Return (values MAP NAMES).

MAP is a vector indexed by the BLUEPRINT's i3 -> the POSE's i3, or NIL where there is no counterpart.
NAMES is a hash POSE-I3 -> (RESIDUE-KEY . ATOM-NAME), so a term key can be printed as chemistry
rather than as three integers.

NIL IS THE INFORMATIVE CASE.  Every rotamer of a monomer carries identical atom NAMES, so a purely
name-based map would happily translate an UNSELECTED rotamer's CB to the pose's CB and produce a
valid-looking, wrong answer.  Only the rotamers MCSTATE selected are mapped; the other ~760 slots'
atoms stay NIL, so a term reaching one of them translates to NIL and is reported as a term that
should never have been in a selected bucket.

Length is the blueprint's full coordinate-vector size, indexed by i3 directly - the same
three-times-larger, no-division-needed convention the scope bit-vectors use."
  (let* ((assembler (design:assembler pose))
         (bp-table (chem:atom-table (topology:energy-function blueprint)))
         (pose-table (chem:atom-table (topology:energy-function assembler)))
         (pose-residues (pose-monomer-index-to-residue assembler))
         (slot-residues (topology::blueprint-slot-residue-set blueprint))
         (bp-agg (topology::aggregate blueprint))
         (pose-agg (topology:aggregate assembler))
         (map (make-array (* 3 (chem:get-number-of-atoms bp-table)) :initial-element nil))
         (names (make-hash-table)))
    ;; ---- the ONE selected rotamer at each locus, FIRST ----
    ;;
    ;; Before the positional walk, because it decides which of the pose's residues that walk must
    ;; SKIP.  MATERIALIZE-BLUEPRINT-BASE declines to build a base residue for exactly the loci that
    ;; own mrkindexes, so the blueprint represents those loci ONLY as slot residues - while the pose
    ;; carries one ordinary residue at every locus.  A ligand with three slot-owning loci therefore
    ;; has 5 non-slot residues in the blueprint against 8 in the pose, and pairing the two sequences
    ;; positionally would misalign every residue after the first slot-owning one.
    (let ((claimed (make-hash-table :test 'eq)))
      (loop for mrkindex across mcstate
            for info = (topology::mrkindex-info blueprint mrkindex)
            for pose-res = (gethash (cons (getf info :shape-index) (getf info :locus-in-shape))
                                    pose-residues)
            do (unless pose-res
                 (error "No pose residue for shape ~s monomer-index ~s (mrkindex ~d, locus ~s)"
                        (getf info :shape-index) (getf info :locus-in-shape)
                        mrkindex (getf info :blueprint-locus)))
               (setf (gethash pose-res claimed) t)
               (pair-residue-atoms
                (topology::residue (aref (topology::rotamer-scans blueprint) mrkindex))
                pose-res bp-table pose-table map names
                (list :locus (getf info :blueprint-locus) (getf info :monomer-name))
                (format nil "locus ~s slot mrkindex ~d" (getf info :blueprint-locus) mrkindex)))
      ;; ---- everything else, paired POSITIONALLY ----
      ;;
      ;; Both aggregates are built by MAKE-ASSEMBLER from the same oligomer-shapes, so once each
      ;; side's slot-owning loci are removed - slot residues on the blueprint, CLAIMED residues on
      ;; the pose - the two sequences describe the same monomers in the same order.  Positional
      ;; rather than by monomer-index because not every fixed residue has an atresidue to read one
      ;; from.
      (dotimes (m (chem:content-size pose-agg))
        (let* ((bp-mol (chem:content-at bp-agg m))
               (pose-mol (chem:content-at pose-agg m))
               (bp-base (loop for i below (chem:content-size bp-mol)
                              for r = (chem:content-at bp-mol i)
                              unless (gethash r slot-residues) collect r))
               (pose-base (loop for i below (chem:content-size pose-mol)
                                for r = (chem:content-at pose-mol i)
                                unless (gethash r claimed) collect r)))
          (unless (= (length bp-base) (length pose-base))
            (error "Molecule ~d has ~d unclaimed residues in the blueprint and ~d in the pose.~%~
                    After removing the blueprint's slot residues and the pose's residues at ~
                    slot-owning loci, the two sequences should describe the same monomers in the ~
                    same order.  They do not, so no pairing below would be meaningful."
                   m (length bp-base) (length pose-base)))
          (loop for bp-res in bp-base
                for pose-res in pose-base
                for n from 0
                do (pair-residue-atoms bp-res pose-res bp-table pose-table map names
                                       (list :backbone m n)
                                       (format nil "molecule ~d backbone residue ~d" m n))))))
    (values map names)))

;;; ------------------------------------------------------------------
;;; Bonded term sets
;;; ------------------------------------------------------------------

(defparameter +bonded-component-classes+
  '(chem:energy-stretch chem:energy-angle chem:energy-linear-angle
    chem:energy-dihedral chem:energy-nonbond14)
  "The classes whose term sets are pure topology, and so must match across two structures whatever
their coordinates.  The nonbond components are deliberately absent - see the file header.

ENERGY-LINEAR-ANGLE is here for a reason worth knowing: the REAL SMIRNOFF path creates linear angles
(smirnoff.lisp:255) and COMMIT-BONDED-FROM-CACHE creates NONE.  So a molecule with a linear group
gets a different bonded term set depending on whether the parameter cache hit - and if the pose took
one path and the blueprint the other, the difference is exactly these terms.  Leaving the class out
of this list would hide precisely the discrepancy the comparison exists to find.")

(defun bonded-component-p (component)
  (some (lambda (c) (typep component c)) +bonded-component-classes+))

(defun collect-bonded-terms (components table &key translate label)
  "Add every bonded term of COMPONENTS to TABLE, keyed by (CLASS . sorted i3 list).

TRANSLATE, when given, maps each i3 through it before keying - that is how the blueprint's indexes
are expressed in the pose's space so the two sets are comparable.  A term with any untranslatable
i3 is counted as UNMAPPED rather than keyed, because on the blueprint side that means the term
reaches a rotamer this mcstate did not select: a real defect, and one silent translation would hide.

Returns (values TABLE UNMAPPED)."
  (let ((unmapped 0))
    (dolist (component components)
      (when (bonded-component-p component)
        (let ((class (class-name (class-of component))))
          (chem:atoms-for-each-term
           component
           (lambda (&rest args)
             ;; ATOMS first, then i3 - the convention EnergyComponent_O::atomsForEachTerm documents.
             ;; Taking the trailing fixnums makes one reader work for every bonded class.
             (let ((i3s (remove-if-not #'integerp args)))
               (let ((mapped (if translate (mapcar translate i3s) i3s)))
                 (if (some #'null mapped)
                     (incf unmapped)
                     (push (or label t)
                           (gethash (cons class (sort mapped #'<)) table nil))))))))))
    (values table unmapped)))

(defun blueprint-selected-bonded-terms (blueprint mcstate backbone-group map)
  "The blueprint's bonded term set for MCSTATE, in POSE index space.

BACKBONE-GROUP plus each selected slot's intramolecular and intermolecular groups.  The pair-scan
groups are not consulted: they are generated :BONDED NIL and hold no bonded terms at all.

Returns (values TABLE UNMAPPED)."
  (let ((table (make-hash-table :test 'equal))
        (total-unmapped 0)
        (translate (lambda (i3) (aref map i3))))
    (flet ((add (group label)
             (when group
               (multiple-value-bind (tbl n)
                   (collect-bonded-terms (chem:energy-component-group/components group)
                                         table :translate translate :label label)
                 (declare (ignore tbl))
                 (incf total-unmapped n)))))
      (add backbone-group :backbone)
      (loop for mrkindex across mcstate
            for scan = (aref (topology::rotamer-scans blueprint) mrkindex)
            do (add (topology::intramolecular-energy-components scan) (list :slot mrkindex))
               (add (topology::intermolecular-energy-components scan) (list :slot-inter mrkindex))))
    (values table total-unmapped)))

(defun pose-bonded-terms (pose)
  "The pose energy function's bonded term set, in its own index space.

Also returns the class names of the components this did NOT read.  Nonbond components are expected
there - they are stage 3, not stage 1 - but so is anything else the pose's factory installed that a
blueprint group has no counterpart for, RESTRAINTS above all: GENERATE-INTO-GROUP never calls
GENERATE-RESTRAINT-ENERGY-FUNCTION-TABLES, so a restraint's energy is in the pose's total with
nothing on the blueprint side to match it.

Reported by CLASS rather than filtered by RESTRAINTP, which is pure-virtual with no Lisp binding -
and reporting every unread class is the more useful answer anyway, since it names whatever is
actually there instead of only what was anticipated."
  (let* ((ef (topology:energy-function (design:assembler pose)))
         (components (chem:all-components ef))
         (unread (remove-duplicates
                  (mapcar (lambda (c) (class-name (class-of c)))
                          (remove-if #'bonded-component-p components)))))
    (multiple-value-bind (table unmapped)
        (collect-bonded-terms components (make-hash-table :test 'equal) :label :pose)
      (declare (ignore unmapped))
      (values table unread))))

;;; ------------------------------------------------------------------
;;; The stage-1 driver
;;; ------------------------------------------------------------------

(defun check-bonded-terms-against-pose (blueprint pose mcstate backbone-group
                                        &key (stream *standard-output*) (limit 20))
  "STAGE 1: compare the blueprint's bonded term set at MCSTATE against the pose's.  Returns T on
an exact match.

An exact match means the partition is right, the copy-remap is right, and the atom correspondence is
right - which together say the blueprint holds the same MOLECULE the pose does.  It says nothing
about where the atoms are; that is stage 2.

MULTISET, NOT SET.  A key is (CLASS . sorted atom list) and that does NOT identify a term uniquely:
one proper dihedral with three Fourier components produces THREE EnergyDihedral terms over the same
four atoms, and an improper over those same four atoms sorts to the same key again.  Comparing
presence alone would call three terms one - so a blueprint with three components and a pose with two
would agree - and would report every multi-component dihedral as a duplicate.  The count per key is
therefore part of the comparison, and a count mismatch is its own category."
  (multiple-value-bind (map names) (blueprint-to-pose-atom-index-map blueprint pose mcstate)
    (multiple-value-bind (bp-table unmapped)
        (blueprint-selected-bonded-terms blueprint mcstate backbone-group map)
      (multiple-value-bind (pose-table unread-classes) (pose-bonded-terms pose)
        (let ((only-bp nil) (only-pose nil) (count-mismatch nil)
              (bp-terms 0) (pose-terms 0))
          (maphash (lambda (k v)
                     (incf bp-terms (length v))
                     (let ((p (gethash k pose-table)))
                       (cond ((null p) (push k only-bp))
                             ((/= (length v) (length p))
                              (push (list k (length v) (length p)) count-mismatch)))))
                   bp-table)
          (maphash (lambda (k v)
                     (incf pose-terms (length v))
                     (unless (gethash k bp-table) (push k only-pose)))
                   pose-table)
          (format stream "~&~%!!!! BONDED TERM SETS: DOES THE BLUEPRINT HOLD THE SAME MOLECULE ~
                            AS THE POSE?~%~
                          ~&  Bonded terms are pure topology - they exist because of the bond ~
                            graph, not because of where the atoms sit - so these two sets must ~
                            match EXACTLY however differently the~%~
                          ~&  two structures are folded.  A difference here explains any energy ~
                            disagreement below; a match means the gap is geometry.~%~
                          ~&  blueprint ~6d terms over ~6d atom-tuples~%~
                          ~&  pose      ~6d terms over ~6d atom-tuples~%~
                          ~&  only in blueprint ~6d     only in pose ~6d     count mismatch ~6d~%"
                  bp-terms (hash-table-count bp-table)
                  pose-terms (hash-table-count pose-table)
                  (length only-bp) (length only-pose) (length count-mismatch))
          ;; PER CLASS, because the totals cannot show which classes were actually exercised.  A
          ;; class with zero terms on BOTH sides matches vacuously and is indistinguishable in the
          ;; totals from one that matched on thousands of terms - and "were there any linear angles
          ;; at all?" is precisely the question this comparison was extended to answer.
          (let ((by-class (make-hash-table :test 'eq)))
            (flet ((tally (table slot)
                     (maphash (lambda (k v)
                                (let ((row (or (gethash (car k) by-class)
                                               (setf (gethash (car k) by-class)
                                                     (list 0 0 0 0)))))
                                  (incf (nth slot row) (length v))
                                  (incf (nth (1+ slot) row))))
                              table)))
              (tally bp-table 0)     ; slots 0,1 = blueprint terms, tuples
              (tally pose-table 2))  ; slots 2,3 = pose terms, tuples
            (format stream "~&  ~32a ~10@a ~10@a~%" "  class" "blueprint" "pose")
            (dolist (class +bonded-component-classes+)
              (let ((row (gethash class by-class)))
                (format stream "~&  ~50a ~10d ~10d~a~%"
                        class (if row (first row) 0) (if row (third row) 0)
                        (if row "" "   <- NO TERMS ON EITHER SIDE: matched vacuously")))))
          (when (plusp unmapped)
            (format stream "~&  *** ~d blueprint terms reach an atom this mcstate did not select. ~
                              A selected bucket must never name an unselected rotamer. ***~%"
                    unmapped))
          (when unread-classes
            (format stream "~&  pose component classes NOT read by this stage:~%~{    ~a~%~}~
                            ~&  (nonbond is expected - it is stage 3.  A RESTRAINT class here is ~
                               the finding: the blueprint generates none, so its energy is in the ~
                               pose total with nothing to match it.)~%"
                    unread-classes))
          ;; Rendered as chemistry, not as indexes.  A key is (CLASS . sorted pose-i3 list), and
          ;; printing it raw - (ENERGY-ANGLE 1023 1026 1029) - names no atom and no residue, so
          ;; acting on it would mean walking the atom table by hand.  NAMES was built for exactly
          ;; this moment.
          (flet ((render-atom (i3)
                   (let ((n (gethash i3 names)))
                     (if n
                         (format nil "~a:~a" (car n) (cdr n))
                         ;; Only reachable for a pose atom outside every paired residue - i.e. the
                         ;; pose holds a residue the correspondence never covered.
                         (format nil "i3=~d[UNNAMED]" i3)))))
            (flet ((show (label keys)
                     (dolist (k (subseq keys 0 (min limit (length keys))))
                       (format stream "~&    ~a ~a~{  ~a~}~%"
                               label (car k) (mapcar #'render-atom (cdr k))))
                     (when (> (length keys) limit)
                       (format stream "~&    ... ~d more~%" (- (length keys) limit)))))
              (show "only-blueprint" only-bp)
              (show "only-pose     " only-pose))
            ;; A count mismatch on a dihedral is almost always the Fourier expansion differing -
            ;; the same torsion parameterized with a different number of components - which is a
            ;; force-field difference rather than a partition one, and reads quite differently
            ;; from a term the blueprint simply failed to produce.
            (dolist (m (subseq count-mismatch 0 (min limit (length count-mismatch))))
              (destructuring-bind (k bp-n pose-n) m
                (format stream "~&    count      ~a~{  ~a~}   blueprint ~d / pose ~d~%"
                        (car k) (mapcar #'render-atom (cdr k)) bp-n pose-n)))
            (when (> (length count-mismatch) limit)
              (format stream "~&    ... ~d more~%" (- (length count-mismatch) limit))))
          (let ((ok (and (null only-bp) (null only-pose)
                         (null count-mismatch) (zerop unmapped))))
            (format stream "~&bonded terms: ~a~%" (if ok "MATCH" "*** DIFFER ***"))
            ok))))))

;;; ==================================================================
;;; STAGE 2 - are the blueprint's atoms WHERE the pose's atoms are?
;;;
;;; Stage 1 established that the two hold the same MOLECULE - same bonded terms over the same atoms.
;;; It said nothing about placement, and it could not: a bonded term exists because of the bond
;;; graph however far apart the two structures are folded.
;;;
;;; This is the question fa_rep actually turns on.  The nonbond components carry essentially all of
;;; it, a term exists there only if a pair falls within RPAIRLIST, and the energy then varies
;;; steeply with distance - so a sub-Angstrom placement difference moves fa_rep by orders of
;;; magnitude while leaving the bonded comparison completely untouched.
;;;
;;; IT ALSO DISAMBIGUATES THE NONBOND TERM SETS, which stage 1 never compared.  Those cannot be
;;; diffed naively - a pair crossing the cutoff changes the set legitimately - so:
;;;
;;;   coordinates MATCH  -> the nonbond sets must match too, and any remaining fa_rep difference is
;;;                         a genuine enumeration bug in the query sets / neighbour list / partition
;;;   coordinates DIFFER -> geometry is the story, and the nonbond difference is downstream of it
;;;
;;; PER RESIDUE, not just globally.  One locus folded differently is the expected failure, and a
;;; global RMSD over ~7000 atoms would average it into invisibility.
;;; ==================================================================

(defun check-coordinates-against-pose (blueprint pose mcstate bp-coords
                                       &key (stream *standard-output*) (limit 20)
                                            (tolerance 1d-6)
                                            ;; The GENERATED conformation, not (DESIGN:COORDS POSE).
                                            ;; The stored vector predates the rotamer labels its own
                                            ;; assembler now carries, so comparing against it
                                            ;; compares two different structures.
                                            (pose-coords (refold-pose-coordinates pose)))
  "STAGE 2: compare the blueprint's coordinates against the pose's, atom by atom.  Returns T when
every mapped atom agrees to within TOLERANCE.

Uses the same BLUEPRINT-TO-POSE-ATOM-INDEX-MAP stage 1 does, so exactly the atoms that were compared
as terms are compared as positions - no separate correspondence to drift out of step."
  (multiple-value-bind (map names) (blueprint-to-pose-atom-index-map blueprint pose mcstate)
    (let* ((by-residue (make-hash-table :test 'equal))
           (n 0)
           (sum-sq 0d0)
           (worst 0d0)
           (worst-key nil))
      (unless pose-coords
        (error "Pose has no COORDS - nothing to compare against"))
      (loop for bp-i3 from 0 below (length map) by 3
            for pose-i3 = (aref map bp-i3)
            when pose-i3
              do (let* ((dx (- (aref bp-coords bp-i3) (aref pose-coords pose-i3)))
                        (dy (- (aref bp-coords (+ bp-i3 1)) (aref pose-coords (+ pose-i3 1))))
                        (dz (- (aref bp-coords (+ bp-i3 2)) (aref pose-coords (+ pose-i3 2))))
                        (d2 (+ (* dx dx) (* dy dy) (* dz dz)))
                        (key (car (gethash pose-i3 names))))
                   (incf n)
                   (incf sum-sq d2)
                   (when (> d2 worst)
                     (setf worst d2 worst-key (gethash pose-i3 names)))
                   ;; (count sum-sq worst) per residue
                   (let ((row (or (gethash key by-residue)
                                  (setf (gethash key by-residue) (list 0 0d0 0d0)))))
                     (incf (first row))
                     (incf (second row) d2)
                     (when (> d2 (third row)) (setf (third row) d2)))))
      (format stream "~&~%!!!! COORDINATES: ARE THE BLUEPRINT'S ATOMS WHERE THE POSE'S ATOMS ARE?~%~
                      ~&  Stage 1 proved the two hold the same MOLECULE.  This asks whether it is ~
                        in the same PLACE - which is~%~
                      ~&  what fa_rep turns on, and what decides whether the uncompared NONBOND ~
                        term sets could differ.~%")
      (if (zerop n)
          (progn (format stream "~&  *** no atoms mapped - nothing compared ***~%") nil)
          (let ((rmsd (sqrt (/ sum-sq n))))
            (format stream "~&  ~d atoms compared   RMSD ~,6f A   worst ~,6f A~@[ at ~a:~a~]~%"
                    n rmsd (sqrt worst) (car worst-key) (cdr worst-key))
            ;; Per residue, worst first.  A single mis-folded locus is the expected failure and it
            ;; would vanish into a global average over ~7000 atoms.
            (let ((rows nil))
              (maphash (lambda (k row)
                         (push (list k (first row) (sqrt (/ (second row) (first row)))
                                     (sqrt (third row)))
                               rows))
                       by-residue)
              (setf rows (sort rows #'> :key #'third))
              (when (> (third (first rows)) tolerance)
                (format stream "~&  worst residues by RMSD:~%")
                (dolist (r (subseq rows 0 (min limit (length rows))))
                  (destructuring-bind (key cnt res-rmsd res-worst) r
                    (when (> res-rmsd tolerance)
                      (format stream "~&    ~34a ~4d atoms   rmsd ~10,4f   worst ~10,4f~%"
                              key cnt res-rmsd res-worst))))))
            (let ((ok (< rmsd tolerance)))
              ;; Two FORMAT calls, not one with a ~A: a ~ line-continuation is processed only in a
              ;; CONTROL string, so a multi-line message passed as an ARGUMENT prints its tildes
              ;; and indentation verbatim.
              (if ok
                  (format stream "~&coordinates: MATCH - so the nonbond term sets must match too, ~
                                  and any remaining fa_rep difference is an ENUMERATION bug~%")
                  (format stream "~&coordinates: *** DIFFER - geometry, not term sets, and the ~
                                  nonbond difference is downstream of it ***~%"))
              ok))))))

;;; ------------------------------------------------------------------
;;; WHICH slot should the mcstate have picked?
;;;
;;; A sidechain in the WRONG PLACE is the same atoms as a sidechain in the right place - so STAGE 1
;;; CANNOT SEE A ROTAMER SELECTION ERROR.  Every rotamer of a monomer carries identical bonded terms,
;;; identical parameters, over identically-named atoms; only the coordinates differ.  A blueprint
;;; that selects slot 53 where the pose holds rotamer 17 passes the bonded comparison perfectly and
;;; is wrong about every atom past CB.
;;;
;;; So when stage 2 reports a locus misplaced, the useful question is not "how far" but "was the
;;; wrong slot selected, and which one was right".  Every rotamer of the monomer is ALREADY
;;; materialized in the blueprint's coordinates - that is what the fan-out is - so the answer is a
;;; scan over that monomer's slots, and it costs nothing.
;;; ------------------------------------------------------------------

(defun slot-rmsd-against-residue (blueprint mrkindex pose-residue bp-coords pose-coords
                                  bp-table pose-table)
  "RMSD between the residue materialized at MRKINDEX and POSE-RESIDUE, matched by atom name.
NIL when no atom names correspond."
  (let ((scan (aref (topology::rotamer-scans blueprint) mrkindex))
        (n 0) (sum 0d0))
    (when scan
      (let ((pose-by-name (make-hash-table :test 'equal)))
        (chem:do-atoms (a pose-residue)
          (setf (gethash (string (chem:get-name a)) pose-by-name) a))
        (chem:do-atoms (a (topology::residue scan))
          (let ((pa (gethash (string (chem:get-name a)) pose-by-name)))
            (when pa
              (let* ((bi (topology::blueprint-atom-i3 bp-table a))
                     (pi3 (chem:get-coordinate-index-times3 pose-table pa))
                     (dx (- (aref bp-coords bi) (aref pose-coords pi3)))
                     (dy (- (aref bp-coords (+ bi 1)) (aref pose-coords (+ pi3 1))))
                     (dz (- (aref bp-coords (+ bi 2)) (aref pose-coords (+ pi3 2)))))
                (incf n)
                (incf sum (+ (* dx dx) (* dy dy) (* dz dz)))))))))
    (when (plusp n) (sqrt (/ sum n)))))

(defun check-rotamer-selection-against-pose (blueprint pose mcstate bp-coords
                                             &key (stream *standard-output*)
                                                  (pose-coords
                                                   (refold-pose-coordinates pose)))
  "For every locus, scan ALL of its monomer's slots and report which one actually matches the pose.

Answers the question stage 2 raises and stage 1 structurally cannot: was the right ROTAMER selected?
Returns T when the mcstate's slot is the best match at every locus."
  (let* ((assembler (design:assembler pose))
         (bp-table (chem:atom-table (topology:energy-function blueprint)))
         (pose-table (chem:atom-table (topology:energy-function assembler)))
         (pose-residues (pose-monomer-index-to-residue assembler))
         (all-best t))
    (format stream "~&~%!!!! ROTAMER SELECTION: DID THE MCSTATE PICK THE SLOT THE POSE ACTUALLY ~
                      HOLDS?~%~
                    ~&  Stage 1 cannot answer this - every rotamer of a monomer has the same bonded ~
                      terms over the same atom names,~%~
                    ~&  so a wrong selection passes it perfectly.  Each row scans every slot of the ~
                      monomer and names the best fit.~%")
    (format stream "~&  ~26a ~6@a ~10@a ~6@a ~10@a ~8@a ~9@a ~8@a ~7@a ~9@a ~9@a~%"
            "  locus" "picked" "its rmsd" "best" "best rmsd" "slots" "db-index" "n-active"
            "pos" "idx-min" "idx-max")
    (loop for mrkindex across mcstate
          for info = (topology::mrkindex-info blueprint mrkindex)
          for pose-res = (gethash (cons (getf info :shape-index) (getf info :locus-in-shape))
                                  pose-residues)
          do (multiple-value-bind (bp-locus bp-monomer picked-slot)
                 (topology::blueprint-locus-of-mrkindex blueprint mrkindex)
               (let ((picked-rmsd (slot-rmsd-against-residue blueprint mrkindex pose-res
                                                             bp-coords pose-coords
                                                             bp-table pose-table))
                     (best nil) (best-rmsd nil)
                     (count (topology::rotamer-slot-count bp-monomer)))
                 (dotimes (slot count)
                   (let ((r (slot-rmsd-against-residue
                             blueprint (topology::blueprint-mrkindex bp-locus bp-monomer slot)
                             pose-res bp-coords pose-coords bp-table pose-table)))
                     (when (and r (or (null best-rmsd) (< r best-rmsd)))
                       (setf best slot best-rmsd r))))
                 (unless (eql best picked-slot) (setf all-best nil))
                 ;; DB-INDEX is what BLUEPRINT-SLOT-FOR-ROTAMER-INDEX was handed - the rotamers
                 ;; DATABASE index off the monomer-shape - and N-ACTIVE is how many rotamers this
                 ;; monomer actually has in its active set.  Printed because PICKED alone cannot
                 ;; distinguish "the conversion mapped it wrongly" from "the conversion was handed
                 ;; a number that already meant something else": if DB-INDEX equals PICKED, the
                 ;; conversion was a no-op and the database index was used as a slot number
                 ;; directly, which is the bug it was written to fix.
                 ;; POS, IDX-MIN and IDX-MAX describe the INDEXES vector the conversion inverts.
                 ;; DB-INDEX alone cannot distinguish the two ways this can go wrong: a database
                 ;; index that is simply absent from the active set (POS "-", and the locus should
                 ;; have been reported unresolved), from one that is present but at the last
                 ;; position - which is what PICKED = N-ACTIVE-1 would mean, and which IDX-MAX
                 ;; explains if DB-INDEX equals it.
                 (let* ((indexes (nth-value 2 (topology::active-rotamers-for-locus-monomer
                                               bp-locus bp-monomer)))
                        (db (getf info :database-index))
                        (pos (and indexes db (position db indexes))))
                   (format stream "~&  ~26a ~6d ~10,4f ~6d ~10,4f ~8d ~9a ~8a ~7a ~9a ~9a~a~%"
                           (list (getf info :blueprint-locus) (getf info :monomer-name))
                           picked-slot (or picked-rmsd -1d0)
                           (or best -1) (or best-rmsd -1d0) count
                           (or db "-")
                           (if indexes (length indexes) "nil")
                           (or pos "-")
                           (if (and indexes (plusp (length indexes)))
                               (reduce #'min indexes) "-")
                           (if (and indexes (plusp (length indexes)))
                               (reduce #'max indexes) "-")
                           (if (eql best picked-slot) "" "   <- WRONG SLOT"))))))
    (format stream "~&rotamer selection: ~a~%"
            (if all-best
                "the mcstate picked the best-matching slot at every locus"
                "*** the mcstate did NOT pick the best-matching slot - see the WRONG SLOT rows ***"))
    all-best))

;;; ------------------------------------------------------------------
;;; Write the two structures out side by side
;;;
;;; The numbers say the sidechains are 3-6 A out of place; a viewer says WHAT is out of place, and
;;; that is usually the faster read.  Two files rather than one because they are separate molecules
;;; as far as any viewer is concerned - load both and they superimpose, since the backbone matches.
;;;
;;; THE BLUEPRINT FILE IS TRIMMED to the structure MCSTATE describes: the backbone plus the ONE
;;; selected rotamer per locus.  Writing its whole aggregate would emit all ~774 rotamers
;;; simultaneously - every sidechain of every monomer superimposed on every locus - which is what
;;; /tmp/blueprint-fanout.mol2 already shows and which is unreadable for this purpose.
;;; ------------------------------------------------------------------

(defun write-blueprint-and-pose-mol2 (blueprint pose mcstate bp-coords
                                      &key (blueprint-path "~/public_html/blueprint-selected.mol2")
                                           (pose-path "~/public_html/pose-selected.mol2")
                                           (stream *standard-output*)
                                           ;; The GENERATED conformation, so the file shows the
                                           ;; structure the comparison actually scored rather than
                                           ;; the receptor's loaded sidechains.
                                           (pose-coords (refold-pose-coordinates pose)))
  "Write the pose and the blueprint's MCSTATE-selected structure to two mol2 files.

Returns (values BLUEPRINT-PATH POSE-PATH)."
  (let* ((bp-agg (topology::aggregate blueprint))
         (slot-residues (topology::blueprint-slot-residue-set blueprint))
         (selected (make-hash-table :test 'eq))
         (trimmed (chem:make-aggregate :blueprint-selected)))
    ;; The atoms already carry coordinates - BUILD-BLUEPRINT-COORDINATES ends with
    ;; COPY-ALL-JOINT-POSITIONS-INTO-ATOMS - but re-applying makes this independent of anything the
    ;; bench did in between, and SAVE-MOL2 writes atom positions, not a coordinate vector.
    (chem:matter/apply-coordinates bp-agg bp-coords)
    (loop for mrkindex across mcstate
          for scan = (aref (topology::rotamer-scans blueprint) mrkindex)
          do (setf (gethash (topology::residue scan) selected) t))
    ;; Walked in the blueprint's own residue order so the output reads like the structure rather
    ;; than like the mcstate.  Residues are SHARED, not copied - ADD-MATTER-DONT-UPDATE-ID leaves
    ;; the child untouched and sets no parent, so the real aggregate is unaffected.
    (dotimes (m (chem:content-size bp-agg))
      (let* ((bp-mol (chem:content-at bp-agg m))
             (mol (chem:make-molecule (chem:get-name bp-mol)))
             (kept 0))
        (loop for i below (chem:content-size bp-mol)
              for res = (chem:content-at bp-mol i)
              when (or (not (gethash res slot-residues)) (gethash res selected))
                do (chem:add-matter-dont-update-id mol res) (incf kept))
        (chem:add-matter trimmed mol)
        (format stream "~&  blueprint molecule ~d: ~d residues written~%" m kept)))
    ;; BREAK THE UNSELECTED FAN-OUT BONDS AROUND THE WRITE.
    ;;
    ;; The detach scan ends with (MAPC #'RESTORE-FAN-OUT-BOND RECORDS), so every backbone attachment
    ;; atom carries a bond to EVERY rotamer's first atom - all ~774 of them.  TRIMMED holds those
    ;; backbone atoms but only 13 of the rotamers, so SAVE-MOL2 walks bonds whose partner is not in
    ;; the file and tries to resolve an index that was never assigned.  That is a GC pointer-tag
    ;; abort, not a Lisp error - it took the process down with no backtrace.
    ;;
    ;; /tmp/blueprint-fanout.mol2 is unaffected because it writes the WHOLE aggregate, where every
    ;; bond partner is present.  Trimming is what makes the dangling bonds reachable.
    ;;
    ;; RECORDS is indexed by mrkindex, one per slot - the same indexing the detach scan uses.
    (let ((records (topology::blueprint-fan-out-bonds blueprint))
          (keep (make-hash-table))
          (broken nil))
      (loop for mrkindex across mcstate do (setf (gethash mrkindex keep) t))
      (unwind-protect
           (progn
             (loop for record in records
                   for mrkindex from 0
                   unless (gethash mrkindex keep)
                     do (topology::break-fan-out-bond record)
                        (push record broken))
             (chem:save-mol2 trimmed blueprint-path t)
             ;; MUTATES the pose's aggregate atom positions - SAVE-MOL2 writes atom positions, not
             ;; a coordinate vector, so the refold has to be pushed into the atoms to appear in the
             ;; file.  The pose's COORDS slot is untouched; only the aggregate the viewer reads is.
             (chem:matter/apply-coordinates (design:pose-aggregate pose) pose-coords)
             (chem:save-mol2 (design:pose-aggregate pose) pose-path t))
        ;; EXACTLY what was broken, not every record.  Restoring all of them would BOND-TO the 13
        ;; selected slots a second time - their bonds were never removed - leaving those CAs
        ;; carrying a duplicate.  RESTORE-FAN-OUT-BOND passes :ERROR-IF-EXCEED-VALENCE NIL, which is
        ;; correct for putting 774 legitimate fan-out bonds back on one CA and is exactly why a
        ;; duplicate would NOT be caught here - it would surface later, somewhere unrelated.
        ;;
        ;; Restored on a non-local exit too: leaving the fan detached would silently change every
        ;; later energy evaluation, and nothing downstream re-attaches it.
        (mapc #'topology::restore-fan-out-bond broken)))
    (format stream "~&  wrote ~a~&  wrote ~a~%" blueprint-path pose-path)
    (values blueprint-path pose-path)))

;;; ------------------------------------------------------------------
;;; What rotamer does each slot actually hold?
;;;
;;; Every rotamer index in this system is an index into ONE vector: the ROTAMER-VECTOR of a
;;; monomer-context's rotamers container.  A slot number is not that index - it is a position in the
;;; ACTIVE subset - and a MONOMER-SHAPE's ROTAMER-INDEX is that index.  Every confusion in this area
;;; has come from treating one as the other.
;;;
;;; So: write out, for every slot, WHICH database rotamer it was built from.  Then the pose's
;;; monomer-shape index can simply be looked up in that table, and "which slot should the mcstate
;;; have picked" stops being an inference about POSITION and INDEXES and becomes a search.
;;; ------------------------------------------------------------------

(defun report-slot-rotamer-indexes (blueprint persona pose mcstate
                                    &key (stream *standard-output*) (limit 24))
  "For each locus the mcstate names, write every slot's database rotamer index and locate the pose's.

Returns T when the pose's rotamer is found in a slot AND that slot is the one the mcstate picked."
  (let ((all-ok t)
        (rows nil))
    (format stream "~&~%!!!! SLOT -> DATABASE ROTAMER INDEX, AND WHERE THE POSE'S ROTAMER LIVES~%~
                    ~&  A slot number is a position in the monomer's ACTIVE subset; a monomer-shape's~
                      ROTAMER-INDEX is an index into the~%~
                    ~&  context's full ROTAMER-VECTOR.  This writes the mapping between them, so the ~
                      pose's index can be looked up~%~
                    ~&  rather than inferred.~%")
    (loop for mrkindex across mcstate
          for info = (topology::mrkindex-info blueprint mrkindex)
          do (multiple-value-bind (bp-locus bp-monomer picked-slot)
                 (topology::blueprint-locus-of-mrkindex blueprint mrkindex)
               (let* ((count (topology::rotamer-slot-count bp-monomer))
                      (pose-index (getf info :database-index))
                      ;; Read off the SCAN, which records what the fold actually applied, rather
                      ;; than recomputed from the persona.  A recomputation would agree with the
                      ;; persona by construction and could not disagree with the fold - and the
                      ;; fold is the thing whose output the coordinates came from.
                      (table (loop for slot below count
                                   collect (topology::rotamer-catalogue-index
                                            (aref (topology::rotamer-scans blueprint)
                                                  (topology::blueprint-mrkindex
                                                   bp-locus bp-monomer slot)))))
                      (where (position pose-index table)))
                 (unless (and where (eql where picked-slot)) (setf all-ok nil))
                 (push (list (list (getf info :blueprint-locus) (getf info :monomer-name))
                             picked-slot
                             (when (and picked-slot (< picked-slot count)) (nth picked-slot table))
                             pose-index)
                       rows)
                 (format stream "~&~%  locus ~a ~a   ~d slots   pose monomer-shape rotamer-index ~a~%"
                         (getf info :blueprint-locus) (getf info :monomer-name) count pose-index)
                 ;; WINDOWED ON THE PICKED SLOT, not on slot 0.  The monomers whose selection is in
                 ;; question carry 55-76 slots and pick slot 53-74, so a table that always starts at
                 ;; the beginning truncates away the single entry the comparison is about.
                 (let* ((start (if (<= count limit)
                                   0
                                   (max 0 (min (- (or picked-slot 0) (floor limit 2))
                                               (- count limit)))))
                        (end (min count (+ start limit))))
                   (format stream "~&    slot -> database index~a:~%      ~{~a~^ ~}~%"
                           (if (> count limit)
                               (format nil " (slots ~d-~d of ~d)" start (1- end) count)
                               "")
                           (loop for slot from start below end
                                 collect (format nil "~d:~a~a" slot (nth slot table)
                                                 (if (eql slot picked-slot) "*" "")))))
                 (cond ((null where)
                        (format stream "~&    *** the pose's rotamer index ~a is in NO slot - the ~
                                          blueprint never materialized it ***~%" pose-index))
                       ((eql where picked-slot)
                        (format stream "~&    pose index ~a is at slot ~d; mcstate picked ~d  ~
                                          <- AGREE~%" pose-index where picked-slot))
                       (t
                        (format stream "~&    pose index ~a is at slot ~d; mcstate picked ~d  ~
                                          *** MCSTATE PICKED THE WRONG SLOT ***~%"
                                pose-index where picked-slot)))))
          finally
             ;; The two numbers, side by side, with nothing between them.  Both are catalogue
             ;; indexes - positions in a monomer-context's ROTAMER-VECTOR - and they are derived
             ;; INDEPENDENTLY: the blueprint's is stamped onto the scan from the rotamer object the
             ;; fan-out actually installed, the pose's is read off its ROTAMER-SHAPE.  So this is a
             ;; real comparison rather than a number checked against itself.
             ;;
             ;; What it CANNOT settle is whether the two indexes point into the same vector: the
             ;; catalogue is context-to-rotamers[CONTEXT] and only THEN rotamer-vector[INDEX].  Equal
             ;; numbers under different contexts name different conformations, which is what
             ;; CHECK-MONOMER-CONTEXT-AGAINST-POSE exists to rule in or out.
             (progn
               (format stream "~&~%  catalogue index at the matching slot, blueprint vs pose:~%")
               (format stream "~&  ~26a ~12@a ~24@a ~20@a~%"
                       "  locus" "picked slot" "blueprint scan index" "pose rotamer-index")
               (dolist (r (nreverse rows))
                 (destructuring-bind (label slot bp-index pose-index) r
                   (format stream "~&  ~26a ~12a ~24a ~20a~a~%"
                           label (or slot "-") (or bp-index "-") (or pose-index "-")
                           (cond ((null bp-index) "   <- no scan at that slot")
                                 ((eql bp-index pose-index) "")
                                 (t "   <- DIFFERENT ROTAMER")))))
               (format stream "~&~%slot/rotamer index audit: ~a~%"
                       (if all-ok
                           "every locus resolved to the slot holding the pose's rotamer"
                           "*** see the flagged loci above ***"))))
    all-ok))

;;; ------------------------------------------------------------------
;;; THE PRECONDITION FOR COMPARING A POSE TO A BLUEPRINT AT ALL
;;;
;;; A blueprint can only reproduce a conformation it can NAME, and the only names it has are
;;; catalogue indices - positions in a monomer-context's ROTAMER-VECTOR.  So the comparison is
;;; meaningful only when every mobile monomer-shape in the pose is a ROTAMER-SHAPE whose
;;; ROTAMER-INDEX is a real catalogue index.
;;;
;;; A RESIDUE-SHAPE breaks that, and does so SILENTLY.  It holds ATOM-COORDINATES directly - a
;;; conformation that came from geometry optimization or from experiment, not from the library - and
;;; shape.lisp:139 defines
;;;
;;;     (defmethod rotamer-index ((shape residue-shape)) (closest-rotamer-index shape))
;;;
;;; so asking it for a rotamer index returns the NEAREST catalogue entry rather than what it holds.
;;; Nothing errors.  The blueprint faithfully builds that nearest entry, and every atom differs by
;;; however far the real conformation sits from the closest thing the library can express - which for
;;; a sidechain is a few tenths of an Angstrom near CB and several Angstroms at the tip.
;;;
;;; That failure is indistinguishable, from the energies alone, from a blueprint bug.  Hence a gate.
;;; ------------------------------------------------------------------

(defun check-pose-shapes-are-rotamers (blueprint pose &key (stream *standard-output*))
  "Verify every slot-owning locus of POSE holds a ROTAMER-SHAPE.  Returns T when the pose is
comparable to a blueprint at all.

POSE is accepted for symmetry with the other checks; the shapes are read through the blueprint's
loci, which share the assembler's oligomer-shapes."
  (declare (ignorable pose))
  (let ((owning (blueprint-slot-owning-loci blueprint))
        (bad nil))
    (format stream "~&~%!!!! CAN THIS POSE BE COMPARED TO A BLUEPRINT AT ALL?~%~
                    ~&  A blueprint can only reproduce a conformation it can NAME, and its only ~
                      names are catalogue indices.  A~%~
                    ~&  RESIDUE-SHAPE holds coordinates instead, and answers ROTAMER-INDEX with the ~
                      CLOSEST catalogue entry - so~%~
                    ~&  the comparison silently becomes 'blueprint vs the nearest library rotamer ~
                      to the pose'.~%")
    (loop for i below (length owning)
          for bp-locus = (aref owning i)
          for shape = (topology::oligomer-shape bp-locus)
          for monomer-shape = (gethash (topology::original-monomer bp-locus)
                                       (topology::monomer-shape-map shape))
          do (let ((class (and monomer-shape (class-name (class-of monomer-shape)))))
               (unless (typep monomer-shape 'topology:rotamer-shape)
                 (push (list (topology::locus bp-locus) class) bad))
               (format stream "~&  locus ~3d  ~a~a~%"
                       (topology::locus bp-locus)
                       (or class "NO MONOMER-SHAPE")
                       (cond ((null monomer-shape) "   <- nothing to read")
                             ((typep monomer-shape 'topology:rotamer-shape) "")
                             (t "   <- NOT a rotamer-shape")))))
    (if bad
        (format stream "~&~%*** ~d locus/loci do not hold a ROTAMER-SHAPE.  Every comparison below ~
                          is measuring the distance from the pose's~%~
                        ~&    real conformation to the nearest catalogue rotamer, NOT a blueprint ~
                          defect.  A pose is comparable only~%~
                        ~&    when it was built FROM catalogue rotamers and has not been ~
                          geometry-optimized since. ***~%"
                (length bad))
        (format stream "~&  every slot-owning locus holds a ROTAMER-SHAPE - the pose is comparable~%"))
    (null bad)))

;;; ------------------------------------------------------------------
;;; Same rotamer, different geometry - so compare the JOINTS
;;;
;;; With selection confirmed correct and the pose holding genuine ROTAMER-SHAPEs, the blueprint's
;;; slot and the pose's residue are built from THE SAME rotamer object, hence the same
;;; INTERNALS-VALUES.  APPLY-FRAGMENT-INTERNALS-TO-ATRESIDUE consumes that vector POSITIONALLY:
;;;
;;;   (loop for joint across (joints atresidue)
;;;         for index3 from 0 by 3
;;;         do (fill-joint-internals-from-vector ... joint (internals-values rotamer) index3))
;;;
;;; so identical internals produce identical geometry ONLY IF the two atresidues present the same
;;; joints in the same order.  The blueprint builds its slot atresidues in
;;; BUILD-BLUEPRINT-FAN-OUT-JOINTS; the pose's come from MAKE-ASSEMBLER.  Two construction paths,
;;; never checked against each other.
;;;
;;; INDEX3 STEPS BY 3 FOR EVERY JOINT, including the ones that consume nothing:
;;; FILL-JOINT-INTERNALS-FROM-VECTOR is a no-op for XYZ-JOINT and JUMP-JOINT (internals.lisp:970-971)
;;; and does real work only for BONDED-JOINT.  So a difference in joint TYPE shifts the vector just
;;; as surely as a difference in order does, and a single extra or missing joint at the head puts
;;; every subsequent internal on the wrong joint.
;;; ------------------------------------------------------------------

(defun pose-monomer-index-to-atresidue (assembler)
  "Hash (SHAPE-INDEX . MONOMER-INDEX) -> ATRESIDUE - the atresidue counterpart of
POSE-MONOMER-INDEX-TO-RESIDUE, for reaching the pose's joints."
  (let ((map (make-hash-table :test 'equal))
        (atmolecules (topology::atmolecules (topology::ataggregate assembler))))
    (dotimes (shape-index (length atmolecules))
      (loop for atres across (topology::atresidues (aref atmolecules shape-index))
            when atres
              do (setf (gethash (cons shape-index (topology::monomer-index atres)) map) atres)))
    map))

(defun joint-signature (atresidue)
  "(NAME . TYPE) per joint, in the order APPLY-FRAGMENT-INTERNALS-TO-ATRESIDUE walks them.

TYPE is included because it decides whether a joint CONSUMES its three internals or ignores them -
and a joint that ignores them still advances INDEX3."
  (loop for joint across (topology:joints atresidue)
        collect (cons (kin:joint/name joint) (class-name (class-of joint)))))

(defun check-joint-order-against-pose (blueprint pose mcstate &key (stream *standard-output*))
  "Compare each selected slot's atresidue joints against the pose's, position by position.

Returns T when every locus matches in count, order and joint type."
  (let* ((assembler (design:assembler pose))
         (pose-atres (pose-monomer-index-to-atresidue assembler))
         (all-ok t))
    (format stream "~&~%!!!! JOINT ORDER: DO THE TWO ATRESIDUES CONSUME THE INTERNALS THE SAME WAY?~%~
                    ~&  The same rotamer's INTERNALS-VALUES are applied POSITIONALLY, one joint at a ~
                      time, three values each -~%~
                    ~&  and a no-op joint still advances the index.  So identical internals give ~
                      identical geometry only if these~%~
                    ~&  two lists agree exactly.~%")
    (loop for mrkindex across mcstate
          for info = (topology::mrkindex-info blueprint mrkindex)
          for scan = (aref (topology::rotamer-scans blueprint) mrkindex)
          for pose-at = (gethash (cons (getf info :shape-index) (getf info :locus-in-shape))
                                 pose-atres)
          do (if (null pose-at)
                 (format stream "~&  locus ~a ~a: no pose atresidue~%"
                         (getf info :blueprint-locus) (getf info :monomer-name))
                 (let* ((bp-sig (joint-signature (topology::atresidue scan)))
                        (pose-sig (joint-signature pose-at))
                        (first-diff (loop for a in bp-sig for b in pose-sig for i from 0
                                          unless (equal a b) return i)))
                   (cond ((and (= (length bp-sig) (length pose-sig)) (null first-diff))
                          (format stream "~&  locus ~a ~a: ~d joints, identical~%"
                                  (getf info :blueprint-locus) (getf info :monomer-name)
                                  (length bp-sig)))
                         (t
                          (setf all-ok nil)
                          (format stream "~&  locus ~a ~a: blueprint ~d joints, pose ~d joints~
                                          ~@[, first difference at position ~d~]~%"
                                  (getf info :blueprint-locus) (getf info :monomer-name)
                                  (length bp-sig) (length pose-sig) first-diff)
                          ;; From the first divergence onward, because everything before it is
                          ;; already known to agree and everything after is downstream of the shift.
                          (let ((from (or first-diff 0)))
                            (loop for i from from below (min (+ from 8)
                                                             (max (length bp-sig) (length pose-sig)))
                                  do (format stream "~&      ~3d  blueprint ~28a  pose ~28a~a~%"
                                             i (or (nth i bp-sig) "-") (or (nth i pose-sig) "-")
                                             (if (equal (nth i bp-sig) (nth i pose-sig))
                                                 "" "  <-")))))))))
    (format stream "~&joint order: ~a~%"
            (if all-ok
                "every selected slot consumes the internals exactly as the pose does"
                "*** the two atresidues do NOT agree - identical internals cannot give identical ~
                 geometry ***"))
    all-ok))

;;; ------------------------------------------------------------------
;;; The joints INSIDE the residue agree - so compare what they hang OFF
;;;
;;; JOINT-SIGNATURE compares the atresidue's own joints, and they match everywhere.  But the first
;;; joint's geometry is not determined by its internals alone: Joint_O::updateXyzCoords builds an
;;; input STUB from the already-placed parent chain (joint.cc:568) and places the atom relative to
;;; THAT.  The chain runs OUT of the residue, into the backbone, where JOINT-SIGNATURE never looked.
;;;
;;; So CB can carry the identical bond/angle/dihedral and still land somewhere else, if the three
;;; atoms its dihedral is measured against are not the same three - and every atom past it inherits
;;; the difference as a rotation about the CA-CB axis.  That is consistent with sidechains that are
;;; chemically reasonable, correctly attached, and progressively wronger toward the tip.
;;; ------------------------------------------------------------------

(defun joint-ancestry (joint depth)
  "(NAME ...) for JOINT and its first DEPTH ancestors, root-ward.  NIL entries where the chain ends."
  (loop repeat (1+ depth)
        for j = joint then (and j (kin:parent j))
        collect (and j (kin:joint/name j))))

(defun check-root-joint-ancestry-against-pose (blueprint pose mcstate
                                               &key (stream *standard-output*) (depth 3))
  "Compare each selected slot's ROOT joint ancestry against the pose's.

The root joint's dihedral is measured against its parent, grandparent and great-grandparent.  If
those differ between the two structures, identical internals place the sidechain differently."
  (let* ((assembler (design:assembler pose))
         (pose-atres (pose-monomer-index-to-atresidue assembler))
         (all-ok t))
    (format stream "~&~%!!!! ROOT JOINT ANCESTRY: IS THE SIDECHAIN MEASURED AGAINST THE SAME ~
                      ATOMS?~%~
                    ~&  The joints INSIDE each residue already match.  But the first one is placed ~
                      relative to a stub built from its~%~
                    ~&  parent chain, which runs out into the backbone - so the same dihedral ~
                      against a different reference puts~%~
                    ~&  the whole sidechain somewhere else.~%")
    (loop for mrkindex across mcstate
          for info = (topology::mrkindex-info blueprint mrkindex)
          for scan = (aref (topology::rotamer-scans blueprint) mrkindex)
          for pose-at = (gethash (cons (getf info :shape-index) (getf info :locus-in-shape))
                                 pose-atres)
          do (if (null pose-at)
                 (format stream "~&  locus ~a ~a: no pose atresidue~%"
                         (getf info :blueprint-locus) (getf info :monomer-name))
                 (let* ((bp-root (aref (topology::joints (topology::atresidue scan)) 0))
                        (pose-root (aref (topology::joints pose-at) 0))
                        (bp-chain (joint-ancestry bp-root depth))
                        (pose-chain (joint-ancestry pose-root depth))
                        (same (equal bp-chain pose-chain)))
                   (unless same (setf all-ok nil))
                   (format stream "~&  locus ~5a ~10a blueprint ~{~a~^<-~}~%~
                                   ~&  ~17a pose      ~{~a~^<-~}~a~%"
                           (getf info :blueprint-locus) (getf info :monomer-name)
                           bp-chain "" pose-chain (if same "" "   <- DIFFERS")))))
    (format stream "~&root joint ancestry: ~a~%"
            (if all-ok
                "every selected slot is measured against the same atoms the pose uses"
                "*** the reference frames DIFFER - identical internals cannot place the sidechain ~
                 identically ***"))
    all-ok))

;;; ------------------------------------------------------------------
;;; The index is right - but is it an index into the SAME VECTOR?
;;;
;;; APPLY-MONOMER-SHAPE-TO-ATRESIDUE-INTERNALS (shape.lisp:97-107) resolves a rotamer in two steps:
;;;
;;;     container = context-to-rotamers[MONOMER-CONTEXT]
;;;     rotamer   = (rotamer-vector container)[ROTAMER-INDEX]
;;;
;;; so ROTAMER-INDEX names a conformation only once the CONTEXT is fixed.  Two sides that agree on
;;; the NUMBER and disagree on the CONTEXT agree on nothing at all.
;;;
;;; Every check written so far is structurally blind to this.  A monomer's topology, atom names,
;;; joints, joint order, parent chain and bonded terms are the SAME in every context it appears in -
;;; contexts differ only in which conformations were catalogued.  So a context mismatch passes stage
;;; 1, passes the joint checks, passes the ancestry check, and passes the slot/database audit (both
;;; sides really do say 5605), while changing every sidechain coordinate.  That is the whole symptom.
;;;
;;; The two sides get their context by different routes, which is what makes this worth testing:
;;; the pose reads a table precomputed at assembler build time (assembler.lisp:1168,
;;; MONOMER-CONTEXTS), the blueprint recomputes it per persona against the oligomer's CURRENT
;;; selections (PERSONA-MONOMER-CONTEXT).  Nothing forces those to agree.
;;;
;;; ROTAMER-VECTOR lengths are printed beside the contexts because they are an independent witness:
;;; two different containers usually hold different numbers of rotamers, so a length mismatch
;;; confirms a context mismatch without trusting the symbols to print distinguishably.
;;; ------------------------------------------------------------------

(defun pose-monomer-context-map (assembler)
  "Hash (SHAPE-INDEX . MONOMER-INDEX) -> MONOMER-CONTEXT, keyed the way
POSE-MONOMER-INDEX-TO-ATRESIDUE keys, so a locus can be looked up with the same key stage 1 uses.

Built by inverting MONOMER-POSITIONS, since the context table is keyed by MONOMER and the pose's
monomer objects are not the blueprint's."
  (let ((map (make-hash-table :test 'equal))
        (atmolecules (topology::atmolecules (topology::ataggregate assembler))))
    (maphash (lambda (monomer context)
               (let ((pos (gethash monomer (topology:monomer-positions assembler))))
                 (when pos
                   (let* ((mi (topology:molecule-index pos))
                          (ri (topology:residue-index pos))
                          (atmol (when (< mi (length atmolecules)) (aref atmolecules mi)))
                          (atres (when (and atmol (< ri (length (topology::atresidues atmol))))
                                   (aref (topology::atresidues atmol) ri))))
                     (when atres
                       (setf (gethash (cons mi (topology::monomer-index atres)) map) context))))))
             (topology:monomer-contexts assembler))
    map))

(defun context-rotamer-count (context-to-rotamers context)
  "How many rotamers CONTEXT's container holds, or NIL when it resolves to no container."
  (let ((container (and context-to-rotamers context (gethash context context-to-rotamers))))
    (when container (length (topology:rotamer-vector container)))))

(defun check-monomer-context-against-pose (blueprint pose mcstate
                                           &key (stream *standard-output*))
  "Compare the MONOMER-CONTEXT each side resolves at every selected locus.

Returns T when the contexts agree everywhere.  A disagreement means the two sides indexed different
ROTAMER-VECTORs with the same number, which reproduces every symptom seen so far."
  (let* ((assembler (design:assembler pose))
         (pose-contexts (pose-monomer-context-map assembler))
         (all-ok t))
    (format stream "~&~%!!!! MONOMER CONTEXT: IS THE ROTAMER INDEX AN INDEX INTO THE SAME VECTOR?~%~
                    ~&  A rotamer index is meaningless without its context - the catalogue is ~
                      context-to-rotamers[CONTEXT] and~%~
                    ~&  THEN rotamer-vector[INDEX].  Both sides saying \"5605\" proves nothing if ~
                      they are indexing different vectors,~%~
                    ~&  and a context mismatch is invisible to every earlier check because topology ~
                      does not vary with context.~%")
    (format stream "~&  ~26a ~38a ~38a ~8@a ~8@a~%"
            "  locus" "blueprint context" "pose context" "bp rots" "pose rots")
    (loop for mrkindex across mcstate
          for info = (topology::mrkindex-info blueprint mrkindex)
          do (multiple-value-bind (bp-locus bp-monomer)
                 (topology::blueprint-locus-of-mrkindex blueprint mrkindex)
               (let* ((bp-context (topology::persona-monomer-context bp-monomer bp-locus))
                      (pose-context (gethash (cons (getf info :shape-index)
                                                   (getf info :locus-in-shape))
                                             pose-contexts))
                      (db (topology:rotamers-database (topology:oligomer-shape bp-locus)))
                      (c2r (when db (topology:context-to-rotamers db)))
                      (bp-n (context-rotamer-count c2r bp-context))
                      (pose-n (context-rotamer-count c2r pose-context))
                      (same (equal bp-context pose-context)))
                 (unless same (setf all-ok nil))
                 (format stream "~&  ~26a ~38a ~38a ~8a ~8a~a~%"
                         (list (getf info :blueprint-locus) (getf info :monomer-name))
                         (or bp-context "-") (or pose-context "-")
                         (or bp-n "-") (or pose-n "-")
                         (if same "" "   <- DIFFERENT CONTEXT")))))
    ;; Written as two FORMAT calls rather than one with a ~A: a ~ line-continuation is processed
    ;; only in a CONTROL string, so a multi-line message passed as an ARGUMENT prints its own
    ;; tildes and indentation verbatim.
    (if all-ok
        (format stream "~&monomer context: both sides index the same ROTAMER-VECTOR - equal ~
                        indexes really do mean equal rotamers, so the coordinate difference is ~
                        NOT a catalogue mismatch~%")
        (format stream "~&monomer context: *** the two sides index DIFFERENT ROTAMER-VECTORs - ~
                        the same index names a different conformation on each side, which ~
                        explains identical terms, identical joints and different ~
                        coordinates ***~%"))
    all-ok))

;;; ------------------------------------------------------------------
;;; THE NUMBERS THAT ACTUALLY GET APPLIED, IN DEGREES
;;;
;;; Every check so far compares a NAME for a conformation - a slot, a catalogue index, a context.
;;; This compares the conformation itself: the bond, angle and dihedral each side hands to each
;;; joint.  It is the last comparison upstream of the coordinates, so it splits the remaining
;;; possibilities cleanly and with no room for a third answer:
;;;
;;;   dihedrals DIFFER  -> the two sides resolved different rotamers.  The lookup is at fault, and
;;;                        the size of the difference says which kind: ~120 deg at chi1 is a
;;;                        different rotamer WELL (a different catalogue entry), a few degrees
;;;                        everywhere is the same entry perturbed.
;;;   dihedrals AGREE   -> the same numbers produced different coordinates.  The lookup is
;;;                        exonerated and the fault is in the FOLD - which, with joint order and
;;;                        parent chains already proven identical, leaves the internals being
;;;                        modified after they are filled (adjustments) as the candidate.
;;;
;;; Both sides are read the way FILL-JOINT-INTERNALS-FROM-VECTOR reads them (internals.lisp:973):
;;; three values per joint at 3i, 3i+1, 3i+2, radians for the two angles, and a joint that ignores
;;; its internals STILL ADVANCES the index.  Reading them any other way would silently shear the
;;; two lists past each other and manufacture a difference.
;;;
;;; The two rotamer objects are obtained by the routes each side genuinely uses - the blueprint's
;;; from the persona (SLOT-ROTAMER-INTERNALS, which holds objects and does no context lookup at
;;; all), the pose's by replaying shape.lisp:97-107 through the pose's own monomer-context.  Not
;;; from each other, or this would compare a number with itself.
;;; ------------------------------------------------------------------

(defun rad-to-degrees (radians)
  (when radians (* radians (/ 180d0 pi))))

(defun degrees-difference (a b)
  "A - B wrapped into (-180, 180].  NIL if either is missing.

Wrapped because these are angles: -179 and 179 are two degrees apart, and a raw subtraction would
call them 358 and flag every locus."
  (when (and a b)
    (let ((d (mod (- a b) 360d0)))
      (if (> d 180d0) (- d 360d0) d))))

(defun rotamer-joint-internals (rotamer joints)
  "Per joint, (NAME BONDEDP BOND ANGLE-DEG DIHEDRAL-DEG) out of ROTAMER's INTERNALS-VALUES.

JOINTS supplies the walk, and the same vector is used for both sides deliberately: CHECK-JOINT-
ORDER-AGAINST-POSE proved the two atresidues carry identical joints in identical order, so one walk
keeps the two readouts index-aligned by construction.  BOND and the angles come along because
reporting 'the dihedrals agree' while the bonds differ would be true and misleading."
  (let ((values (when rotamer (topology:internals-values rotamer))))
    (loop for joint across joints
          for index3 from 0 by 3
          for bondedp = (typep joint 'kin:bonded-joint)
          for avail = (and values bondedp (< (+ index3 2) (length values)))
          collect (list (kin:joint/name joint)
                        bondedp
                        (when avail (aref values index3))
                        (when avail (rad-to-degrees (aref values (+ index3 1))))
                        (when avail (rad-to-degrees (aref values (+ index3 2))))))))

(defun pose-applied-rotamer (bp-locus assembler context-map info)
  "The ROTAMER object the pose applies at BP-LOCUS, resolved exactly as
APPLY-MONOMER-SHAPE-TO-ATRESIDUE-INTERNALS resolves it (shape.lisp:97-107):

    container = context-to-rotamers[MONOMER-CONTEXT]
    rotamer   = (rotamer-vector container)[ROTAMER-INDEX]

The oligomer-shape is SHARED with the pose's assembler, so the monomer-shape and the rotamers
database are the pose's own, not copies.  Returns (values ROTAMER ROTAMER-INDEX CONTEXT)."
  (let* ((shape (topology:oligomer-shape bp-locus))
         (monomer (topology::original-monomer bp-locus))
         (monomer-shape (gethash monomer (topology:monomer-shape-map shape)))
         (index (when (typep monomer-shape 'topology:rotamer-shape)
                  (topology:rotamer-index monomer-shape)))
         ;; Keyed by MONOMER first because the blueprint's locus and the assembler share the
         ;; oligomer's monomer objects; the positional map is the fallback if that ever stops
         ;; being true.
         (context (or (gethash monomer (topology:monomer-contexts assembler))
                      (gethash (cons (getf info :shape-index) (getf info :locus-in-shape))
                               context-map)))
         (db (topology:rotamers-database shape))
         (container (when (and db context)
                      (gethash context (topology:context-to-rotamers db))))
         (vector (when container (topology:rotamer-vector container))))
    (values (when (and vector index (< index (length vector))) (aref vector index))
            index context)))

(defun check-internal-dihedrals-against-pose (blueprint persona pose mcstate
                                              &key (stream *standard-output*)
                                                   (tolerance 1d-3)
                                                   (show-all nil))
  "Compare, in DEGREES, the internal dihedral each side applies to every joint of every mobile
sidechain.  Returns T when every dihedral agrees to within TOLERANCE degrees.

With SHOW-ALL nil only joints that differ are listed per locus, plus a per-locus summary; the
agreeing joints are counted, not printed."
  (let* ((assembler (design:assembler pose))
         (context-map (pose-monomer-context-map assembler))
         (all-ok t)
         (rows nil))
    (format stream "~&~%!!!! INTERNAL DIHEDRALS: DO THE TWO SIDES APPLY THE SAME ANGLES?~%~
                    ~&  Not a name for a conformation - the conformation.  Each side's rotamer is ~
                      fetched the way that side really~%~
                    ~&  fetches it, and its INTERNALS-VALUES are read the way ~
                      FILL-JOINT-INTERNALS-FROM-VECTOR reads them: three per~%~
                    ~&  joint, angles in radians, no-op joints still advancing the index.  ~
                      Differences near 120 deg mean a different~%~
                    ~&  rotamer well; a few degrees everywhere means the same entry perturbed.~%")
    (loop for mrkindex across mcstate
          for info = (topology::mrkindex-info blueprint mrkindex)
          do (multiple-value-bind (bp-locus bp-monomer picked-slot)
                 (topology::blueprint-locus-of-mrkindex blueprint mrkindex)
               (declare (ignore bp-monomer))
               (let* ((scan (aref (topology::rotamer-scans blueprint) mrkindex))
                      (joints (topology:joints (topology::atresidue scan)))
                      (bp-rotamer (topology::slot-rotamer-internals persona mrkindex))
                      (bp-rows (rotamer-joint-internals bp-rotamer joints)))
                 (multiple-value-bind (pose-rotamer pose-index pose-context)
                     (pose-applied-rotamer bp-locus assembler context-map info)
                   (let* ((pose-rows (rotamer-joint-internals pose-rotamer joints))
                          (compared 0) (differing 0)
                          (worst 0d0) (worst-name nil)
                          (worst-bond 0d0) (worst-angle 0d0))
                     (format stream "~&~%  locus ~a ~a   slot ~a   blueprint catalogue index ~a   ~
                                     pose rotamer-index ~a~%"
                             (getf info :blueprint-locus) (getf info :monomer-name)
                             picked-slot
                             (if scan (topology::rotamer-catalogue-index scan) "-")
                             (or pose-index "-"))
                     (cond
                       ((null bp-rotamer)
                        (setf all-ok nil)
                        (format stream "~&    *** the persona has NO rotamer at this slot - ~
                                          nothing was applied ***~%"))
                       ((null pose-rotamer)
                        (setf all-ok nil)
                        (format stream "~&    *** could not resolve the pose's rotamer ~
                                          (context ~a) - nothing to compare ***~%"
                                (or pose-context "unknown")))
                       (t
                        (when show-all
                          (format stream "~&    ~8a ~12@a ~12@a ~12@a~%"
                                  "joint" "pose deg" "blueprint" "delta"))
                        (loop for b in bp-rows
                              for p in pose-rows
                              do (destructuring-bind (name bondedp bond angle dihedral) b
                                   (declare (ignore bond angle))
                                   (destructuring-bind (pname pbondedp pbond pangle pdihedral) p
                                     (declare (ignore pname pbondedp))
                                     (when (and bondedp dihedral pdihedral)
                                       (incf compared)
                                       ;; BOND and ANGLE are tracked alongside so the verdict
                                       ;; cannot say "the dihedrals agree" while the bonds do not.
                                       (let ((delta (degrees-difference pdihedral dihedral))
                                             (dbond (abs (- pbond (third b))))
                                             (dangle (abs (or (degrees-difference pangle (fourth b))
                                                              0d0))))
                                         (when (> dbond worst-bond) (setf worst-bond dbond))
                                         (when (> dangle worst-angle) (setf worst-angle dangle))
                                         (when (> (abs delta) (abs worst))
                                           (setf worst delta worst-name name))
                                         (when (> (abs delta) tolerance) (incf differing))
                                         (when (or show-all (> (abs delta) tolerance))
                                           (format stream "~&    ~8a ~12,3f ~12,3f ~12,3f~a~%"
                                                   name pdihedral dihedral delta
                                                   (if (> (abs delta) tolerance) "  <-" ""))))))))
                        (when (plusp differing) (setf all-ok nil))
                        (format stream "~&    ~d dihedrals compared, ~d differ by more than ~
                                        ~,3f deg~@[, worst ~,3f deg at ~a~]~%"
                                compared differing tolerance
                                (when worst-name worst) worst-name)
                        (when (or (> worst-bond 1d-4) (> worst-angle 1d-3))
                          (format stream "~&    (bonds differ by up to ~,4f A, bond angles by up ~
                                          to ~,3f deg - so this is not a dihedral-only ~
                                          difference)~%"
                                  worst-bond worst-angle))
                        (push (list (list (getf info :blueprint-locus) (getf info :monomer-name))
                                    compared differing worst worst-name)
                              rows)))))))
          finally
             (progn
               (format stream "~&~%  summary, worst dihedral difference per locus:~%")
               (format stream "~&  ~26a ~10@a ~10@a ~14@a ~8a~%"
                       "  locus" "compared" "differing" "worst deg" "  at")
               (dolist (r (nreverse rows))
                 (destructuring-bind (label compared differing worst worst-name) r
                   (format stream "~&  ~26a ~10d ~10d ~14,3f ~8a~a~%"
                           label compared differing worst (or worst-name "-")
                           (if (plusp differing) "   <-" ""))))
               ;; Two calls, not one with a ~A - see the note in CHECK-MONOMER-CONTEXT-AGAINST-POSE.
               (if all-ok
                   (format stream "~&~%internal dihedrals: both sides apply the SAME angles - the ~
                                   rotamer lookup is exonerated, and identical internals over ~
                                   identical joints producing different coordinates puts the ~
                                   fault in the FOLD, not the catalogue~%")
                   (format stream "~&~%internal dihedrals: *** the two sides apply DIFFERENT ~
                                   angles - they resolved different rotamers, so the coordinate ~
                                   difference is a LOOKUP fault and the fold is not ~
                                   implicated ***~%"))))
    all-ok))

;;; ------------------------------------------------------------------
;;; WHICH INTERNALS VECTOR IS THE GROUND TRUTH?
;;;
;;; The rotamers agree exactly - 0 dihedrals differ at all 13 loci - and the coordinates do not.
;;; So something between the catalogue and the coordinates rewrites the angles, and the internals
;;; VECTOR is the only thing in between.  Both sides fill one and then fold it:
;;;
;;;   blueprint   UPDATE-BLUEPRINT-INTERNALS = base fill + slot fill + ADJUST-ALL-INTERNALS
;;;   pose        UPDATE-INTERNALS           = FILL-INTERNALS-FROM-OLIGOMER-SHAPE + ADJUST-INTERNALS
;;;               then BUILD-...-AND-ADJUST also runs ADJUST-ATOM-TREE-EXTERNAL-COORDINATES
;;;
;;; Both adjustment tables - INTERNAL-ADJUSTMENTS and EXTERNAL-ADJUSTMENTS - are keyed by ATRESIDUE
;;; (assembler.lisp:1065,1101).  The blueprint's slot atresidues are objects the FAN-OUT created,
;;; so an adjustment registered against the pose's atresidue does not reach them, and the two sides
;;; can diverge here while agreeing on everything upstream.
;;;
;;; A joint's three internals live at POSITION-INDEX-X3 + 0/1/2 of the internals vector
;;; (bondedJoint.h:95-100) - the SAME index that addresses the coordinate vector.  So the i3x for a
;;; dihedral is POSITION-INDEX-X3 + 2, and it is printed for both sides because the two vectors are
;;; indexed independently and a difference there is itself worth seeing.
;;;
;;; Five comparisons per locus, which between them name the culprit without further guessing:
;;;
;;;   rotamer  vs bp-int     did the BLUEPRINT's adjustments move it off the catalogue value?
;;;   rotamer  vs pose-int   did the POSE's adjustments?
;;;   bp-int   vs bp-meas    does the blueprint's fold honour its own internals?
;;;   pose-int vs pose-meas  does the pose's fold honour its own internals?
;;;   bp-meas  vs pose-meas  the discrepancy being chased, per dihedral
;;;
;;; MEASURED means recomputed from the built coordinates over the joint's own parent chain, so it
;;; is the ground truth by construction: whatever the vectors say, this is what the structure has.
;;; A vector that disagrees with its OWN side's measurement did not survive the fold.  A constant
;;; offset appearing in both MEAS columns equally is a sign convention, not a defect - which is why
;;; both sides are measured the same way rather than one being trusted.
;;; ------------------------------------------------------------------

(defun coords-dihedral-degrees (coords i3 j3 k3 l3)
  "Dihedral in degrees over four i3 offsets into COORDS, via GEOM:CALCULATE-DIHEDRAL-ARRAY.

Cando's own routine, not a local one.  A hand-rolled version here produced every angle with the
opposite sign to BONDED-JOINT/GET-PHI - a whole column of exact negations - and a sign convention
guessed wrong in the other direction would have been just as invisible.

The argument order is the one ASSEMBLER.LISP:667 uses for a joint's phi (joint, parent,
grandparent, great-grandparent), which is also the order EXTRACT-DIHEDRAL-LINE-SEGMENT-FROM-JOINT
returns, so the result is directly comparable to the value in the internals vector.

COORDS is an EXTERNAL COORDINATE NVECTOR indexed by POSITION-INDEX-X3 - not ATOM_O positions.  The
atoms only carry what COPY-ALL-JOINT-POSITIONS-INTO-ATOMS last wrote into them, so they can lag the
vector the fold actually produced, and the internals being compared against are defined against
this indexing anyway."
  (when (and i3 j3 k3 l3 coords)
    (rad-to-degrees (geom:calculate-dihedral-array i3 j3 k3 l3 coords))))

(defun joint-phi-degrees (joint internals)
  "JOINT's dihedral out of INTERNALS in degrees, or NIL when it has none.

Goes through EXTRACT-DIHEDRAL-RAD-FROM-JOINT rather than BONDED-JOINT/GET-PHI directly: a
COMPLEX-BONDED-JOINT's phi can be undefined, and reading the raw slot would report whatever happens
to sit at that index as though it were an angle.  NIL for a JUMP- or XYZ-JOINT, which carry no phi
at all - the two XYZ CB joints in this structure are placed absolutely and read NIL here."
  (when (and internals (typep joint 'kin:bonded-joint))
    (multiple-value-bind (rad defined)
        (topology::extract-dihedral-rad-from-joint joint internals)
      (when defined (rad-to-degrees rad)))))

(defun joint-dihedral-readout (joint internals coords)
  "(values I3X FROM-INTERNALS-DEG MEASURED-DEG) for JOINT.

I3X is where this joint's dihedral lives in the internals vector - POSITION-INDEX-X3 + 2 - and the
measurement uses the same parent chain EXTRACT-DIHEDRAL-LINE-SEGMENT-FROM-JOINT reports, so the two
numbers describe the same angle and are directly comparable."
  (let* ((bondedp (typep joint 'kin:bonded-joint))
         (i3x (when bondedp (+ (kin:joint/position-index-x3 joint) 2)))
         (from-internals (joint-phi-degrees joint internals))
         (segment (when bondedp (topology::extract-dihedral-line-segment-from-joint joint)))
         (measured (when segment
                     (apply #'coords-dihedral-degrees coords segment))))
    (values i3x from-internals measured)))

(defun pose-regenerated-internals (assembler)
  "Rebuild the pose's internals vector the way MOPT-SIDECHAIN does (montecarlo.lisp:42-44):
MAKE-INTERNALS-FOR-ASSEMBLER then UPDATE-INTERNALS per oligomer-shape, which is
FILL-INTERNALS-FROM-OLIGOMER-SHAPE followed by ADJUST-INTERNALS.

Regenerated rather than read because a POSE keeps COORDS and its ASSEMBLER but no internals vector.
That makes this the one number here that is reconstructed, so if it disagrees with POSE-MEAS the
reconstruction is a suspect alongside the fold."
  (let ((internals (topology:make-internals-for-assembler assembler)))
    (loop for oligomer-shape in (topology:oligomer-shapes assembler)
          do (topology:update-internals assembler internals :oligomer-shape oligomer-shape))
    internals))

;;; ------------------------------------------------------------------
;;; The coordinates behind the angles
;;;
;;; POSE-INT matches BP-INT exactly while POSE-MEAS does not match BP-MEAS.  Identical internals
;;; over identical joints cannot fold to different geometry, so one of the two measurements is not
;;; reading what it claims to read.
;;;
;;; The measurement is the suspect part.  It reaches coordinates through the POSE atresidue found
;;; by (SHAPE-INDEX . MONOMER-INDEX), and every check that has validated that mapping so far is
;;; blind to picking the wrong residue OF THE SAME KIND: joint names, joint count, joint types and
;;; the CB<-CA<-N<-C parent chain are identical for every GLU in the structure.  Land on a
;;; different GLU and all of it still passes while the coordinates come from somewhere else.
;;;
;;; So print the positions.  Stage 2 already proved the two backbones are identical, which makes
;;; this decisive: if the pose's CA, N and C are not sitting on the blueprint's CA, N and C, the
;;; pose-side lookup is reading the wrong residue and POSE-MEAS was never comparable.
;;;
;;; NOTE the index used here.  A joint's dihedral is at POSITION-INDEX-X3 + 2 of the INTERNALS
;;; vector - that is the i3x in the table above - but the joint's COORDINATE is at
;;; POSITION-INDEX-X3 itself.  Printing the +2 index against the coordinate vector would read the
;;; z of the previous atom and look almost plausible.
;;; ------------------------------------------------------------------

(defun joint-xyz (joint coords)
  "(values X Y Z I3) for JOINT out of COORDS, or NIL when out of range.

I3 is POSITION-INDEX-X3 with NO +2 - that offset addresses the dihedral in the internals vector,
not the position in the coordinate vector."
  (let ((i3 (kin:joint/position-index-x3 joint)))
    (when (and coords i3 (>= i3 0) (< (+ i3 2) (length coords)))
      (values (aref coords i3) (aref coords (+ i3 1)) (aref coords (+ i3 2)) i3))))

(defun ancestry-joints (joint depth)
  "JOINT and its DEPTH ancestors, root-ward.  NIL entries once the chain ends."
  (loop repeat (1+ depth)
        for j = joint then (when (and j (kin:joint/parent-bound-p j)) (kin:parent j))
        collect j))

(defun print-joint-pair-row (label bp-joint pose-joint bp-coords pose-coords stream)
  "One row: the same joint's position on each side, and the distance between them."
  (multiple-value-bind (bx by bz bi)
      (if bp-joint (joint-xyz bp-joint bp-coords) (values nil nil nil nil))
    (multiple-value-bind (px py pz pidx)
        (if pose-joint (joint-xyz pose-joint pose-coords) (values nil nil nil nil))
      (let ((d (when (and bx px)
                 (sqrt (+ (expt (- bx px) 2) (expt (- by py) 2) (expt (- bz pz) 2))))))
        (format stream "~&      ~16a ~7a ~28a ~7a ~28a ~9a~a~%"
                label (or bi "-")
                (if bx (format nil "(~8,3f ~8,3f ~8,3f)" bx by bz) "-")
                (or pidx "-")
                (if px (format nil "(~8,3f ~8,3f ~8,3f)" px py pz) "-")
                (if d (format nil "~,4f" d) "-")
                (if (and d (> d 1d-4)) "  <-" ""))))))

(defun report-root-frame-and-joints (bp-joints pose-joints bp-coords pose-coords
                                     &key (stream *standard-output*) (limit 10))
  "Print the root joint's parent chain and every joint, as positions, on both sides.

The parent chain runs OUT of the residue into the backbone, which stage 2 proved identical - so
those four rows are the ones that decide whether the two sides are looking at the same place."
  (let ((bp-root (when (plusp (length bp-joints)) (aref bp-joints 0)))
        (pose-root (when (plusp (length pose-joints)) (aref pose-joints 0))))
    (format stream "~&      root frame - the stub the sidechain is placed against ~
                    (backbone, proven identical in stage 2):~%")
    (format stream "~&      ~16a ~7@a ~28@a ~7@a ~28@a ~9@a~%"
            "atom" "bp i3" "blueprint xyz" "pose i3" "pose xyz" "dist")
    (loop for tag in '("root" "parent" "grandparent" "great-grandpa")
          for bj in (ancestry-joints bp-root 3)
          for pj in (ancestry-joints pose-root 3)
          do (print-joint-pair-row
              (format nil "~a ~a" tag (if bj (kin:joint/name bj) "-"))
              bj pj bp-coords pose-coords stream))
    (format stream "~&      sidechain joints:~%")
    (loop for i below (min (length bp-joints) (length pose-joints) limit)
          for bj = (aref bp-joints i)
          do (print-joint-pair-row (format nil "~a" (kin:joint/name bj))
                                   bj (aref pose-joints i) bp-coords pose-coords stream))))

;;; ------------------------------------------------------------------
;;; Are the pose's stored coordinates the output of its own pipeline?
;;;
;;; Every input to CB agrees: same rotamer index, same internals, same BETA-CARBON adjustment, same
;;; CA/N/C stub.  The blueprint's coordinates honour those internals and the pose's do not.
;;;
;;; CB's phi is set by BETA-CARBON from the backbone template alone (c-dih - 120), overwriting the
;;; rotamer's own value - so it does not depend on which rotamer is selected, and no
;;; rotamer-selection or catalogue-index story can explain a CB mismatch.  Two possibilities are
;;; left, and they need opposite fixes:
;;;
;;;   the stored coordinates are not this pipeline's output (built or minimized elsewhere), or
;;;   POSE-REGENERATED-INTERNALS is not what the pose actually folded - a pose keeps no internals
;;;   vector, so that reconstruction is the one number here that is inferred rather than read.
;;;
;;; Folding the pose's own assembler from its own internals separates them.  If the fresh fold
;;; reproduces POSE-INT, the pipeline is faithful and the stored coordinates did not come from it.
;;; If the fresh fold reproduces the STORED coordinates instead, the reconstruction is what is
;;; wrong and the pose is self-consistent.
;;;
;;; Compared by DIHEDRAL, not by position: a fresh fold leaves the ligand in its own local frame -
;;; TRANSFORM-EXTERNALS-TO-GLOBAL-FRAME is a separate step the bench performs with the POSE's
;;; orientation - so positions would differ by a rigid motion that says nothing about the fold.
;;; Dihedrals are invariant under it, so ligand and receptor loci read the same way.
;;; ------------------------------------------------------------------

(defun pose-refolded-coordinates (assembler internals ligand-orientation)
  "Fold ASSEMBLER from INTERNALS into a fresh vector via UPDATE-EXTERNALS.

NOT through BUILD-ALL-ATOM-TREE-EXTERNAL-COORDINATES-AND-ADJUST, the obvious-looking entry point:
it passes OLIGOMER-SHAPE twice to a four-parameter callee (assembler.lisp:1075) and so signals
WRONG-NUMBER-OF-ARGUMENTS on every call.  Nothing has ever folded through it, which also retires
its callers at montecarlo.lisp:27,44,63,82 as dead paths.

UPDATE-EXTERNALS reaches the same build-then-adjust pair by way of
UPDATE-EXTERNALS-FOR-{LIGAND,RECEPTOR}-OLIGOMER-SHAPE (assembler.lisp:1490-1497), and additionally
applies the ligand's orientation transform - the ligand folds in its own local frame, so without it
the ligand lands at the origin ~150 A from the receptor."
  (let ((coords (topology:make-coordinates-for-assembler assembler)))
    (topology:update-externals assembler internals
                               :ligand-orientation ligand-orientation
                               :coords coords)
    coords))

(defun refold-pose-coordinates (pose &key stream)
  "Fresh coordinates for POSE's assembler, folded from the conformation its monomer-shapes NAME.

Silent unless STREAM is supplied, because this is the DEFAULT source of pose coordinates for every
check below and a banner per caller would bury the output.  The bench prints it once.

The pose's stored COORDS were folded before MAKE-ASSEMBLER-FOR-COMPLEX-SCANS replaced the moveable
receptor monomers' shapes with fresh ROTAMER-SHAPEs (protein-binder.lisp:5186), and nothing
re-folded afterwards - so those receptor sidechains sit at their loaded conformation while the
labels say otherwise.  This builds what the labels describe.

Returns a NEW vector; POSE is not modified.  Scoring the blueprint against this compares two
structures that are supposed to be the same one, which is the only comparison that means anything
until the pipeline itself re-folds."
  (let* ((assembler (design:assembler pose))
         (coords (pose-refolded-coordinates assembler
                                            (pose-regenerated-internals assembler)
                                            (design:orientation pose))))
    ;; WHEN, not a NIL stream: (FORMAT NIL ...) builds and discards a string rather than printing
    ;; nothing, which is silent but wasteful on every defaulted call.
    (when stream
      (format stream "~&~%!!!! POSE RE-FOLD FOR SCORING~%~
                      ~&  Every check below reads THIS, not (DESIGN:COORDS POSE): the pose's ~
                        assembler rebuilt from the~%~
                      ~&  rotamers its own monomer-shapes name, so both sides describe the same ~
                        structure.~%"))
    coords))

(defun check-pose-refold (blueprint pose mcstate &key (stream *standard-output*)
                                                      (tolerance 1d-2)
                                                      (per-locus 3))
  "Fold the pose's assembler from its own internals and compare the result against both the pose's
INTERNALS and its STORED coordinates, per joint, in degrees.

Returns T when the fresh fold reproduces the internals."
  (let* ((assembler (design:assembler pose))
         (pose-coords (design:coords pose))
         (pose-atres (pose-monomer-index-to-atresidue assembler))
         (pose-internals (pose-regenerated-internals assembler))
         (refold (pose-refolded-coordinates assembler pose-internals
                                            (design:orientation pose)))
         (rows nil))
    (format stream "~&~%!!!! POSE RE-FOLD: ARE THE STORED COORDINATES THIS PIPELINE'S OUTPUT?~%~
                    ~&  CB's phi comes from BETA-CARBON via the backbone alone, so it cannot ~
                      depend on the rotamer - which rules~%~
                    ~&  out a stale rotamer or a wrong catalogue index as the explanation.  ~
                      Folding the pose's own assembler from~%~
                    ~&  its own internals decides between 'the stored coordinates did not come ~
                      from this pipeline' and 'the~%~
                    ~&  internals I reconstructed are not what it folded'.  Read by dihedral, so ~
                      the ligand's local frame is moot.~%")
    (format stream "~&  ~-24a ~-8a ~10@a ~10@a ~10@a ~11@a ~13@a~%"
            "  locus" "joint" "pose-int" "refold" "stored" "refold-int" "refold-stored")
    (loop for mrkindex across mcstate
          for info = (topology::mrkindex-info blueprint mrkindex)
          do (let ((pose-at (gethash (cons (getf info :shape-index)
                                           (getf info :locus-in-shape))
                                     pose-atres)))
               (when pose-at
                 (let ((worst-ri 0d0) (worst-rs 0d0) (worst-name nil) (shown 0))
                   (loop for joint across (topology:joints pose-at)
                         when (typep joint 'kin:bonded-joint)
                           do (let* ((int (joint-phi-degrees joint pose-internals))
                                     (segment
                                       (topology::extract-dihedral-line-segment-from-joint joint))
                                     (re (apply #'coords-dihedral-degrees refold segment))
                                     (st (apply #'coords-dihedral-degrees pose-coords segment))
                                     (d-ri (degrees-difference re int))
                                     (d-rs (degrees-difference re st)))
                                (when (and d-ri (> (abs d-ri) (abs worst-ri)))
                                  (setf worst-ri d-ri))
                                (when (and d-rs (> (abs d-rs) (abs worst-rs)))
                                  (setf worst-rs d-rs worst-name (kin:joint/name joint)))
                                ;; A few joints per locus is enough to read the pattern - the
                                ;; per-locus worst lines below carry the verdict.
                                (when (< shown per-locus)
                                  (incf shown)
                                  (format stream "~&  ~-24a ~-8a ~10a ~10a ~10a ~11a ~13a~%"
                                          (if (= shown 1)
                                              (list (getf info :blueprint-locus)
                                                    (getf info :monomer-name))
                                              "")
                                          (kin:joint/name joint)
                                          (if int (format nil "~,2f" int) "-")
                                          (if re (format nil "~,2f" re) "-")
                                          (if st (format nil "~,2f" st) "-")
                                          (if d-ri (format nil "~,2f" d-ri) "-")
                                          (if d-rs (format nil "~,2f" d-rs) "-")))))
                   (push (list (list (getf info :blueprint-locus) (getf info :monomer-name))
                               worst-ri worst-rs worst-name)
                         rows)))))
    (setf rows (nreverse rows))
    (format stream "~&~%  worst per locus:~%")
    (format stream "~&  ~-24a ~14@a ~17@a ~-8a~%"
            "  locus" "refold vs int" "refold vs stored" "  at")
    (dolist (r rows)
      (destructuring-bind (label ri rs name) r
        (format stream "~&  ~-24a ~14,2f ~17,2f ~-8a~%" label ri rs (or name "-"))))
    (let ((faithful (notany (lambda (r) (> (abs (second r)) tolerance)) rows))
          (matches-stored (notany (lambda (r) (> (abs (third r)) tolerance)) rows)))
      (cond ((and faithful matches-stored)
             (format stream "~&~%pose re-fold: the fold reproduces BOTH the internals and the ~
                             stored coordinates - there is nothing left to explain on the pose ~
                             side~%"))
            (faithful
             (format stream "~&~%pose re-fold: the fold reproduces the pose's INTERNALS but NOT ~
                             its stored coordinates.  The pipeline is faithful and the internals ~
                             reconstruction is sound, so (DESIGN:COORDS POSE) did not come from ~
                             it - those coordinates carry geometry from somewhere else~%"))
            (matches-stored
             (format stream "~&~%pose re-fold: the fold reproduces the STORED coordinates but not ~
                             the internals I regenerated - the pose is self-consistent and ~
                             POSE-REGENERATED-INTERNALS is the wrong reconstruction, so every ~
                             POSE-INT column above should be re-read~%"))
            (t
             (format stream "~&~%pose re-fold: the fold matches NEITHER - it disagrees with the ~
                             reconstructed internals and with the stored coordinates ~
                             separately~%")))
      faithful)))

;;; ------------------------------------------------------------------
;;; Is the adjustment registered, initialized, and live - on each side?
;;;
;;; The angles say the pose adjusts the spiro root dihedral and the blueprint does not.  That is an
;;; observation about outcomes; this asks the mechanism directly, because there are three distinct
;;; ways for an adjustment to do nothing and they need different fixes:
;;;
;;;   1. NOT REGISTERED    no entry under this atresidue in INTERNAL-ADJUSTMENTS at all.
;;;   2. REGISTERED, DEAD  an entry exists but INITIALIZE-ADJUSTMENT never bound OTHER, and
;;;                        INTERNAL-ADJUST is a no-op unless (SLOT-BOUNDP adjustment 'OTHER)
;;;                        (joint-templates.lisp:72).
;;;   3. LIVE              bound, and the phi really is rewritten.
;;;
;;; Case 2 is the one to expect here.  INITIALIZE-ADJUSTMENT binds OTHER only when
;;; (= (KIN:NUMBER-OF-CHILDREN JPARENT) 2) - the adjustment means \"sit 180 degrees from my only
;;; sibling\", which is meaningless without exactly one sibling.  The blueprint's fan-out hangs the
;;; root joint of EVERY SLOT off the same backbone parent, so that parent has two children plus one
;;; per slot, the test fails, OTHER stays unbound, and the adjustment silently does nothing.  The
;;; pose has one residue per locus, so its parent keeps two children and the adjustment lives.
;;;
;;; PARENT CHILDREN is therefore the column that matters, and it is why this is worth printing
;;; rather than reasoning about: it distinguishes \"the blueprint never got the adjustment\" from
;;; \"the blueprint got it and the fan-out disabled it\", which are different bugs.
;;; ------------------------------------------------------------------

(defun adjustment-parent-note (adjustment)
  "How many children the adjusted joint's parent has.

Not part of the adjustment's own description - it is JOINT TREE state, and it is precisely what the
fan-out changes, so it belongs beside the adjustment rather than inside it.  Reaches only JOINT,
the one slot every ADJUSTMENT has by definition."
  (let* ((joint (topology::joint adjustment))
         (parent (when (and joint (kin:joint/parent-bound-p joint)) (kin:parent joint))))
    (if parent
        (format nil "   parent ~a has ~d children"
                (kin:joint/name parent) (kin:number-of-children parent))
        "")))

(defun adjustment-signature (adjustment)
  "(CLASS JOINT-NAME OTHER) for comparing one side's adjustment against the other's.

Reaches into slots, where REPORT-ATRESIDUE-ADJUSTMENTS deliberately does not.  The distinction is
DESCRIBING versus COMPARING: printing is the class's own job and PRINT-OBJECT does it, but its
output embeds joint object addresses, which differ between two separately-built trees and so cannot
be compared.  This reduces to names.

OTHER is :N/A for a class without the slot (BETA-CARBON binds JC and JHA instead and is not gated
on OTHER), :UNBOUND when INITIALIZE-ADJUSTMENT declined to bind it - which is the state that made
the fan-out bug silent."
  (list (class-name (class-of adjustment))
        (string (kin:joint/name (topology::joint adjustment)))
        (cond ((not (slot-exists-p adjustment 'topology::other)) :n/a)
              ((not (slot-boundp adjustment 'topology::other)) :unbound)
              (t (let ((o (topology::other adjustment)))
                   (if o (string (kin:joint/name o)) :nil))))))

(defun atresidue-adjustment-signatures (atresidue adjustments)
  "Sorted signatures of every adjustment registered against ATRESIDUE, internal and external.

Sorted so the comparison does not depend on PUSH order - ADD-TO-ADJUSTMENTS pushes, so two sides
that registered the same adjustments in a different order would otherwise look different."
  (when atresidue
    (sort (mapcar #'adjustment-signature
                  (append (gethash atresidue (topology::internal-adjustments adjustments))
                          (gethash atresidue (topology::external-adjustments adjustments))))
          #'string< :key (lambda (sig) (format nil "~s" sig)))))

(defun report-atresidue-adjustments (label atresidue adjustments stream)
  "Print every adjustment registered against ATRESIDUE, letting each one describe itself.

PRINT-OBJECT rather than slot reads: only the class knows which of its slots decides whether it
fires.  Reading OTHER from here reported BETA-CARBON - which has no OTHER and is gated on nothing -
as DEAD, and that mislabelling sent me looking for a second broken adjustment that did not exist."
  (let ((internal (when atresidue
                    (gethash atresidue (topology::internal-adjustments adjustments))))
        (external (when atresidue
                    (gethash atresidue (topology::external-adjustments adjustments)))))
    (if (and (null internal) (null external))
        (format stream "~&      ~10a none registered~%" label)
        (progn
          (dolist (a internal)
            (format stream "~&      ~10a internal ~a~a~%" label a (adjustment-parent-note a)))
          (dolist (a external)
            (format stream "~&      ~10a external ~a~a~%" label a (adjustment-parent-note a)))))))

(defun check-slot-adjustments-against-pose (blueprint pose mcstate
                                            &key (stream *standard-output*))
  "For every selected locus, report the adjustments registered on the blueprint's slot atresidue
and on the pose's atresidue, with enough state to say whether each one actually fires.

Returns T when both sides resolve the same adjustments to the same partner joints."
  (let* ((assembler (design:assembler pose))
         (pose-atres (pose-monomer-index-to-atresidue assembler))
         (bp-adjustments (topology::adjustments blueprint))
         (pose-adjustments (topology::adjustments assembler))
         (all-match t))
    (format stream "~&~%!!!! ADJUSTMENTS: REGISTERED, INITIALIZED, OR DEAD - ON EACH SIDE?~%~
                    ~&  The angles say the pose adjusts the root dihedral and the blueprint does ~
                      not.  This asks WHY, and there~%~
                    ~&  are three different answers with three different fixes: not registered, ~
                      registered but never~%~
                    ~&  initialized (OTHER unbound, so INTERNAL-ADJUST returns immediately), or ~
                      live.~%~
                    ~&  PARENT CHILDREN is the column to watch: INITIALIZE-ADJUSTMENT binds OTHER ~
                      only when the parent has~%~
                    ~&  exactly 2 children, and the fan-out hangs every slot's root joint off one ~
                      shared backbone parent.~%")
    (loop for mrkindex across mcstate
          for info = (topology::mrkindex-info blueprint mrkindex)
          do (let* ((scan (aref (topology::rotamer-scans blueprint) mrkindex))
                    (bp-at (when scan (topology::atresidue scan)))
                    (pose-at (gethash (cons (getf info :shape-index)
                                            (getf info :locus-in-shape))
                                      pose-atres)))
               (format stream "~&~%    locus ~a ~a~%"
                       (getf info :blueprint-locus) (getf info :monomer-name))
               (report-atresidue-adjustments "blueprint" bp-at bp-adjustments stream)
               (report-atresidue-adjustments "pose" pose-at pose-adjustments stream)
               ;; THE VERDICT, separate from the printing above.  The two sides must resolve the
               ;; same adjustments to the same partners; the fan-out bug showed up exactly here, as
               ;; blueprint OTHER unbound against pose OTHER = CAC.
               (let ((bp-sig (atresidue-adjustment-signatures bp-at bp-adjustments))
                     (pose-sig (atresidue-adjustment-signatures pose-at pose-adjustments)))
                 (unless (equal bp-sig pose-sig)
                   (setf all-match nil)
                   (format stream "~&      *** the two sides resolve DIFFERENT adjustments ***~%~
                                   ~&        blueprint ~s~%~
                                   ~&        pose      ~s~%"
                           bp-sig pose-sig)))))
    (if all-match
        (format stream "~&~%adjustments: both sides resolve the same adjustments to the same ~
                        partner joints~%")
        (format stream "~&~%adjustments: *** the two sides DIFFER - an adjustment that resolves on ~
                        one side and not the other silently changes that joint's geometry ***~%"))
    all-match))

(defun check-internals-vectors-against-pose (blueprint persona pose mcstate bp-coords bp-internals
                                             &key (stream *standard-output*)
                                                  (tolerance 1d-2)
                                                  (limit 20)
                                                  (pose-coords
                                                   (refold-pose-coordinates pose)))
  "Compare, per joint, the dihedral in each side's INTERNALS VECTOR against the catalogue value and
against the angle measured from that side's own built coordinates.

Returns T when every internals vector agrees with its own coordinates AND the two sides agree."
  (let* ((assembler (design:assembler pose))
         (pose-atres (pose-monomer-index-to-atresidue assembler))
         (pose-internals (pose-regenerated-internals assembler))
         (rows nil))
    (format stream "~&~%!!!! INTERNALS VECTORS: WHICH ONE IS THE GROUND TRUTH?~%~
                    ~&  The rotamers agree exactly and the coordinates do not, so the rewrite ~
                      happens in the internals VECTOR - the~%~
                    ~&  only stage in between.  A joint's dihedral lives at POSITION-INDEX-X3 + 2 ~
                      of that vector, the same index that~%~
                    ~&  addresses the coordinates, and MEAS is recomputed from the built structure ~
                      with GEOM:CALCULATE-DIHEDRAL-ARRAY,~%~
                    ~&  so it is true by construction and in the same convention as ~
                      BONDED-JOINT/GET-PHI.  A vector disagreeing with~%~
                    ~&  its OWN side's MEAS did not survive the fold - those coordinates are not ~
                      what that vector describes.~%")
    (loop for mrkindex across mcstate
          for info = (topology::mrkindex-info blueprint mrkindex)
          do (multiple-value-bind (bp-locus bp-monomer picked-slot)
                 (topology::blueprint-locus-of-mrkindex blueprint mrkindex)
               (declare (ignore bp-monomer))
               (let* ((scan (aref (topology::rotamer-scans blueprint) mrkindex))
                      (bp-joints (topology:joints (topology::atresidue scan)))
                      (pose-at (gethash (cons (getf info :shape-index)
                                              (getf info :locus-in-shape))
                                        pose-atres))
                      (bp-rotamer (topology::slot-rotamer-internals persona mrkindex))
                      (cat (rotamer-joint-internals bp-rotamer bp-joints))
                      (shown 0)
                      (n 0)
                      ;; The ROOT joint carries the sidechain's attachment dihedral and is the one
                      ;; joint ADJUST-INTERNALS rewrites - every other joint reads back its
                      ;; catalogue value exactly.  Kept apart because a 112 deg adjustment and a
                      ;; fold error in one column hide each other: the adjustment is EXPECTED and
                      ;; only interesting when the two sides disagree about it, while any fold
                      ;; difference at all is a defect.
                      (adjusted nil)
                      (n-rest 0)
                      (worst-bp-fold 0d0) (worst-pose-fold 0d0)
                      (worst-meas 0d0) (worst-meas-name nil))
                 (format stream "~&~%  locus ~a ~a   slot ~a~%"
                         (getf info :blueprint-locus) (getf info :monomer-name) picked-slot)
                 (if (null pose-at)
                     (format stream "~&    no pose atresidue - nothing to compare~%")
                     (let ((pose-joints (topology:joints pose-at)))
                       (format stream "~&    ~8a ~7@a ~7@a ~9@a ~9@a ~9@a ~9@a ~9@a~%"
                               "joint" "bp i3x" "pos i3x" "rotamer" "bp-int" "pose-int"
                               "bp-meas" "pose-meas")
                       (loop for i below (min (length bp-joints) (length pose-joints))
                             for bp-joint = (aref bp-joints i)
                             for pose-joint = (aref pose-joints i)
                             for catrow = (nth i cat)
                             do (multiple-value-bind (bp-i3x bp-int bp-meas)
                                    (joint-dihedral-readout bp-joint bp-internals bp-coords)
                                  (multiple-value-bind (pose-i3x pose-int pose-meas)
                                      (joint-dihedral-readout pose-joint pose-internals pose-coords)
                                    (when (and catrow (second catrow))
                                      (let* ((cat-deg (fifth catrow))
                                             (d-bp-adjust (degrees-difference bp-int cat-deg))
                                             (d-pose-adjust (degrees-difference pose-int cat-deg))
                                             (d-bp-fold (degrees-difference bp-meas bp-int))
                                             (d-pose-fold (degrees-difference pose-meas pose-int))
                                             (d-meas (degrees-difference bp-meas pose-meas))
                                             ;; DETECTED, not assumed to be joint 0: a joint whose
                                             ;; internals no longer equal its catalogue value is
                                             ;; one something rewrote.  Assuming the root would be
                                             ;; wrong the moment a second adjustment is registered.
                                             (adjustedp
                                               (or (and d-bp-adjust
                                                        (> (abs d-bp-adjust) tolerance))
                                                   (and d-pose-adjust
                                                        (> (abs d-pose-adjust) tolerance)))))
                                        (incf n)
                                        (if adjustedp
                                            (push (list (first catrow) d-bp-adjust d-pose-adjust
                                                        (degrees-difference bp-int pose-int)
                                                        d-bp-fold d-pose-fold d-meas)
                                                  adjusted)
                                            (progn
                                              (incf n-rest)
                                              (macrolet ((track (place value)
                                                           `(when (and ,value
                                                                       (> (abs ,value)
                                                                          (abs ,place)))
                                                              (setf ,place ,value))))
                                                (track worst-bp-fold d-bp-fold)
                                                (track worst-pose-fold d-pose-fold))
                                              (when (and d-meas
                                                         (> (abs d-meas) (abs worst-meas)))
                                                (setf worst-meas d-meas
                                                      worst-meas-name (first catrow)))))
                                        ;; EVERY joint prints, up to LIMIT.  Printing only the
                                        ;; disagreeing ones made a locus where nothing downstream
                                        ;; differs look like a locus whose downstream joints had
                                        ;; gone missing - the agreement is the result, so it has to
                                        ;; be on the page.
                                        (when (or adjustedp (< shown limit))
                                          (incf shown)
                                          (format stream "~&    ~8a ~7a ~7a ~9,2f ~9@a ~9@a ~
                                                          ~9@a ~9@a~a~%"
                                                  (first catrow) (or bp-i3x "-") (or pose-i3x "-")
                                                  cat-deg
                                                  (if bp-int (format nil "~,2f" bp-int) "-")
                                                  (if pose-int (format nil "~,2f" pose-int) "-")
                                                  (if bp-meas (format nil "~,2f" bp-meas) "-")
                                                  (if pose-meas
                                                      (format nil "~,2f" pose-meas) "-")
                                                  (cond (adjustedp "  * adjusted")
                                                        ((and d-meas
                                                              (> (abs d-meas) tolerance))
                                                         "  <- differs")
                                                        (t "")))))))))
                       (format stream "~&    ~d dihedrals: ~d adjusted, ~d unadjusted, ~d listed~%"
                               n (length adjusted) n-rest shown)
                       ;; The adjusted joints, one line each.  BP-VS-POSE on the INTERNALS is the
                       ;; question for these: both sides are supposed to apply the same adjustment,
                       ;; and the size of the adjustment itself is not a defect.
                       (dolist (a (reverse adjusted))
                         (destructuring-bind (name bpa posea intdiff bpf posef meas) a
                           (declare (ignore bpf posef))
                           (format stream "~&      adjusted ~6a bp ~8,2f  pose ~8,2f  ~
                                           internals differ ~8,2f~a   (measured differ ~,2f)~%"
                                   name (or bpa 0d0) (or posea 0d0) (or intdiff 0d0)
                                   (if (and intdiff (> (abs intdiff) tolerance))
                                       "  <- DIFFERENT ADJUSTMENT" "")
                                   (or meas 0d0))))
                       (format stream "~&      unadjusted: worst fold bp ~,2f / pose ~,2f   ~
                                       bp-vs-pose ~,2f~@[ at ~a~]~%"
                               worst-bp-fold worst-pose-fold worst-meas worst-meas-name)
                       ;; The positions the angles above were computed from.  Identical internals
                       ;; over identical joints cannot fold differently, so if the two sides really
                       ;; do disagree, one of them is reading coordinates from somewhere else -
                       ;; and the backbone rows say which.
                       (report-root-frame-and-joints bp-joints pose-joints bp-coords pose-coords
                                                     :stream stream)
                       (flet ((worst-abs (key)
                                (reduce #'max
                                        (mapcar (lambda (a) (abs (or (nth key a) 0d0))) adjusted)
                                        :initial-value 0d0)))
                         (push (list (list (getf info :blueprint-locus)
                                           (getf info :monomer-name))
                                     (length adjusted)
                                     (worst-abs 1) (worst-abs 2) (worst-abs 3)
                                     worst-bp-fold worst-pose-fold worst-meas)
                               rows)))))))
    ;; The summary is the point: five columns, and which of them is non-zero names the stage.
    (setf rows (nreverse rows))
    (format stream "~&~%  ADJUSTED joints and the REST, kept apart - worst per locus, degrees:~%")
    (format stream "~&  ~24a ~5@a ~8@a ~9@a ~11@a ~9@a ~10@a ~12@a~%"
            "  locus" "nadj" "bp adj" "pose adj" "adj differ" "fold bp" "fold pose" "rest bp-pose")
    (dolist (r rows)
      (destructuring-bind (label nadj bpa posea adjdiff bpf posef meas) r
        (format stream "~&  ~24a ~5d ~8,2f ~9,2f ~11,2f ~9,2f ~10,2f ~12,2f~a~%"
                label nadj bpa posea adjdiff bpf posef meas
                (cond ((> adjdiff tolerance) "  <- adjustment differs")
                      ((or (> (abs bpf) tolerance) (> (abs posef) tolerance)) "  <- fold")
                      ((> (abs meas) tolerance) "  <- rest")
                      (t "")))))
    (format stream "~&~%  how to read it:~%~
                    ~&    BP ADJ / POSE ADJ  how far ADJUST-INTERNALS moved the root dihedral off ~
                      its catalogue value.~%~
                    ~&                       Large is normal and not a defect on its own - the ~
                      adjustment is the point.~%~
                    ~&    ADJ DIFFER         the two sides' ADJUSTED internals against each other. ~
                      This is the one that~%~
                    ~&                       matters for the root joint: both sides should apply ~
                      the SAME adjustment,~%~
                    ~&                       and the adjustment tables being keyed by ATRESIDUE ~
                      is why they might not.~%~
                    ~&    FOLD BP / POSE     that side's coordinates against its OWN internals ~
                      vector.  Non-zero means~%~
                    ~&                       the vector is not what got built.  Note POSE-INT is ~
                      regenerated, so a pose~%~
                    ~&                       fold difference indicts the regeneration as much as ~
                      the fold.~%~
                    ~&    REST BP-POSE       the two structures against each other over the ~
                      UNADJUSTED joints - the~%~
                    ~&                       joints both sides agree the catalogue fully ~
                      determines.~%")
    (notany (lambda (r) (or (> (abs (fifth r)) tolerance)      ; adjustments disagree
                            (> (abs (sixth r)) tolerance)      ; blueprint fold
                            (> (abs (seventh r)) tolerance)    ; pose fold
                            (> (abs (eighth r)) tolerance)))   ; unadjusted joints disagree
            rows)))
