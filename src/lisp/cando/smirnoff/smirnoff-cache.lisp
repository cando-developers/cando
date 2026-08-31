;;;; smirnoff-cache.lisp
;;;;
;;;; Build a lookup cache of SMIRNOFF force-field parameters keyed by
;;;; (atom-name . constitution-context) tuples, so that bonded-parameter
;;;; assignment can become a table lookup instead of per-SMIRKS-rule VF2
;;;; subgraph-isomorphism matching.
;;;;
;;;; Each of the four interaction types (bond, angle, dihedral, improper) gets
;;;; its own hash-table.  A key is a canonically ordered list of
;;;; (atom-name . constitution-context) conses - one per atom in the
;;;; interaction - and the value is a parameter object.
;;;;
;;;; Proper and improper torsions share the dihedral energy component and are
;;;; separated by the term's PROPER boolean flag.
;;;;
;;;; C++ PREREQUISITES (the bond/nonbond side already exists; add the rest as
;;;; small mirrors of chem:walk-stretch-terms):
;;;;   - chem:walk-stretch-terms                     [EXISTS]
;;;;       callback: (index atom1 atom2 i1 i2 kb r0)
;;;;   - chem:energy-function/get-stretch-component   (mirror of the existing
;;;;     chem:energy-function/get-nonbond-component)
;;;;   - chem:energy-function/get-angle-component
;;;;   - chem:energy-function/get-dihedral-component
;;;;   - chem:walk-angle-terms                        [ADD]
;;;;       callback: (index atom1 atom2 atom3 i1 i2 i3 kt t0)
;;;;   - chem:walk-dihedral-terms                     [ADD]
;;;;       callback: (index atom1 atom2 atom3 atom4 i1 i2 i3 i4 v n phase proper)

(in-package :smirnoff)

;;; ---------------------------------------------------------------------------
;;; Parameter objects - one class per interaction type
;;; ---------------------------------------------------------------------------

(defclass bond-parameters ()
  ((kb :initarg :kb :reader kb :documentation "Bond force constant.")
   (r0 :initarg :r0 :reader r0 :documentation "Equilibrium bond length."))
  (:documentation "Force-field parameters for a bond (stretch) interaction."))

(defclass angle-parameters ()
  ((kt :initarg :kt :reader kt :documentation "Angle force constant.")
   (t0 :initarg :t0 :reader t0 :documentation "Equilibrium angle (radians)."))
  (:documentation "Force-field parameters for an angle interaction."))

(defclass fourier-term ()
  ((v           :initarg :v           :reader v           :documentation "Barrier amplitude.")
   (periodicity :initarg :periodicity :reader periodicity :documentation "Periodicity / multiplicity n.")
   (phase-rad   :initarg :phase-rad   :reader phase-rad   :documentation "Phase in radians."))
  (:documentation "One Fourier component of a torsion potential."))

(defclass torsion-parameters ()
  ((fourier-terms :initarg :fourier-terms :initform nil :accessor fourier-terms
                  :documentation "List of FOURIER-TERM - a torsion may have several."))
  (:documentation "Base class for proper and improper torsion parameters."))

(defclass dihedral-parameters (torsion-parameters) ()
  (:documentation "Force-field parameters for a proper dihedral interaction."))

(defclass improper-parameters (torsion-parameters) ()
  (:documentation "Force-field parameters for an improper torsion interaction."))

(defclass nonbond-parameters ()
  ((vdw-type :initarg :vdw-type :reader vdw-type
             :documentation "The vdw/nonbond atom-type symbol assigned by ASSIGN-FORCE-FIELD-TYPES.")
   (charge   :initarg :charge   :reader charge
             :documentation "The partial atomic charge in elementary-charge units, read from the atom."))
  (:documentation "Cached nonbonded parameters for a single atom: its vdw atom type (which
indexes the ffnonbond-db for radius/epsilon at energy-evaluation time) and its partial charge.
Keyed by one (atom-name . constitution-context) part rather than a tuple, since nonbond
interactions are per-atom."))

;;; ---------------------------------------------------------------------------
;;; The cache - four hash-tables, one per interaction type
;;; ---------------------------------------------------------------------------

(defclass parameter-cache ()
  ((bond-table     :initarg :bond-table     :initform (make-hash-table :test 'equal) :accessor bond-table)
   (angle-table    :initarg :angle-table    :initform (make-hash-table :test 'equal) :accessor angle-table)
   (dihedral-table :initarg :dihedral-table :initform (make-hash-table :test 'equal) :accessor dihedral-table)
   (improper-table :initarg :improper-table :initform (make-hash-table :test 'equal) :accessor improper-table)
   (nonbond-table  :initarg :nonbond-table  :initform (make-hash-table :test 'equal) :accessor nonbond-table)
   (trained-context-table
    :initarg :trained-context-table
    :initform (make-hash-table :test 'equal)
    :accessor trained-context-table
    :documentation "Constitution contexts whose chemically real trainer completed successfully.")
   (trained-oligomer-table
    :initarg :trained-oligomer-table
    :initform (make-hash-table :test 'equal)
    :accessor trained-oligomer-table
    :documentation "Chemical/focus signatures of completely harvested training oligomers.")
   (parameter-provenance-table
    :initarg :parameter-provenance-table
    :initform (make-hash-table :test 'equal)
    :accessor parameter-provenance-table
    :documentation "Maps (KIND . KEY) to the trainer that first supplied that parameter." )
   )
  (:documentation
   "A cache of SMIRNOFF force-field parameters.  The four bonded hash-tables
map a canonically ordered key - a list of (atom-name . constitution-context)
conses - to a parameter object (BOND-PARAMETERS / ANGLE-PARAMETERS /
DIHEDRAL-PARAMETERS / IMPROPER-PARAMETERS).  The NONBOND-TABLE maps a single
(atom-name . constitution-context) part to a NONBOND-PARAMETERS (vdw type +
charge).  All tables use an EQUAL test so the cons/list keys compare
structurally."))

;;; ---------------------------------------------------------------------------
;;; Canonical key ordering
;;;
;;; A key part is a (atom-name . constitution-context) cons.  The key for an
;;; interaction is the list of its atoms' parts, ordered canonically so the
;;; first part is <= the last.  We take the lexicographically smaller of the
;;; part list and its reverse; for symmetric interactions (bond, angle, proper
;;; dihedral - each equal to its own reverse) this yields ONE key regardless of
;;; the order the atoms were presented in.
;;; ---------------------------------------------------------------------------

(defun key-part< (a b)
  "Total order on (atom-name . constitution-context) key parts: compare the
  atom-name, then (on a tie) the constitution-context, as string designators.
  No allocation - the old version concatenated a fresh comparison string per call,
  which dominated the cached-path runtime (concatenate-into-sequence)."
  (let ((name-a (car a))
        (name-b (car b)))
    (if (eq name-a name-b) 
        (and (string< (cdr a) (cdr b)) t) ; same name -> break the tie on context
        (and (string< name-a name-b) t))))   ; string< accepts symbols/keywords directly


(defun key-list< (a b)
  "Lexicographic < on two equal-length lists of key parts."
  (loop for pa in a
        for pb in b
        do (cond ((key-part< pa pb) (return t))
                 ((key-part< pb pa) (return nil)))
        finally (return nil)))

(defun canonicalize-key (parts)
  "Return PARTS in canonical order: the lexicographically smaller of PARTS and its
  reverse, so the first key part is <= the last.  The common 2- and 3-part keys (bonds,
  angles) are decided in place - nothing is consed unless the reversed order is the one
  returned - while 4+ part keys (dihedrals, impropers) take the general path."
  (let ((rest (cdr parts)))
    (cond
      ;; 2 parts (bonds): reverse = (b a); smaller iff b < a
      ((null (cdr rest))
       (let ((a (car parts)) (b (car rest)))
         (if (key-part< b a) (list b a) parts)))
      ;; 3 parts (angles): reverse = (c b a); the center is shared, so it's decided by the ends
      ((null (cddr rest))
       (let ((a (car parts)) (c (cadr rest)))
         (if (key-part< c a) (list c (car rest) a) parts)))
      ;; 4+ parts (dihedrals/impropers): general path
      (t (let ((reversed (reverse parts)))
           (if (key-list< reversed parts) reversed parts))))))

;;; ---------------------------------------------------------------------------
;;; Fourier-term accumulation for torsions
;;;
;;; The dihedral energy component stores one term per Fourier component, so a
;;; single torsion is several terms sharing the same four atoms (hence the same
;;; key).  We accumulate them into the torsion's FOURIER-TERMS list, deduping by
;;; value so that recurring instances of the same torsion type don't duplicate.
;;; ---------------------------------------------------------------------------

(defun fourier-term-equal (a b)
  (and (= (v a) (v b))
       (eql (periodicity a) (periodicity b))
       (= (phase-rad a) (phase-rad b))))

(defun add-fourier-term (torsion-parameters fourier-term)
  "Add FOURIER-TERM to TORSION-PARAMETERS unless an equal one is already present.
Returns TORSION-PARAMETERS."
  (unless (member fourier-term (fourier-terms torsion-parameters)
                  :test #'fourier-term-equal)
    (push fourier-term (fourier-terms torsion-parameters)))
  torsion-parameters)

(defun torsion-parameters-equal (a b)
  (let ((aterms (fourier-terms a))
        (bterms (fourier-terms b)))
    (and (= (length aterms) (length bterms))
         (every (lambda (term)
                  (member term bterms :test #'fourier-term-equal))
                aterms))))

(defun parameter-entry-equal (a b)
  (typecase a
    (bond-parameters
     (and (typep b 'bond-parameters)
          (= (kb a) (kb b))
          (= (r0 a) (r0 b))))
    (angle-parameters
     (and (typep b 'angle-parameters)
          (= (kt a) (kt b))
          (= (t0 a) (t0 b))))
    (torsion-parameters
     (and (typep b (class-name (class-of a)))
          (torsion-parameters-equal a b)))
    (nonbond-parameters
     (and (typep b 'nonbond-parameters)
          (equal (vdw-type a) (vdw-type b))
          (= (charge a) (charge b))))
    (t (equalp a b))))

(defun parameter-entry-description (entry)
  (typecase entry
    (bond-parameters
     (list :kb (kb entry) :r0 (r0 entry)))
    (angle-parameters
     (list :kt (kt entry) :t0 (t0 entry)))
    (torsion-parameters
     (sort (mapcar (lambda (term)
                     (list :v (v term)
                           :periodicity (periodicity term)
                           :phase-rad (phase-rad term)))
                   (fourier-terms entry))
           #'string< :key #'prin1-to-string))
    (nonbond-parameters
     (list :vdw-type (vdw-type entry) :charge (charge entry)))
    (t entry)))

(defun install-parameter-cache-entry (cache kind table key value provenance)
  "Install VALUE unless KEY already has an identical value; signal on a real conflict."
  (multiple-value-bind (existing foundp) (gethash key table)
    (cond
      ((not foundp)
       (setf (gethash key table) value)
       (when provenance
         (setf (gethash (cons kind key) (parameter-provenance-table cache)) provenance)))
      ((parameter-entry-equal existing value)
       existing)
      (t
       (error "Conflicting SMIRNOFF parameters for ~s key ~s.~%~
               Existing trainer: ~s~%New trainer: ~s~%~
               Existing parameters: ~s~%New parameters: ~s"
              kind key
              (gethash (cons kind key) (parameter-provenance-table cache))
              provenance
              (parameter-entry-description existing)
              (parameter-entry-description value))))))

;;; ---------------------------------------------------------------------------
;;; Build the cache from an energy function
;;; ---------------------------------------------------------------------------

(defun expand-parameter-cache (energy-function molecule atom-to-constitution-context atom-types
                               &optional cache &key provenance)
  "Expand an existing parameter-cache or build a new one and return a PARAMETER-CACHE from ENERGY-FUNCTION.

ATOM-TYPES is the atom -> vdw-type-symbol hash-table produced by
CHEM:ASSIGN-FORCE-FIELD-TYPES (the same hash handed to
ENERGY-FUNCTION/DEFINE-FOR-MATTER-WITH-ATOM-TYPES).  It is required: the vdw type
symbol lives ONLY in this hash - CHEM:GET-TYPE is just a gethash into it and the
atom table stores only the integer type index - so without it the nonbond
parameters cannot be recovered.  The nonbond-table is populated with one
NONBOND-PARAMETERS per typed atom, keyed by that atom's single
(atom-name . constitution-context) part; the vdw type is the ATOM-TYPES value and
the charge is read from the atom with CHEM:GET-CHARGE (the same value
EnergyAtom::defineForAtom stores).

ATOM-TO-CONSTITUTION-CONTEXT is a hash-table mapping each atom to its
constitution-context.  For every bond, angle, dihedral and improper term in
ENERGY-FUNCTION we form a key - a canonically ordered list of
(atom-name . constitution-context) conses, one per atom in the interaction -
and store the parameters for that interaction under it.  Proper and improper
torsions come from the same dihedral component and are split by the term's
PROPER flag.
MOLECULE is the single molecule this harvest describes.  Both ENERGY-FUNCTION's components
and ATOM-TYPES span the WHOLE aggregate, so every walk below has to be restricted to
MOLECULE's own atoms.  This matters because defineForAggregate types EVERY molecule (pass 1)
before it generates ANY bonded terms (pass 2): by the time the first molecule is harvested,
ATOM-TYPES already holds the receptor's atoms, and a receptor atom has no :given-atom-type
and therefore no constitution-context.  An atom that IS in MOLECULE but has no context is
still an error - that is the check that catches missing :given-atom-type.

Provide an existing cache or NIL if you don't have one yet."
  (let ((cache (if cache
                   cache
                   (make-instance 'parameter-cache)))
        (mine (make-hash-table :test #'eq)))
    (chem:map-atoms nil (lambda (a) (setf (gethash a mine) t)) molecule)
    (flet ((minep (&rest atoms) (every (lambda (a) (gethash a mine)) atoms))
           (part (atom)
             (cons (chem:get-name atom)
                   (or (gethash atom atom-to-constitution-context)
                       (error "No constitution-context for atom ~s" atom)))))
      ;; --- bonds ---
      (let ((stretch (chem:energy-function/get-stretch-component energy-function)))
        (when stretch
          (chem:walk-stretch-terms
           stretch
           (lambda (index a1 a2 i1 i2 kb r0)
             (declare (ignore index i1 i2))
             (when (minep a1 a2)
               (install-parameter-cache-entry
                cache :bond (bond-table cache)
                (canonicalize-key (list (part a1) (part a2)))
                (make-instance 'bond-parameters :kb kb :r0 r0)
                provenance))))))
      ;; --- angles ---
      (let ((angle (chem:energy-function/get-angle-component energy-function)))
        (when angle
          (chem:walk-angle-terms
           angle
           (lambda (index a1 a2 a3 i1 i2 i3 kt t0)
             (declare (ignore index i1 i2 i3))
             (when (minep a1 a2 a3)
               (install-parameter-cache-entry
                cache :angle (angle-table cache)
                (canonicalize-key (list (part a1) (part a2) (part a3)))
                (make-instance 'angle-parameters :kt kt :t0 t0)
                provenance))))))
      ;; --- dihedrals + impropers (one component, split by PROPER) ---
      (let ((dihedral (chem:energy-function/get-dihedral-component energy-function)))
        (when dihedral
          ;; Accumulate every Fourier component for one physical torsion before comparing its
          ;; context key.  Comparing term-by-term would mistake the second legitimate Fourier
          ;; component for a conflict with the first.
          (let ((physical-torsions (make-hash-table :test 'equal)))
            (chem:walk-dihedral-terms
             dihedral
             (lambda (index a1 a2 a3 a4 i1 i2 i3 i4 v n phase proper)
               (declare (ignore index i1 i2 i3 i4))
               (when (minep a1 a2 a3 a4)
                 (let* ((physical-key (list proper a1 a2 a3 a4))
                        (cache-key
                          (canonicalize-key
                           (list (part a1) (part a2) (part a3) (part a4))))
                        (entry (gethash physical-key physical-torsions))
                        (parameters
                          (or (cdr entry)
                              (make-instance
                               (if proper 'dihedral-parameters 'improper-parameters)))))
                   (when (and entry (not (equal cache-key (car entry))))
                     (error "One physical torsion produced two cache keys: ~s and ~s"
                            (car entry) cache-key))
                   (add-fourier-term
                    parameters
                    (make-instance 'fourier-term
                                   :v v :periodicity n :phase-rad phase))
                   (setf (gethash physical-key physical-torsions)
                         (cons cache-key parameters))))))
            (maphash
             (lambda (physical-key entry)
               (let* ((proper (first physical-key))
                      (kind (if proper :dihedral :improper))
                      (table (if proper (dihedral-table cache) (improper-table cache))))
                 (install-parameter-cache-entry
                  cache kind table (car entry) (cdr entry) provenance)))
             physical-torsions))))
      ;; --- nonbonded (vdw type + partial charge), one entry per typed atom ---
      ;; Iterate ATOM-TYPES itself: its keys are the typed atoms and its values
      ;; are the vdw type symbols, so the type is always present.  Charge is read
      ;; straight off the atom, matching EnergyAtom::defineForAtom.  PART supplies
      ;; the (atom-name . constitution-context) key and errors if an atom has no
      ;; constitution-context, the same contract as the bonded walks above.
      ;; ATOM-TYPES covers the whole AGGREGATE, so MINEP is what keeps this to MOLECULE.
      (maphash (lambda (atom vdw-type)
                 (when (minep atom)
                   (install-parameter-cache-entry
                    cache :nonbond (nonbond-table cache) (part atom)
                    (make-instance 'nonbond-parameters
                                   :vdw-type vdw-type
                                   :charge (chem:get-charge atom))
                    provenance)))
               atom-types))
    cache))

;;;; ==========================================================================
  ;;;; cached-smirnoff-force-field  (predictive)
  ;;;;
  ;;;; Composition, not inheritance: holds the real combined-smirnoff force field
  ;;;; plus a parameter-cache.  On a per-molecule basis it checks cached bonded and nonbonded
  ;;;; terms, otherwise it delegates the whole flow to real SMIRNOFF and harvests.
  ;;;; ==========================================================================

(defclass cached-smirnoff-force-field ()
  ((smirnoff      :initarg :smirnoff      :accessor smirnoff-force-field)   ; the real combined-smirnoff
   (smirnoff-name :initarg :smirnoff-name :initform :smirnoff :accessor smirnoff-name)
   (cache         :initform (make-instance 'parameter-cache) :initarg :cache :accessor cache)))

(define-condition grouped-parameter-cache-miss (chem:force-field-cache-miss)
  ((force-field :initarg :force-field :reader cache-miss-force-field)
   (molecule :initarg :molecule :reader cache-miss-molecule)
   (kind :initarg :kind :reader cache-miss-kind)
   (key :initarg :key :reader cache-miss-key)
   (atom :initarg :atom :reader cache-miss-atom))
  (:report
   (lambda (condition stream)
     (format stream
             "Grouped SMIRNOFF generation requires a parameter that is not cached.~%~
              Molecule: ~a~%Reason: ~a~@[~%Missing key: ~s~]~@[~%First atom: ~a~]"
             (chem:get-name (cache-miss-molecule condition))
             (cache-miss-kind condition)
             (cache-miss-key condition)
             (let ((atom (cache-miss-atom condition)))
               (and atom (chem:get-name atom)))))))

(defun grouped-parameter-cache-miss-resolved-p (condition)
  "True when CONDITION's exact atom or bond key is now in its SMIRNOFF cache."
  (let* ((cache (cache (cache-miss-force-field condition)))
         (kind (cache-miss-kind condition))
         (key (cache-miss-key condition)))
    (case kind
      (:atom-not-in-cache
       (nth-value 1 (gethash key (nonbond-table cache))))
      (:bond-not-in-cache
       (nth-value 1 (gethash key (bond-table cache))))
      (otherwise nil))))

(defmethod chem:retry-force-field-cache-miss
    ((condition grouped-parameter-cache-miss))
  "Verify this exact SMIRNOFF key before invoking the private grouped-lookup restart."
  (unless (grouped-parameter-cache-miss-resolved-p condition)
    (error "The force-field cache-miss handler returned success, but uncached SMIRNOFF did not ~
            produce the requested ~s entry ~s"
           (cache-miss-kind condition) (cache-miss-key condition)))
  (invoke-restart 'retry-grouped-cache-lookup))

(defmethod chem:force-fields-as-list ((ff cached-smirnoff-force-field))
    ;; nonbond merging (leap.core:nonbond-force-field-component) needs the real
    ;; SMIRNOFF parts: their ffnonbond-db resolves the $-vdw types the cache stored.
  (chem:force-fields-as-list (smirnoff-force-field ff)))

(defmethod chem:nonbond-force-field-name ((ff cached-smirnoff-force-field) force-field-name)
  (declare (ignore force-field-name))
  ;; share the wrapped SMIRNOFF's nonbond db: bonded terms come from the cache, so
  ;; the nonbond merge should treat this molecule as :smirnoff and dedup with it.
  (smirnoff-name ff))


  ;;; --- cache key + coverage -------------------------------------------------
 (defun atom-cache-key (atom &optional (errorp t))
    "An atom's cache key is its :given-atom-type = (atom-name . constitution-context).

  ERRORP NIL returns NIL rather than signalling.  A COVERAGE check wants that: an atom with no
  :given-atom-type is exactly an atom the cache cannot hold a key for, so the answer is 'not
  covered' - fall back to real SMIRNOFF - not an error.  Signalling there turns a cache miss into
  a crash for any molecule carrying atoms this cache has never seen.

  Harvesting still wants the error: an untyped atom at that point means something upstream
  failed to type it."
    (or (chem:matter-get-property-or-default atom :given-atom-type nil)
        (and errorp
             (error "cached-smirnoff-force-field: atom ~s has no :given-atom-type" atom))))

(defun bonded-cache-covers-molecule-p (molecule cache)
  "T iff every per-atom nonbond key and every bond key of MOLECULE is already cached.
  Bonds and nonbonds are the reliable coverage signal - SMIRNOFF parameterizes every
  bond and every atom (no skips), unlike angles/dihedrals which are legitimately absent
  for linear groups (e.g. -CN).  All bonds present => every atom's constitution-context
  has been harvested => every angle/dihedral SMIRNOFF would assign is present too, and
  any absent angle/dihedral is a genuine SMIRNOFF skip (safe to skip in commit)."
  (block covered 
    ;; ATOM-CACHE-KEY with ERRORP NIL - an atom with no :given-atom-type is not covered, which
    ;; is an answer this predicate is entitled to give.  Signalling made an untyped atom fatal
    ;; instead of merely uncached.
    (chem:map-atoms nil
                    (lambda (atom)
                      (let ((key (atom-cache-key atom nil)))
                        (unless (and key
                                     (nth-value 1 (gethash key (nonbond-table cache))))
                          (return-from covered nil))))
                    molecule)
    (chem:map-bonds nil
                    (lambda (a1 a2 o b) (declare (ignore o b))
                      (let ((k1 (atom-cache-key a1 nil))
                            (k2 (atom-cache-key a2 nil)))
                        (unless (and k1 k2
                                     (nth-value 1 (gethash (canonicalize-key (list k1 k2))
                                                           (bond-table cache))))
                          (return-from covered nil))))
                    molecule)
    t))

(defun uncovered-cache-reason (molecule cache)
  "Why BONDED-CACHE-COVERS-MOLECULE-P said no: (values KIND KEY ATOM), or NIL when it said yes.

  The predicate answers T or NIL and a NIL is expensive - it drops the whole molecule onto real
  SMIRKS, which is vf2 over every atom for every term in every force field.  On a blueprint's
  fanned-out molecule that is tens of thousands of vertices, and coverage is ALL OR NOTHING: one
  untyped atom out of 27000 costs the entire molecule.  So the one thing worth knowing on a miss is
  WHICH atom, and the predicate is not allowed to tell you."
  (block reason
    (chem:map-atoms nil
                    (lambda (atom)
                      (let ((key (atom-cache-key atom nil)))
                        (cond ((null key)
                               (return-from reason (values :no-given-atom-type nil atom)))
                              ((not (nth-value 1 (gethash key (nonbond-table cache))))
                               (return-from reason (values :atom-not-in-cache key atom))))))
                    molecule)
    (chem:map-bonds nil
                    (lambda (a1 a2 o b) (declare (ignore o b))
                      (let ((k1 (atom-cache-key a1 nil))
                            (k2 (atom-cache-key a2 nil)))
                        (when (and k1 k2
                                   (not (nth-value 1 (gethash (canonicalize-key (list k1 k2))
                                                              (bond-table cache)))))
                          (return-from reason
                            (values :bond-not-in-cache (canonicalize-key (list k1 k2)) a1)))))
                    molecule)
    nil))

(defvar *reported-coverage-failures* (make-hash-table :test #'equal)
  "Coverage failures already reported, keyed by (molecule-name kind key).

  Throttling, not memoization.  A blueprint detach scan re-parameterizes the same molecule once per
  rotamer slot, so an unthrottled report would print the same line ~2000 times and bury the run.")

(defun report-coverage-failure (molecule cache)
  "Warn once about why MOLECULE misses the cache.  Called only on the slow branch."
  (multiple-value-bind (kind key atom) (uncovered-cache-reason molecule cache)
    (when kind
      (let ((id (list (chem:get-name molecule) kind key)))
        (unless (gethash id *reported-coverage-failures*)
          (setf (gethash id *reported-coverage-failures*) t)
          #+(or)(format *error-output*
                  "~&;;; SMIRNOFF CACHE MISS on molecule ~a - falling back to real SMIRKS.~%~
                   ;;;   reason: ~a~@[  key: ~s~]~@[  first atom: ~a~]~%~
                   ;;;   Real SMIRKS is vf2 over the WHOLE molecule per force-field term.  On a~%~
                   ;;;   blueprint's fanned-out molecule that dominates everything else, and the~%~
                   ;;;   result is not harvested (a scoped pass must never teach the cache), so~%~
                   ;;;   it is paid again for every later pass.  Train this context first.~%"
                  (chem:get-name molecule) kind key
                  (and atom (chem:get-name atom)))
          )))))

  ;;; --- atom typing on a HIT: memoized vdw $types, no SMIRKS re-run ----------
(defmethod chem:assign-force-field-types ((ff cached-smirnoff-force-field) molecule atom-types)
  ;; Only reached on the covered branch, so every key is present.  Charge stays on
  ;; the atom (construct reads it), so we only need to restore the vdw type.
  (let ((nb (nonbond-table (cache ff))))
    (chem:map-atoms nil
                    (lambda (atom)
                      (let ((p (gethash (atom-cache-key atom) nb)))
                        (unless p (error "cached vdw type missing for ~s" atom))  ; defensive
                        (setf (gethash atom atom-types) (vdw-type p))))
                    molecule)))

 ;;; --- the per-molecule entry points, one per pass ---------------------------
 ;;;
 ;;; Coverage is tested SEPARATELY in each pass rather than memoized from pass 1, because an
 ;;; EARLIER molecule's harvest during pass 2 can cover us in between.  The cache only ever
 ;;; grows, so a pass-1 hit implies a pass-2 hit; the only new combination is pass-1 miss /
 ;;; pass-2 hit, which is a win - bonded terms come from the cache and no harvest is needed,
 ;;; since whoever covered us stored the same keys.
 ;;;
 ;;; One cost is inherent to the two-pass order: every coverage test in pass 1 happens before
 ;;; any harvest, so N identical molecules (a solvent box) each pay real-SMIRNOFF vdw TYPING
 ;;; instead of only the first.  Bonded generation - the expensive half, vf2 over the whole
 ;;; molecule per term - is still cached from the second molecule on.

(defmethod chem:construct-atom-table-for-molecule
    (ef molecule (ff cached-smirnoff-force-field) ff-name atom-types nonbond-force-field keep)
  (declare (ignore ff-name))
  (if (bonded-cache-covers-molecule-p molecule (cache ff))
      (chem:assign-force-field-types ff molecule atom-types)                        ; memoized vdw types
      (chem:assign-force-field-types (smirnoff-force-field ff) molecule atom-types)) ; real SMIRKS
  (chem:construct-from-molecule (chem:atom-table ef) molecule nonbond-force-field keep atom-types))

(defmethod chem:generate-for-molecule-using-force-field
    (ef molecule (ff cached-smirnoff-force-field) ff-name atom-types nonbond-force-field keep group)
  (declare (ignore ff-name atom-types nonbond-force-field))
  (cond
    (group
     ;; A grouped pass is deliberately forbidden from teaching the cache from its scoped/fanned
     ;; molecule.  Signal a structured miss instead.  A caller that knows how to construct a
     ;; chemically real trainer can fill the cache and invoke RETRY-GROUPED-CACHE-LOOKUP; callers
     ;; without that knowledge get the condition as an ordinary error.
     (loop
       (when (bonded-cache-covers-molecule-p molecule (cache ff))
         (return (commit-bonded-from-cache ef molecule (cache ff) keep group)))
       (multiple-value-bind (kind key atom)
           (uncovered-cache-reason molecule (cache ff))
         (restart-case
             (error 'grouped-parameter-cache-miss
                    :force-field ff :molecule molecule
                    :kind kind :key key :atom atom)
           (retry-grouped-cache-lookup ()
             :report "Retry the grouped SMIRNOFF cache lookup after filling the cache.")))))
    ((bonded-cache-covers-molecule-p molecule (cache ff))
     (commit-bonded-from-cache ef molecule (cache ff) keep group))
    (t
     (progn
        ;; Say WHY, once.  This branch is orders of magnitude slower than the other and there is
        ;; nothing in the output to distinguish "the run is big" from "the cache is missing one
        ;; atom's context and every molecule is being re-matched from scratch".
        (report-coverage-failure molecule (cache ff))
        (chem:generate-molecule-energy-function-tables
         ef molecule (smirnoff-force-field ff) keep group)
        ;; Harvest only from an UNFILTERED parameterization.  KEEP is T for a normal build; NIL
        ;; means no components were generated at all (an energy function built for its atom table
        ;; alone), and a FUNCTION means the terms present are a deliberate subset.  Feeding either
        ;; into the cache would store an incomplete parameterization that a later full build then
        ;; treats as covered - a silently wrong energy function rather than a slow one.
        ;;
        ;; A GROUP likewise means a scoped pass (one blueprint rotamer slot), never the whole
        ;; molecule - so it must not be harvested either.
        (when (and (eq keep t) (null group))
          (harvest-into-cache ef molecule (cache ff)))))))

(defun commit-bonded-from-cache (ef molecule cache keep-interaction-factory group)
  "Coverage guaranteed all bonds+nonbonds present.  Angles/dihedrals/impropers apply
  skip-if-absent: an absent key is a legitimate SMIRNOFF skip (e.g. a linear-group
  angle) that real SMIRNOFF would also omit, so skipping reproduces its energy function.
  GROUP is NIL or the ENERGY-COMPONENT-GROUP this pass's components belong to."
  (when keep-interaction-factory
    (let ((at       (chem:atom-table ef))
          (stretch-test  (or (eq keep-interaction-factory t)
                             (funcall keep-interaction-factory
                                      (load-time-value (find-class 'chem:energy-stretch)))))
          (angle-test  (or (eq keep-interaction-factory t)
                           (funcall keep-interaction-factory
                                    (load-time-value (find-class 'chem:energy-angle)))))
          (dihedral-test  (or (eq keep-interaction-factory t)
                              (funcall keep-interaction-factory
                                       (load-time-value (find-class 'chem:energy-dihedral)))))
          )
      (when stretch-test
        (let ((stretch (chem:make-energy-stretch ef group)))
          (chem:map-bonds
           nil
           (lambda (a1 a2 o b) (declare (ignore o b))
             (when (or (eq stretch-test t) (funcall stretch-test a1 a2))
               (let ((p (gethash (canonicalize-key (list (atom-cache-key a1) (atom-cache-key a2)))
                                 (bond-table cache))))
                 (when p (chem:add-stretch-term stretch at a1 a2 (kb p) (r0 p))))))
           molecule)))
      (when angle-test
        (let ((angle (chem:make-energy-angle ef group)))
          (chem:map-angles
           nil
           (lambda (a1 a2 a3)
             (when (or (eq angle-test t) (funcall angle-test a1 a2 a3))
               (let ((p (gethash (canonicalize-key (list (atom-cache-key a1) (atom-cache-key a2)
                                                         (atom-cache-key a3)))
                                 (angle-table cache))))
                 (when p (chem:add-angle-term angle at a1 a2 a3 (kt p) (t0 p))))))
           molecule)))
      (when dihedral-test
        (let ((dihedral (chem:make-energy-dihedral ef group)))
          (chem:map-dihedrals
           nil
           (lambda (a1 a2 a3 a4)
             (when (or (eq dihedral-test t) (funcall dihedral-test a1 a2 a3 a4))
               (let ((p (gethash (canonicalize-key (list (atom-cache-key a1) (atom-cache-key a2)
                                                         (atom-cache-key a3) (atom-cache-key a4)))
                                 (dihedral-table cache))))
                 (when p (loop for fc in (fourier-terms p)
                               do (chem:add-dihedral-term dihedral at a1 a2 a3 a4
                                                          (phase-rad fc) t (v fc) (periodicity fc)))))))
           molecule)
          (chem:map-impropers
           nil
           (lambda (a1 a2 a3 a4)
             (when (or (eq dihedral-test t) (funcall dihedral-test a1 a2 a3 a4))
               (let ((p (gethash (canonicalize-key (list (atom-cache-key a1) (atom-cache-key a2)
                                                         (atom-cache-key a3) (atom-cache-key a4)))
                                 (improper-table cache))))
                 (when p (loop for fc in (fourier-terms p)
                               do (chem:add-dihedral-term dihedral at a1 a2 a3 a4
                                                          (phase-rad fc) nil (v fc) (periodicity fc)))))))
           molecule))))))


;;; --- harvest after a miss: fill tables + record the signature -------------
(defun harvest-into-cache (ef molecule cache)
  "After a real-SMIRNOFF fallback, copy the parameterization into CACHE (bonded + the
  vdw $type/charge nonbond memo, both via expand-parameter-cache)"
  (let ((atom-to-context (make-hash-table :test 'eq)))
    (chem:map-atoms nil
                    (lambda (atom)
                      (let ((given (chem:matter-get-property-or-default atom :given-atom-type nil)))
                        (when given (setf (gethash atom atom-to-context) (cdr given)))))
                    molecule) 
    (expand-parameter-cache ef molecule atom-to-context (chem:atom-types ef) cache)
    ))
