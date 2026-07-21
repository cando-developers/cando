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

;;; ---------------------------------------------------------------------------
;;; Build the cache from an energy function
;;; ---------------------------------------------------------------------------

(defun expand-parameter-cache (energy-function atom-to-constitution-context atom-types &optional cache)
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
Provide an existing cache or NIL if you don't have one yet."
  (let ((cache (if cache
                   cache
                   (make-instance 'parameter-cache))))
    (flet ((part (atom)
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
             (setf (gethash (canonicalize-key (list (part a1) (part a2)))
                            (bond-table cache))
                   (make-instance 'bond-parameters :kb kb :r0 r0))))))
      ;; --- angles ---
      (let ((angle (chem:energy-function/get-angle-component energy-function)))
        (when angle
          (chem:walk-angle-terms
           angle
           (lambda (index a1 a2 a3 i1 i2 i3 kt t0)
             (declare (ignore index i1 i2 i3))
             (setf (gethash (canonicalize-key (list (part a1) (part a2) (part a3)))
                            (angle-table cache))
                   (make-instance 'angle-parameters :kt kt :t0 t0))))))
      ;; --- dihedrals + impropers (one component, split by PROPER) ---
      (let ((dihedral (chem:energy-function/get-dihedral-component energy-function)))
        (when dihedral
          (chem:walk-dihedral-terms
           dihedral
           (lambda (index a1 a2 a3 a4 i1 i2 i3 i4 v n phase proper)
             (declare (ignore index i1 i2 i3 i4))
             (let* ((key (canonicalize-key (list (part a1) (part a2) (part a3) (part a4))))
                    (table (if proper (dihedral-table cache) (improper-table cache)))
                    (fourier (make-instance 'fourier-term :v v :periodicity n :phase-rad phase))
                    (existing (gethash key table)))
               (if existing
                   (add-fourier-term existing fourier)
                   (setf (gethash key table)
                         (add-fourier-term
                          (make-instance (if proper 'dihedral-parameters 'improper-parameters))
                          fourier))))))))
      ;; --- nonbonded (vdw type + partial charge), one entry per typed atom ---
      ;; Iterate ATOM-TYPES itself: its keys are the typed atoms and its values
      ;; are the vdw type symbols, so the type is always present.  Charge is read
      ;; straight off the atom, matching EnergyAtom::defineForAtom.  PART supplies
      ;; the (atom-name . constitution-context) key and errors if an atom has no
      ;; constitution-context, the same contract as the bonded walks above.
      (maphash (lambda (atom vdw-type)
                 (setf (gethash (part atom) (nonbond-table cache))
                       (make-instance 'nonbond-parameters
                                      :vdw-type vdw-type
                                      :charge (chem:get-charge atom))))
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
 (defun atom-cache-key (atom)
    "An atom's cache key is its :given-atom-type = (atom-name . constitution-context)."
    (or (chem:matter-get-property-or-default atom :given-atom-type nil)
        (error "cached-smirnoff-force-field: atom ~s has no :given-atom-type" atom)))

(defun bonded-cache-covers-molecule-p (molecule cache)
  "T iff every per-atom nonbond key and every bond key of MOLECULE is already cached.
  Bonds and nonbonds are the reliable coverage signal - SMIRNOFF parameterizes every
  bond and every atom (no skips), unlike angles/dihedrals which are legitimately absent
  for linear groups (e.g. -CN).  All bonds present => every atom's constitution-context
  has been harvested => every angle/dihedral SMIRNOFF would assign is present too, and
  any absent angle/dihedral is a genuine SMIRNOFF skip (safe to skip in commit)."
  (block covered 
    (chem:map-atoms nil
                    (lambda (atom)
                      (unless (nth-value 1 (gethash (atom-cache-key atom) (nonbond-table cache)))
                        (return-from covered nil)))
                    molecule)
    (chem:map-bonds nil
                    (lambda (a1 a2 o b) (declare (ignore o b))
                      (unless (nth-value 1 (gethash (canonicalize-key
                                                     (list (atom-cache-key a1) (atom-cache-key a2)))
                                                    (bond-table cache)))
                        (return-from covered nil)))
                    molecule)
    t))

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

 ;;; --- the per-molecule entry point -----------------------------------------
(defmethod chem:define-for-molecule-using-force-field
    (ef molecule (ff cached-smirnoff-force-field) ff-name atom-types nonbond-force-field keep)
  (declare (ignore ff-name))
  (let ((smirnoff (smirnoff-force-field ff))
        (cache    (cache ff)))
    (cond
      ((bonded-cache-covers-molecule-p molecule cache)
       ;; COVERED: memoized vdw typing (this class's assign), construct, bonded FROM CACHE
       (chem:assign-force-field-types ff molecule atom-types)
       (chem:construct-from-molecule (chem:atom-table ef) molecule nonbond-force-field keep atom-types)
       (commit-bonded-from-cache ef molecule cache))
      (t       
       ;; MISS: full real SMIRNOFF (assign vdw, construct, generate w/ perception) - pass :smirnoff
       (chem:define-for-molecule-using-force-field
           ef molecule smirnoff (smirnoff-name ff) atom-types nonbond-force-field keep)
       (harvest-into-cache ef molecule cache)))))

(defun commit-bonded-from-cache (ef molecule cache)
  "Coverage guaranteed all bonds+nonbonds present.  Angles/dihedrals/impropers apply
  skip-if-absent: an absent key is a legitimate SMIRNOFF skip (e.g. a linear-group
  angle) that real SMIRNOFF would also omit, so skipping reproduces its energy function."
  (let ((at       (chem:atom-table ef))
        (stretch  (chem:make-energy-stretch ef))
        (angle    (chem:make-energy-angle ef))
        (dihedral (chem:make-energy-dihedral ef)))
    (chem:map-bonds
     nil
     (lambda (a1 a2 o b) (declare (ignore o b))
       (let ((p (gethash (canonicalize-key (list (atom-cache-key a1) (atom-cache-key a2)))
                         (bond-table cache))))
         (when p (chem:add-stretch-term stretch at a1 a2 (kb p) (r0 p)))))
     molecule)
    (chem:map-angles
     nil
     (lambda (a1 a2 a3)
       (let ((p (gethash (canonicalize-key (list (atom-cache-key a1) (atom-cache-key a2)
                                                 (atom-cache-key a3)))
                         (angle-table cache))))
         (when p (chem:add-angle-term angle at a1 a2 a3 (kt p) (t0 p)))))
     molecule)
    (chem:map-dihedrals
     nil
     (lambda (a1 a2 a3 a4)
       (let ((p (gethash (canonicalize-key (list (atom-cache-key a1) (atom-cache-key a2)
                                                 (atom-cache-key a3) (atom-cache-key a4)))
                         (dihedral-table cache))))
         (when p (loop for fc in (fourier-terms p)
                       do (chem:add-dihedral-term dihedral at a1 a2 a3 a4
                                                  (phase-rad fc) t (v fc) (periodicity fc))))))
     molecule)
    (chem:map-impropers
     nil
     (lambda (a1 a2 a3 a4)
       (let ((p (gethash (canonicalize-key (list (atom-cache-key a1) (atom-cache-key a2)
                                                 (atom-cache-key a3) (atom-cache-key a4)))
                         (improper-table cache))))
         (when p (loop for fc in (fourier-terms p)
                       do (chem:add-dihedral-term dihedral at a1 a2 a3 a4
                                                  (phase-rad fc) nil (v fc) (periodicity fc))))))
     molecule)))


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
    (expand-parameter-cache ef atom-to-context (chem:atom-types ef) cache)
    ))

