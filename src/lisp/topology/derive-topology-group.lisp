(in-package #:topology)

(defun %derivation-function (value option)
  (unless (or (functionp value)
              (and (symbolp value) (fboundp value) (not (macro-function value))))
    (error "~s must be a function designator, got ~s" option value))
  (coerce value 'function))

(defun %derivation-symbol (function original)
  (let ((result (funcall function original)))
    (unless (and result (symbolp result))
      (error "Name transformation of ~s must return a non-NIL symbol, got ~s"
             original result))
    result))

(defun %derivation-plug-function (function)
  "Memoize canonical plug transformations and preserve .0/.1 bond suffixes."
  (let ((cache (make-hash-table :test #'eq)))
    (lambda (name)
      (unless (symbolp name)
        (error "Plug names must be symbols, got ~s" name))
      (let* ((text (symbol-name name))
             (dot (position #\. text))
             (base (subseq text 0 dot)))
        (unless (and (> (length base) 1)
                     (find (char base 0) "+-"))
          (error "Invalid plug name ~s" name))
        (let* ((canonical (intern base :keyword))
               (replacement
                 (or (gethash canonical cache)
                     (let* ((result (%derivation-symbol function canonical))
                            (result-text (symbol-name result)))
                       (unless (and (> (length result-text) 1)
                                    (char= (char result-text 0) (char base 0))
                                    (not (position #\. result-text)))
                         (error "Plug transformation ~s -> ~s must preserve direction and return an unindexed plug name"
                                canonical result))
                       (setf (gethash canonical cache)
                             (intern result-text :keyword))))))
          (if dot
              (intern (concatenate 'string (symbol-name replacement)
                                   (subseq text dot))
                      :keyword)
              replacement))))))

(defun %derivation-graph (tree plug-function)
  "Copy a definition graph, transforming only :PLUGS property values."
  (cond
    ((atom tree) tree)
    ((eq (car tree) :plugs)
     (list* :plugs (mapcar plug-function (second tree))
            (%derivation-graph (cddr tree) plug-function)))
    (t (cons (%derivation-graph (car tree) plug-function)
             (%derivation-graph (cdr tree) plug-function)))))

(defun %derivation-ring-info (tree plug-function)
  (cond
    ((consp tree)
     (cons (%derivation-ring-info (car tree) plug-function)
           (%derivation-ring-info (cdr tree) plug-function)))
    ((and (symbolp tree)
          (> (length (symbol-name tree)) 1)
          (find (char (symbol-name tree) 0) "+-"))
     (funcall plug-function tree))
    (t tree)))

(defun %derivation-options (definition plug-function omit-dihedrals)
  (let* ((options (copy-tree (getf definition :options)))
         (dihedrals (getf options :dihedrals))
         (kept-dihedrals
           (loop for entry in dihedrals
                 unless (member (second entry) omit-dihedrals :test #'eq)
                   collect
                   (if (consp (first entry))
                       (destructuring-bind ((path external-name) local-name) entry
                         (list (list (mapcar plug-function path) external-name)
                               local-name))
                       entry))))
    (setf (getf options :dihedrals) kept-dihedrals
          (getf options :cluster-dihedrals)
          (remove-if (lambda (name) (member name omit-dihedrals :test #'eq))
                     (getf options :cluster-dihedrals))
          (getf options :plug-names)
          (mapcar plug-function (getf options :plug-names))
          (getf options :rings)
          (%derivation-ring-info (getf options :rings) plug-function))
    (values options (or (getf definition :%no-dihedrals)
                        (and dihedrals (null kept-dihedrals))))))

(defun %derivation-dihedral-data (dihedrals)
  (loop for info in dihedrals
        collect
        (etypecase info
          (dihedral-info-atom (list (atom-name info) (name info)))
          (dihedral-info-external
           (list (list (copy-list (plug-path info)) (external-dihedral-name info))
                 (name info))))))

(defun %check-derived-group-name (group owner planned-owners)
  (multiple-value-bind (previous presentp) (gethash group planned-owners)
    (when (and presentp (not (equal previous owner)))
      (error "Different definitions transform to the same name ~s" group))
    (setf (gethash group planned-owners) owner))
  (dolist (member (gethash group *topology-groups*))
    (let ((topology (chem:find-topology member nil)))
      (unless (and topology
                   (equal (getf (property-list topology) :%derived-definition)
                          owner))
        (error "Derived name ~s conflicts with an existing topology group" group)))))

(defun derive-topology-group
    (destination-group &key from
                            (name-transform nil name-transform-p)
                            (plug-name-transform nil plug-name-transform-p)
                            (name-prefix nil name-prefix-p)
                            (rename-plugs nil rename-plugs-p)
                            omit-dihedrals)
  "Derive registered topologies in DESTINATION-GROUP from the group FROM.

NAME-TRANSFORM maps a constitution's base name and each concrete topology name
to a new non-NIL symbol. PLUG-NAME-TRANSFORM receives canonical keyword plug
names (without .0/.1 bond suffixes) and returns a replacement of the same
direction. Both callbacks run during derivation; they are not stored.

NAME-PREFIX is a string convenience alternative to NAME-TRANSFORM. RENAME-PLUGS
is an alist convenience alternative to PLUG-NAME-TRANSFORM; unlisted plugs keep
their names. Supplying both alternatives for either transformation is an error.
OMIT-DIHEDRALS names annotations to omit only from the derived definitions.

The original objects and groups are unchanged. Atom graphs, stereochemistry,
restraints, local dihedrals, ring information, synthesis codes, and rotamer
limits are rebuilt from recorded DEFINE-TOPOLOGY data. Plug changes also apply
to external dihedral paths and ring-start information. Source group subsets
produce only their selected stereoisomers.

All variants are built and checked before registration. Unrelated name
collisions signal errors; repeated derivations replace their own variants.
Return the list of derived topology objects. Topologies defined before source
recording was installed must have their DEFINE-TOPOLOGY forms reloaded first."
  (unless (and destination-group (symbolp destination-group)
               from (symbolp from) (not (eq destination-group from)))
    (error "FROM and DESTINATION-GROUP must be distinct non-NIL symbols"))
  (when (and name-transform-p name-prefix-p)
    (error "Supply only one of :NAME-TRANSFORM and :NAME-PREFIX"))
  (when (and plug-name-transform-p rename-plugs-p)
    (error "Supply only one of :PLUG-NAME-TRANSFORM and :RENAME-PLUGS"))
  (unless (or name-transform-p name-prefix-p)
    (error "Supply :NAME-TRANSFORM or :NAME-PREFIX"))
  (when (and name-prefix-p (not (stringp name-prefix)))
    (error ":NAME-PREFIX must be a string"))
  (dolist (entry rename-plugs)
    (unless (and (consp entry) (symbolp (car entry)) (symbolp (cdr entry))
                 (car entry) (cdr entry))
      (error ":RENAME-PLUGS must be an alist of non-NIL symbols, got ~s" entry)))
  (let* ((name-function
           (if name-transform-p
               (%derivation-function name-transform :name-transform)
               (lambda (name)
                 (intern (concatenate 'string name-prefix (symbol-name name))
                         (or (symbol-package name) *package*)))))
         (plug-function
           (%derivation-plug-function
            (if plug-name-transform-p
                (%derivation-function plug-name-transform :plug-name-transform)
                (lambda (name) (or (cdr (assoc name rename-plugs :test #'eq))
                                  name)))))
         (names (copy-list
                 (multiple-value-bind (members presentp)
                     (gethash from *topology-groups*)
                   (unless presentp (error "Unknown source topology group ~s" from))
                   members)))
         (families (make-hash-table :test #'eq))
         (ordered-definitions nil)
         (planned-owners (make-hash-table :test #'eq))
         (planned-topologies (make-hash-table :test #'eq))
         (entries nil))
    ;; Group stereoisomers by their shared recorded definition.
    (dolist (name names)
      (let* ((source (chem:find-topology name))
             (definition (getf (property-list source) :%topology-definition)))
        (unless definition
          (error "Topology ~s has no recorded definition; reload its DEFINE-TOPOLOGY form before deriving"
                 name))
        (when (member destination-group
                      (let ((groups (getf definition :name)))
                        (if (consp groups) groups (list groups)))
                      :test #'eq)
          (error "Destination ~s is also an original definition group" destination-group))
        (unless (gethash definition families)
          (push definition ordered-definitions))
        (push source (gethash definition families))))
    (dolist (definition (nreverse ordered-definitions))
      (let* ((source-name (getf definition :name))
             (source-base (if (consp source-name) (first source-name) source-name))
             (derived-base (%derivation-symbol name-function source-base))
             (owner (list destination-group source-base))
             (sources (nreverse (gethash definition families)))
             (plug-targets (make-hash-table :test #'eq)))
        (%check-derived-group-name derived-base owner planned-owners)
        ;; Distinct plugs within a definition must not collapse onto one name.
        (maphash (lambda (name plug)
                   (declare (ignore plug))
                   (let ((target (funcall plug-function name)))
                     (when (gethash target plug-targets)
                       (error "Distinct plugs in ~s transform to ~s" source-base target))
                     (setf (gethash target plug-targets) name)))
                 (plugs (first sources)))
        (multiple-value-bind (options omit-all-p)
            (%derivation-options definition plug-function omit-dihedrals)
          (let* ((*register-defined-topologies* nil)
                 (variants
                   (apply #'do-define-topology (list derived-base destination-group)
                          (%derivation-graph (getf definition :sexp) plug-function)
                          options)))
            ;; An explicit omission of all annotations must not turn into the
            ;; builder's default automatically generated dihedrals.
            (when omit-all-p
              (setf (getf (residue-properties (constitution (first variants)))
                          :dihedrals)
                    nil))
            (let* ((properties (residue-properties (constitution (first variants))))
                   (dihedrals
                     (remove-if (lambda (info)
                                  (member (name info) omit-dihedrals :test #'eq))
                                (getf properties :dihedrals)))
                   (derived-definition
                     (getf (property-list (first variants)) :%topology-definition)))
              (setf (getf (residue-properties (constitution (first variants))) :dihedrals)
                    dihedrals
                    (getf (getf derived-definition :options) :dihedrals)
                    (%derivation-dihedral-data dihedrals)
                    (getf derived-definition :%no-dihedrals) (null dihedrals)))
            (dolist (source sources)
              (let* ((index (getf (property-list source) :%definition-isomer-index))
                     (variant (and (integerp index) (nth index variants)))
                     (target-name
                       (if (eq (name source) source-base)
                           derived-base
                           (%derivation-symbol name-function (name source))))
                     (origin (list destination-group (name source)))
                     (existing (chem:find-topology target-name nil)))
                (unless variant
                  (error "Invalid recorded stereoisomer index for ~s" (name source)))
                (when (gethash target-name planned-topologies)
                  (error "Different source topologies transform to ~s" target-name))
                (when (and existing
                           (not (equal (getf (property-list existing) :%derived-from)
                                       origin)))
                  (error "Derived topology ~s collides with an existing topology" target-name))
                (%check-derived-group-name target-name owner planned-owners)
                (setf (gethash target-name planned-topologies) t
                      (name variant) target-name
                      (name (stereoisomer variant)) target-name
                      (getf (property-list variant) :%derived-from) origin
                      (getf (property-list variant) :%derived-definition) owner)
                (push (list variant derived-base) entries)))))))
    ;; Publish only after every definition and callback has passed its checks.
    (loop for (variant derived-base) in (nreverse entries)
          for name = (name variant)
          do (chem:register-topology variant name)
             (dolist (group (remove-duplicates
                            (list name derived-base destination-group) :test #'eq))
               (pushnew name (gethash group *topology-groups*) :test #'eq))
          collect variant)))
