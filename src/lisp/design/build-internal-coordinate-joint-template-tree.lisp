(in-package :design.joint-tree)


(defun interpret-builder-info-instruction (instr residue)
  (cond
    ((eq (car instr) 'chem:in-plug)
     (destructuring-bind (cmd
                          atom-name
                          plug-name
                          &key
                            frame
                            (bond-order :single-bond)
                            atom2-name
                            (bond-order2 :single-bond))
         instr
       (declare (ignore cmd))
       (let* ((atom (when (chem:has-atom-with-name residue atom-name)
                      (chem:first-atom-with-name residue atom-name))))
         (when atom
           (chem:make-in-plug plug-name nil atom-name bond-order atom2-name bond-order2)))))
    ((eq (car instr) 'chem:out-plug)
     (destructuring-bind (cmd
                          atom-name
                          plug-name
                          &key
                            (bond-order :single-bond)
                            atom2-name
                            (bond-order2 :single-bond))
         instr
       (declare (ignore cmd))
       (let* ((atom (when (chem:has-atom-with-name residue atom-name)
                      (chem:first-atom-with-name residue atom-name))))
         (when atom
           (chem:make-out-plug plug-name nil nil atom-name bond-order atom2-name bond-order2)))))
    (t (error "Unknown instruction ~a" instr))))
     
            
(defun interpret-plugs (fragment codes)
  (unless (consp codes)
    (error "code must be a list"))
  (loop for code in codes
        do (unless (and (consp code) (eq (car code) 'chem::builder-info))
             (let ((demo-code-block "(builder-info
 (out-plug :fg1 :fg1)
 (out-plug :fg2 :fg2)
 (out-plug :fg3 :fg3)
 (out-plug :+default :+default)
 (in-plug :-default :-default :frame :amine)
 (in-plug :fg :-side :amide-dkp))"))
               (error "interpret-plugs needs a code block that looks like the following~%~a~%" demo-code-block)))
        append (let ((instructions (cdr code)))
                 (loop for instr in instructions
                       collect (interpret-builder-info-instruction instr fragment)))))

(defclass prepare-topology ()
  ((name :initarg :name :reader name)
   (constitution-atoms :initarg :constitution-atoms :reader constitution-atoms)
   (plugs :initarg :plugs :reader plugs)
   (fragment :initarg :fragment :reader fragment))
  (:documentation "Store information needed to create a new topology"))

(defmethod print-object ((obj prepare-topology) stream)
  (print-unreadable-object (obj stream)
    (format stream "~a ~a" (class-name (class-of obj)) (name obj))))

(defun extract-topologys (chem-draw-object)
  "Loop through all of the fragments in the **chem-draw-object** and then loop through all of the residues
in each fragment. Those residues that have the property :topology-name are topology's that we want to extract.
Return a list of prepare-topology objects"
  (let ((code (chem:chem-draw-code chem-draw-object))
        (extracted-topologys nil))
    (loop for fragment in (chem:all-fragments-as-list chem-draw-object)
          do (cando:do-residues (residue (chem:get-molecule fragment))
               (when (chem:has-property residue :topology-name)
                 (push (make-instance 'prepare-topology
                                      :name (chem:matter-get-property residue :topology-name)
                                      :constitution-atoms (chem:make-constitution-atoms-from-residue residue)
                                      :plugs (remove nil (interpret-plugs residue code))
                                      :fragment fragment)
                       extracted-topologys))))
    extracted-topologys))

(defun find-in-plug (plugs)
  (loop for p in plugs
        when (chem:get-is-in p)
          do (return-from find-in-plug p)))

(defun find-out-plugs (plugs)
  (let (outs)
    (loop for p in plugs
          when (not (chem:get-is-in p))
            collect p)))



(defun get-all-out-plug-bond0atoms-as-set (plugs)
  (let ((res nil))
    (loop for p in plugs
          do (let ((atom-name (chem:get-b0 p)))
               (unless (not (contains res atomName))
                 (error "The atom name[~s] is already in the set[~s]" atom-name res))
               (push atom-name res)))
    res))

(defun get-all-out-plug-bond1atoms-as-set (plugs out-plug-bond0-atoms-set )
  (let ((res nil))
    (loop for p in plugs
          do (when (chem:get-b1 p)
               (let ((chem:atom-name (chem:get-b1 p)))
                 (unless (not (member atom-name res))
                   (error "The atom name[~s] is already in the listset[~s]" atom-name res))
                 (unless (not (member atom-name out-plug-bond0-atoms-set))
                   (error "An atom can't be both a bond0 and bond1 atom!!! Something is wrong with the plugs")))
                 (push atomName res)))
    res))

(defun get-all-stub-pivot-atoms-as-set (plugs out-plug-bond0-atoms-set )
  (let ((res (ObjectSet))
        (plugCounter 0))
    (loop for p in plugs
          do (when (chem:has-stub-pivot-atom p)
               (let ((atomName (getB1 p)))
                 (unless (not (member atom-name res))
                   (error "The atom name[~s] is already in the set[~s]" atom-name res))
                 (unless (not (member atom-name out-plug-bond0-atoms-set))
                   (error "An atom can't be both a bond0 and bond1 atom!!! Something is wrong with the plugs")))
               (push atom-name res))
             (incf plug-counter))
    res))

;;
;; Define the stub atoms in terms of
;; 1) parent older-sibling even-older-sibling
;; or if no even-older-sibling 2) parent older-sibling grand-parent
;; or if no older-sibling      3) parent grand-parent great-grand-parent
;; unless your parent is a OddIcoorAtom then just use
;; parent grand-parent
;;
;; (defun setStubAtoms (buildAtom parent)
;;   (IMPLEMENT_ME))
;; 
;; 
;; (defun addAtomToICoorTableRecursively (icoorTable parent me)
;;   (let ((buildAtom (cond
;;                  ((getProperty me :bond1Atom) (OddICoorAtom (atomName me)))
;;                  ( true (NormalICoorAtom (atomName me))))))
;;     (setProperty me :buildAtom buildAtom)
;; ;    (setStubAtoms buildAtom parent)
;;     (appendEntry icoorTable buildAtom)
;;     (for child (getProperty me :children)
;;       (addAtomToIcoorTableRecursively icoorTable me child))))
;; 

(defun one-atom-as-string (me constitution-atoms)
  (let ((target (make-string-output-stream))
        (result (make-string-output-stream)))
    (format target "~a" (kin:joint-template-atom-name me constitution-atoms))
    (format target "/~a" (class-name (class-of me)))
    (format result "~30a: " (get-output-stream-string target))
    (loop for child in (kin:children me)
          do (let ((child-name (kin:joint-template-atom-name child constitution-atoms)))
               (format result "~a/~a" child-name (kin:joint-template-comment child))))
    (get-output-stream-string result)))

(defun dump-build-order-recursively (parent-template cur-template constitution-atoms sout)
  (format sout "~a~%" (one-atom-as-string cur-template constitution-atoms))
  (dolist (child (kin:children cur-template))
    (dump-build-order-recursively cur-template child constitution-atoms sout)))


(defun joint-template-factory (parent-template atom constitution-atoms constitution-name topology-name fragment)
  (let* ((out-plug-atom-prop (chem:matter-get-property-or-default atom :out-plug nil))
         (entity-to-delay-children-for (chem:matter-get-property-or-default atom :entity-to-delay-children-for nil))
         (root-atom-prop (chem:matter-get-property-or-default atom :root-atom nil))
         (atom-index (chem:constitution-atoms-index constitution-atoms (chem:atom-name atom)))
         (atom-name (chem:atom-name (chem:atom-with-id constitution-atoms atom-index)))
         (comment (format nil "~s [~a]" atom-name (chem:matter-get-property-or-default atom :weight -1))) )
    (cond
      (entity-to-delay-children-for
       (when root-atom-prop
         (error "joint-template-factory problem - atom has entity-to-delay-children-for and root-atom-prop - atom[~s] constitutionName[~s] topologyName[~s]"  atom constitution-name topology-name))
       (let ((checkpoint (if (type-of entity-to-delay-children-for 'chem:plug)
                             (core:make-cxx-object 'kin:checkpoint-out-plug-joint
                                                   :constitution-name constitution-name
                                                   :topology-name topology-name
                                                   :out-plug entity-to-delay-children-for)
                             (core:make-cxx-object 'kin:checkpoint-joint
                                                   :constitution-name constitution-name
                                                   :topology-name topology-name
                                                   :atom-name (chem:atom-name entity-to-delay-children-for)))))
         (core:make-cxx-object 'kin:delayed-bonded-joint-template :id atom-index
                                                                  :parent parent-template
                                                                  :children #()
                                                                  :checkpoint checkpoint
                                                                  :comment comment
                                                                  :out-plug out-plug-atom-prop )))
      (root-atom-prop
       (core:make-cxx-object 'kin:root-bonded-joint-template :id atom-index
                                                             :parent parent-template
                                                             :children #()
                                                             :constitution-name constitution-name
                                                             :topology-name topology-name
                                                             :in-plug root-atom-prop
                                                             :comment comment
                                                             :out-plug out-plug-atom-prop))
      (t (core:make-cxx-object 'kin:bonded-joint-template :id atom-index
                                                          :parent parent-template
                                                          :children #()
                                                          :comment comment
                                                          :out-plug out-plug-atom-prop)))))


;;
;; Build an AtomTreeTemplate recursively using the properties defined
;; for each atom
;;
(defun build-atom-tree-template-recursively (parent fragment root constitution-atoms constitution-name topology-name )
  (let ((root-template (joint-template-factory parent root constitution-atoms constitution-name topology-name fragment))
        (children (progn
                    (chem:matter-get-property root :children))))
    (loop for child in children
          do (let ((child-template (build-atom-tree-template-recursively root-template fragment child constitution-atoms constitution-name topology-name )))
               (kin:add-child root-template child-template)))
    root-template))


(defun build-internal-coordinate-joint-template-tree (prepare-topology)
;;;chemdraw-fragment constitution-atoms plugs constitution-name topology-name)
  (with-accessors ((name name) (constitution-atoms constitution-atoms)
                   (plugs plugs) (fragment fragment))
      prepare-topology
    (let* ((residue (chem:build-residue constitution-atoms))
           (in-plug (let ((in-plug (find-in-plug plugs)))
                      (if in-plug
                          in-plug
                          (error "There has to be an in-plug in topology ~a" name))))
           (outplugs (find-out-plugs plugs))
           (root-atom-name (chem:root-atom-name in-plug))
           (root-atom (chem:atom-with-name residue root-atom-name))
           (spanning-loop (chem:make-spanning-loop root-atom))
           (all-spanning-atoms (chem:all-atoms spanning-loop)))
      (chem:set-property root-atom :root-atom in-plug)
      (loop for atom in all-spanning-atoms
            do (chem:set-property atom :out-plug nil)
               (chem:set-property atom :bond1-atom nil)
               (chem:set-property atom :weight 1)
               (chem:set-property atom :children nil))
      ;;
      ;; Now iterate over the plugs and assign atom properties that will 
      ;; define the build order and identify for OutPlugs
      ;; 1) :highestPriority atoms (atoms that bond0 atoms depend on)
      ;; 2) :highPriority atoms (bond0 atoms)
      ;; 3) atoms that need to delay building their children until another atom is defined
      ;;     identify by setting :entityToDelayChildrenFor XXXX where
      ;;           XXXX == :nextResidueBond1
      ;;        or XXXX == the stubPivotAtom that you need to build your children
      (loop for plug in plugs
            do (cond
                 ((typep plug 'chem:out-plug)
                  (let ((bond0-atom (chem:atom-with-name residue (chem:get-b0 plug)))
                        (bond1-atom (when (chem:get-b1 plug)
                                      (chem:atom-with-name residue (chem:get-b1 plug))))
                        (stub-pivot-atom (when (chem:has-stub-pivot-atom plug)
                                           (chem:atom-with-name residue (chem:get-stub-pivot-atom plug)))))
                    (unless (not (and (or bond1-atom stub-pivot-atom) (eq bond1-atom stub-pivot-atom)))
                      (error "The same plug cannot define a bond1Atom[~a] and a stubPivotAtom[~a] at the same time" bond1-atom stub-pivot-atom))
                    (chem:set-property bond0-atom :out-plug plug)
                    (chem:set-property bond0-atom :high-priority t)
                    (cond
                      (bond1-atom
                       (chem:set-property bond1-atom :bond1-atom t)
                       (chem:set-property bond1-atom :highest-priority t)
                       (chem:set-property bond1-atom :build-before bond0-atom)
                       (chem:set-property bond1-atom :entity-to-delay-children-for plug ))
                      (stub-pivot-atom
                       (let ((bond0-depends-on-atom (chem:matter-get-property-or-default bond0-atom :build-before nil)))
                         (if (eq bond0-depends-on-atom stub-pivot-atom)
                             (progn
                               ;; (% "Setting atom[%s] :entity-to-delay-children-for to-> %s" (description bond0-atom) stub-pivot-atom)
                               (chem:set-property bond0-atom :entity-to-delay-children-for stub-pivot-atom))
                             (progn
                               ;; (% "Setting stubPivotAtom[%s] :highestPriority=t and :buildBefore to-> %s" (description stubPivotAtom) (description bond0-atom))
                               (chem:set-property stub-pivot-atom :build-before bond0-atom)
                               (chem:set-property stub-pivot-atom :highest-priority t))))))))))
      ;; Now iterate over the atoms from farthest to the nearest to the root
      ;; and add yourself to your backSpanning atoms :children
      (dolist (atom all-spanning-atoms)
        (let ((parent-atom (chem:get-back-span atom)))
          (when parent-atom
            (let ((parent-children (chem:matter-get-property parent-atom :children)))
              (setf parent-children (cons atom parent-children))
              (chem:set-property parent-atom :children parent-children)))))
      ;; Now iterate over the atoms from farthest to the nearest and assign a branch weight
      ;; to each atom
      ;; 1) Get the largest weight of this atoms children
      ;; 2) If the current atom is a bond0Atom then add 1000 to its weight
      ;; 3) otherwise if the current atom is a bond1Atom then add 10000 to its weight
      ;; 4) otherwise add 1 to its weight
      ;;   This is to try and ensure that branches that involve Plugs that make two bonds are
      ;;   built first and that the bond1Atom of these plugs gets built before the bond0Atom.
      ;;   This is so that the bond1Atom will define the dihedral to the root atom of the next residue
      ;;   that is connected through the bond0Atom.
      ;;   This will help to isolate residues from each other in terms of internal dihedral coordinates
      ;;   so that bis-amino acids are built from the amino acid atoms of the previous residue rather
      ;;   than building from the inside of the previous bis-amino acid. 
      ;;   Here is a concrete example - If we had a pro4(SS) residue X modeled on another pro4(SS) residue Y and
      ;;   the dihedral for X's root were defined relative to the CB of Y then if we built X onto a pro4(SR) 
      ;;   residue Z then that dihedral would be wrong.  If instead we define X's root atom dihedral in terms
      ;;   of the CGC of Y then it won't matter what stereoisomer of Z we put on.
      ;;
      (dolist (atom all-spanning-atoms)
        (let ((atom-weight (cond
                             ((chem:matter-get-property-or-default atom :highest-priority nil) 10000)
                             ((chem:matter-get-property-or-default atom :high-priority nil)     1000)
                             (t                                                                    1))))
          (let ((max-child-weight 0))
            (loop for child in (chem:matter-get-property atom :children)
                  do (setf max-child-weight (max max-child-weight (chem:matter-get-property child :weight))))
            (setf atom-weight (+ atom-weight max-child-weight))
            (chem:set-property atom :weight atom-weight))))
      ;;
      ;; Now sort the children of each atom in descending order of their weight
      ;; and if there are ties (like methylene or methyl hydrogens) sort them
      ;; in increasing alphabetical order by their name
      ;;
      (dolist (atom all-spanning-atoms)
        (let* ((children (chem:matter-get-property atom :children))
               (sorted-children (sort children (lambda (a b)
                                                 (let ((wa (chem:matter-get-property a :weight))
                                                       (wb (chem:matter-get-property b :weight)))
                                                   (cond
                                                     ((> wa wb) t)
                                                     ((< wa wb) nil)
                                                     (t (string< (string (chem:atom-name a)) (string (chem:atom-name b))))))))))
          (chem:set-property atom :children sorted-children)))
      ;;
      ;; Ok, now every atom has :children defined and they are sorted in order of how they
      ;; should be built - now we are ready to build the InternalCoordinateTable for the residue
      ;;
      (let ((tree-template (build-atom-tree-template-recursively nil fragment root-atom constitution-atoms name name)))
        (let ((build-order (with-output-to-string (sout)
                             (dump-build-order-recursively nil tree-template constitution-atoms sout))))
          (format t "Build order:~%~a~%" build-order))
        tree-template))))





