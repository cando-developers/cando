
(in-package :pas)

(defclass build-context ()
  ((match :initarg :match :accessor match)
   (molecule :initarg :molecule :accessor molecule)
   (name-counter :initarg :name-counter :accessor name-counter)
   (ring-bonds :initform (make-hash-table) :accessor ring-bonds)))

(defclass ring-bond ()
  ((source-atom :initarg :source-atom :accessor source-atom)
   (target-atom :initarg :target-atom :accessor target-atom)
   (bond-order :initform :single-bond :initarg :bond-order :accessor bond-order)))

(defclass ring-start-close ()
  ((ring-id :initform nil :initarg :ring-id :accessor ring-id)))

(defclass ring-start (ring-start-close) ())
(defclass ring-close (ring-start-close) ())

(defclass atom-node ()
  ((children :initform nil :initarg :children :accessor children)
   (parent :initarg :parent :accessor parent)
   (bond-order :initarg :bond-order :initform nil :accessor bond-order)
   (element :initarg :element :accessor element)
   (hydrogens :initform 0 :initarg :hydrogens :accessor hydrogens)
   (tag :initform nil :initarg :tag :accessor tag)
   (ring-closing :initform nil :initarg :ring-closing :accessor ring-closing)))

(defmethod print-object ((object atom-node) stream)
  (print-unreadable-object (object stream :type t)
    (format stream ":element ~a :hydrogens ~a :tag ~a :ring-closing ~a :parent ~a :children-number ~a"
            (element object)
            (hydrogens object)
            (tag object)
            (ring-closing object)
            (parent object)
            (length (children object)))))

(defgeneric walk-for-atom-node (node parent node-map atom-node))

(defmethod walk-for-atom-node ((node chem:logical) parent node-map atom-node)
  (let ((children (chem:chem-info-node-children node)))
    (loop for child in children
          do (walk-for-atom-node child parent node-map atom-node))))

(defmethod walk-for-atom-node ((node chem:atom-test) parent node-map atom-node)
  (case (chem:atom-test-type node)
    (:sapelement
     (setf (element atom-node) (chem:get-symbol-arg node)))
    (:sapatom-map
     (setf (tag atom-node) (chem:get-int-arg node)))
    (:saptotal-hcount
     (setf (hydrogens atom-node) (chem:get-int-arg node)))
    (otherwise
     (format t "What do I do with atom-test: ~a test: ~a symbol: ~a~%"
             node
             (chem:atom-test-type node)
             (chem:get-symbol-arg node))))
  (case (chem:get-ring-test node)
    (:sarring-set
     (setf (ring-closing atom-node) (make-instance 'ring-start
                                                    :ring-id (chem:get-ring-id node))))
    (:sarring-test
     (setf (ring-closing atom-node) (make-instance 'ring-close
                                                    :ring-id (chem:get-ring-id node))))
    (otherwise nil)))

(defmethod walk-for-atom-node ((node chem:bond-to-atom-test) parent node-map atom-node)
  (setf (bond-order atom-node) (chem:bond-type node))
  (loop for child in (chem:chem-info-node-children node)
        do (walk-for-atom-node child parent node-map atom-node)))

(defgeneric walk-chem-info-tree (node parent node-map atom-stack))

(defmethod walk-chem-info-tree ((node chem:smarts-root) parent node-map atom-stack)
  (loop for child in (chem:chem-info-node-children node)
        do (walk-chem-info-tree child node node-map atom-stack)))

(defun register-node-map (node node-map atom-node)
  (setf (gethash node node-map) atom-node))

(defun lookup-node-map (node node-map)
  (let ((val (gethash node node-map)))
    val))

(defmethod walk-chem-info-tree ((node chem:bond-match-node) parent node-map atom-stack)
  (let ((atom-node (make-instance 'atom-node :parent (first atom-stack))))
    (walk-for-atom-node node parent node-map atom-node)
    (register-node-map node node-map atom-node)
    atom-node))

(defmethod walk-chem-info-tree ((node chem:chain) parent node-map atom-stack)
  (let* ((children (chem:chem-info-node-children node))
         (first (first children))
         (rest (second children)))
    (let ((an (walk-chem-info-tree first node node-map atom-stack)))
      (push an atom-stack))
    (when rest
      (walk-chem-info-tree rest node node-map atom-stack))))

(defmethod walk-chem-info-tree ((node chem:branch) parent node-map atom-stack)
  (let* ((children (chem:chem-info-node-children node))
         (left (first children))
         (right (second children)))
    (let ((atm-node (lookup-node-map parent node-map)))
      (register-node-map node node-map atm-node))
    (when left
      (walk-chem-info-tree left node node-map atom-stack))
    (when right
      (walk-chem-info-tree right node node-map atom-stack))))

(defun translate-chem-info-tree (root)
  (let ((node-map (make-hash-table)))
    (walk-chem-info-tree root nil node-map nil)
    (let ((seen-atoms (make-hash-table)))
      (maphash (lambda (key atm)
                 (unless (gethash atm seen-atoms)
                   (setf (gethash atm seen-atoms) t)
                   (when (and atm (parent atm))
                     (push atm (children (parent atm))))))
               node-map)
      (values seen-atoms node-map))))

(defun find-matches (molecule smarts)
  (let* ((smarts-graph (chem:make-chem-info-graph smarts))
         (mol-graph (chem:make-molecule-graph-from-molecule molecule)))
    (chem:boost-graph-vf2 smarts-graph mol-graph)))

(defun name-atom (context element)
  (let ((name (format nil "~a~a" (string element) (incf (name-counter context)))))
    (intern name :keyword)))

(defun find-residue (molecule cur-atm)
  (let* ((residue (block find-residue
                    (cando:do-residues (res molecule)
                      (cando:do-atoms (atm res)
                        (when (eq atm cur-atm)
                          (return-from find-residue res))))
                    nil)))
    (unless residue
      (error "Could not find residue containing ~s" cur-atm))
    residue))

(defun remove-atom (molecule atm)
  (let ((residue (find-residue molecule atm)))
    (chem:remove-atom-delete-bonds residue atm)))
        
(defun join-atom (context cur-atm new-atm bond-order)
  (unless cur-atm
    (error "The cur-atm is nil - it must be defined to join-atom"))
  (let* ((molecule (molecule context))
         (residue (find-residue molecule cur-atm)))
    (chem:add-matter residue new-atm)
    (let ((bo (case bond-order
                ((:sabsingle-bond :sabsingle-or-aromatic-bond) :single-bond)
                ((:sabdouble-bond :sabdouble-or-aromatic-bond) :double-bond)
                ((:sabtriple-bond :sabtriple-or-aromatic-bond) :triple-bond)
                (otherwise :single-bond))))
      (chem:bond-to cur-atm new-atm bo))))

(defun walk-build-instructions (build context cur-atm)
  (if (tag build)
      (let ((tag-atm (elt (match context) (tag build))))
        (setf cur-atm tag-atm)
        (chem:set-element tag-atm (element build)))
      (let* ((element (element build))
             (new-atm (chem:make-atom (name-atom context element) element))
             (bond-order (bond-order build)))
        (join-atom context cur-atm new-atm bond-order)
        (setf cur-atm new-atm)))
  ;; Add any ring closing bond information
  (when (ring-closing build)
    (if (typep (ring-closing build) 'ring-start)
        (setf (gethash (ring-id (ring-closing build)) (ring-bonds context))
              (make-instance 'ring-bond
                             :source-atom cur-atm))
        (let ((ring-bond (gethash (ring-id (ring-closing build)) (ring-bonds context))))
          (setf (target-atom ring-bond) cur-atm))))
  ;; Now walk the children
  (loop for build-child in (children build)
        do (walk-build-instructions build-child context cur-atm)))

(defun calculate-name-counter (molecule)
  (let ((max-counter -1))
    (cando:do-atoms (atm molecule)
      (let* ((name (string (chem:get-name atm)))
             (first-digit-pos (1+ (position-if-not #'digit-char-p name :from-end t))))
        (when (< first-digit-pos (length name))
          (let ((numeric-part (parse-integer name :start first-digit-pos)))
            (if (> numeric-part max-counter)
                (setf max-counter numeric-part))))))
    (1+ max-counter)))

(defun finish-build (context)
  (maphash (lambda (ring-id ring-bond)
             (chem:bond-to (source-atom ring-bond) (target-atom ring-bond) :single-bond))
           (ring-bonds context)))
    
(defun apply-match (old-match old-molecule product-smarts)
  (let* ((new-to-old (make-hash-table))
         (new-mol (chem:copy old-molecule new-to-old))
         (build (translate-chem-info-tree product-smarts))
         (old-to-new (make-hash-table))
         (tag-to-atm-node (make-hash-table :test #'eql)))
    (maphash (lambda (atm-node dummy)
               (declare (ignore dummy))
               (when (tag atm-node)
                 (setf (gethash (tag atm-node) tag-to-atm-node) atm-node)))
             build)
    (maphash (lambda (k v)
               (setf (gethash v old-to-new) k))
             new-to-old)
    ;; Delete atoms that have tags in the reactant but not the product
    ;; And generate a new vector of matching atoms in the new molecule
    (let ((new-match (make-array (length old-match) :initial-element nil)))
      (loop for index from 0 below (length old-match)
            for old-atm = (elt old-match index)
            when old-atm
              do (let ((new-atm (gethash old-atm old-to-new)))
                   (if (null (gethash index tag-to-atm-node))
                       (progn
                         (remove-atom new-mol new-atm))
                       (progn
                         (setf (elt new-match index) new-atm)))))
      ;; Find the root build instruction
      (let (root)
        (maphash (lambda (build dummy)
                   (unless (parent build)
                     (setf root build)))
                 build)
        ;; Now walk the build instructions and apply them to the molecule
        (let ((context (make-instance 'build-context
                                      :match new-match
                                      :molecule new-mol
                                      :name-counter (calculate-name-counter new-mol))))
          (walk-build-instructions root context nil)
          (finish-build context)
          (let ((mol (molecule context)))
            (cando:simple-build-unbuilt-atoms mol)
            mol))))))

(defun pas (pattern molecule)
  (let* ((smirks (chem:compile-smirks pattern))
         (reactant-smarts (chem:reactant smirks))
         (product-smarts (chem:product smirks))
         (matches (find-matches molecule reactant-smarts)))
    ;; Rename the pas-molecules so that they all have unique names
    (format t "matches -> ~s~%" matches)
    (let ((pas-molecules (loop for match in matches
                               collect (apply-match match molecule product-smarts))))
      (loop for counter from 1
            for mol in pas-molecules
            for name = (chem:get-name mol)
            for name-str = (string name)
            for new-name-str = (format nil "~a-~d" name-str counter)
            do (chem:set-name mol (intern new-name-str :keyword)))
      pas-molecules)))

