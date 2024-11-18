(in-package :topology)

(defclass joint-template (cando.serialize:serializable)
  ((parent :initform nil :initarg :parent :accessor parent)
   (atom-name :initarg :atom-name :accessor atom-name)
   (constitution-atoms-index :initarg :constitution-atoms-index :accessor constitution-atoms-index)
   (children :initform nil :initarg :children :accessor children)))


(defmethod print-object ((obj joint-template) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (format stream "~a" (atom-name obj)))))

(defclass bonded-joint-template (joint-template) ())

(defclass adjustable-bonded-joint-mixin ()
  ((adjustment :initarg :adjustment :accessor adjustment)))

(defclass adjustment ()
  ((joint :initarg :joint :accessor joint)))

(defclass internal-adjustment (adjustment) ()
  (:documentation "Adjust the internal coordinates of the joint when you can use internal coordinates. This is the best way to adjust coordinates.  If you need to adjust an internal coordinate to be planar relative to the other child of the parent of joing use topology:internal-planar-adjustment"))

(defclass external-adjustment (adjustment) ()
  (:documentation "Adjust the external coordinates of the joint when you cannot use internal coordinates but can only use external coordinates"))

(defclass adjustments ()
  ((internal-adjustments :initform (make-hash-table) :initarg :internal-adjustments :accessor internal-adjustments)
   (external-adjustments :initform (make-hash-table) :initarg :external-adjustments :accessor external-adjustments)))

(defgeneric add-to-adjustments (adjustments adjustment atresidue))

(defmethod add-to-adjustments (adjustments (adjustment internal-adjustment) atresidue)
  (push adjustment (gethash atresidue (internal-adjustments adjustments) nil)))

(defmethod add-to-adjustments (adjustments (adjustment external-adjustment) atresidue)
  (push adjustment (gethash atresidue (external-adjustments adjustments) nil)))

(defgeneric initialize-adjustment (adjustment assembler))

(defgeneric external-adjust (adjustment assembler coordinates))
(defmethod external-adjust ((adjustment internal-adjustment) assembler coordinates)
  "default - don't do anything"
  nil)

(defgeneric internal-adjust (adjustment assembler))

(defmethod internal-adjust ((adjustment external-adjustment) assembler)
  "default - don't do anything"
  nil)

(defclass internal-planar-adjustment (internal-adjustment)
  ((other :initarg :other :accessor other))
  (:documentation "Adjust joint to be rotated 180 degrees from other"))

(defmethod initialize-adjustment ((adjustment internal-planar-adjustment) assembler)
  (let* ((joint (joint adjustment))
         (jparent (kin:parent joint))
         (jother (when (= (kin:number-of-children jparent) 2)
                   (kin:joint/only-other-child jparent joint))))
    ;; Only when there is an other joint do we set (other adjustment)
    ;; for ring closing connections we need to do something else
    (when jother
      (setf (other adjustment) jother))))

(defmethod internal-adjust ((adjustment internal-planar-adjustment) assembler)
  ;; Only when the "other" slot is bound do we make adjustment
  (when (slot-boundp adjustment 'other)
    (let* ((joint (joint adjustment))
           (other (other adjustment))
           (phi-original (kin:bonded-joint/get-phi joint))
           (phi-other (kin:bonded-joint/get-phi other))
           (phi-adjust (radians-add phi-other PI)))
      #+(or)(format t "internal-adjust ~s phi-original ~8,3f  phi-adjust ~8,3f~%"
              joint
              (topology:rad-to-deg phi-original) (topology:rad-to-deg phi-adjust))
      #+(or)(kin:bonded-joint/set-distance joint 1.47)
      #+(or)(kin:bonded-joint/set-theta joint (topology:deg-to-rad 120.0))
      (fill-joint-phi joint phi-adjust)
      )))

(defun make-bonded-joint-template (constitution-atoms-index &key atom-name parent)
  (make-instance 'bonded-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent))

(defun make-adjustable-bonded-joint-template (constitution-atoms-index &key atom-name parent adjustment)
  (make-instance 'adjustable-bonded-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent
                 :adjustment adjustment))
(defclass in-plug-bonded-joint-template (bonded-joint-template)
  ((in-plug :initarg :in-plug :accessor in-plug)))

(defclass adjustable-in-plug-bonded-joint-template (in-plug-bonded-joint-template
                                                    adjustable-bonded-joint-mixin)
  ())

(defclass adjustable-bonded-joint-template (bonded-joint-template
                                            adjustable-bonded-joint-mixin)
  ())

(defun make-in-plug-bonded-joint-template (constitution-atoms-index &key atom-name parent in-plug)
  (make-instance 'in-plug-bonded-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent
                 :in-plug in-plug))

(defun make-adjustable-in-plug-bonded-joint-template (constitution-atoms-index &key atom-name parent in-plug adjustment)
  (make-instance 'adjustable-in-plug-bonded-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent
                 :in-plug in-plug
                 :adjustment adjustment))

(defclass complex-bonded-joint-template (bonded-joint-template)
  ((input-stub-joints :initform (make-array 2) :initarg :input-stub-joints :accessor input-stub-joints)))

(defmethod print-object ((obj complex-bonded-joint-template) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (let ((stubs (loop for stub across (input-stub-joints obj)
                           when (not (null stub))
                             count 1)))
          (format stream "~a :stubs ~a" (atom-name obj) (1+ stubs))))))

(defun make-complex-bonded-joint-template (constitution-atoms-index &key atom-name stub-joints)
  (cond
    ((eq (length stub-joints) 1)
     (make-instance 'complex-bonded-joint-template
                    :constitution-atoms-index constitution-atoms-index
                    :atom-name atom-name
                    :parent (first stub-joints)
                    :input-stub-joints (make-array 2 :initial-element nil)))
    ((eq (length stub-joints) 2)
     (make-instance 'complex-bonded-joint-template
                    :constitution-atoms-index constitution-atoms-index
                    :atom-name atom-name
                    :parent (first stub-joints)
                    :input-stub-joints (make-array 2 :initial-contents (list (second stub-joints)
                                                                             nil))))
    ((eq (length stub-joints) 3)
     (make-instance 'complex-bonded-joint-template
                    :constitution-atoms-index constitution-atoms-index
                    :atom-name atom-name
                    :parent (first stub-joints)
                    :input-stub-joints (make-array 2 :initial-contents (list (second stub-joints)
                                                                             (third stub-joints)))))
    (t (error "Illegal number of stub-joints: ~s" stub-joints))))

(defclass jump-joint-template (joint-template) ())

(defclass xyz-joint-template (joint-template) ())

(defun make-jump-joint-template (constitution-atoms-index &key atom-name)
  (make-instance 'jump-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name))

(defun make-xyz-joint-template (constitution-atoms-index &key atom-name parent)
  (make-instance 'xyz-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent))

(defun add-child (joint-template child-template)
  (setf (parent child-template) joint-template)
  (setf (children joint-template) (append (children joint-template) (list child-template))))

(defun sibling (joint-template index)
  (nth index (children joint-template)))


(defun walk-joint-template (joint-template callback)
  "Walk a joint-template tree"
  (funcall callback joint-template)
  (loop for child in (children joint-template)
        do (walk-joint-template child callback)))

(defun new-joint-template-factory (parent-template node in-plug xyz-joints)
  (let* ((atom-name (topology:name node))
         (child-indexes (mapcar #'topology:constitution-atom-index (topology:children node)))
         (constitution-atoms-index (topology:constitution-atom-index node))
         (gparent-template (if parent-template
                               (parent parent-template)
                               nil))
         (ggparent-template (if gparent-template
                                (parent gparent-template)
                                nil)))
    (cond
      ((and xyz-joints
            (or (eq xyz-joints :all)
                (member atom-name xyz-joints)))
       (make-xyz-joint-template constitution-atoms-index
                                :atom-name atom-name
                                :parent parent-template))
      ((and (null parent-template) (typep in-plug 'topology:in-plug))
       (let ((adjust (getf (property-list node) :adjust)))
         (if adjust
             (make-adjustable-in-plug-bonded-joint-template
              constitution-atoms-index
              :atom-name atom-name
              :parent nil
              :in-plug in-plug
              :adjustment adjust)
             (make-in-plug-bonded-joint-template constitution-atoms-index
                                                 :atom-name atom-name
                                                 :parent nil
                                                 :in-plug in-plug))))
      ((typep in-plug 'topology:in-plug)
       (let ((adjust (getf (property-list node) :adjust)))
         (if adjust
             (make-adjustable-bonded-joint-template
              constitution-atoms-index
              :atom-name atom-name
              :parent nil
              :adjustment adjust)
             (make-bonded-joint-template constitution-atoms-index
                                         :atom-name atom-name
                                         :parent parent-template))))
      ((null parent-template)
       (make-jump-joint-template constitution-atoms-index
                                 :atom-name atom-name
                                 ))
      ((null gparent-template)
       (let ((stub-joints (cond
                            ((= 1 constitution-atoms-index)
                             (list parent-template))
                            ((= 2 constitution-atoms-index)
                             (list parent-template (sibling parent-template 0)))
                            (t
                             (list parent-template
                                   (sibling parent-template 0)
                                   (sibling parent-template 1))))))
         (make-complex-bonded-joint-template constitution-atoms-index
                                             :atom-name atom-name
                                             :stub-joints stub-joints)))
      ((null ggparent-template)
       (let* ((sibling (when (> (length (children parent-template)) 0)
                         (elt (children parent-template) 0)))
              (aunt (when (> (length (children gparent-template)) 0)
                      (elt (children gparent-template) 0)))
              (stub-joints (cond
                            ((= constitution-atoms-index 2)
                             (list parent-template
                                   gparent-template))
                            (sibling
                             (list parent-template
                                   sibling
                                   gparent-template))
                            (aunt
                             (list parent-template
                                   gparent-template
                                   aunt))
                            (t
                             (error "One of the previous cases must be true")))))
         (make-complex-bonded-joint-template constitution-atoms-index
                                             :atom-name atom-name
                                             :stub-joints stub-joints)))
      (t
       (let ((adjust (getf (property-list node) :adjust)))
         (if adjust
             (make-adjustable-bonded-joint-template
              constitution-atoms-index
              :atom-name atom-name
              :parent parent-template
              :adjustment adjust)
             (make-bonded-joint-template constitution-atoms-index
                                         :atom-name atom-name
                                         :parent parent-template)))))))


(defun build-joint-template-recursively (parent root in-plug xyz-joints)
  (let ((root-template (new-joint-template-factory parent root in-plug xyz-joints))
        (children (topology:children root)))
    (loop for child in children
          for sub-child-index from 0
          do (let ((child-template (build-joint-template-recursively root-template
                                                                     child
                                                                     in-plug
                                                                     xyz-joints)))
               (add-child root-template child-template)))
    root-template))

(defun build-joint-template (graph xyz-joints)
  (let ((root-node (topology:root-node graph)))
    (build-joint-template-recursively nil root-node (topology:in-plug graph) xyz-joints)))

(defun kin:joint-calculate-position-index (joint atom-table)
  (let ((atom-id (kin:joint/id joint)))
    (destructuring-bind (molecule-id residue-id atom-id)
        atom-id
      (let* ((mol (if (< molecule-id (length (chem:atom-table/atom-table-molecules atom-table)))
                      (aref (chem:atom-table/atom-table-molecules atom-table) molecule-id)
                      (error "The atom-table/atom-table-molecules is wrong")))
             (res (chem:content-at mol residue-id))
             (atm (chem:content-at res atom-id)))
        #+(or)(format t "joint-calculate-position-index number-of-atoms: ~a~%" (chem:atom-table/get-number-of-atoms atom-table))
        (unless (eq (chem:get-name atm) (kin:joint/name joint))
          (break "The joint ~a does not match the atom-table ~a" joint atom-table))
        (let* ((index*3 (let ((res (gethash atm (chem:atom-table/get-atom-table-indexes atom-table))))
                          (unless res
                            (error "Could not find atm ~a in atom-table ~a for joint ~a" atm atom-table joint))
                          res)))
          index*3)))))

(defun residue-shape-atom (nil-or-monomer-shape atom-name)
  (let* ((residue (residue nil-or-monomer-shape))
         (atm (chem:matter/first-atom-with-name residue atom-name)))
    (unless atm (error "Could not find atom with name ~s in ~s" atom-name residue))
    atm))

(defgeneric write-into-joint-tree (joint-template
                                   parent-joint
                                   atresidue
                                   atmolecule-index
                                   atresidue-index
                                   atom-table
                                   adjustments
                                   monomer
                                   nil-or-monomer-shape
                                   monomer-subset))

(defmethod write-into-joint-tree ((joint-template t)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset)
  (error "write-into-joint-tree - handle joint-template ~a" joint-template))

(defmethod write-into-joint-tree ((joint-template jump-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset)
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (etypecase nil-or-monomer-shape
                  (null (kin:make-jump-joint atomid atom-name atom-table))
                  (rotamer-shape (kin:make-jump-joint atomid atom-name atom-table))
                  (residue-shape (kin:make-xyz-joint atomid atom-name atom-table (residue-shape-atom nil-or-monomer-shape atom-name)))))
         )
    (put-joint atresidue joint constitution-atoms-index)
    (when parent-joint (kin:joint/add-child parent-joint joint))
    joint))

(defmethod write-into-joint-tree ((joint-template complex-bonded-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset
                                  )
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (input-stub-joints (input-stub-joints joint-template)))
    (flet ((make-complex-bonded-joint ()
             (let ((joint (kin:make-complex-bonded-joint atomid atom-name atom-table)))
               (put-joint atresidue joint constitution-atoms-index)
               (let ((input-stub0-template (aref input-stub-joints 0))
                     (input-stub1-template (aref input-stub-joints 1)))
                 (cond
                   ((and (null input-stub0-template) (null input-stub1-template))
                    ;; Do nothing
                    )
                   ((and input-stub0-template (null input-stub1-template))
                    (let* ((input-stub0-index (constitution-atoms-index input-stub0-template))
                           (input-stub0 (aref (joints atresidue) input-stub0-index)))
                      (kin:set-input-stub-joint1 joint input-stub0)))
                   ((and (null input-stub0-template) input-stub1-template)
                    (let* ((input-stub1-index (constitution-atoms-index input-stub1-template))
                           (input-stub1 (aref (joints atresidue) input-stub1-index)))
                      (kin:set-input-stub-joint1 joint input-stub1)))
                   (t
                    (let* ((input-stub0-index (constitution-atoms-index input-stub0-template))
                           (input-stub0 (aref (joints atresidue) input-stub0-index))
                           (input-stub1-index (constitution-atoms-index input-stub1-template))
                           (input-stub1 (aref (joints atresidue) input-stub1-index)))
                      (kin:set-input-stub-joint1 joint input-stub0)
                      (kin:set-input-stub-joint2 joint input-stub1))))
                 (when parent-joint (kin:joint/add-child parent-joint joint))
                 joint))))
      (etypecase nil-or-monomer-shape
        (null (make-complex-bonded-joint))
        (rotamer-shape (make-complex-bonded-joint))
        (residue-shape
         (let ((joint (kin:make-xyz-joint atomid atom-name atom-table (residue-shape-atom nil-or-monomer-shape atom-name))))
           (put-joint atresidue joint constitution-atoms-index)
           (when parent-joint (kin:joint/add-child parent-joint joint))
           joint))))))

(defmethod write-into-joint-tree ((joint-template bonded-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset
                                  )
  (when (and monomer-subset (null parent-joint))
    (error "Handle ~s joint if we have a monomer-subset" joint-template))
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint
           (etypecase nil-or-monomer-shape
             (null (kin:make-bonded-joint atomid atom-name atom-table))
             (rotamer-shape (kin:make-bonded-joint atomid atom-name atom-table))
             (residue-shape (kin:make-xyz-joint atomid atom-name atom-table (residue-shape-atom nil-or-monomer-shape atom-name))))))
    (put-joint atresidue joint constitution-atoms-index)
    (when parent-joint (kin:joint/add-child parent-joint joint))
    joint))

(defun create-stub-joint (monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
  (let* ((base-assembler (backbone-assembler monomer-subset))
         (base-coordinates (coordinates monomer-subset))
         (base-ataggregate (ataggregate base-assembler))
         (base-monomer-pos (gethash monomer (monomer-positions base-assembler)))
         (base-atresidue (at-position base-ataggregate base-monomer-pos))
         (base-joint0 (aref (joints base-atresidue) 0))
         (base-parent-joint (kin:parent base-joint0))
         (base-grand-parent-joint (kin:parent base-parent-joint))
         (base-great-grand-parent-joint (kin:parent base-grand-parent-joint))
         (base-j-atom (find-atom-for-joint base-assembler base-joint0))
         (base-pj-atom (find-atom-for-joint base-assembler base-parent-joint))
         (base-gpj-atom (find-atom-for-joint base-assembler base-grand-parent-joint))
         (base-ggpj-atom (find-atom-for-joint base-assembler base-great-grand-parent-joint))
         (constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (kin:make-stub-joint atomid atom-name atom-table base-j-atom base-pj-atom base-gpj-atom base-ggpj-atom))
         (backbone-atom-tree-spliced-parent (gethash monomer (monomers monomer-subset))))
    (put-joint atresidue joint constitution-atoms-index)
    (kin:joint/set-parent joint backbone-atom-tree-spliced-parent)
    joint))

(defun create-xyz-joint-connect-to-backbone (monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (kin:make-xyz-joint atomid atom-name atom-table))
         (backbone-atom-tree-spliced-parent (gethash monomer (monomers monomer-subset))))
    (put-joint atresidue joint constitution-atoms-index)
    (kin:joint/set-parent joint backbone-atom-tree-spliced-parent)
    joint))


(defun deal-with-pair-scan-sidechain (monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
  #+(or) (create-xyz-joint-connect-to-backbone monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
  (create-stub-joint monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
  )

(defmethod write-into-joint-tree ((joint-template xyz-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset
                                  )
  (if (and monomer-subset (null parent-joint))
      (deal-with-pair-scan-sidechain monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
      (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
             (atom-name (atom-name joint-template))
             (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
             (new-joint (kin:make-xyz-joint atomid atom-name atom-table)))
        (put-joint atresidue new-joint constitution-atoms-index)
        (when parent-joint (kin:joint/add-child parent-joint new-joint))
        new-joint)))

(defmethod write-into-joint-tree ((joint-template adjustable-bonded-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset
                                  )
  (when (and monomer-subset (null parent-joint))
    (error "Handle ~s joint if we have a monomer-subset" joint-template))
  (let ((joint (call-next-method)))
    (when adjustments
      (add-to-adjustments adjustments
                          (make-instance (adjustment joint-template) :joint joint)
                          atresidue))
    joint))

(defmethod write-into-joint-tree ((joint-template in-plug-bonded-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset
                                  )
  (if (and monomer-subset (null parent-joint))
      (deal-with-pair-scan-sidechain monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
      (let* ((joint (call-next-method))
             (in-plug (in-plug joint-template))
             (in-plug-name (name in-plug))
             (coupling-name (coupling-name in-plug-name))
             (out-plug-name (out-plug-name coupling-name)))
        (if (= (length (plug-bonds in-plug)) 1)
            (kin:joint/set-property joint :out-plug-name out-plug-name)
            (let ((indexed-out-plug-name (intern (format nil "~a.0" out-plug-name) :keyword)))
              (kin:joint/set-property joint :out-plug-name indexed-out-plug-name)))
        joint)))

(defmethod write-into-joint-tree ((joint-template adjustable-in-plug-bonded-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset
                                  )
  (if (and monomer-subset (null parent-joint))
      (deal-with-pair-scan-sidechain monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
      (let ((joint (call-next-method)))
        (when adjustments
          (add-to-adjustments adjustments
                              (make-instance (adjustment joint-template) :joint joint)
                              atresidue))
        joint)))

(defclass topology-graph ()
  ((topology :initarg :topology :reader topology)
   (nodes :initarg :nodes :reader nodes)
   (nodes-to-ids :initarg :nodes-to-ids :reader nodes-to-ids)
   (names-to-nodes :initarg :names-to-nodes :reader names-to-nodes)
   (undirected-edges :initarg :undirected-edges :reader undirected-edges)))

(defclass topology-edge ()
  ((from :initarg :from :reader from)
   (to :initarg :to :reader to)
   (label :initarg :label :initform nil :reader label)))

(defmethod topology.graphviz:make-graph ((obj topology:topology))
  (let* ((nodes (labels ((accumulate-nodes (joint-template)
                           (list* joint-template
                                  (loop for child in (topology:children joint-template)
                                        append (accumulate-nodes child)))))
                  (accumulate-nodes (topology::joint-template obj))))
         (nodes-to-ids (make-hash-table))
         (names-to-nodes (make-hash-table)))
    (loop for node in nodes
          do (setf (gethash node nodes-to-ids) (gensym))
          do (setf (gethash (atom-name node) names-to-nodes) node))
    (let* ((constitution (topology:constitution obj))
           (constitution-atoms (topology:constitution-atoms constitution))
           (constitution-bond-edges nil))
      (loop for catom across constitution-atoms
            for cbond-from-index from 0
            for from-atom-name = (atom-name catom)
            do (loop for cbond in (topology:bonds catom)
                     for cbond-to-index = (topology:to-atom-index cbond)
                     for catom-to = (aref constitution-atoms cbond-to-index)
                     for to-atom-name = (atom-name catom-to)
                     when (< cbond-from-index cbond-to-index)
                       do (let ((edge (make-instance 'topology-edge
                                                     :from from-atom-name
                                                     :to to-atom-name
                                                     :label (case (order cbond)
                                                              (:single-bond nil)
                                                              (:double-bond "=")
                                                              (otherwise (string (order cbond)))))))
                            (push edge constitution-bond-edges))))
      (make-instance 'topology-graph
                     :topology obj
                     :nodes nodes
                     :nodes-to-ids nodes-to-ids
                     :names-to-nodes names-to-nodes
                     :undirected-edges constitution-bond-edges
                     ))))

(defmethod topology.graphviz:graph-label ((graph topology-graph))
  (string (name (topology graph))))

(defmethod topology.graphviz:graph-name ((graph topology-graph))
  (string (name (topology graph))))

(defmethod topology.graphviz:nodes ((obj topology-graph))
  (nodes obj))

(defmethod topology.graphviz:node-label ((node topology:joint-template) graph)
  (declare (ignore graph))
  (format nil "~a" (symbol-name (topology::atom-name node))))

(defmethod topology.graphviz:node-id ((node topology:joint-template) graph)
  (format nil "~a" (symbol-name (gethash node (nodes-to-ids graph)))))


(defmethod topology.graphviz:directed-edges ((obj topology-graph))
  (labels ((accumulate-directed-edges (node)
             (loop for child in (children node)
                   append (list*
                            (make-instance 'topology-edge
                                           :from node
                                           :to child)
                            (accumulate-directed-edges child)))))
    (accumulate-directed-edges (joint-template (topology obj)))))

(defmethod topology.graphviz:directed-edge-from (edge) (from edge))
(defmethod topology.graphviz:directed-edge-to (edge) (to edge))
(defmethod topology.graphviz:directed-edge-label (edge) (declare (ignore edge)) nil)

(defmethod topology.graphviz:undirected-edges ((obj topology-graph))
  (undirected-edges obj))

(defmethod topology.graphviz:undirected-edge-from (edge (graph topology-graph))
  (let ((node (gethash (from edge) (names-to-nodes graph))))
    node))
(defmethod topology.graphviz:undirected-edge-to (edge (graph topology-graph))
  (let ((node (gethash (to edge) (names-to-nodes graph))))
    node))

(defmethod topology.graphviz:undirected-edge-label (edge (graph topology-graph))
  (declare (ignore graph))
  (label edge))



