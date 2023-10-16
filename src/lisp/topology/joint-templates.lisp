(in-package :topology)

(defclass joint-template ()
  ((parent :initform nil :initarg :parent :accessor parent)
   (atom-name :initarg :atom-name :accessor atom-name)
   (constitution-atoms-index :initarg :constitution-atoms-index :accessor constitution-atoms-index)
   ))

(cando.serialize:make-class-save-load joint-template
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (atom-name obj)))))

(defclass bonded-joint-template (joint-template)
  ((children :initform nil :initarg :children :accessor children)
   ))

(cando.serialize:make-class-save-load
 bonded-joint-template
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (atom-name obj)))))

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
         (jother (kin:joint/only-other-child jparent joint)))
    (setf (other adjustment) jother)))

(defmethod internal-adjust ((adjustment internal-planar-adjustment) assembler)
  (let* ((joint (joint adjustment))
         (other (other adjustment))
         (phi-other (kin:bonded-joint/get-phi other))
         (adjust-phi (radians-add phi-other PI)))
    (kin:bonded-joint/set-phi joint adjust-phi)))




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

(cando.serialize:make-class-save-load
 in-plug-bonded-joint-template
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (atom-name obj)))))

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

(cando.serialize:make-class-save-load
 complex-bonded-joint-template
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (atom-name obj)))) )

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

(defclass jump-joint-template (joint-template)
  ((children :initform nil :initarg :children :accessor children)))

(cando.serialize:make-class-save-load jump-joint-template)

(defun make-jump-joint-template (constitution-atoms-index &key atom-name)
  (make-instance 'jump-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name))

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

(defun new-joint-template-factory (parent-template node in-plug)
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
                            ((eql 0 (first child-indexes))
                             (list parent-template))
                            ((eql 1 (first child-indexes))
                             (list parent-template (sibling parent-template 0)))
                            (t
                             (list parent-template
                                   (sibling parent-template 1)
                                   (sibling parent-template 0))))))
         (make-complex-bonded-joint-template constitution-atoms-index
                                             :atom-name atom-name
                                             :stub-joints stub-joints)))
      ((null ggparent-template)
       (let ((stub-joints (cond
                            ((and (>= (length child-indexes) 2)
                                  (eql 0 (first child-indexes))
                                  (eql 0 (second child-indexes)))
                             (list parent-template
                                   gparent-template))
                            ((and (>= (length child-indexes) 2)
                                  (> 0 (first child-indexes))
                                  (eql 0 (second child-indexes)))
                             (list parent-template
                                   gparent-template
                                   (sibling parent-template 0)))
                            (t
                             (list parent-template
                                   gparent-template
                                   (sibling gparent-template 0))))))
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


(defun build-joint-template-recursively (parent root in-plug)
  (let ((root-template (new-joint-template-factory parent root in-plug))
        (children (topology:children root)))
    (loop for child in children
          for sub-child-index from 0
          do (let ((child-template (build-joint-template-recursively root-template
                                                                     child
                                                                     in-plug)))
               (add-child root-template child-template)))
    root-template))

(defun build-joint-template (graph)
  (let ((root-node (topology:root-node graph)))
    (build-joint-template-recursively nil root-node (topology:in-plug graph))))

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
        (let* ((index*3 (let ((res (gethash atm (chem:atom-table/get-atom-table-indices atom-table))))
                          (unless res
                            (error "Could not find atm ~a in atom-table ~a for joint ~a" atm atom-table joint))
                          res)))
          index*3)))))

(defgeneric write-into-joint-tree (joint-template parent-joint atresidue atmolecule-index atresidue-index atom-table adjustments one-orientation))

(defmethod write-into-joint-tree ((joint-template t) parent-joint atresidue atmolecule-index atresidue-index atom-table adjustments one-orientation)
  (error "write-into-joint-tree - handle joint-template ~a" joint-template))


(defmethod write-into-joint-tree ((joint-template jump-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  one-orientation)
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (if one-orientation
                    (kin:make-jump-joint atomid atom-name atom-table
                                         (parent-to-relative-transform one-orientation)
                                         (lab-frame one-orientation)
                                         )
                    (kin:make-jump-joint atomid atom-name atom-table nil nil))))
    (put-joint atresidue joint constitution-atoms-index)
    (when parent-joint (kin:joint/add-child parent-joint joint))
    joint))

(defmethod write-into-joint-tree ((joint-template complex-bonded-joint-template) parent-joint atresidue atmolecule-index atresidue-index atom-table adjustments one-orientation)
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (kin:make-complex-bonded-joint atomid atom-name atom-table))
         (input-stub-joints (input-stub-joints joint-template)))
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
      joint)))

(defmethod write-into-joint-tree ((joint-template bonded-joint-template) parent-joint atresidue atmolecule-index atresidue-index atom-table adjustments one-orientation)
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (kin:make-bonded-joint atomid atom-name atom-table)))
    (put-joint atresidue joint constitution-atoms-index)
    (when parent-joint (kin:joint/add-child parent-joint joint))
    joint))

(defmethod write-into-joint-tree ((joint-template adjustable-bonded-joint-template) parent-joint atresidue atmolecule-index atresidue-index atom-table adjustments one-orientation)
  (let ((joint (call-next-method)))
    (add-to-adjustments adjustments
                        (make-instance (adjustment joint-template) :joint joint)
                        atresidue)
    joint))

(defmethod write-into-joint-tree ((joint-template in-plug-bonded-joint-template) parent-joint atresidue atmolecule-index atresidue-index atom-table adjustments one-orientation)
  (let* ((joint (call-next-method))
         (in-plug (in-plug joint-template))
         (in-plug-name (name in-plug))
         (coupling-name (coupling-name in-plug-name))
         (out-plug-name (out-plug-name coupling-name)))
    (if (= (length (plug-bonds in-plug)) 1)
        (kin:joint/set-property joint :out-plug-name out-plug-name)
        (let ((indexed-out-plug-name (intern (format nil "~a.0" out-plug-name) :keyword)))
          (kin:joint/set-property joint :out-plug-name indexed-out-plug-name)))
    joint))

(defmethod write-into-joint-tree ((joint-template adjustable-in-plug-bonded-joint-template) parent-joint atresidue atmolecule-index atresidue-index atom-table adjustments one-orientation)
  (let ((joint (call-next-method)))
    (add-to-adjustments adjustments
                        (make-instance (adjustment joint-template) :joint joint)
                        atresidue)
    joint))

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



