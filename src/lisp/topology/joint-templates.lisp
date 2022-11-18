(in-package :topology)

(defclass joint-template ()
  ((parent :initform nil :initarg :parent :accessor parent)
   (atom-name :initarg :atom-name :accessor atom-name)
   (constitution-atoms-index :initarg :constitution-atoms-index :accessor constitution-atoms-index)
   ))

(defmethod print-object ((obj joint-template) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (atom-name obj))))

(defclass bonded-joint-template (joint-template)
  ((children :initform nil :initarg :children :accessor children)
   ))

(defun make-bonded-joint-template (constitution-atoms-index &key atom-name parent)
  (make-instance 'bonded-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent))

(defclass in-plug-bonded-joint-template (bonded-joint-template)
  ((in-plug :initarg :in-plug :accessor in-plug)))

(defun make-in-plug-bonded-joint-template (constitution-atoms-index &key atom-name parent in-plug)
  (make-instance 'in-plug-bonded-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent
                 :in-plug in-plug))

(defclass complex-bonded-joint-template (bonded-joint-template)
  ((input-stub-joints :initform (make-array 2) :initarg :input-stub-joints :accessor input-stub-joints)))

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

(defun make-jump-joint-template (constitution-atoms-index &key atom-name)
  (make-instance 'jump-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name))

(defun add-child (joint-template child-template)
  (setf (parent child-template) joint-template)
  (setf (children joint-template) (append (children joint-template) (list child-template))))

(defun sibling (joint-template index)
  (nth index (children joint-template)))

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
       (make-in-plug-bonded-joint-template constitution-atoms-index
                                           :atom-name atom-name
                                           :parent nil
                                           :in-plug in-plug))
      ((typep in-plug 'topology:in-plug)
       (make-bonded-joint-template constitution-atoms-index
                                   :atom-name atom-name
                                   :parent parent-template))
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
      (t (make-bonded-joint-template constitution-atoms-index
                                                       :atom-name atom-name
                                                       :parent parent-template)))))


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
    (format t "root-node = ~a~%" root-node)
    (build-joint-template-recursively nil root-node (topology:in-plug graph))))

(defun topologies-from-graph (graph)
  (let* ((constitution (topology:constitution-from-graph graph))
         (plugs (topology:plugs constitution))
         (stereo-information (topology:stereo-information constitution))
         (tops (loop for stereoisomer in stereo-information
                     for name = (topology:name stereoisomer)
                     for configurations = (topology:configurations stereoisomer)
                     for joint-template = (build-joint-template graph)
                     for topology = (make-instance 'topology:topology
                                                   :name name
                                                   :constitution constitution
                                                   :plugs plugs
                                                   :joint-template joint-template
                                                   :stereoisomer-atoms configurations)
                     do (push topology (topology:topology-list constitution))
                     do (format t "stereoisomer ~a ~s~%" name configurations)
                     do (setf (topology:property-list topology) (list* :joint-template joint-template (topology:property-list topology)))
                     do (cando:register-topology topology name)
                     collect topology)))
    (setf (topology:topology-list constitution) tops)
    tops))
