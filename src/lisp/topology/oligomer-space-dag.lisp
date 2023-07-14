(in-package :topology.dag)


(define-condition no-matching-context ()
  ((context :initarg :context :accessor context))
  (:report (lambda (condition stream)
             (format stream "no-matching-context ~a" (context condition)))))

(defclass node ()
  ((id :initform (gensym) :initarg :id :accessor id)
   (name :initarg :name :accessor name)
   (spanning-depth :initform nil :initarg :spanning-depth :initform nil :accessor spanning-depth)
   (label :initarg :label :accessor label)))

(cando.serialize:make-class-save-load node
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t :identity t)
     (format stream "~a :spanning-depth ~a :label ~a" (name obj) (spanning-depth obj) (label obj)))))

(defclass cap-node (node)
  ())

(cando.serialize:make-class-save-load cap-node
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t :identity t)
     (format stream "~a :spanning-depth ~a :label ~a" (name obj) (spanning-depth obj) (label obj)))))

(defclass edge ()
  ((from-node :initarg :from-node :accessor from-node)
   (to-node :initarg :to-node :accessor to-node)
   (name :initarg :name :accessor name)
   (raw-name :initarg :raw-name :accessor raw-name)))

(cando.serialize:make-class-save-load edge
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a ~a" (from-node obj) (raw-name obj) (to-node obj)))))

(defun other-node (edge node)
  "If the edge contains node then return the other node, otherwise nil"
  (cond
    ((eq (to-node edge) node)
     (from-node edge))
    ((eq (from-node edge) node)
     (to-node edge))
    (t nil)))

(defclass dag ()
  ((label :initarg :label :accessor label)
   (root :initarg :root :accessor root)
   (nodes :initform nil :accessor nodes)
   (edges :initform nil :accessor edges)))

(cando.serialize:make-class-save-load dag
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~s" (label obj)))))

(defun recursive-spanning-tree (root depth dag seen)
  (if (gethash root seen)
      nil
      (progn
        (setf (spanning-depth root) depth
              (gethash root seen) t)
        (let ((neighbors (loop for edge in (edges dag)
                               for other-node = (other-node edge root)
                               when other-node
                                 collect other-node)))
          (append neighbors
                  (loop for neighbor in neighbors
                        append (recursive-spanning-tree neighbor (1+ depth) dag seen)))))))

#|
(defclass done-info ()
  ((number-of-atoms :initarg :number-of-atoms :accessor number-of-atoms)
   (build-info-msg :initarg :build-info-msg :accessor build-info-msg)
   (finished-steps :initarg

|#                   
(defun walk-spanning-tree (dag)
  (let ((root (root dag))
        (seen (make-hash-table)))
    (cons root (recursive-spanning-tree root 0 dag seen))))

(defun ensure-monomer (node map)
  (let ((monomer (gethash node map)))
    (unless monomer
      (error "Could not find monomer for node ~a" node))
    monomer))

(defun oligomer-space-from-dag (foldamer dag topology-groups)
  (let ((focus-node (root dag))
        (node-to-monomer (make-hash-table))
        (oligomer-space (make-instance 'topology:oligomer-space
                                       :foldamer foldamer)))
    (loop for node in (nodes dag)
          for name = (name node)
          for names = (gethash name topology-groups)
          for monomer = (make-instance 'topology:monomer
                                       :monomers names
                                       :id name)
          do (setf (gethash node node-to-monomer) monomer)
          do (vector-push-extend monomer (topology:monomers oligomer-space)))
    (loop for edge in (edges dag)
          for source-node = (from-node edge)
          for target-node = (to-node edge)
          for source-monomer = (ensure-monomer source-node node-to-monomer)
          for target-monomer = (ensure-monomer target-node node-to-monomer)
          for coupling-name = (name edge)
          for source-plug-name = (topology:out-plug-name coupling-name)
          for target-plug-name = (topology:in-plug-name coupling-name)
          for directional-coupling = (make-instance 'topology:directional-coupling
                                                    :name coupling-name
                                                    :source-monomer source-monomer
                                                    :target-monomer target-monomer
                                                    :source-plug-name source-plug-name
                                                    :target-plug-name target-plug-name)
          do (setf (gethash source-plug-name (topology:couplings source-monomer))
                   directional-coupling
                   (gethash target-plug-name (topology:couplings target-monomer))
                   directional-coupling)
          do (vector-push-extend directional-coupling (topology:couplings oligomer-space)))
    (values oligomer-space (gethash focus-node node-to-monomer))))

(defun node-from-name (foldamer-name tree-name maybe-name label)
  (cond
    ((and (consp maybe-name) (eq (car maybe-name) :cap))
     (when (keywordp (cadr maybe-name))
       (error "Parsing ~s ~s - the cap name ~s must be a non-keyword symbol" foldamer-name tree-name (cadr maybe-name)))
     (make-instance 'cap-node :name (cadr maybe-name) :label label))
    ((symbolp maybe-name)
     (make-instance 'node :name maybe-name :label label))
    (t (error "Illegal name in node-from-name ~s" maybe-name))))

(defun parse-recursive (name tree-name sub-tree prev-node dag label)
  (cond
    ((null sub-tree))
    ((consp (car sub-tree))
     (parse-recursive name tree-name (car sub-tree) prev-node dag :head-cons )
     (parse-recursive name tree-name (cdr sub-tree) prev-node dag :tail-cons))
    ((symbolp (car sub-tree))
     (let* ((plug-name (car sub-tree))
            (node (node-from-name name tree-name (cadr sub-tree) label)))
       (cond
         ((topology:in-plug-name-p plug-name)
          (let ((edge (make-instance 'edge
                                     :raw-name plug-name
                                     :name (topology:coupling-name plug-name)
                                     :from-node node
                                     :to-node prev-node)))
            (push node (nodes dag))
            (push edge (edges dag))))
         ((topology:out-plug-name-p plug-name)
          (let ((edge (make-instance 'edge
                                     :raw-name plug-name
                                     :name (topology:coupling-name plug-name)
                                     :from-node prev-node
                                     :to-node node)))
            (push node (nodes dag))
            (push edge (edges dag))))
         (t (error "Illegal plug-name ~s" plug-name)))
       (parse-recursive name tree-name (cddr sub-tree) node dag :symbol)))
    (t (error "Illegal sub-tree ~s" sub-tree))))

(defun parse-label (head)
  head)

(defun parse-dag-for-oligomer-space (name tree &key label)
  (let* ((node (node-from-name name label (car tree) :top))
         (dag (make-instance 'dag :root node :label label)))
    (push node (nodes dag))
    (parse-recursive name label (cdr tree) node dag :top)
    (walk-spanning-tree dag)
    dag))

(defun parse-labeled-dag-for-oligomer-space (name labeled-tree)
  (let ((label (parse-label (car labeled-tree)))
        (tree (cdr labeled-tree)))
    (parse-dag-for-oligomer-space name tree :label label)))

(defun parse-oligomer-space-dag (foldamer context &key (topology-groups topology:*topology-groups*))
  (let ((dag (parse-dag-for-oligomer-space context)))
    (multiple-value-bind (oligomer-space focus-monomer)
        (oligomer-space-from-dag foldamer dag topology-groups)
      (values oligomer-space focus-monomer dag))))

(defun parse-oligomer-space-labeled-dag (name foldamer labeled-context &key (topology-groups topology:*topology-groups*))
  (let ((dag (parse-labeled-dag-for-oligomer-space name labeled-context)))
    (multiple-value-bind (oligomer-space focus-monomer)
        (oligomer-space-from-dag foldamer dag topology-groups)
      (values oligomer-space focus-monomer dag))))

(defun validate-dag (dag)
  (let ((label (label dag)))
    (loop for node in (nodes dag)
          for node-name = (name node)
          for depth = (spanning-depth node)
          do (cond
               ((= 1 depth)
                (when (typep node 'cap-node)
                  (warn "dag ~s has a cap-node ~s directly connected to the root" label node-name)))
               ((> depth 1)
                (when (not (typep node 'cap-node))
                  (warn "dag ~s has node ~s at level ~a that should be a cap node or eliminated" label node-name depth)))))))




