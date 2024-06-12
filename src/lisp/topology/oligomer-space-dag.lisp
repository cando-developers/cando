(in-package :topology.dag)

#|

Example of oligomer-space-dag


(topology:is-oligomer-space-supported
 spiros:*spiros*
 `(spiros::cala
   :amide spiros:apro4 (:side spiros:side)
   (:sideamide (spiros::sideamide :label :ringc))
   :dkp spiros:pro4 (:side spiros:side)
   :dkp spiros:pro4 (:side spiros:side)
   :dkp spiros::ampro
   :amide spiros::ala
   :amide spiros::apro4 (:side spiros:side)
   (:sideamide spiros:sideamide)
   :dkp spiros:pro4 (:side spiros:side)
   :dkp spiros:pro4 (:side spiros:side)
   :dkp spiros::ampro-ring ((ring :+ring :+ring :ringc))
   )
   )


|#

(define-condition no-matching-context ()
  ((context :initarg :context :accessor context))
  (:report (lambda (condition stream)
             (format stream "no-matching-context ~a" (context condition)))))

(defclass node (cando.serialize:serializable)
  ((id :initform (gensym) :initarg :id :accessor id)
   (index :initarg :index :reader index)
   (name :initarg :name :accessor name)
   (spanning-depth :initform nil :initarg :spanning-depth :initform nil :accessor spanning-depth)
   (label :initarg :label :accessor label)))

(defmethod print-object ((obj node) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t :identity t)
        (format stream "#~d/~a :spanning-depth ~a :label ~a" (index obj) (name obj) (spanning-depth obj) (label obj)))))

(defclass cap-node (node)
  ())

(defclass ring-cap-node (node)
  ())

(defclass edge (cando.serialize:serializable)
  ((from-node :initarg :from-node :accessor from-node)
   (to-node :initarg :to-node :accessor to-node)
   (name :initarg :name :accessor name)
   (raw-name :initarg :raw-name :accessor raw-name)))

(defmethod print-object ((obj edge) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (format stream "#~d/~a ~a #~d/~a"
                (index (from-node obj)) (name (from-node obj))
                (raw-name obj)
                (index (to-node obj)) (name (to-node obj))
                ))))

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
   (nodes :initform (make-array 16 :adjustable t :fill-pointer 0) :accessor nodes)
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

(defun topology:verify-oligomer-space (oligomer-space &key print)
  "Check that every monomer in the oligomer space has a monomer-context within the foldamer"
  (let ((used-contexts-set (make-hash-table :test 'equal))
        (foldamer (topology:foldamer oligomer-space)))
    (loop for monomer across (topology:monomers oligomer-space)
          for monomer-context = (topology:foldamer-monomer-context monomer oligomer-space foldamer)
          do (when print (format t "monomer-context = ~s~%" monomer-context))
          do (setf (gethash monomer-context used-contexts-set) t)
          if monomer-context
            do (when print
                 (format t "monomer-context ~a~%   is matched by ~a~%" (recursive-dump-local-monomer-context monomer nil 1) monomer-context))
          else
            do (error "Foldamer does not describe oligomer space"))
    (let (used-contexts unused-contexts)
      (maphash (lambda (key value)
                 (declare (ignore value))
                 (if (gethash key used-contexts-set)
                     (push key used-contexts)
                     (push key unused-contexts)))
               used-contexts-set)
      (values t used-contexts unused-contexts))))

(defun oligomer-space-from-dag (foldamer dag topology-groups)
  (let ((focus-node (root dag))
        (node-to-monomer (make-hash-table))
        (oligomer-space (make-instance 'topology:oligomer-space
                                       :foldamer foldamer)))
    (loop for node across (nodes dag)
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

(defun node-from-name (foldamer-name tree-name maybe-name label index)
  (cond
    ((and (consp maybe-name) (eq (car maybe-name) :cap))
     (when (keywordp (cadr maybe-name))
       (error "Parsing ~s ~s - the cap name ~s must be a non-keyword symbol" foldamer-name tree-name (cadr maybe-name)))
     (make-instance 'cap-node :name (cadr maybe-name) :label label :index index))
    ((and (consp maybe-name) (eq (car maybe-name) :ring-cap))
     (when (keywordp (cadr maybe-name))
       (error "Parsing ~s ~s - the ring-cap name ~s must be a non-keyword symbol" foldamer-name tree-name (cadr maybe-name)))
     (make-instance 'ring-cap-node :name (cadr maybe-name) :label label :index index))
    ((symbolp maybe-name)
     (make-instance 'node :name maybe-name :label label :index index))
    (t (error "Illegal name in node-from-name ~s" maybe-name))))

(defun parse-recursive (name tree-name sub-tree prev-node dag label)
  (cond
    ((null sub-tree))
    ((and (consp (car sub-tree))
          (consp (car (car sub-tree)))
          (string= (string (car (car (car sub-tree)))) "RING"))
     (break "break subtree starts with a ring")
     (parse-recursive name tree-name (car sub-tree) prev-node dag :head-cons )
     (parse-recursive name tree-name (cdr sub-tree) prev-node dag :tail-cons))
    ((consp (car sub-tree))
     (parse-recursive name tree-name (car sub-tree) prev-node dag :head-cons )
     (parse-recursive name tree-name (cdr sub-tree) prev-node dag :tail-cons))
    ((symbolp (car sub-tree))
     (let* ((plug-name (car sub-tree))
            (node (node-from-name name tree-name (cadr sub-tree) label (length (nodes dag)))))
       (cond
         ((topology:in-plug-name-p plug-name)
          (let ((edge (make-instance 'edge
                                     :raw-name plug-name
                                     :name (topology:coupling-name plug-name)
                                     :from-node node
                                     :to-node prev-node)))
            (vector-push-extend node (nodes dag))
            (push edge (edges dag))))
         ((topology:out-plug-name-p plug-name)
          (let ((edge (make-instance 'edge
                                     :raw-name plug-name
                                     :name (topology:coupling-name plug-name)
                                     :from-node prev-node
                                     :to-node node)))
            (vector-push-extend node (nodes dag))
            (push edge (edges dag))))
         (t (error "Illegal plug-name ~s" plug-name)))
       (parse-recursive name tree-name (cddr sub-tree) node dag :symbol)))
    (t (error "Illegal sub-tree ~s" sub-tree))))

(defun parse-label (head)
  head)

(defun parse-dag-for-oligomer-space (name tree &key label)
  (let* ((node (node-from-name name label (car tree) :top 0))
         (dag (make-instance 'dag :root node :label label)))
    (vector-push-extend node (nodes dag))
    (parse-recursive name label (cdr tree) node dag :top)
    (walk-spanning-tree dag)
    dag))

(defun parse-labeled-dag-for-oligomer-space (name labeled-tree)
  (error "Deprecated")
  (let ((label (parse-label (car labeled-tree)))
        (tree (cdr labeled-tree)))
    (parse-dag-for-oligomer-space name tree :label label)))

(defun parse-oligomer-space-dag (foldamer context &key (topology-groups topology:*topology-groups*))
  (let ((dag (parse-dag-for-oligomer-space (car context) context)))
    (multiple-value-bind (oligomer-space focus-monomer)
        (oligomer-space-from-dag foldamer dag topology-groups)
      (values oligomer-space focus-monomer dag))))

(defun parse-oligomer-space-labeled-dag (name foldamer labeled-context &key (topology-groups topology:*topology-groups*))
  (destructuring-bind (&key description pattern shape-kind)
      labeled-context
    (let ((dag (parse-dag-for-oligomer-space name pattern :label description)))
      (multiple-value-bind (oligomer-space focus-monomer)
          (oligomer-space-from-dag foldamer dag topology-groups)
        (values oligomer-space focus-monomer dag description shape-kind)))))



(defclass check-node ()
  ((node :initarg :node :accessor node)
   (node-index :initarg :node-index :reader node-index)
   (node-name :initarg :node-name :accessor node-name)
   (in-plug-names :initarg :in-plug-names :initform nil :accessor in-plug-names)
   (out-plug-names :initarg :out-plug-names :initform nil :accessor out-plug-names)))

(defmethod print-object ((obj check-node) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream ":node-name ~s :in-plug-names ~s  :out-plug-names ~s~%"
            (node-name obj)
            (in-plug-names obj)
            (out-plug-names obj))))

(defun interpret-node-name (name)
  (cond
    ((symbolp name)
     (gethash name topology:*topology-groups*))
    (t
     (error "interpret-node-name: What do we do with ~s" name))))

(defun compare-check-node-to-topology (check-node topology label)
  (when (in-plug-names check-node)
    (let ((in-plug-name (first (in-plug-names check-node))))
      (unless (gethash in-plug-name (topology:plugs topology))
        (error "Pattern ~s node ~a the in-plug name ~s was not found in the topology ~s"
               label
               (node-id check-node)
               in-plug-name topology))))
  (loop for out-plug-name in (out-plug-names check-node)
        for top-plug = (gethash out-plug-name (topology:plugs topology))
        unless top-plug
          do (error "Pattern ~s node ~a the out-plug name ~s was not found in the topology ~s"
                    label
                    (node-id check-node)
                    out-plug-name topology))
  (maphash (lambda (plug-name plug)
             (declare (ignore plug))
             (unless (topology:in-plug-name-p plug-name)
               (unless (member plug-name (out-plug-names check-node))
                 (error "Pattern ~s node ~a the topology ~s out-plug ~s was not found"
                        label
                        (node-id check-node)
                        topology plug-name))))
           (topology:plugs topology)))

(defun node-id (check-node)
  (format nil "#~d/~s" (node-index check-node) (node-name check-node)))

(defun add-unique-out-plug-name (plug-name check-node label)
  (when (member plug-name (out-plug-names check-node))
    (error "In the context ~s and node named ~s there is already an out-plug named ~s"
           label (node-id check-node) plug-name))
  (push plug-name (out-plug-names check-node)))

(defun add-unique-in-plug-name (plug-name check-node label)
  (when (> (length (in-plug-names check-node)) 0)
    (error "In the context ~s and node named ~s there is already an in-plug named ~s and you are trying to add ~s"
           label check-node (first (in-plug-names check-node)) plug-name))
  (push plug-name (in-plug-names check-node)))

(defun add-unique-plug-name (plug-name check-node label)
  (if (topology:in-plug-name-p plug-name)
      (add-unique-in-plug-name plug-name check-node label)
      (add-unique-out-plug-name plug-name check-node label)))


(defun validate-dag (dag)
  (let ((label (label dag))
        (node-ht (make-hash-table)))
    (loop for node across (nodes dag)
          for node-index from 0
          for node-name = (name node)
          for depth = (spanning-depth node)
          for check-node = (make-instance 'check-node
                                          :node-index node-index
                                          :node node
                                          :node-name node-name)
          do (setf (gethash node node-ht) check-node)
          do (cond
               ((= 1 depth)
                (when (typep node 'cap-node)
                  (warn "dag ~s has a cap-node ~s directly connected to the root" label node-name)))
               ((> depth 1)
                (when (not (typep node 'cap-node))
                  (warn "dag ~s has node ~s at level ~a that should be a cap node or eliminated" label node-name depth)))))
    (loop for edge in (edges dag)
          for from-node = (from-node edge)
          for from-check-node = (gethash from-node node-ht)
          for to-node = (to-node edge)
          for to-check-node = (gethash to-node node-ht)
          for from-edge-name = (raw-name edge)
          for to-edge-name = (topology:other-plug-name from-edge-name)
          if (topology:is-in-plug-name from-edge-name)
            do (progn
                 (add-unique-plug-name to-edge-name from-check-node label)
                 (add-unique-plug-name from-edge-name to-check-node label))
          else
            do (progn
                 (add-unique-plug-name from-edge-name from-check-node label)
                 (add-unique-plug-name to-edge-name to-check-node label)))
    ;; Debugging print check-nodes
    #+(or)
    (maphash (lambda (node check-node)
               (declare (ignore node))
               (format t "DEBUG: ~a :in-plug-names ~s  :out-plug-names ~s~%"
                       (node-id check-node)
                       (in-plug-names check-node)
                       (out-plug-names check-node)))
             node-ht)
    ;; Check that all check-nodes have zero or one in-plug
    ;; and that all out-plug-names are unique
    (maphash (lambda (node check-node)
               (when (> (length (in-plug-names check-node)) 1)
                 (error "The pattern ~s has node ~s that has more than 1 in-edge-names: ~s"
                        label (node-name check-node) (in-plug-names check-node)))
               ;; Sort the names and then compare adjacent pairs
               (let ((sorted-out-plug-names (sort (copy-seq (out-plug-names check-node)) #'string< :key #'string)))
                 (mapc (lambda (first next)
                         (when (eq first next)
                           (error "There are duplicate out-plug names for context: ~s node: ~s dup out-plug name: ~s" label (node-id check-node) first)))
                       sorted-out-plug-names (rest sorted-out-plug-names))))
             node-ht)
    ;; Check that all the check-nodes match a topology
    (maphash (lambda (node check-node)
               (let* ((node-name (node-name check-node))
                      (names (interpret-node-name node-name))
                      (in-plug-names (in-plug-names check-node))
                      (out-plug-names (out-plug-names check-node)))
                 (loop for name in names
                       for topology = (chem:find-topology name t)
                       do (compare-check-node-to-topology check-node topology label))))
             node-ht)))


