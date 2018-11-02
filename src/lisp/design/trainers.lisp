(in-package :design)


(defun monomer-context (monomer)
  (if (chem:has-in-coupling monomer)
      (let* ((in-coupling (chem:get-in-coupling monomer))
             (in-coupling-name (chem:get-name in-coupling))
             (previous-monomer (chem:get-in-monomer in-coupling))
             (previous-monomer-name (chem:current-stereoisomer-name previous-monomer))
             (current-monomer-name (chem:current-stereoisomer-name monomer)))
        (list previous-monomer-name in-coupling-name current-monomer-name))
      (list nil nil (chem:current-stereoisomer-name monomer))))

(defclass trainers ()
  ((oligomer :initarg :oligomer :accessor oligomer)
   (focus-monomer-sequence-number :initarg :focus-monomer-sequence-number :accessor focus-monomer-sequence-number)
   (context :initarg :context :accessor context)
   (conformation :initarg :conformation :accessor conformation)
   (aggregate :initarg :aggregate :accessor aggregate)
   (atom-id-map :initarg :atom-id-map :accessor atom-id-map)
   (superposable-conformation-collection :accessor superposable-conformation-collection)))

(defun collect-focus-joints (conformation focus-monomer-index)
   (let ((focus-joints (make-hash-table)))
     ;; First put all of the joints in the focus-monomer into the focus-joints hash-table
    (kin:walk (kin:get-atom-tree conformation)
              (lambda (joint)
                (let ((atom-id (kin:atom-id joint)))
                  (when (= (second atom-id) focus-monomer-index)
                    (setf (gethash joint focus-joints) t)))))
     ;; Then find the root joint
     (let ((root-joint
             (block root-find
               (alexandria:maphash-keys
                (lambda (joint)
                  (let ((parent-joint (kin:get-parent joint)))
                    (when (null (gethash parent-joint focus-joints))
                      (return-from root-find joint))))
                focus-joints))))
       ;; Add the parent, grandparent and great-grandparent if they exist and correspond to atoms
       (let* ((parent (kin:get-parent root-joint))
              (grandparent (and parent (kin:get-parent parent)))
              (great-grandparent (and grandparent (and grandparent (kin:get-parent grandparent)))))
         (when (and parent (kin:corresponds-to-atom parent))
           (setf (gethash parent focus-joints) t))
         (when (and grandparent (kin:corresponds-to-atom grandparent))
           (setf (gethash grandparent focus-joints) t))
         (when (and great-grandparent (kin:corresponds-to-atom great-grandparent))
           (setf (gethash great-grandparent focus-joints) t))))
     ;; Return the focus-joints as a list
     (alexandria:hash-table-keys focus-joints)))

(defun augment-trainer-with-superposable-conformation-collection (trainer &key generate-atom-tree-dot)
  (let* ((oligomer (oligomer trainer))
         (mol (chem:get-molecule oligomer))
         (agg (let ((agg (chem:make-aggregate)))
                (chem:add-molecule agg mol)
                agg))
         (conformation (kin:make-conformation (list oligomer)))
         (focus-joints (collect-focus-joints conformation (focus-monomer-sequence-number trainer)))
         (atom-id-map (chem:build-atom-id-map agg))
         (superpose-atoms (loop for focus-joint in focus-joints
                                for atom = (chem:lookup-atom atom-id-map (kin:atom-id focus-joint))
                                when (not (eq (chem:get-element atom) :H))
                                  collect atom))
         ;; Now collect the atoms in the aggregate that correspond to the focus
         (superposable-conformation-collection (let ((scc (core:make-cxx-object 'chem:superposable-conformation-collection)))
                                                 (chem:set-rms-cut-off scc 1.0)
                                                 (chem:set-matter scc agg)
                                                 (loop for superpose-atom in superpose-atoms
                                                       do (chem:add-superpose-atom scc superpose-atom))
                                                 scc)))
    (when generate-atom-tree-dot
      (let* ((context (design::context trainer))
             (filename (format nil "/tmp/atom-tree-~a-~a-~a.dot"
                               (first context)
                               (second context)
                               (third context))))
        (design.graphviz-draw-joint-tree:draw-joint-tree (kin:get-atom-tree conformation)
                                                         filename)))
    (setf (conformation trainer) conformation
          (aggregate trainer) agg
          (atom-id-map trainer) atom-id-map
          (superposable-conformation-collection trainer) superposable-conformation-collection)))


(defun build-trainers (oligomers)
  "Given a list of training oligomers expand them into a list of trainers"
  (let ((result (make-hash-table :test #'equalp))
        (all-oligomers (loop for oligomer in oligomers
                             append (loop for sequence-id below (chem:number-of-sequences oligomer)
                                          collect (progn
                                                    (chem:goto-sequence oligomer sequence-id)
                                                    (chem:copy oligomer))))))
    (loop for oligomer in all-oligomers
          do (loop for monomer in (chem:monomers-as-list oligomer)
                       for monomer-sequence-number = (chem:get-sequence-number monomer)
                   for context = (monomer-context monomer)
                   do (format t "context: ~a~%" context)
                   do (setf (gethash context result) (make-instance 'trainers
                                                                    :oligomer oligomer
                                                                    :focus-monomer-sequence-number monomer-sequence-number
                                                                    :context (monomer-context monomer)))))
    (alexandria:hash-table-values result)))

    


(defun jostle-trainer (trainer &key test)
  (let ((aggregate (aggregate trainer))
        (atom-id-map (atom-id-map trainer))
        (superposable-conformation-collection (superposable-conformation-collection trainer)))
    (cando:jostle aggregate)
    (unless test (energy:minimize aggregate :max-tn-steps 500))
    (chem:create-entry-if-conformation-is-new superposable-conformation-collection aggregate)))


(defun combine-all-coordinates (trainer)
  (let ((conf-col (superposable-conformation-collection trainer))
        (conf-col-agg (aggregate trainer))
        (agg (chem:make-aggregate)))
    (loop for index below (chem:number-of-entries conf-col)
          for entry = (chem:get-entry conf-col index)
          do (chem:write-coordinates-to-matter entry conf-col-agg)
          do (let ((new-mol (chem:matter-copy (chem:content-at conf-col-agg 0))))
               (chem:add-matter agg new-mol)))
    agg))

(defun extract-internal-coordinates (trainer)
  (let ((conf-col (superposable-conformation-collection trainer))
        (conf-col-agg (aggregate trainer))
        (conformation (conformation trainer))
        (atom-id-map (atom-id-map trainer))
        (focus-monomer-sequence-number (focus-monomer-sequence-number trainer)))
    (loop for index below (chem:number-of-entries conf-col)
          for entry = (chem:get-entry conf-col index)
          do (chem:write-coordinates-to-matter entry conf-col-agg)
             (kin:walk (kin:get-atom-tree conformation) 
                       (lambda (o) 
                         (let* ((atom (chem:lookup-atom atom-id-map (kin:atom-id o)))
                                (pos (chem:get-position atom)))
                           (kin:set-position o pos))))
             (kin:update-internal-coords (kin:get-atom-tree conformation))
          collect (let* ((internals (make-array 16 :adjustable t :fill-pointer 0))
                         (monomer (kin:lookup-monomer-id (kin:get-fold-tree conformation)
                                                         (list 0 focus-monomer-sequence-number))))
                    (format t "Extracting internals for monomer: ~s parentPlugName: ~s~%"
                            monomer
                            (kin:parent-plug-name monomer))
                    (kin:walk-joints monomer
                                     (lambda (index joint)
                                       (format t "Extracting internal for ~a index: ~a atom-id: ~a~%"
                                               (kin:name joint)
                                               index
                                               (third (kin:atom-id joint)))
                                       (cond
                                         ((typep joint 'kin:bonded-atom)
                                          (vector-push-extend (third (kin:atom-id joint)) internals)
                                          (vector-push-extend (kin:name joint) internals)
                                          (vector-push-extend (kin:get-distance joint) internals)
                                          (vector-push-extend (kin:get-theta joint) internals)
                                          (vector-push-extend (kin:get-phi joint) internals)))))
                    internals))))


(defun describe-internal-coordinates (coordinates)
  (let (context-coords)
    (maphash (lambda (key coords)
               (push (cons key coords) context-coords))
             coordinates)
    (let ((sorted (sort context-coords
                        (lambda (x y)
                          (string< (format nil "~a" (reverse key))
                                   (format nil "~a" (reverse key)))))))
      (loop for (context . coords ) in sorted
            do (format t "context: ~s~%" context)
               (loop with onec = (first coords)
                     for index = 0 then (+ index 5)
                     while (< index (length onec))
                     do (let ((atom-index (elt onec index))
                              (atom-name (elt onec (+ index 1))))
                          (format t "  ~a ~a~%" atom-index atom-name)))))))
