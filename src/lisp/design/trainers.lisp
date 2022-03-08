(in-package :kin)


(defun monomer-context (monomer)
  "A context is a list of three names - the previous monomer name, the in coupling name and the current monomer name"
  (if (chem:has-in-coupling monomer)
      (let* ((in-coupling (chem:get-in-coupling monomer))
             (in-coupling-name (chem:get-name in-coupling))
             (previous-monomer (chem:get-source-monomer in-coupling))
             (previous-monomer-name (chem:current-stereoisomer-name previous-monomer))
             (current-monomer-name (chem:current-stereoisomer-name monomer)))
        (list previous-monomer-name in-coupling-name current-monomer-name))
      (list nil nil (chem:current-stereoisomer-name monomer))))

(defclass trainers ()
  ((oligomer :initarg :oligomer :accessor oligomer)
   (focus-residue-sequence-number :initarg :focus-residue-sequence-number :accessor focus-residue-sequence-number)
   (context :initarg :context :accessor context)
   (conformation :initarg :conformation :accessor conformation)
   (aggregate :initarg :aggregate :accessor aggregate)
   (atom-id-map :initarg :atom-id-map :accessor atom-id-map)
   (superposable-conformation-collection :initarg :superposable-conformation-collection
                                         :accessor superposable-conformation-collection)))

(cando:make-class-save-load trainers)

(defun collect-focus-joints (conformation focus-residue-index)
  (let ((focus-joints (make-hash-table)))
    ;; First put all of the joints in the focus-residue into the focus-joints hash-table
    (kin:walk (kin:joint-tree conformation)
              (lambda (joint)
                (let ((atom-id (kin:atom-id joint)))
                  (when (= (second atom-id) focus-residue-index)
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
      (let* ((parent (let ((relative (kin:get-parent root-joint)))
                       (if (typep relative 'kin:jump-joint)
                           nil
                           relative)))
             (grandparent (and parent
                               (let ((relative (kin:get-parent parent)))
                                 (if (typep relative 'kin:jump-joint)
                                     nil
                                     relative))))
             (great-grandparent (and grandparent
                                     (let ((relative (kin:get-parent grandparent)))
                                       (if (typep relative 'kin:jump-joint)
                                           nil
                                           relative)))))
        (when (and parent (kin:corresponds-to-joint parent))
          (setf (gethash parent focus-joints) t))
        (when (and grandparent (kin:corresponds-to-joint grandparent))
          (setf (gethash grandparent focus-joints) t))
        (when (and great-grandparent (kin:corresponds-to-joint great-grandparent))
          (setf (gethash great-grandparent focus-joints) t))))
    ;; Return the focus-joints as a list
    (alexandria:hash-table-keys focus-joints)))

(defun generate-superposable-conformation-collection (oligomer focus-residue-sequence-number context &key generate-atom-tree-dot)
  (let* ((conformation (kin:make-conformation oligomer))
         (aggregate (aggregate conformation))
         (_ (break "Check the aggregate ~a and conformation ~a" aggregate conformation))
         (focus-joints (collect-focus-joints conformation focus-residue-sequence-number))
         (atom-id-map (chem:build-atom-id-map aggregate))
         (superpose-atoms (loop for focus-joint in focus-joints
                                for atom-id = (kin:atom-id focus-joint)
                                for atom = (chem:lookup-atom atom-id-map atom-id)
                                when (not (eq (chem:get-element atom) :H))
                                  collect atom))
         ;; Now collect the atoms in the aggregate that correspond to the focus
         (superposable-conformation-collection (let ((scc (core:make-cxx-object 'chem:superposable-conformation-collection)))
                                                 (chem:set-rms-cut-off scc 1.0)
                                                 (chem:set-matter scc aggregate)
                                                 (loop for superpose-atom in superpose-atoms
                                                       do (chem:add-superpose-atom scc superpose-atom))
                                                 scc)))
    (when generate-atom-tree-dot
      (let* ((filename (format nil "/tmp/atom-tree-~a-~a-~a.dot"
                               (first context)
                               (second context)
                               (third context))))
        (format t "Wrote joint tree to ~a~%" filename)
        (design.graphviz-draw-joint-tree:draw-joint-tree (joint-tree conformation) filename)))
    (values conformation aggregate atom-id-map superposable-conformation-collection)))


(defun build-trainers-from-oligomers (oligomers)
  "Given a list of training oligomers expand them into a list of trainers"
  (let ((result (make-hash-table :test #'equalp))
        (all-oligomers (loop for oligomer in oligomers
                             append (loop for sequence-id below (chem:|Oligomer_O::numberOfSequences| oligomer)
                                          collect (progn
                                                    (chem:|Oligomer_O::gotoSequence| oligomer sequence-id)
                                                    (chem:|Oligomer_O::deepCopy| oligomer))))))
    (loop for oligomer in all-oligomers
          do (loop for monomer in (chem:|Oligomer_O::monomersAsList| oligomer)
                   for monomer-sequence-number = (chem:get-sequence-number monomer)
                   for context = (monomer-context monomer)
                   do (multiple-value-bind (conformation aggregate atom-id-map superposable-conformation-collection)
                          (generate-superposable-conformation-collection oligomer monomer-sequence-number (monomer-context monomer) :generate-atom-tree-dot nil)
                        (leap:assign-atom-types aggregate)
                        (setf (gethash context result) (make-instance 'trainers
                                                                      :oligomer oligomer
                                                                      :focus-residue-sequence-number monomer-sequence-number
                                                                      :context (monomer-context monomer)
                                                                      :conformation conformation
                                                                      :aggregate aggregate
                                                                      :atom-id-map atom-id-map
                                                                      :superposable-conformation-collection superposable-conformation-collection)))))
    (let ((trainers (alexandria:hash-table-values result)))
      trainers)))

(defun build-trainers (grammar)
  "Return a new design that includes trainers"
  (let* ((oligomers (design:build-training-oligomers grammar))
         (trainers (build-trainers-from-oligomers oligomers)))
    (make-instance 'grammar
                   :topologys (topologys grammar)
                   :cap-name-map (cap-name-map grammar)
                   :trainers trainers)))


(defun jostle-trainer (trainer &key test verbose)
  (let ((aggregate (aggregate trainer))
        (atom-id-map (atom-id-map trainer))
        (superposable-conformation-collection (superposable-conformation-collection trainer)))
    (cando:jostle aggregate 5.0)
    (unless test
      (when verbose (format t "Restraints off~%"))
      (energy:minimize aggregate :max-tn-steps 0
                                 :restraints-on nil
                                 :tn-tolerance 100.0
                                 :print-intermediate-results verbose)
      (when verbose (format t "Restraints on~%"))
      (energy:minimize aggregate :max-tn-steps 1000
                       :tn-tolerance 1.0e-12
                       :print-intermediate-results verbose))
    (prog1
        (chem:create-entry-if-conformation-is-new superposable-conformation-collection aggregate)
      (when verbose (format t "Number of entries: ~d~%" (chem:number-of-entries superposable-conformation-collection))))))


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
  "Extract internal coordinates for the focus monomer of the trainer"
  (let ((conformation-collection (superposable-conformation-collection trainer))
        (conformation-collection-agg (aggregate trainer))
        (conformation (conformation trainer))
        (atom-id-map (atom-id-map trainer))
        (focus-residue-sequence-number (focus-residue-sequence-number trainer)))
    (format t "extract-internal-coordinates ~%")
    (loop for index below (chem:number-of-entries conformation-collection)
          for entry = (chem:get-entry conformation-collection index)
          do (format t "index = ~a~%" index)
          do (chem:write-coordinates-to-matter entry conformation-collection-agg)
             (kin:walk-children (joint-tree conformation) 
                                (lambda (o) 
                                  (let* ((atom (chem:lookup-atom atom-id-map (kin:atom-id o)))
                                         (pos (chem:get-position atom)))
                                    (kin:set-position o pos))))
             (kin:update-internal-coords (joint-tree conformation))
          collect (let* ((internals (make-array 16 :adjustable t :fill-pointer 0))
                         (residue (lookup-residue-id (fold-tree conformation)
                                                     (list 0 focus-residue-sequence-number))))
                    (walk-joints residue
                                 (lambda (index joint)
                                   (format t "Extracting internal for ~a index: ~a atom-id: ~a~%"
                                           (kin:name joint)
                                           index
                                           (third (kin:atom-id joint)))
                                   (cond
                                     ((typep joint 'kin:bonded-joint)
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
                          (string< (format nil "~a" (reverse (car x)))
                                   (format nil "~a" (reverse (car y))))))))
      (loop for (context . coords ) in sorted
            do (format t "context: ~s~%" context)
               (loop with onec = (first coords)
                     for index = 0 then (+ index 5)
                     while (< index (length onec))
                     do (let ((atom-index (elt onec index))
                              (atom-name (elt onec (+ index 1))))
                          (format t "  ~a ~a~%" atom-index atom-name)))))))
