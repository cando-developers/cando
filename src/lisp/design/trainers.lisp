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
              (great-grandparent (and grandparent (and grandparent (kin:get-parent grand-parent)))))
         (when (and parent (kin:corresponds-to-atom parent))
           (setf (gethash parent focus-joints) t))
         (when (and grandparent (kin:corresponds-to-atom grandparent))
           (setf (gethash grandparent focus-joints) t))
         (when (and great-grandparent (kin:corresponds-to-atom great-grandparent))
           (setf (gethash great-grandparent focus-joints) t))))
     ;; Return the focus-joints as a list
     (alexandria:hash-table-keys focus-joints)))

(defun augment-trainer-with-superposable-conformation-collection (trainer)
  (loop for oligomer = (oligomer trainer)
        for mol = (chem:get-molecule oligomer)
        for agg = (let ((agg (chem:make-aggregate)))
                    (chem:add-molecule agg mol)
                    agg)
        for conformation = (kin:make-conformation (list oligomer))
        for focus-joints = (collect-focus-joints conformation (focus-monomer-sequence-number trainer))
        for atom-id-map = (chem:build-atom-id-map agg)
        for superpose-atoms = (loop for focus-joint in focus-joints
                                    for atom = (chem:lookup-atom atom-id-map (kin:atom-id focus-joint))
                                    when (not (eq (chem:get-element atom) :H))
                                      collect atom)
        ;; Now collect the atoms in the aggregate that correspond to the focus
        for superposable-conformation-collection = (let ((scc (core:make-cxx-object 'chem:superposable-conformation-collection)))
                                                     (chem:set-rms-cut-off scc 1.0)
                                                     (chem:set-matter scc agg)
                                                     (loop for superpose-atom in superpose-atoms
                                                           do (chem:add-superpose-atom scc superpose-atom))
                                                     scc)
        do (setf (superposable-conformation-collection trainer) superposable-conformation-collection)
        ))


(defun build-trainers (oligomers)
  "Given a list of training oligomers expand them into a list of trainers"
  (let ((result (make-hash-table :test #'equalp))
        (all-oligomers (loop for oligomer in oligomers
                             append (loop for sequence-id below (chem:number-of-sequences oligomer)
                                          collect (progn
                                                    (chem:goto-sequence oligomer sequence-id)
                                                    (chem:copy oligomer))))))
    (loop for oligomer in all-oligomers
          append (loop for monomer in (chem:monomers-as-list oligomer)
                       for monomer-sequence-number = (chem:get-sequence-number monomer)
                       for context = (monomer-context monomer)
                       collect (make-instance 'trainers
                                               :oligomer oligomer
                                               :focus-monomer-sequence-number monomer-sequence-number
                                               :context (monomer-context monomer))))))

    

