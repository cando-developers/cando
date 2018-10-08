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
   (conformation :initarg :conformation :accessor conformation)))

(defun collect-focus-joints (conformation focus-monomer-index)
  (let ((focus-atoms (make-hash-table)))
    (kin:walk (kin:get-atom-tree conformation)
              (lambda (joint)
                (let ((atom-id (kin:atom-id joint)))
                  (when (= (second atom-id) focus-monomer-index)
                    (setf (gethash joint focus-atoms) t)
                    (let ((parent-joint (kin:get-parent joint)))
                      (if (kin:corr
                      (format t "Walking joints atom-id: ~a  joint: ~a parent: ~a~%" atom-id joint parent-joint))))))))

(defun build-subset-molecules (trainers &optional (start 0) (end (length trainers)))
  (loop for index from start below end
        for trainer = (elt trainers index)
        for oligomer = (oligomer trainer)
        for mol = (chem:get-molecule oligomer)
        for agg = (let ((agg (chem:make-aggregate)))
                    (chem:add-molecule agg mol)
                    agg)
        for conformation = (kin:make-conformation (list oligomer))
        for focus-joints = (collect-focus-joints conformation (focus-monomer-sequence-number trainer))
        for atom-id-map = (chem:build-atom-id-map agg)
        ;; Now collect the atoms in the aggregate that correspond to the focus
        for superposable-conformation-collection = (let ((scc (core:make-cxx-object 'chem:superposable-conformation-collection)))
                                                     (chem:set-rms-cut-off scc 1.0)
                                                     (chem:set-matter scc agg)
                                                     scc)
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

    

