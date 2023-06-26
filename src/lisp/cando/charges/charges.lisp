(in-package :charges)


(defun constrain-residue-charges (matter)
  "Adjust charges so that residues have the correct net charges.
Do this by calculating the deviation of the total charge for each residue from
its net-charge and then subtract a fraction of the deviation from each atomic
charge."
  (chem:do-residues (res matter)
    (let ((residue-charge 0.0)
          (atoms 0))
      (chem:do-atoms (atm res)
        (incf atoms)
        (incf residue-charge (chem:atom/get-charge atm)))
      (let ((adjust (/ (- residue-charge (float (chem:residue/get-net-charge res))) atoms)))
        (chem:do-atoms (atm res)
          (chem:atom/set-charge atm (- (chem:atom/get-charge atm) adjust)))))))


(defclass name-charge ()
  ((name :initarg :name :reader name)
   (charge :initarg :charge :reader charge)))

(defmethod print-object ((obj name-charge) stream)
  (format stream "(~s . ~,4f) " (name obj) (charge obj)))

(defclass residue-charge ()
  ((residue-name :initarg :residue-name :accessor residue-name)
   (name-charges :initarg :name-charges :accessor name-charges)))

(defmethod print-object ((obj residue-charge) stream)
  (format stream "(~a ~s)" (residue-name obj) (name-charges obj)))

(defun residue-charges (residue)
  (let ((atom-charges nil))
    (chem:do-atoms (atm residue)
      (push (make-instance 'name-charge
                           :name (chem:get-name atm)
                           :charge (chem:get-charge atm))
            atom-charges))
    (make-instance 'residue-charge
                   :residue-name (chem:get-name residue)
                   :name-charges (nreverse atom-charges))))

(defun dump-residue-charges (residue-charges &optional (stream t))
  (loop for residue-charge in residue-charges
        do (format stream "~a~%" residue-charge)))

(defun calculate-charge-trainers (foldamer charge-trainers)
  (let ((residue-charges nil))
    (loop for os-desc in charge-trainers
          do (multiple-value-bind (oligomer-space focus-monomer)
                 (topology.dag:parse-oligomer-space-dag foldamer os-desc)
               (multiple-value-bind (aggs monomer-to-mol-res-ids)
                   (topology:build-all-molecules oligomer-space)
                 (lparallel:pmapcar #'cando:starting-geometry aggs)
                 (lparallel:pmapcar #'charges:calculate-am1-bcc-charges aggs)
                 (mapc #'constrain-residue-charges aggs)
                 (loop for agg in aggs
                       for mol-res = (gethash focus-monomer monomer-to-mol-res-ids)
                       for residue-index = (topology:monomer-position-residue-index mol-res)
                       for residue = (chem:content-at (chem:content-at agg 0) residue-index)
                       do (push (residue-charges residue) residue-charges)))))
    residue-charges))

