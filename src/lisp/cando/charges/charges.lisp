(in-package :charges)


(defun constrain-residue-charges (matter &key verbose)
  "Adjust charges so that residues have the correct net charges.
Do this by calculating the deviation of the total charge for each residue from
its net-charge and then subtract a fraction of the deviation from each atomic
charge."
  (chem:do-residues (res matter)
    (let ((target-residue-charge (chem:residue/get-net-charge res))
          (residue-charge 0.0)
          (atoms 0))
      (when verbose
        (format t "Constraining residue charge for ~s to ~d~%" res target-residue-charge))
      (chem:do-atoms (atm res)
        (incf atoms)
        (incf residue-charge (chem:atom/get-charge atm)))
      (when verbose
        (format t "Starting residue charge: ~f for ~d atoms~%" residue-charge atoms))
      (let ((adjust (/ (- residue-charge (float target-residue-charge)) atoms)))
        (when verbose
          (format t "Charge adjustment ~f~%" adjust))
        (chem:do-atoms (atm res)
          (chem:atom/set-charge atm (- (chem:atom/get-charge atm) adjust)))))))


(defclass name-charge (cando.serialize:serializable)
  ((name :initarg :name :reader name)
   (charge :initarg :charge :reader charge)))

(defmethod print-object ((obj name-charge) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (format stream "(~s . ~,4f) " (name obj) (charge obj)))))

(defclass residue-charge (cando.serialize:serializable)
  ((residue-name :initarg :residue-name :accessor residue-name)
   (name-charges :initarg :name-charges :accessor name-charges)))

(defmethod print-object ((obj residue-charge) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (format stream "(~a ~s)" (residue-name obj) (name-charges obj)))))

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

(defun calculate-charge-trainers (foldamer charge-trainers &key intermediate-file-directory)
  (let ((residue-charges nil))
    (loop for os-desc in charge-trainers
          do (multiple-value-bind (oligomer-space focus-monomer)
                 (topology.dag:parse-oligomer-space-dag foldamer os-desc)
               (format t "Building molecules for ~s~%" os-desc)
               (multiple-value-bind (aggs monomer-to-mol-res-ids)
                   (topology:build-all-molecules oligomer-space)
                 (let ((mol-res-id (gethash focus-monomer monomer-to-mol-res-ids)))
                   (format t "starting-geometry for ~s~%" os-desc)
                   (lparallel:pmapcar #'cando:starting-geometry aggs)
                   (when intermediate-file-directory
                     (loop for agg in aggs
                           for mol = (cando:mol agg 0)
                           for focus-res = (cando:res mol (topology:residue-index mol-res-id))
                           for focus-res-name = (string (chem:matter/get-name focus-res))
                           for filename = (merge-pathnames (make-pathname :name focus-res-name :type "mol2" :directory (list :relative focus-res-name)) intermediate-file-directory)
                           do (format t "Writing agg ~s to ~s~%" focus-res-name filename)
                           do (ensure-directories-exist filename)
                           do (cando:save-mol2 agg filename :use-sybyl-types t)
                           ))
                   (format t "calculate-am1-bcc-charges for ~s~%" os-desc)
                   (lparallel:pmapcar #'charges:calculate-am1-bcc-charges aggs)
                   (format t "Finishing ~s~%" os-desc)
                   (mapc #'constrain-residue-charges aggs)
                   (loop for agg in aggs
                         for mol-res = (gethash focus-monomer monomer-to-mol-res-ids)
                         for residue-index = (topology:residue-index mol-res)
                         for residue = (chem:content-at (chem:content-at agg 0) residue-index)
                         do (push (residue-charges residue) residue-charges))))))
    residue-charges))

