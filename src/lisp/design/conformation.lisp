
(in-package :kin)

(defclass conformation ()
  ((aggregate :initarg :aggregate :accessor aggregate)
   (ataggregate :initarg :ataggregate :accessor ataggregate)
   (joint-tree :initarg :joint-tree :accessor joint-tree)))

(defun make-conformation (oligomer)
  (let* ((aggregate (chem:make-aggregate))
         (ataggregate (make-instance 'ataggregate :aggregate aggregate))
         (joint-tree (make-joint-tree)))
    (resize-atmolecules ataggregate 1)
    (multiple-value-bind (molecule monomer-positions)
        (design:build-molecule oligomer)
      (format t "monomer-positions: ~a~%" (alexandria:hash-table-alist monomer-positions))
      (multiple-value-bind (molecule molecule-index)
          (chem:add-matter aggregate molecule)
        ;; This is where I would invoke Conformation_O::buildMoleculeUsingOligomer
        ;; Use the monomers-to-topologys
        (let* ((atmolecule (build-atmolecule-using-oligomer oligomer molecule molecule-index monomer-positions joint-tree)))
          (put-atmolecule ataggregate atmolecule molecule-index)
          (let ((conf (make-instance 'conformation
                                     :aggregate aggregate
                                     :ataggregate ataggregate
                                     :joint-tree joint-tree)))
            conf
            ))))))

(defun walk-atoms-joints (conformation callback)
  (let ((aggregate (aggregate conformation))
        (ataggregate (ataggregate conformation)))
    (loop for molecule-index below (chem:content-size aggregate)
          for molecule = (chem:content-at aggregate molecule-index)
          for atmolecule = (aref (atmolecules ataggregate) molecule-index)
          do (loop for residue-index below (chem:content-size molecule)
                   for residue = (chem:content-at molecule residue-index)
                   for atresidue = (aref (atresidues atmolecule) residue-index)
                   do (loop for atom-index below (chem:content-size residue)
                           for atom = (chem:content-at residue atom-index)
                           for joint = (aref (joints atresidue) atom-index)
                           do (funcall callback atom joint (list molecule-index residue-index atom-index)))))))

(defun copy-atom-positions-into-joints (conformation)
  (walk-atoms-joints conformation
                     (lambda (atm jnt atomid)
                       (declare (ignore atomid))
                       (kin:set-position jnt (chem:get-position atm)))))

(defun copy-joint-positions-into-atoms (conformation)
  (walk-atoms-joints conformation
                     (lambda (atm jnt atomid)
                       (declare (ignore atomid))
                       (chem:set-position atm (kin:get-position jnt)))))

(defun update-joint-tree-internal-coordinates (conformation)
  (let ((ataggregate (ataggregate conformation)))
    (walk-ataggregate-joints ataggregate
                             (lambda (joint atom-id)
                               (declare (ignore atom-id))
                               (kin:update-internal-coord joint)))))
  
(defun build-all-atom-tree-external-coordinates (conformation)
  (let ((joint (kin:root (joint-tree conformation))))
    (kin:update-xyz-coords joint)))


(defun zero-all-atom-tree-external-coordinates (conf)
  "Set the external coordinates for each joint to the origin"
  (let ((ataggregate (ataggregate conf)))
    (walk-ataggregate-joints ataggregate
                             (lambda (joint atom-id)
                               (declare (ignore atom-id))
                               (kin:set-position joint (geom:vec 0.0 0.0 0.0))))))

