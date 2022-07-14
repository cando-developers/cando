
(in-package :kin)

(defclass conformation ()
  ((oligomer :initarg :oligomer :accessor oligomer)
   (aggregate :initarg :aggregate :accessor aggregate)
   (ataggregate :initarg :ataggregate :accessor ataggregate)
   (joint-tree :initarg :joint-tree :accessor joint-tree)))

(cando:make-class-save-load conformation)

(defun make-conformation (oligomer)
  (let* ((aggregate (chem:make-aggregate))
         (ataggregate (make-instance 'ataggregate :aggregate aggregate))
         (joint-tree (make-joint-tree)))
    (resize-atmolecules ataggregate 1)
    (multiple-value-bind (molecule monomer-positions)
        (build-molecule oligomer)
      (when (chem:verbose 1)
        (format t "monomer-positions: ~a~%" (alexandria:hash-table-alist monomer-positions)))
      (multiple-value-bind (molecule molecule-index)
          (chem:add-matter aggregate molecule)
        ;; This is where I would invoke Conformation_O::buildMoleculeUsingOligomer
        ;; Use the monomers-to-topologys
        (let* ((atmolecule (build-atmolecule-using-oligomer oligomer molecule molecule-index monomer-positions joint-tree)))
          (put-atmolecule ataggregate atmolecule molecule-index)
          (let ((conf (make-instance 'conformation
                                     :oligomer oligomer
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
                       (set-position jnt (chem:get-position atm)))))

(defun copy-joint-positions-into-atoms (conformation)
  (walk-atoms-joints conformation
                     (lambda (atm jnt atomid)
                       (declare (ignore atomid))
                       (chem:set-position atm (get-position jnt)))))

(defun update-joint-tree-internal-coordinates (conformation)
  (let ((ataggregate (ataggregate conformation)))
    (walk-ataggregate-joints ataggregate
                             (lambda (joint atom-id)
                               (declare (ignore atom-id))
                               (update-internal-coord joint)))))
  
(defun build-all-atom-tree-external-coordinates (conformation)
  (let ((joint (root (joint-tree conformation))))
    (update-xyz-coords joint)))


(defun zero-all-joint-tree-external-coordinates (conf)
  "Set the external coordinates for each joint to the origin"
  (let ((ataggregate (ataggregate conf)))
    (walk-ataggregate-joints ataggregate
                             (lambda (joint atom-id)
                               (declare (ignore atom-id))
                               (set-position joint (geom:vec 0.0 0.0 0.0))))))




(defun build-structure (conf rotamer-library)
  (let* ((oligomer (oligomer conf))
         (monomers (chem:oligomer/monomers-as-list
                    oligomer))
         (atmolecule (kin:atmolecule-aref (kin:ataggregate conf) 0))
         )
    (loop for monomer in monomers
          for index from 0
          for context = (kin::monomer-context monomer)
          for rotamer-list = (gethash context (context-to-rotamer-list rotamer-library))
          for rotamer-index = (chem:get-monomer-index monomer)
          for rotamer = (elt (entries rotamer-list) rotamer-index)
          for atresidue = (kin::atresidue-aref atmolecule index)
          if (null rotamer-list)
            do (format t "Could not find context: ~s~%" context)
          else
            do (progn
                 #+(or)(progn
                         (format t "~a : ~a ~s~%" index monomer context)
                         (format t "      rotamer-index: ~a~%" rotamer-index)
                         (format t "      ~a~%" (length (entries rotamer-list)))
                         (format t "      ~a~%" atresidue)
                         (format t "   rotamer: ~a~%" (with-output-to-string (sout) (kin::print-internal rotamer sout))))
                 (apply-rotamer-to-atresidue atresidue rotamer))
          )
    (kin:build-all-atom-tree-external-coordinates conf)
    (kin:copy-joint-positions-into-atoms conf)
    monomers))

(defun increment-sequence (conformation)
  (let* ((next-oligomer (chem:oligomer/deep-copy (oligomer conformation)))
         (res (chem:increment-sequence next-oligomer)))
    (if res
        (values (make-conformation next-oligomer) res)
        (values conformation res))))


(defun build-all-sequences (oligomer rotamer-library sdf-filename)
  (chem:oligomer/goto-sequence oligomer 0)
  (with-open-file (fout sdf-filename :direction :output)
    (loop for seq-index below (chem:oligomer/number-of-sequences oligomer)
          for olig-temp = (chem:oligomer/deep-copy oligomer)
          do (chem:oligomer/goto-sequence olig-temp seq-index)
          do (format t "Generating sequence ~a~%" seq-index)
          do (let ((conf (make-conformation olig-temp)))
               (build-structure conf rotamer-library)
               (energy:minimize (kin:aggregate conf))
               (sdf:write-sdf-stream (kin:aggregate conf) fout)))))



