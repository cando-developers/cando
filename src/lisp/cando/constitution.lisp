(in-package :cando)

(list-all-packages)

(defun make-topology-from-residue (residue &optional (pdb-name (chem:get-name residue)))
  (let* ((constitution-atoms (chem:make-constitution-atoms-from-residue residue))
         (stereoisomer (chem:make-stereoisomer (chem:get-name residue) pdb-name nil))
         (stereoinfo (chem:make-stereoinformation (list stereoisomer) nil))
         (constitution (chem:make-constitution (chem:get-name residue)
                                               ""
                                               constitution-atoms
                                               stereoinfo
                                               nil
                                               nil))
         (topology (chem:make-topology-from-residue residue (chem:get-name residue) constitution)))
    (chem:add-topology constitution topology)
    topology))
