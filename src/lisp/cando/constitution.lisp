(in-package :cando)

(defun make-topology-from-residue (residue &key (pdb-name (chem:get-name residue)) verbose)
  (or pdb-name (error "The residue must have a non-nil name"))
  (let* ((constitution-atoms (chem:make-constitution-atoms-from-residue residue verbose))
         (stereoisomer (chem:make-stereoisomer (chem:get-name residue) pdb-name nil))
         (stereoinfo (chem:make-stereoinformation (list stereoisomer) nil))
         (constitution (chem:make-constitution pdb-name
                                               ""
                                               constitution-atoms
                                               stereoinfo
                                               nil
                                               nil))
         (topology (chem:make-topology-from-residue residue pdb-name constitution)))
    (chem:add-topology constitution topology)
    (values topology constitution)))
