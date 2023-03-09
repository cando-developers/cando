(in-package :cando)

(defun make-simple-topology-from-residue (residue &key (pdb-name (chem:get-name residue)) verbose)
  "Build everything necessary to create a simple topology from a residue."
  (or pdb-name (error "The residue must have a non-nil name"))
  (let* ((constitution-atoms (topology:make-constitution-atoms-from-residue residue verbose))
         (stereoisomer (make-instance 'topology:stereoisomer
                                      :name (chem:get-name residue)
                                      :pdb pdb-name
                                      :stereoconfigurations nil))
         (constitution (make-instance 'topology:constitution
                                      :name pdb-name
                                      :constitution-atoms constitution-atoms))
         (topology (topology:make-topology-from-residue residue pdb-name constitution)))
    (values topology constitution)))
