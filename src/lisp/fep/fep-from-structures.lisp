(in-package :fep)

(defun only-residue (molecule)
  (unless (= (chem:content-size molecule) 1)
    (error "The molecule ~a must have only one residue" molecule))
  (chem:content-at molecule 0))

(defun fep-calculation-from-ligands (ligands)
  (let* ((fep-structures (loop for molecule in ligands
                              for residue = (only-residue molecule)
                              collect (make-instance 'simple-fep-structure
                                                     :name (chem:get-name molecule)
                                                     :drawing molecule
                                                     :core-residue residue
                                                     :core-residue-name (chem:get-name residue)
                                                     :core-atoms (chem:all-atoms-as-list molecule)
                                                     :molecule molecule
                                                     :net-charge 0.0 #|Not always!|#)))
         (calc (make-instance 'fep:fep-calculation
                              :ligands fep-structures)))
    calc))
                              
                              
                              
