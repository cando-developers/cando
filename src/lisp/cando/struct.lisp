(in-package :struct)


(defun remove-atom (residue atoms)
  (dolist (a atoms)
    (chem:remove-atom-delete-bonds residue a)))
