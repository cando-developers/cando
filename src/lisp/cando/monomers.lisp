
(in-package :cando)


(defun constitution-from-residue (residue)
  (let ((constitution-atoms (chem:make-constitution-atoms-from-residue residue)))
    (print constitution-atoms)))



