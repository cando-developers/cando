(in-package :cando-build)


(defun make-atom (atom-name element &key (charge 0.0 chargep))
  "Docstring for make-atom"
  (let ((atm (chem:make-atom atom-name element)))
    (chem:set-charge atm charge)
    atm))


