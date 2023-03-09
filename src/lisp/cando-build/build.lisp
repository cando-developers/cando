(in-package :cando-build)


(defun make-atom (atom-name element &rest rest &key &allow-other-keys)
  "Make an atom and set all of its properties, next step needs to bond atom to something."
  (let ((atm (chem:make-atom atom-name element)))
    (alexandria:doplist (key value rest) atm
                        (leap.commands::set-property atm key value))
    atm))

(defun calc-mol-charge (molecule)
  "Calculate the charge on all atoms assigned to -molecule-, can be a single atom up to a whole molecule."
  (let ((charge-a 0.0))
    (chem:do-atoms (atm molecule)
      (incf charge-a (chem:get-charge atm)))
    charge-a))

(defun bond-atoms (atom1 atom2 bond-order)
  "Bond atoms 1 and 2 using specified bond-order. If atom is new, include it in a new or existing residue"
  (chem:bond-to atom1 atom2 bond-order))

(defun add-atom-to-residue (residue atom-name)
  "Add a single atom to an existing residue"
  (chem:add-atom residue atom-name))

(defun make-new-residue (name)
  "Create a new residue assigned -name-"
  (chem:make-residue name))


