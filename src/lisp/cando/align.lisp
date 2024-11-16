(in-package :cando)


(defun align-atoms (&key moveable-atoms fixed-atoms moveable-matter)
  "Align the atoms in MOVEABLE-ATOMS-SEQUENCE to the FIXED-ATOMS-SEQUENCE and apply the transformation to MOVEABLE-MATTER (an aggregate, molecule, or residue)."
  (chem:align-matter moveable-matter moveable-atoms fixed-atoms))

(defun align-transform (&key moveable-atoms fixed-atoms)
  "Align the atoms in MOVEABLE-ATOMS-SEQUENCE to the FIXED-ATOMS-SEQUENCE and return the transform"
  (chem:align-transform moveable-atoms fixed-atoms))
