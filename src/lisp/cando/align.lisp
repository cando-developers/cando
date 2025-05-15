(in-package :cando)


(defun align-atoms (&key moveable-atoms fixed-atoms moveable-matter)
  "Align the atoms in MOVEABLE-ATOMS-SEQUENCE to the FIXED-ATOMS-SEQUENCE and apply the transformation to MOVEABLE-MATTER (an aggregate, molecule, or residue)."
  (chem:align-matter moveable-matter moveable-atoms fixed-atoms))

(defun align-transform (&key moveable-atoms fixed-atoms)
  "Align the atoms in MOVEABLE-ATOMS-SEQUENCE to the FIXED-ATOMS-SEQUENCE and return the transform"
  (chem:align-transform moveable-atoms fixed-atoms))


(defun align-heavy-atoms (aggregate coordinates)
  "Return a vector of coordinates, one for each of COORDINATES.  Use the heavy atoms in AGGREGATE for the alignment"
  (let ((heavy-atoms-index3 (make-array 256 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (index3 0))
    ;; Build the indexes3
    (chem:do-molecules (mol aggregate)
      (chem:do-residues (res mol)
        (chem:do-atoms (atm res)
          (when (chem:element-is-hydrogen (chem:atom/get-element atm))
            (vector-push-extend index3 heavy-atoms-index3)
            (incf index3 3)
            ))))
    (setf heavy-atoms-index3 (copy-seq heavy-atoms-index3))
    (let* ((transforms (loop with fixed-coords = (aref coordinates 0)
                             for index from 0 below (length coordinates)
                             for moveable-coords = (aref coordinates index)
                             collect (chem:align-coordinates heavy-atoms-index3 fixed-coords moveable-coords)))
           (results (loop with results = (make-array (length transforms))
                          for coords across coordinates
                          for transform in transforms
                          for index from 0
                          for transformed = (chem:nvector-transform coords transform)
                          do (setf (aref results index) transformed)
                          finally (return results))))
      results)))
