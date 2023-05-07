(in-package :topology)


(defun make-constitution-atom-from-atom (atom constitution-atom-index)
  (make-instance 'constitution-atom
                 :atom-name (chem:atom-name atom)
                 :element (chem:get-element atom)
                 :properties (chem:matter/properties atom)
                 :index constitution-atom-index))

(defun define-constitution-atom-bonding (constitution-atom atm atom-to-index-map)
  (let ((constitution-bonds (loop for bond-index below (chem:number-of-bonds atm)
                                  for other-atom = (chem:atom/bonded-neighbor atm bond-index)
                                  for order = (chem:atom/bonded-order atm bond-index)
                                  for other-atom-index = (gethash other-atom atom-to-index-map)
                                  when other-atom-index
                                  collect (make-instance 'constitution-bond
                                                         :to-atom-index other-atom-index
                                                         :order order))))
    (setf (bonds constitution-atom) constitution-bonds)))

(defun ensure-all-atom-names-are-unique (residue)
  (let ((ht (make-hash-table)))
    (chem:do-atoms (atm residue)
      (when (gethash (chem:atom/atom-name atm) ht)
        (error "The atom name ~a is not unique" (chem:atom/atom-name atm))))))
      
(defun make-constitution-atoms-from-residue (residue &optional verbose)
  (ensure-all-atom-names-are-unique residue)
  (let ((caindex -1)
        (constitution-atoms (make-array 16 :adjustable t :fill-pointer 0))
        (atoms-to-index-map (make-hash-table)))
    (chem:do-atoms (atm residue)
      (let ((constitution-atom (make-constitution-atom-from-atom atm (incf caindex))))
        (vector-push-extend constitution-atom constitution-atoms)
        (setf (gethash atm atoms-to-index-map) caindex)))
    (chem:do-atoms (atm residue)
      (let* ((cai (gethash atm atoms-to-index-map))
             (constitution-atom (elt constitution-atoms cai)))
        (define-constitution-atom-bonding constitution-atom atm atoms-to-index-map)))
    constitution-atoms))

(defun constitution-atom-names-with-dihedrals (constitution)
  (let ((constitution-atoms (constitution-atoms constitution)))
    (loop for ca across constitution-atoms
          when (getf (properties ca) :dihedral)
            collect (atom-name ca))))
