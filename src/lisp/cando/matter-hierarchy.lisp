;;;
;;; Code to deal with the aggregate/molecule/residue/atom hierarchy


(in-package :cando)


;;; ------------------------------------------------------------
;;;
;;; build-molecules-from-atom-connectivity
;;;
(defun build-molecules-from-atom-connectivity (atoms residues)
  "* Exported
* Arguments
- atoms :: Hash-table of atom ids to atoms.
- residues :: Hash-table of residue ids to residues.
* Description
Return a hash table of molecule-ids to molecules that contains molecules identified by walking the
spanning trees of the atoms and identifying the residues that belong to
those atoms."
  (let ((atoms-seen (make-hash-table :test #'eq))
        (atoms-to-residues (create-atoms-to-residues residues))
        (residues-to-molecules (make-hash-table :test #'eq))
        (molecules (make-hash-table :test #'eql))
        (molecule-id 0))
    (maphash (lambda (id atom)
               (unless (gethash atom atoms-seen)
                 ;; Start a new molecule
                 (let ((mol (chem:make-molecule)))
                   (walk-atom-spanning-tree atom mol
                                            atoms-seen
                                            residues-to-molecules
                                            atoms-to-residues)
                   (setf (gethash (incf molecule-id) molecules) mol))))
             atoms)
    molecules))
             


(defun create-atoms-to-residues (residues)
  "* Arguments
- residues :: A hashtable of residue ids to residues.
* Description
Return a hash table that maps atoms to residues that contain them."
  (let ((a2r (make-hash-table :test #'eq)))
    (maphash (lambda (ri res)
               (chem:map-atoms 'vector (lambda (a)
                                         (setf (gethash a a2r) res))
                               res))
             residues)
    a2r))


(defun walk-atom-spanning-tree (atom molecule
                                atoms-seen
                                residues-to-molecules
                                atoms-to-residues)
  "* Arguments
- atom :: An atom.
- molecule :: A molecule.
- atoms-seen :: A hash-table of atoms to T.
- residues-to-molecules :: A hash-table of residues to molecules.
- atoms-to-residues :: A hash-table of atoms to the residues that contain them.
* Description
Walk the spanning tree of the ATOM and add the ATOM to ATOMS-SEEN.
For each atom walked look up the residue in ATOMS-TO-RESIDUES
and if the residue is not in RESIDUES-TO-MOLECULES add it."
  (let ((spanning-loop (chem:make-spanning-loop atom)))
    (loop for a = (chem:next spanning-loop (lambda (a b) t))
       while a
       do (unless (gethash a atoms-seen)
            (setf (gethash a atoms-seen) t)
            (let ((residue (gethash a atoms-to-residues)))
              (or residue (error "The residue could not be found for atom: ~a" a))
              (unless (gethash residue residues-to-molecules)
                (setf (gethash residue residues-to-molecules) molecule)
                (chem:add-matter molecule residue)))))))


