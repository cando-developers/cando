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





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; maybe-split-molecule
;;;
;;; Walk a spanning tree of the atoms in a molecule and split it
;;; into separate molecules and return them in a list.

(defun maybe-split-molecules-in-aggregate (aggregate)
  "* Arguments
- aggregate :: An aggregate.
* Description
Break up the molecules in the aggregate into a list of molecules using spanning trees and return that list.
"
  (let ((atoms-to-residues (make-hash-table))
        (residues-to-molecules (make-hash-table))
        (atoms-seen (make-hash-table))
        molecules
        atoms)
    (chem:map-residues
     nil
     (lambda (residue)
       (chem:map-atoms
        nil
        (lambda (atom)
          (setf (gethash atom atoms-to-residues) residue)
          (push atom atoms))
        residue))
     aggregate)
    (loop for atom in atoms
          when (not (gethash atom atoms-seen))
            do (push (span-across-molecule-from-atom atom atoms-seen residues-to-molecules atoms-to-residues) molecules))
;;; Remove the molecules in the aggregate
    (let ((mols (chem:map-molecules 'list #'identity aggregate)))
      (loop for mol in mols
            do (chem:remove-molecule aggregate mol)))
    ;; Add the new molecules
    (loop for mol in molecules
          do (chem:add-matter aggregate mol))
    aggregate))


(defun span-across-molecule-from-atom (atom
                                       atoms-seen
                                       residues-to-molecules
                                       atoms-to-residues)
  (let ((new-molecule (chem:make-molecule))
        (spanning-tree (chem:make-spanning-loop atom))
        residues)
    (loop
      (if (chem:advance-loop-and-process spanning-tree)
          (let* ((atom (chem:get-atom spanning-tree))
                 (residue (gethash atom atoms-to-residues))
                 (residue-molecule (gethash residue residues-to-molecules)))
            (setf (gethash atom atoms-seen) atom)
            (unless residue-molecule
              (push residue residues)
              (setf (gethash residue residues-to-molecules) new-molecule)))
          (return nil)))
    (let ((residue-list (sort residues #'< :key #'chem:get-id)))
      (loop for res in residue-list
            do (chem:add-matter new-molecule res)))
    (when (= (chem:content-size new-molecule) 1)
      (chem:set-name new-molecule (chem:get-name (chem:content-at new-molecule 0))))
    new-molecule))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Remove molecules of a particular type
;;;

(defun remove-molecules (aggregate molecule-type)
  (let (molecules)
    (do-molecules (molecule aggregate)
      (when (eq (chem:molecule-type molecule) 'cando:solvent)
        (push molecule molecules)))
    (loop for molecule in molecules
          do (chem:remove-molecule aggregate molecule))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Loops for molecules, residues, atoms
;;;


(defmacro do-molecules ((molecule-var matter) &body body)
  `(chem:map-molecules
    nil
    (lambda (,molecule-var)
      ,@body)
    ,matter))

(defmacro do-residues ((residue-var matter) &body body)
  `(chem:map-residues
    nil
    (lambda (,residue-var)
      ,@body)
    ,matter))

(defmacro do-atoms ((atom-var matter) &body body)
  `(chem:map-atoms
    nil
    (lambda (,atom-var)
      ,@body)
    ,matter))
