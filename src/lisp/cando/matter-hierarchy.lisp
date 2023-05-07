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
  (let ((spanning-loop (chem:spanning-loop/make atom)))
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

(defun merge-molecules (aggregate mol-dest mol-source)
  #+(or)(format *debug-io* "Merging the molecules ~s ~s~%" mol-dest mol-source)
  (let ((source-residues (chem:map-residues 'list #'identity mol-source)))
    (loop for res in source-residues
          do (chem:add-matter mol-dest res)
             (chem:remove-residue mol-source res)))
  (chem:remove-molecule aggregate mol-source))

(defun maybe-join-molecules-in-aggregate (aggregate)
  "* Arguments
- aggregate :: An aggregate.
* Description
The aggregate may have chem:molecules that have bonds between them.
If this is the case then merge the residues of the two bonded molecules into one molecule.
This can happen if there are disulphide-bonds between chains."
  (tagbody
   top
     (let ((atoms-to-molecules (make-hash-table)))
       (chem:map-molecules
        nil
        (lambda (mol)
          (chem:map-atoms
           nil
           (lambda (atom)
             (setf (gethash atom atoms-to-molecules) mol))
           mol))
        aggregate)
       (chem:map-bonds
        nil
        (lambda (a1 a2 order bond)
          (declare (ignore bond))
          (let ((m1 (gethash a1 atoms-to-molecules))
                (m2 (gethash a2 atoms-to-molecules)))
            (unless (eq m1 m2)
              (merge-molecules aggregate m1 m2)
              (go top))))
        aggregate)
       (return-from maybe-join-molecules-in-aggregate nil))))

(defun maybe-split-molecules-in-aggregate (aggregate)
  "* Arguments
- aggregate :: An aggregate.
* Description
Break up the molecules in the aggregate into a list of molecules using spanning trees and return that list.
"
  (let ((atoms-to-residues (make-hash-table))
        (residues-to-molecules (make-hash-table))
        (atoms-seen (make-hash-table))
        molecules)
    (chem:map-molecules
     nil
     (lambda (mol)
       (let (atoms)
         (chem:map-residues
          nil
          (lambda (residue)
            (chem:map-atoms
             nil
             (lambda (atom)
               (setf (gethash atom atoms-to-residues) residue)
               (push atom atoms))
             residue))
          mol)
         (loop for atom in atoms
               for iteration from 0
               when (not (gethash atom atoms-seen))
                 do (multiple-value-bind (new-molecule atoms-spanned)
                        (span-across-molecule-from-atom mol
                                                        atom
                                                        atoms-seen
                                                        residues-to-molecules
                                                        atoms-to-residues)
                      (if (= atoms-spanned (length atoms))
                          (push mol molecules) ; the new-molecule has all the atoms of the old molecule - so keep the oldq
                          (push new-molecule molecules))))))
     aggregate)
;;; Remove the molecules in the aggregate
    (let ((mols (chem:map-molecules 'list #'identity aggregate)))
      (loop for mol in mols
            do (chem:remove-molecule aggregate mol)))
    ;; Add the new molecules
    (loop for mol in (nreverse molecules)
          do (chem:add-matter aggregate mol))
    aggregate))


(defun span-across-molecule-from-atom (old-molecule
                                       atom
                                       atoms-seen
                                       residues-to-molecules
                                       atoms-to-residues)
  (format *debug-io* "span-across-molecule-from-atom old-molecule ~a atom ~a~%" old-molecule atom)
  (let ((new-molecule (chem:make-molecule (chem:get-name old-molecule)))
        (spanning-tree (chem:spanning-loop/make atom))
        (atoms-spanned 0)
        residues)
    (loop
      (if (chem:advance-loop-and-process spanning-tree)
          (progn
            (format *debug-io* "Advancing loop with atom: ~a~%" (chem:get-atom spanning-tree))
            (let* ((atom (chem:get-atom spanning-tree))
                   (residue (gethash atom atoms-to-residues))
                   (_ (unless residue
                        (error "Could not find the residue for atom ~s" atom)))
                   (residue-molecule (gethash residue residues-to-molecules)))
              (incf atoms-spanned)
              (format *debug-io* "Atom ~a is in residue ~a in molecule ~a~%" atom residue residue-molecule)
              (setf (gethash atom atoms-seen) atom)
              (unless residue-molecule
                (push residue residues)
                (setf (gethash residue residues-to-molecules) new-molecule))))
          (progn
            (return nil))))
    (let ((residue-list (sort residues #'< :key #'chem:get-id)))
      (loop for res in residue-list
            do (chem:add-matter new-molecule res)))
    (when (= (chem:content-size new-molecule) 1)
      (chem:set-name new-molecule (chem:get-name (chem:content-at new-molecule 0))))
    (when (= (chem:content-size new-molecule) 0)
      (error "I am about to return an empty molecule!!!!"))
    (values new-molecule atoms-spanned)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Remove molecules of a particular type
;;;

(defun remove-molecules (aggregate molecule-type)
  (let (molecules)
    (chem:do-molecules (molecule aggregate)
      (when (eq (chem:molecule-type molecule) molecule-type)
        (push molecule molecules)))
    (loop for molecule in molecules
          do (chem:remove-molecule aggregate molecule))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; combine molecules into a new aggregate
;;;
(defun combine (&rest matters)
  "Combine the contents of the aggregates into a new aggregate and return that"
  (let ((new-agg (chem:make-aggregate)))
    (loop for matter in matters
          do (cond
               ((typep matter 'chem:molecule)
                (chem:add-matter new-agg matter))
               ((typep matter 'chem:aggregate)
                (chem:do-molecules (mol matter)
                  (chem:add-matter new-agg mol)))
               (t (error "You cannot combine a ~s" matter))))
    new-agg))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Rename atoms in a matter using a smirks string with maps
;;;

(defun rename-atoms-using-smirks (matter smirks-root name-map)
  "Search through the matter for an atom that matches the smirks string and
   then rename the atoms using the smirks labels.
   The atom names become <element><smirks#> and attached hydrogens 
   are named using a scheme that I'll have to come up with."
  (chem:do-atoms (a matter)
    (multiple-value-bind (matched match)
        (smarts:match smirks-root a)
      (format t "Do something with the match ~s~%" match))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Manipulate solvent in matter
;;;

(defun assign-solvent-molecules-using-residue-name (matter solvent-residue-name)
  "Search through the matter for molecules that contain a single residue that has the 
solvent-residue-name. Set the molecule-type of those molecules to :solvent"
  (chem:do-molecules (mol matter)
    (when (and (= (chem:content-size mol) 1)
               (eq (chem:get-name (chem:content-at mol 0)) solvent-residue-name))
      (chem:setf-molecule-type mol :solvent))))

(defun identify-solvent-residue-name (matter)
  "Figure out the residue name of solvent by looking at the name of the single
residue in any molecules that have the molecule-type :solvent."
  (chem:do-molecules (mol matter)
    (when (eq (chem:molecule-type mol) :solvent)
      (if (= (chem:content-size mol) 1)
          (chem:get-name (chem:content-at mol 0))
          (error "A molecule with molecule-type :solvent was found that has ~a residues - it needs to have a single residue" (chem:content-size mol))))))

(defun ensure-unique-hydrogen-names (matter)
  (let ((hydrogen-counter 0))
    (chem:do-atoms (atom1 matter)
      (when (eq (chem:get-element atom1) :H)
        (let ((hydrogen-name (intern (format nil "H~a" (incf hydrogen-counter)) :keyword)))
          (chem:set-name atom1 hydrogen-name))))))
