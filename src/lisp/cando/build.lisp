(in-package :cando)

;;;
;;;  Hydrogen builder tools
;;;
;;; Code to build the coordinates of hydrogens based on the atoms they are connected to
;;;


(defun find-neighbor-excluding (atom excluded-neighbor)
  "* Arguments
- atom :: A chem:atom
- excluded-neighbor :: A chem:atom
* Description 
Return a neighbor of atom that is not excluded-neighbor."
  (let ((neighbors (chem:atom/bonded-atoms-as-list atom)))
    (dolist (neighbor neighbors)
      (when (not (eq neighbor excluded-neighbor))
        (when (geom:is-defined (chem:get-position neighbor))
          (return-from find-neighbor-excluding neighbor)))))
  (error "Could not find an neighbor with a defined position connected to ~a that is not ~a" atom excluded-neighbor))

(defun build-hydrogen (hydrogen a0 a1 a2 distance angle-deg dihedral-deg res)
  "* Arguments
- hydrogen :: A chem::atom
- a0 :: A chem::atom
- a1 :: A chem::atom
- a2 :: A chem::atom
- distance :: Number
- angle-deg :: Number
- dihedral-deg :: Number
* Description
Build a coordinate for atom hydrogen using the other angles and a
distance/angle-deg/dihedral-deg.  The angle-deg and dihedral-deg are
provided in degrees."
  (let ((angle-rad (* 0.0174533 angle-deg))
        (dihedral-rad (* 0.0174533 dihedral-deg)))
    (unless (geom:is-defined (chem:get-position a0))
      (error "Cannot build hydrogen ~a having ~a with undefined position in residue ~s" hydrogen a0 res))
    (unless (geom:is-defined (chem:get-position a1))
      (error "Cannot build hydrogen ~a having ~a with undefined position in residue ~s" hydrogen a1 res))
    (unless (geom:is-defined (chem:get-position a2))
      (error "Cannot build hydrogen ~a having ~a with undefined position in residue ~s" hydrogen a2 res))
    (let ((pos (geom:build-using-bond-angle-dihedral distance (chem:get-position a0)
                                                     angle-rad (chem:get-position a1)
                                                     dihedral-rad (chem:get-position a2))))
      (chem:set-position hydrogen pos)
      (chem:setf-needs-build hydrogen nil))))

(defun sort-atoms-by-name (atoms)
  "* Arguments
- atoms :: List
* Description
Sort the atoms by name in increasing alphabetical order."
  (sort atoms
        (lambda (x y)
          (string< (string (chem:get-name x))
                   (string (chem:get-name y))))))

(defparameter *rad-to-deg* (/ 1.0 0.0174533))

(defparameter *one-over-tet-vec-length* (/ 1.0 (geom:vlength (geom:vec 1.0 1.0 1.0))))
(defparameter *corners-of-tetrahedron*
  (list (geom:v* (geom:vec  1.0  1.0  1.0) *one-over-tet-vec-length*)
        (geom:v* (geom:vec -1.0 -1.0  1.0) *one-over-tet-vec-length*)
        (geom:v* (geom:vec -1.0  1.0 -1.0) *one-over-tet-vec-length*)
        (geom:v* (geom:vec  1.0 -1.0 -1.0) *one-over-tet-vec-length*)))
        
(defun build-unbuilt-hydrogens-on (central res)
  (let ((neighbors (chem:atom/bonded-atoms-as-list central))
        unbuilt-hydrogens built-atoms)
    (dolist (neighbor neighbors)
      (if (and (eq (chem:get-element neighbor) :H)
               (chem:needs-build neighbor))
          (push neighbor unbuilt-hydrogens)
          (push neighbor built-atoms)))
    ;; Sort the unbuilt hydrogens in alphabetical order
    (setf unbuilt-hydrogens (sort-atoms-by-name unbuilt-hydrogens))
    (let ((num-hydrogens (length unbuilt-hydrogens))
          (num-built (length built-atoms))
          (hybridization (calculate-hybridization central)))
      (case hybridization
        ;; Build sp3 -CH3, -NH3, -NH2, -OH groups
        (:sp3
         (case num-built
           (1 (let* ((angle-atom (car built-atoms))
                     (dihedral-atom (find-neighbor-excluding angle-atom central)))
                (loop for dihedral in '(180.0 -60.0 60.0)
                   for h in unbuilt-hydrogens
                   do (build-hydrogen h central angle-atom dihedral-atom 1.0 109.5 dihedral res))))
           (2
            (let ((sorted-built (sort-atoms-by-name built-atoms)))
              (loop for improper in '(120.0 -120.0)
                 for h in unbuilt-hydrogens
                 do (build-hydrogen h central (first sorted-built) (second sorted-built) 1.0 109.5 improper res))))
           (3
            (let* ((sorted-built (sort-atoms-by-name built-atoms)))
              (destructuring-bind (n0 n1 n2)
                  sorted-built
                (let* ((used-improper (geom:calculate-dihedral (chem:get-position n0)
                                                               (chem:get-position central)
                                                               (chem:get-position n1)
                                                               (chem:get-position n2)))
                       (other-improper (- (* used-improper *rad-to-deg*))))
                  (build-hydrogen (car unbuilt-hydrogens) central n1 n2 1.0 109.5 other-improper res)))))
           (0
            (let ((central-pos (chem:get-position central)))
              (loop for h in unbuilt-hydrogens
                    for index from 0
                    do (chem:set-position h (geom:v+ central-pos (elt *corners-of-tetrahedron* index))))))
           (otherwise (error "Handle building hydrogens on sp3 centers with ~a built neighbors - central atom ~a - unbuilt hydrogens: ~a" num-built central unbuilt-hydrogens))))
        (:sp2
         (case num-built
           (1 (let* ((angle-atom (car built-atoms))
                     (dihedral-atom (find-neighbor-excluding angle-atom central)))
                (loop for dihedral in '(180.0 0.0)
                   for h in unbuilt-hydrogens
                   do (build-hydrogen h central angle-atom dihedral-atom 1.0 120.0 dihedral res))))
           (2 (let ((sorted-built (sort-atoms-by-name built-atoms)))
                (build-hydrogen (first unbuilt-hydrogens)
                                central
                                (first sorted-built)
                                (second sorted-built)
                                1.0 120.0 180.0
                                res)))
           (otherwise (error "Handle building hydrogens on sp2 centers with ~a built neighbors" num-built))))
        (otherwise (warn "Handle building hydrogens on atom with hybridization ~a" hybridization)))
      num-hydrogens)))

(defun build-unbuilt-hydrogens (matter &key progress)
  "* Arguments
- matter :: chem:matter
- progress :: boolean
* Description
Search the matter for hydrogens that need to be build and build good geometry 
for them.  Show progress if progress is T."
  (let* ((atoms-to-residues
           (cond
             ((typep matter 'chem:aggregate)
              (let ((ator (make-hash-table)))
                (chem:do-molecules (mol matter)
                  (chem:do-residues (res mol)
                    (chem:do-atoms (atm res)
                      (setf (gethash atm ator) res))))
                ator))
             ((typep matter 'chem:residue)
              (let ((ator (make-hash-table)))
                (chem:do-atoms (atm matter)
                  (setf (gethash atm ator) matter))
                ator))
             (t (error "Handle matter ~s" matter))
            ))
         (number-of-unbuilt-hydrogens
           (let ((num 0))
             (chem:map-atoms
              nil
              (lambda (atom)
                (when (and (eq (chem:get-element atom) :H)
                           (chem:needs-build atom))
                  (incf num)))
              matter)
             num)))
    ;; Now build the hydrogens
    (let ((bar (when progress
                 (make-progress-bar :total number-of-unbuilt-hydrogens
                                    :message "Hydrogens built")))
          (built-hydrogens 0))
      (chem:map-atoms
       nil
       (lambda (atom)
         (when (and (eq (chem:get-element atom) :H)
                    (chem:needs-build atom))
           (let* ((res (gethash atom atoms-to-residues))
                  (built (build-unbuilt-hydrogens-on (chem:bonded-neighbor atom 0) res)))
             (incf built-hydrogens built))
           (when progress (progress-advance bar built-hydrogens))))
       matter)
      (when progress (progress-done bar))
      built-hydrogens)))
