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
  (let ((neighbors (chem:bonded-atoms-as-list atom)))
    (dolist (neighbor neighbors)
      (when (not (eq neighbor excluded-neighbor))
        (return-from find-neighbor-excluding neighbor)))))

(defun build-hydrogen (hydrogen a0 a1 a2 distance angle-deg dihedral-deg)
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
(defun build-unbuilt-hydrogens-on (central)
  (let ((neighbors (chem:bonded-atoms-as-list central))
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
                   do (build-hydrogen h central angle-atom dihedral-atom 1.0 109.5 dihedral))))
           (2
            (let ((sorted-built (sort-atoms-by-name built-atoms)))
              (loop for improper in '(120.0 -120.0)
                 for h in unbuilt-hydrogens
                 do (build-hydrogen h central (first sorted-built) (second sorted-built) 1.0 109.5 improper))))
           (3
            (let* ((sorted-built (sort-atoms-by-name built-atoms)))
              (destructuring-bind (n0 n1 n2)
                  sorted-built
                (format t "Building from sorted-built: ~a~%" sorted-built)
                (let* ((used-improper (geom:calculate-dihedral (chem:get-position n0)
                                                               (chem:get-position central)
                                                               (chem:get-position n1)
                                                               (chem:get-position n2)))
                       (other-improper (- (* used-improper *rad-to-deg*))))
                  (build-hydrogen (car unbuilt-hydrogens) central n1 n2 1.0 109.5 other-improper)))))
           (otherwise (error "Handle building hydrogens on sp3 centers with ~a built neighbors - unbuilt hydrogens: ~a" num-built unbuilt-hydrogens))))
        (:sp2
         (case num-built
           (1 (let* ((angle-atom (car built-atoms))
                     (dihedral-atom (find-neighbor-excluding angle-atom central)))
                (loop for dihedral in '(180.0 0.0)
                   for h in unbuilt-hydrogens
                   do (build-hydrogen h central angle-atom dihedral-atom 1.0 120.0 dihedral))))
           (2 (let ((sorted-built (sort-atoms-by-name built-atoms)))
                (build-hydrogen (first unbuilt-hydrogens)
                                central
                                (first sorted-built)
                                (second sorted-built)
                                1.0 120.0 180.0)))
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
  (let ((number-of-unbuilt-hydrogens
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
    (let ((bar (make-progress-bar :total number-of-unbuilt-hydrogens
                                  :message "Hydrogens built"))
          (built-hydrogens 0))
      (chem:map-atoms
       nil
       (lambda (atom)
         (when (and (eq (chem:get-element atom) :H)
                    (chem:needs-build atom))
           (let ((built (build-unbuilt-hydrogens-on (chem:bonded-neighbor atom 0))))
             (incf built-hydrogens built))
           (progress-advance bar built-hydrogens)))
       matter)
      (progress-done bar)
      built-hydrogens)))
