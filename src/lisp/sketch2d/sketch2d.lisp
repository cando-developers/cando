(in-package :sketch2d)

#| 
Implement the algorithm described in:
Simulation-Based Algorithm for two-Dimensional Chemical Structure Diagram Generation of
Complex Molecules and Ligand-Protein Interactions
Tomasz Fraczek, 
J. Chem. Inf. Model. 2016, 56, 2320-2335
DOI; 10.1021/acs.jcim.6b00391

I'll try a few changes and improvements.
I'll use SMARTS to identify features.
I'll use an angle term instead of a bond term.
|#


(defclass sketch2d ()
  ((molecule :initarg :molecule :accessor molecule)
   (rings :initarg :rings :accessor rings)
   (aromaticity-info :initarg :aromaticity-info :accessor aromaticity-info)
   (dynamics :initarg :dynamics :accessor dynamics)
   (sketch-atoms-to-original :initarg :sketch-atoms-to-original :accessor sketch-atoms-to-original)
   (original-molecule :initarg :original-molecule :accessor original-molecule)))


(defconstant +oozp-scale+ 1.0)
(defconstant +150-degrees-rad+ (* 150.0 0.0174533))
(defparameter *stage1-flatten-force-components* (list 0.85 0.85 0.85))
(defparameter *stage3-flatten-force-components* (list 0.85 0.85 0.25))
(defconstant +stage1-nonbond-constant+ 0.04)
(defconstant +stage2-nonbond-constant+ 1.0)
(defconstant +first-bond-force+ 0.1)
(defconstant +stage1-bond-length+ 5.0)
(defconstant +stage1-scale-sketch-nonbond+ 0.04)
(defconstant +stage1-sketch-nonbond-force+ 1.0)

(defparameter *stage3-flatten-force* 1.0)
(defconstant +stage4-bond-length+ 1.5)
(defconstant +stage4-bond-force+ (* 2.0 +first-bond-force+))
(defconstant +angle-120+ (* 0.0174533 120.0))
(defconstant +angle-180+ (* 0.0174533 180.0))

(defparameter *methyl-smarts* (chem:compile-smarts "[#6H3:0]"))
(defparameter *methyl* (chem:make-chem-info-graph *methyl-smarts*))
(defparameter *hydroxyl* (chem:make-chem-info-graph (chem:compile-smarts "[#8H1:0]")))
(defparameter *thiol* (chem:make-chem-info-graph (chem:compile-smarts "[#16H1:0]")))
(defparameter *primary-ammonium* (chem:make-chem-info-graph (chem:compile-smarts "[#7H3:0]")))
(defparameter *primary-amine* (chem:make-chem-info-graph (chem:compile-smarts "[#7H2:0]")))
(defparameter *secondary-carbon* (chem:make-chem-info-graph (chem:compile-smarts "[#6X4H2:0]")))
(defparameter *secondary-carbon-missing-h* (chem:make-chem-info-graph (chem:compile-smarts "[#6X2H0:0](-*)-*")))
(defparameter *tertiary-carbon* (chem:make-chem-info-graph (chem:compile-smarts "[#6X4H1:0]")))
(defparameter *tertiary-amine-with-h* (chem:make-chem-info-graph (chem:compile-smarts "[#7X4H1:0]")))
(defparameter *lp-oxygen* (chem:make-chem-info-graph (chem:compile-smarts "[#8X2:0]")))
(defparameter *lp-nitrogen* (chem:make-chem-info-graph (chem:compile-smarts "[#7X2:0]")))
(defparameter *lp-sulfur* (chem:make-chem-info-graph (chem:compile-smarts "[#16X2:0]")))
(defparameter *bad-angle* (chem:make-chem-info-graph (chem:compile-smarts "[*X3:0](~[*X1:1])(~[*:2])~[*:3]")))

(defclass sketch-nonbond-force-field () ())

(defmethod chem:find-atom-type-position ((force-field sketch-nonbond-force-field) type) 0)


(defparameter *quaternary-center* (chem:make-chem-info-graph (chem:compile-smarts "[*X4:0](~[*:1])(~[*:2])(~[*:3])~[*:4]")))
(defparameter *three-bond-center* (chem:make-chem-info-graph (chem:compile-smarts "[*X3:0](~[*:1])(~[*:2])~[*:3]")))
(defparameter *linear-bond-center* (chem:make-chem-info-graph (chem:compile-smarts "[*X2:0](#[*:1])-[*:2]")))
(defparameter *cumulated-double-bond-center* (chem:make-chem-info-graph (chem:compile-smarts "[*X2:0](=[*:1])=[*:2]")))

(defmethod prepare-stage1-sketch-function (sketch-function &key (bond-length +stage1-bond-length+))
  "Generate bond and angle energy terms for a 2D sketch."
  (let* ((molecule (chem:get-graph sketch-function))
         (atom-table (chem:node-table sketch-function))
         (nonbond-energy (chem:get-sketch-nonbond-component sketch-function))
         (bond-energy (chem:get-stretch-component sketch-function)))
    (chem:set-scale-sketch-nonbond nonbond-energy +stage1-scale-sketch-nonbond+)
    (loop for ia1 from 0 below (1- (chem:get-number-of-atoms atom-table))
          for atom1 = (chem:elt-atom atom-table ia1)
          for atom1-coordinate-index = (chem:get-coordinate-index atom-table atom1)
          for freeze1 = (or (= (chem:get-atomic-number atom1) 1)
                            (eq (chem:get-element atom1) :lp))
          do (loop for ia2 from (1+ ia1) below (chem:get-number-of-atoms atom-table)
                   for atom2 = (chem:elt-atom atom-table ia2)
                   for atom2-coordinate-index = (chem:get-coordinate-index atom-table atom2)
                   for freeze2 = (or (= (chem:get-atomic-number atom2) 1)
                                     (eq (chem:get-element atom2) :lp))
                   for freeze-flag = (if (or freeze1 freeze2) 1 0)
                   do (chem:add-sketch-nonbond-term nonbond-energy
                                                    atom1-coordinate-index
                                                    atom2-coordinate-index
                                                    freeze-flag
                                                    +stage1-sketch-nonbond-force+)))
    (chem:map-bonds 'nil
                    (lambda (a1 a2 bond-order)
                      (let ((atom1-coord-index (chem:get-coordinate-index atom-table a1))
                            (atom2-coord-index (chem:get-coordinate-index atom-table a2)))
                        (chem:add-sketch-stretch-term bond-energy
                                                      atom1-coord-index
                                                      atom2-coord-index
                                                      +first-bond-force+ bond-length)))
                    molecule)
    sketch-function))

(defun ordered-key (ia1 ia2)
  (if (< ia1 ia2)
      (cons ia1 ia2)
      (cons ia2 ia1)))

(defun add-angle-terms (sketch-function bond-force bond-length)
  (let* ((molecule (chem:get-graph sketch-function))
         (angle120-bond-length (* (sqrt 3.0) bond-length))
         (angle-quat-bond-length (* (sqrt 2.0) bond-length))
         (atom-table (chem:node-table sketch-function))
         (bond-energy (chem:get-stretch-component sketch-function))
         (sketch-nonbond-component (chem:get-sketch-nonbond-component sketch-function))
         (mol-graph (chem:make-molecule-graph-from-molecule molecule))
         (linear-bond-centers (chem:boost-graph-vf2 *linear-bond-center* mol-graph))
         (cumulated-double-bond-centers (chem:boost-graph-vf2 *cumulated-double-bond-center* mol-graph))
         (three-bond-centers (chem:boost-graph-vf2 *three-bond-center* mol-graph))
         (unique-three-bond-centers (make-hash-table :test #'equal))
         (quat-centers (chem:boost-graph-vf2 *quaternary-center* mol-graph))
         (unique-quat-centers (make-hash-table :test #'equal))
         double-nonbond-pairs)
    (progn
      (loop for tbc in three-bond-centers
            for center = (elt tbc 0)
            for atom1 = (elt tbc 1)
            for atom2 = (elt tbc 2)
            for atom3 = (elt tbc 3)
            for center-name = (chem:get-name center)
            for atom1-name = (chem:get-name atom1)
            for atom2-name = (chem:get-name atom2)
            for atom3-name = (chem:get-name atom3)
            for sorted-names = (sort (list (string atom1-name) (string atom2-name) (string atom3-name)) #'string<)
            do (setf (gethash (cons center-name sorted-names) unique-three-bond-centers) (list center atom1 atom2 atom3)))
      (loop for tbc being the hash-values of unique-three-bond-centers
            for center = (elt tbc 0)
            for atom1 = (elt tbc 1)
            for atom2 = (elt tbc 2)
            for atom3 = (elt tbc 3)
            for atom1-coord-index = (chem:get-coordinate-index atom-table atom1)
            for atom2-coord-index = (chem:get-coordinate-index atom-table atom2)
            for atom3-coord-index = (chem:get-coordinate-index atom-table atom3)
            do (chem:add-sketch-stretch-term bond-energy atom1-coord-index atom2-coord-index bond-force angle120-bond-length)
               (chem:add-sketch-stretch-term bond-energy atom2-coord-index atom3-coord-index bond-force angle120-bond-length)
               (chem:add-sketch-stretch-term bond-energy atom3-coord-index atom1-coord-index bond-force angle120-bond-length)))
    ;; Gather up all the unique quat-centers
    (progn
      (loop for tbc in quat-centers
            for center = (elt tbc 0)
            for atom1 = (elt tbc 1)
            for atom2 = (elt tbc 2)
            for atom3 = (elt tbc 3)
            for atom4 = (elt tbc 4)
            for center-name = (chem:get-name center)
            for atom1-name = (chem:get-name atom1)
            for atom2-name = (chem:get-name atom2)
            for atom3-name = (chem:get-name atom3)
            for atom4-name = (chem:get-name atom4)
            for sorted-names = (sort (list (string atom1-name) (string atom2-name) (string atom3-name) (string atom4-name)) #'string<)
            do (setf (gethash (cons center-name sorted-names) unique-quat-centers) (list center atom1 atom2 atom3 atom4)))
      ;; Push the atom pairs that are part of angles centered on quat centers
      (loop for fbc being the hash-values of unique-quat-centers
            for center = (elt fbc 0)
            for atom1 = (elt fbc 1)
            for atom2 = (elt fbc 2)
            for atom3 = (elt fbc 3)
            for atom4 = (elt fbc 4)
            for atom1-index = (chem:get-coordinate-index atom-table atom1)
            for atom2-index = (chem:get-coordinate-index atom-table atom2)
            for atom3-index = (chem:get-coordinate-index atom-table atom3)
            for atom4-index = (chem:get-coordinate-index atom-table atom4)
            do (push (ordered-key atom1-index atom2-index) double-nonbond-pairs)
               (push (ordered-key atom1-index atom3-index) double-nonbond-pairs)
               (push (ordered-key atom1-index atom4-index) double-nonbond-pairs)
               (push (ordered-key atom2-index atom3-index) double-nonbond-pairs)
               (push (ordered-key atom2-index atom4-index) double-nonbond-pairs)
               (push (ordered-key atom3-index atom4-index) double-nonbond-pairs)))
    ;; Double the nonbond constant for quat center angle outer atom pairs
    (chem:walk-sketch-nonbond-terms
     sketch-nonbond-component
     (lambda (index freeze-flags ia1 ia2 constant)
       (declare (ignore freeze-flags))
       (let ((key (ordered-key ia1 ia2)))
         (when (member key double-nonbond-pairs :test #'equal)
           (chem:modify-sketch-nonbond-term-constant sketch-nonbond-component index (* 2.0 constant))))))
    #+(or)(loop for lbc in linear-bond-centers
                for center = (elt lbc 0)
                for atom1 = (elt lbc 1)
                for atom2 = (elt lbc 2)
                do (chem:add-angle-term angle-energy atom-table atom1 center atom2 +initial-linear-angle-force+ +angle-180+))
    #+(or)(loop for lbc in cumulated-double-bond-centers
          for center = (elt lbc 0)
          for atom1 = (elt lbc 1)
          for atom2 = (elt lbc 2)
          do (chem:add-angle-term angle-energy atom-table atom1 center atom2 +initial-linear-angle-force+ +angle-180+)))
  sketch-function)

(defun add-flatten-function (sketch-function molecule &optional (scale +oozp-scale+))
  (let ((atom-table (chem:node-table sketch-function))
        (oozp-component (chem:get-out-of-zplane-component sketch-function)))
    (cando:do-atoms (atom molecule)
      (let ((atom-coord-index (chem:get-coordinate-index atom-table atom)))
      (chem:add-out-of-zplane-term oozp-component atom-coord-index scale 0.0))
      )))

(defun add-atom-annotations (molecule)
  (cando:do-atoms (atom molecule)
    (let ((hydrogens 0)
          (charge 0))
    (cond
      ((/= (chem:get-atomic-number atom) 1)
       (let ((bonds (chem:bonds-as-list atom)))
         (loop for bond in bonds
               for other-atom = (chem:get-other-atom bond atom)
               when (= (chem:get-atomic-number other-atom) 1)
                 do (incf hydrogens))
         (chem:set-property atom :hydrogens hydrogens)))))))

(defun sketch-preparation (molecule)
  (labels ((find-groups (chem-info-graph molecule-graph)
             (loop for match in (chem:boost-graph-vf2 chem-info-graph molecule-graph)
                   collect (elt match 0)))
           (remove-connected-hydrogens (mol atoms &optional number-to-remove)
             (cando:do-residues (res mol)
               (loop for connected-hydrogens = nil
                     for atom in atoms
                     for bonds = (chem:bonds-as-list atom)
                     do (loop for bond in bonds
                              for a1 = (chem:get-atom1 bond)
                              for a2 = (chem:get-atom2 bond)
                              when (eq (chem:get-element a1) :H)
                                do (push a1 connected-hydrogens)
                              when (eq (chem:get-element a2) :H)
                                do (push a2 connected-hydrogens))
                        (loop for hatom in (subseq connected-hydrogens 0 number-to-remove)
                              do (chem:remove-atom-delete-bonds res hatom)))))
           (add-lone-pair (mol atoms idx)
             (cando:do-residues (res mol)
               (loop for atom in atoms
                     for lp-idx from idx
                     for lp-name = (intern (format nil "lp~a" lp-idx) :keyword)
                     for lp = (chem:make-atom lp-name :lp)
                     do (incf idx)
                        (chem:add-matter res lp)
                        (chem:bond-to atom lp :single-bond)
                     ))
             idx)
           (unique-secondary-carbon-missing-h (mol-graph)
             (let ((missing-h (find-groups *secondary-carbon-missing-h* mol-graph))
                   (uniques (make-hash-table)))
               (loop for one in missing-h
                     do (setf (gethash one uniques) one))
               (let ((unique-list nil))
                 (maphash (lambda (key value)
                            (declare (ignore key))
                            (push value unique-list))
                          uniques)
                 unique-list))))
    (let* ((sketch-atoms-to-original (make-hash-table))
           (mol-copy (chem:matter-copy molecule sketch-atoms-to-original))
           (rings (chem:identify-rings mol-copy))
           (chem:*current-rings* rings)
           (aromaticity-info (chem:identify-aromatic-rings mol-copy :mdl)))
      (add-atom-annotations mol-copy)
      (let ((mol-graph (chem:make-molecule-graph-from-molecule mol-copy)))
        (remove-connected-hydrogens mol-copy (find-groups *methyl* mol-graph))
        (remove-connected-hydrogens mol-copy (find-groups *hydroxyl* mol-graph))
        (remove-connected-hydrogens mol-copy (find-groups *thiol* mol-graph))
        (remove-connected-hydrogens mol-copy (find-groups *primary-ammonium* mol-graph))
        (remove-connected-hydrogens mol-copy (find-groups *primary-amine* mol-graph))
        (remove-connected-hydrogens mol-copy (find-groups *secondary-carbon* mol-graph) 1)
        (remove-connected-hydrogens mol-copy (find-groups *tertiary-carbon* mol-graph) 1)
        (remove-connected-hydrogens mol-copy (find-groups *tertiary-amine-with-h* mol-graph) 1)
        (let ((idx 0))
          (let ((missing-h (unique-secondary-carbon-missing-h mol-graph)))
            (format t "Missing-h ~a~%" missing-h)
            (setf idx (add-lone-pair mol-copy missing-h idx))
            (setf idx (add-lone-pair mol-copy (find-groups *lp-oxygen* mol-graph) idx))
            (setf idx (add-lone-pair mol-copy (find-groups *lp-nitrogen* mol-graph) idx))
            (setf idx (add-lone-pair mol-copy (find-groups *lp-sulfur* mol-graph) idx))))
        (chem:setf-force-field-name mol-copy :.hidden.sketch2d-1)
;;;        (cando:jostle mol-copy 20.0 t)
        (make-instance 'sketch2d
                       :molecule mol-copy
                       :rings rings
                       :aromaticity-info aromaticity-info
                       :original-molecule molecule
                       :sketch-atoms-to-original sketch-atoms-to-original)))))


(defun randomize-atoms (atom-table &key frozen from-zero (width 40.0) &aux (half-width (/ width 2.0)))
  "Randomly jostle atoms from their current positions"
  (flet ((jostle-atom (atom)
           (when from-zero
             (chem:set-position atom (geom:vec 0.0 0.0 0.0)))
           (let* ((cp (chem:get-position atom))
                  (pos (geom:vec
                        (+ (- (random width) half-width) (geom:vx cp))
                        (+ (- (random width) half-width) (geom:vy cp))
                        (+ (- (random width) half-width) (geom:vz cp)))))
             (chem:set-position atom pos))))
    (loop for index from 0 below (chem:get-number-of-atoms atom-table)
          for atom = (chem:elt-atom atom-table index)
          do (if frozen
                 (when (= (aref frozen index) 0)
                   (jostle-atom atom))
                 (jostle-atom atom)))))

(defun randomize-coordinates (coordinates &key frozen from-zero (width 40.0) &aux (half-width (/ width 2.0)))
  "Randomly jostle atoms from their current positions"
  (flet ((jostle-atom (index)
           (when from-zero
             (setf (elt coordinates (+ index 0)) 0.0
                   (elt coordinates (+ index 1)) 0.0
                   (elt coordinates (+ index 2)) 0.0))
           (let* ((xp (elt coordinates (+ index 0)))
                  (yp (elt coordinates (+ index 1)))
                  (zp (elt coordinates (+ index 2)))
                  (xn (+ (- (random width) half-width) xp))
                  (yn (+ (- (random width) half-width) yp))
                  (zn (+ (- (random width) half-width) zp)))
             (setf (elt coordinates (+ index 0)) xn
                   (elt coordinates (+ index 1)) yn
                   (elt coordinates (+ index 2)) zn))))
    (loop for coord-index below (length coordinates) by 3
          for index from 0
          do (if frozen
                 (when (= (aref frozen index) 0)
                   (jostle-atom coord-index))
                 (jostle-atom coord-index)))))

(defparameter *edited-mol* nil)
(defun setup-simulation (molecule &key accumulate-coordinates)
  (let* ((sketch (sketch-preparation molecule))
         (edited-mol (molecule sketch)))
    (setf *edited-mol* edited-mol)
    (cando:do-atoms (atom edited-mol)
      (chem:set-type atom :sketch))
    (let* ((dummy-sketch-nonbond-ff (make-instance 'sketch-nonbond-force-field))
           (sketch-function (chem:make-sketch-function edited-mol dummy-sketch-nonbond-ff))
           (atom-table (chem:node-table sketch-function))
           (dynamics (dynamics:make-atomic-simulation sketch-function
                                                      :accumulate-coordinates accumulate-coordinates)))
      (setf (dynamics sketch) dynamics)
      sketch)))

(defun advance-simulation (dynamics &key frozen)
  (dynamics:velocity-verlet-step
   dynamics
   :velocity-verlet-function #'chem:sketch-function-velocity-verlet-step
   :frozen frozen))



#|
Use the approach described at:
https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
to check if two line segments (bonds) overlap/intersect
|#
(defun orientation (p1 p2 p3)
  (let ((orient (- (* (- (geom:vy p2) (geom:vy p1))
                      (- (geom:vx p3) (geom:vx p2)))
                   (* (- (geom:vy p3) (geom:vy p2))
                      (- (geom:vx p2) (geom:vx p1))))))
    (cond
      ((= 0.0 orient) 0)
      ((< 0.0 orient) 1)
      (t 2))))

;;; Given three colinear points p, q, r, the function checks if 
;;; point q lies on line segment 'pr' 
(defun on-segment (p q r)
  (and (<= (geom:vx q) (max (geom:vx p) (geom:vx r)))
       (>= (geom:vx q) (max (geom:vx p) (geom:vx r)))
       (<= (geom:vy q) (max (geom:vy p) (geom:vy r)))
       (>= (geom:vy q) (max (geom:vy p) (geom:vy r)))))

(defun bonds-intersect (p1 q1 p2 q2)
  (let ((o1 (orientation p1 q1 p2))
        (o2 (orientation p1 q1 q2))
        (o3 (orientation p2 q2 p1))
        (o4 (orientation p2 q2 q1)))
    (let ((intersect (cond
                       ((and (not (= o1 o2)) (not (= o3 o4))) t)
                       ((and (= o1 0) (on-segment p1 p2 q1)) t)
                       ((and (= o2 0) (on-segment p1 q2 q1)) t)
                       ((and (= o3 0) (on-segment p2 p1 q2)) t)
                       ((and (= o4 0) (on-segment p2 q1 q2)) t)
                       (t nil))))
      intersect)))


(defun intersection-point (p1 q1 p2 q2)
  (let ((x1 (geom:vx p1))
        (y1 (geom:vy p1))
        (x2 (geom:vx q1))
        (y2 (geom:vy q1))
        (x3 (geom:vx p2))
        (y3 (geom:vy p2))
        (x4 (geom:vx q2))
        (y4 (geom:vy q2)))
    (let* ((tat (infix:infix ((y3 - y4)*(x1 - x3)+(x4 - x3)*(y1 - y3))))
           (tab (infix:infix ((x4 - x3)*(y1 - y2)-(x1 - x2)*(y4 - y3))))
           (ta (infix:infix tat / tab)))
      (geom:v+ (geom:v* (geom:v- q1 p1) ta) p1))))

(defun test-intersect ()
  (let ((p1 (geom:vec 1 1 0))
        (q1 (geom:vec 10 1 0))
        (p2 (geom:vec 1 2 0))
        (q2 (geom:vec 10 2 0)))
    (if (not (bonds-intersect p1 q1 p2 q2))
        (format t "Passed test1~%")
        (error "failed test   (not (bonds-intersect ~a ~a ~a ~a))")))
  (let ((p1 (geom:vec 10 0 0))
        (q1 (geom:vec 0 10 0 )) 
        (p2 (geom:vec 0 0 0 ))
        (q2 (geom:vec 10 10 0)))
    (if (bonds-intersect p1 q1 p2 q2)
        (format t "Passed test2~%")
        (error "failed test   (bonds-intersect ~a ~a ~a ~a)~%" p1 q1 p2 q2))
    (format t "Intersection-point ~a~%" (intersection-point p1 q1 p2 q2)))
  (let ((p1 (geom:vec -5 -5 0))
        (q1 (geom:vec 0 0 0 )) 
        (p2 (geom:vec 1 1 0 ))
        (q2 (geom:vec 10 10 0))) 
    (if (not (bonds-intersect p1 q1 p2 q2))
        (format t "Passed test3~%")
        (error "failed test   (not (bonds-intersect ~a ~a ~a ~a))"))))


(defstruct stretch-term atom1 atom2)

(defun detect-intersecting-bonds (matter)
  (let ((stretches (chem:map-bonds
                    'vector
                    (lambda (a1 a2 order)
                      (make-stretch-term :atom1 a1 :atom2 a2))
                    matter))
        intersecting)
    (loop named outer
          for idx1 from 0 below (1- (length stretches))
          for str1 = (elt stretches idx1)
          for atomp1 = (stretch-term-atom1 str1)
          for atomq1 = (stretch-term-atom2 str1)
          for p1 = (chem:get-position atomp1)
          for q1 = (chem:get-position atomq1)
          do (loop named inner
                   for idx2 from (1+ idx1) below (length stretches)
                   for str2 = (elt stretches idx2)
                   for atomp2 = (stretch-term-atom1 str2)
                   for atomq2 = (stretch-term-atom2 str2)
                   for p2 = (chem:get-position atomp2)
                   for q2 = (chem:get-position atomq2)
                   do (progn
                        (cond
                          ((or (eq atomp1 atomp2)
                               (eq atomp1 atomq2)
                               (eq atomq1 atomp2)
                               (eq atomq1 atomq2)))
                          ((bonds-intersect p1 q1 p2 q2)
                           (push (cons str1 str2) intersecting))))))
    intersecting))

(defstruct distance-atom distance index atom)

(defun generate-problem-area (point cutoff atom-table)
  (let* ((number-of-atoms (chem:get-number-of-atoms atom-table))
         close-atoms)
    (loop for ia from 0 below number-of-atoms
          for atom = (chem:elt-atom atom-table ia)
          for atom-pos = (chem:get-position atom)
          for delta = (geom:v- atom-pos point)
          for distance = (geom:vlength delta)
          for da = (make-distance-atom :distance distance
                                       :index ia
                                       :atom atom)
          do (when (< distance cutoff)
               (push da close-atoms)))
    close-atoms))

(defun identify-problem-areas (dynamics &optional (cutoff (* 2.4 +stage4-bond-length+)))
  (let* ((energy-function (dynamics:scoring-function dynamics))
         (molecule (chem:get-graph energy-function))
         (atom-table (chem:atom-table energy-function))
         (energy-stretch (chem:get-stretch-component energy-function))
         (intersections (detect-intersecting-bonds molecule)))
    (if intersections
      (let (worst-problem-area)
        (loop for intersect in intersections
              for str1 = (car intersect)
              for str2 = (cdr intersect)
              for atomp1 = (stretch-term-atom1 str1)
              for atomq1 = (stretch-term-atom2 str1)
              for atomp2 = (stretch-term-atom1 str2)
              for atomq2 = (stretch-term-atom2 str2)
              for p1 = (chem:get-position atomp1)
              for q1 = (chem:get-position atomq1)
              for p2 = (chem:get-position atomp2)
              for q2 = (chem:get-position atomq2)
              for intersect-point = (intersection-point p1 q1 p2 q2)
              for problem-area = (generate-problem-area intersect-point cutoff atom-table)
              do (cond
                   ((null worst-problem-area)
                    (setf worst-problem-area problem-area))
                   ((> (length problem-area)
                       (length worst-problem-area))
                    (setf worst-problem-area problem-area))))
        (if worst-problem-area
            (let* ((number-of-atoms (chem:get-number-of-atoms atom-table))
                   (frozen (make-array number-of-atoms :element-type 'bit :initial-element 1)))
              (loop for unfrozen in worst-problem-area
                    for index = (distance-atom-index unfrozen)
                    for atom = (distance-atom-atom unfrozen)
                    do (setf (elt frozen index) 0))
              (values (length intersections) frozen))
            (values 0 nil)))
      (values 0 nil))))

(defun randomize-neighbor-near-unfrozen-atom (coordinates neighbor-index*3 index*3 &optional (bond-length +stage4-bond-length+))
  (let ((xpos (elt coordinates index*3))
        (ypos (elt coordinates (+ 1 index*3)))
        (zpos (elt coordinates (+ 2 index*3))))
    (incf xpos (* bond-length (- (random 2.0) 1.0)))
    (incf ypos (* bond-length (- (random 2.0) 1.0)))
    (incf zpos (* bond-length (- (random 2.0) 1.0)))
    (setf (elt coordinates neighbor-index*3) xpos
          (elt coordinates (+ neighbor-index*3 1)) ypos
          (elt coordinates (+ neighbor-index*3 2)) zpos)))
              
(defun unfreeze (atom-table coordinates frozen)
  (let ((new-frozen (copy-seq frozen)))
    (loop for index from 0 below (length frozen)
          for atom = (chem:elt-atom atom-table index)
          do (when (= (elt frozen index) 0)
               (loop named inner
                     for bondi from 0 below (chem:number-of-bonds atom)
                     for neighbor = (chem:bonded-neighbor atom bondi)
                     for neighbor-index*3 = (chem:get-coordinate-index atom-table neighbor)
                     for neighbor-index = (/ neighbor-index*3 3)
                     do (when (= (elt frozen neighbor-index) 1)
                          (progn
                            (setf (elt new-frozen neighbor-index) 0)
                            (randomize-neighbor-near-unfrozen-atom coordinates neighbor-index*3 (* index 3))
                            (return-from inner nil))))))
    new-frozen))
                          
(defun geometric-center-of-unfrozen-atoms (atom-table coordinates frozen)
  (let ((xpos 0.0)
        (ypos 0.0)
        (zpos 0.0)
        (count 0))
    (loop for index from 0 below (length frozen)
          for coord-index = (* 3 index)
          do (when (not (elt frozen index))
               (setf xpos (+ (elt coordinates coord-index) xpos)
                     ypos (+ (elt coordinates (+ coord-index 1)) ypos)
                     zpos (+ (elt coordinates (+ coord-index 2)) zpos))
               (incf count)))
    (geom:vec (/ xpos (float count 1.0d0))
              (/ ypos (float count 1.0d0))
              (/ zpos (float count 1.0d0)))))

(defun transform-unfrozen-atoms-center-to (to-center atom-table coordinates frozen)
  (let* ((from-center (geometric-center-of-unfrozen-atoms atom-table coordinates frozen))
         (offset (geom:v- to-center from-center)))
    (loop with count = 0
          for index from 0 below (length frozen)
          for coord-index = (* 3 index)
          do (when (not (elt frozen index))
               (let ((xpos (+ (geom:vx offset) (elt coordinates coord-index)))
                     (ypos (+ (geom:vy offset) (elt coordinates (+ 1 coord-index))))
                     (zpos (+ (geom:vz offset) (elt coordinates (+ 2 coord-index)))))
                 (setf (elt coordinates coord-index) xpos
                       (elt coordinates (+ coord-index 1)) ypos
                       (elt coordinates (+ coord-index 2)) zpos))))))


(defstruct bad-angle center bisected in-small-ring mono left right)

(defun has-bisected-angle (center mono left right)
  (let ((vcenter (chem:get-position center))
        (vmono (chem:get-position mono))
        (vleft (chem:get-position left))
        (vright (chem:get-position right)))
    (let ((vleft-right (geom:v+ (geom:v- vleft vcenter) (geom:v- vright vcenter)))
          (vdelta-mono (geom:v- vmono vcenter)))
      (> (geom:vdot vleft-right vdelta-mono) 0.0)))) ; mono pointed into V formed by left-center-right

(defun has-wide-angle (center mono left right)
  (let ((vcenter (chem:get-position center))
        (vmono (chem:get-position mono))
        (vleft (chem:get-position left))
        (vright (chem:get-position right)))
    (let ((angle-rad (geom:vangle vleft vcenter vright)))
      (> angle-rad +150-degrees-rad+))))


(defun connected (ba1 ba2)
  (or (eq (bad-angle-center ba1) (bad-angle-left ba2))
      (eq (bad-angle-center ba1) (bad-angle-right ba2))))


(defun add-bad-angle-to-groups (bad-angle bonded-groups)
  (loop for ia1 from 0 below (length bonded-groups)
        for group = (elt bonded-groups ia1)
        when (position bad-angle group :test #'connected)
          do (progn
               (vector-push-extend bad-angle group)
               (return-from add-bad-angle-to-groups nil)))
  (let ((new-group (make-array 16 :fill-pointer 0 :adjustable t)))
    (vector-push-extend bad-angle new-group)
    (vector-push-extend new-group bonded-groups)))

(defun bad-angle-swap (bad-angle)
  (let ((temp-left (bad-angle-left bad-angle)))
    (setf (bad-angle-left bad-angle) (bad-angle-right bad-angle)
          (bad-angle-right bad-angle) temp-left)))

(defun has-two-neighbors (bad-angle bad-angle-vector)
  (let ((left-pos (position (bad-angle-left bad-angle) bad-angle-vector :key #'bad-angle-center))
        (right-pos (position (bad-angle-right bad-angle) bad-angle-vector :key #'bad-angle-center)))
    (cond
      ((and left-pos right-pos)) ; we found neighbors on both sides, keep going
      (left-pos                 ; we found a neighbor on our left side
       (bad-angle-swap bad-angle) ; swap neighbor so that its on right
       nil)                       ; return nil to return
      (right-pos ; we found a neighbor on our right side - and none on our left - return nil to return
       nil)
      (t (error "either we find a neighbor on left or right or both - but not this")))))

(defun find-leftmost-position (bad-angle-vector)
  (if (= (length bad-angle-vector) 1)
      0            ; If there is only one bad angle then it's leftmost
      (progn
        (loop for ia from 0 below (length bad-angle-vector)
              for bad-angle = (elt bad-angle-vector ia)
              for has-two-neighbors = (has-two-neighbors bad-angle bad-angle-vector)
              unless has-two-neighbors ; if pos is nil then we found a leftmost and we return
                do (progn
                     (return-from find-leftmost-position ia)))
        0                       ; default index 0 since we have a ring
        )))

(defun is-neighbor-with-atom (center-atom right-atom bad-angle-or-nil)
  (if bad-angle-or-nil
      (if (eq (bad-angle-center bad-angle-or-nil) right-atom)
        (cond
          ((eq center-atom (bad-angle-left bad-angle-or-nil)))
          ((eq center-atom (bad-angle-right bad-angle-or-nil))
           (bad-angle-swap bad-angle-or-nil)
           t)
          (nil)))
      nil))

(defun find-neighbor-to-right (bad-angle bad-angle-vector)
  (let ((right-atom (bad-angle-right bad-angle))
        (center-atom (bad-angle-center bad-angle)))
    (let ((found (position-if (lambda (ba) (is-neighbor-with-atom center-atom right-atom ba)) bad-angle-vector)))
      found)))

(defun order-group (bad-angle-vector)
  ;; Find the left most bad angle -
  (let ((new-bad-angle-vector (make-array (length bad-angle-vector) :fill-pointer 0 :adjustable t))
        (leftmost (find-leftmost-position bad-angle-vector)))
    (vector-push (elt bad-angle-vector leftmost) new-bad-angle-vector)
    (let ((temp-index nil)
          (left (elt bad-angle-vector leftmost)))
      (setf (elt bad-angle-vector leftmost) nil)
      (loop for ia from 1 below (length bad-angle-vector)
            for next-left = left then next-left
            for right-idx = (find-neighbor-to-right next-left bad-angle-vector)
            when right-idx
              do (let ((right (elt bad-angle-vector right-idx)))
                   (vector-push right new-bad-angle-vector)
                   (setf (elt bad-angle-vector right-idx) nil)
                   (setf next-left right))))
    new-bad-angle-vector))


(defun is-in-small-ring (atom rings)
  (loop for ring in rings
        when (and (< (length ring) 7) (position atom ring))
          do (return-from is-in-small-ring t))
  nil)
             
(defun find-bad-angle-groups (sketch)
  ;; Find the bad angles
  (let* ((molecule (molecule sketch))
         (rings (rings sketch))
         (mol-graph (chem:make-molecule-graph-from-molecule molecule))
         (matches (chem:boost-graph-vf2 *bad-angle* mol-graph))
         (bad-angles-ht (make-hash-table :test #'equalp)))
    (loop for match in matches
          for center = (elt match 0)
          for mono = (elt match 1)
          for left = (elt match 2)
          for right = (elt match 3)
          for in-small-ring = (is-in-small-ring center rings)
          do (let ((bisected (has-bisected-angle center mono left right))
                   (wide-angle (has-wide-angle center mono left right)))
               (when (or bisected wide-angle)
                 (let ((ba23 (make-bad-angle :center (elt match 0)
                                             :in-small-ring in-small-ring
                                             :bisected bisected
                                             :mono (elt match 1)
                                             :left (elt match 2) :right (elt match 3)))
                       (ba32 (make-bad-angle :center (elt match 0)
                                             :in-small-ring in-small-ring
                                             :bisected bisected
                                             :mono (elt match 1)
                                             :left (elt match 3) :right (elt match 2))))
                   (unless (or (gethash ba23 bad-angles-ht)
                               (gethash ba32 bad-angles-ht))
                     (setf (gethash ba23 bad-angles-ht) t))))))
    ;; group together the bad angles that connect to each other
    (let (bad-angles)
      (maphash (lambda (ba dummy)
                 (declare (ignore dummy))
                 (push ba bad-angles))
               bad-angles-ht)
      (let ((bonded-groups (make-array (length bad-angles) :fill-pointer 0)))
        (loop for bad-angle in bad-angles
              do (add-bad-angle-to-groups bad-angle bonded-groups))
        (let ((merged t))
          (loop while merged
                do (setf merged nil)
                do (loop named maybe-merge
                         for ia from 0 below (1- (length bonded-groups))
                         for groupa = (elt bonded-groups ia)
                         do (loop for ib from (1+ ia) below (length bonded-groups)
                                  for groupb = (elt bonded-groups ib)
                                  when (and groupa groupb
                                            (some (lambda (bad-angle)
                                                    (position bad-angle groupa :test #'connected))
                                                  groupb))
                                    do (progn
                                         (setf merged t)
                                         (setf (elt bonded-groups ia) (concatenate 'vector groupa groupb))
                                         (setf (elt bonded-groups ib) nil))))))
        ;; Gather up the groups and order them center -> right
        (let ((groups (loop for ia from 0 below (length bonded-groups)
                            for group = (elt bonded-groups ia)
                            when group
                              collect group)))
          (loop for group in groups
                collect (order-group group)))))))

(defun is-ring (group)
  (when (> (length group) 1)
    (let ((last-bad-angle (elt group (1- (length group))))
          (first-bad-angle (elt group 0)))
      (eq (bad-angle-right last-bad-angle)
          (bad-angle-center first-bad-angle)))))

(defun maybe-bring-heteroatom-to-front (group)
  (let ((hetero-idx (position-if (lambda (a) (not (eq :c (chem:get-element a)))) group)))
    (if (and hetero-idx (> hetero-idx 0))
        (concatenate 'vector (subseq group hetero-idx (length group)) (subseq group 0 (1- hetero-idx)))
        group)))

(defun fix-one-bad-angle (bad-angle force)
  (let* ((center (chem:get-position (bad-angle-center bad-angle)))
         (mono (chem:get-position (bad-angle-mono bad-angle)))
         (left (chem:get-position (bad-angle-left bad-angle)))
         (right (chem:get-position (bad-angle-right bad-angle)))
         (in-small-ring (bad-angle-in-small-ring bad-angle))
         (bisected (bad-angle-bisected bad-angle))
         (delta (geom:v- mono center))
         (between-left-right (geom:v+ (geom:v- left center) (geom:v- right center)))
         (delta-dot-between-left-right (geom:vdot between-left-right delta))
         (would-flip-in (< delta-dot-between-left-right 0.0)))
    (flet ((flip-mono (why)
             (let ((flipped-mono (geom:v+ center (geom:v* delta -1.0))))
               (chem:set-position (bad-angle-mono bad-angle) flipped-mono))))
      (cond
        ((and in-small-ring bisected)
         (flip-mono :in-small-ring-bisected))
        ((and would-flip-in in-small-ring)
         nil)
        (force
         (flip-mono :force))))))

(defun fix-bad-angle-group (group)
  ;; Handle double bonds
  (cond
    ((= (length group) 1)
     (fix-one-bad-angle (elt group 0) t))
    ((= (length group) 2)
     (fix-one-bad-angle (elt group 0) t)
     (fix-one-bad-angle (elt group 1) nil))
    ((= (length group) 3)
     (fix-one-bad-angle (elt group 0) nil)
     (fix-one-bad-angle (elt group 1) t)
     (fix-one-bad-angle (elt group 2) nil))
    (t
     (fix-one-bad-angle (elt group 0) t)
     (fix-one-bad-angle (elt group 1) nil)
     (fix-one-bad-angle (elt group 2) nil)
     (fix-bad-angle-group (subseq group 3 (length group))))))


(defun fix-bad-angles (sketch)
  (let ((groups (find-bad-angle-groups sketch)))
    (loop for group in groups
          do (progn
               (when (is-ring group)
                 (setf group (maybe-bring-heteroatom-to-front group)))
               (fix-bad-angle-group group)))
    groups))

(defun sketch2d-dynamics (dynamics &key accumulate-coordinates (unfreeze t) frozen (bond-length +stage1-bond-length+))
  (let* ((sketch-function (dynamics:scoring-function dynamics))
         (atom-table (chem:node-table sketch-function))
         (energy-stretch (chem:get-stretch-component sketch-function))
         (energy-sketch-nonbond (chem:get-sketch-nonbond-component sketch-function)))
    (randomize-coordinates (dynamics:coordinates dynamics) :from-zero t :frozen frozen)
    (chem:disable (chem:get-point-to-line-restraint-component sketch-function))
    (apply #'chem:setf-velocity-scale sketch-function *stage1-flatten-force-components*)
    ;; stage 1
    (prepare-stage1-sketch-function sketch-function :bond-length bond-length)
    (chem:set-scale-sketch-nonbond energy-sketch-nonbond +stage1-nonbond-constant+)
    (dotimes (i 5000) (advance-simulation dynamics :frozen frozen))
    ;; stage 2
    (chem:set-scale-sketch-nonbond energy-sketch-nonbond +stage2-nonbond-constant+)
    (dotimes (i 5000) (advance-simulation dynamics :frozen frozen))
    ;; stage 3
    ;; Add the flattening force
    (apply #'chem:setf-velocity-scale sketch-function *stage3-flatten-force-components*)
    (add-flatten-function sketch-function (chem:get-graph sketch-function) *stage3-flatten-force*)
    (if (and frozen unfreeze)
        (let ((center (geometric-center-of-unfrozen-atoms atom-table (dynamics:coordinates dynamics) frozen)))
          (loop for stage from 0 below 4
                do (dotimes (i 1000) (advance-simulation dynamics :frozen frozen))
                do (setf frozen (unfreeze atom-table (dynamics:coordinates dynamics) frozen)))
          (dotimes (i 1000) (advance-simulation dynamics :frozen frozen))
          (transform-unfrozen-atoms-center-to center atom-table (dynamics:coordinates dynamics) frozen))
        (dotimes (i 5000) (advance-simulation dynamics :frozen frozen)))
    ;; stage 4
    ;; Double the bond spring force
    ;; Drop the bond lengths to 1.5
    ;; Add the angle terms
    ;; Progressively reduce the nonbond scale
    (chem:walk-sketch-stretch-terms 
     energy-stretch
     (lambda (index atom1-I1 atom2-I2 kb r0)
       (chem:modify-sketch-stretch-term-kb energy-stretch index +stage4-bond-force+)
       (chem:modify-sketch-stretch-term-r0 energy-stretch index +stage4-bond-length+)))
    (add-angle-terms sketch-function +stage4-bond-force+ +stage4-bond-length+)
    (loop for i from 0 below 5000
          for cutoff from 80.0 downto 1.0 by (/ 79.0 5000.0)
          do (chem:set-long-distance-cutoff energy-sketch-nonbond (* cutoff +stage4-bond-length+))
             (advance-simulation dynamics :frozen (unless unfreeze frozen)))
    ;; stage 5
    (chem:enable (chem:get-point-to-line-restraint-component sketch-function))
    (chem:set-freeze-flags energy-sketch-nonbond 1)
    (loop for i from 0 below 5000
          do (advance-simulation dynamics :frozen (unless unfreeze frozen)))
    (dynamics:write-coordinates-back-to-matter dynamics)
    dynamics
    ))

(defun sketch2d-molecule (molecule &key accumulate-coordinates)
  (let* ((sketch (setup-simulation molecule :accumulate-coordinates accumulate-coordinates))
         (dynamics (dynamics sketch))
         (scoring-function (dynamics:scoring-function dynamics)))
    (sketch2d-dynamics dynamics :accumulate-coordinates accumulate-coordinates)
    sketch))

(defun organize-moment-of-inertia-eigen-system (eigen-vector-matrix eigen-values)
  "Reorganize the eigenvalues and eigen-vector-matrix so that the eigen-values
are in the order (low, middle, high) and the column eigen-vectors are in the same order."
  (let* ((values (list (cons (first eigen-values) 0)
                      (cons (second eigen-values) 1)
                      (cons (third eigen-values) 2)))
         (sorted-values (sort values #'< :key #'car))
         (new-eigen-vectors (geom:make-matrix t))
         (col0 (cdr (first sorted-values)))
         (col1 (cdr (second sorted-values)))
         (col2 (cdr (third sorted-values))))
      (geom:at-row-col-put new-eigen-vectors 0 col0 (geom:at-row-col-get eigen-vector-matrix 0 0))
      (geom:at-row-col-put new-eigen-vectors 1 col0 (geom:at-row-col-get eigen-vector-matrix 1 0))
      (geom:at-row-col-put new-eigen-vectors 2 col0 (geom:at-row-col-get eigen-vector-matrix 2 0))
      (geom:at-row-col-put new-eigen-vectors 0 col1 (geom:at-row-col-get eigen-vector-matrix 0 1))
      (geom:at-row-col-put new-eigen-vectors 1 col1 (geom:at-row-col-get eigen-vector-matrix 1 1))
      (geom:at-row-col-put new-eigen-vectors 2 col1 (geom:at-row-col-get eigen-vector-matrix 2 1))
      (geom:at-row-col-put new-eigen-vectors 0 col2 (geom:at-row-col-get eigen-vector-matrix 0 2))
      (geom:at-row-col-put new-eigen-vectors 1 col2 (geom:at-row-col-get eigen-vector-matrix 1 2))
    (geom:at-row-col-put new-eigen-vectors 2 col2 (geom:at-row-col-get eigen-vector-matrix 2 2))
    (values new-eigen-vectors (mapcar #'car sorted-values))))

(defun align-molecule-horizontal (molecule)
  "Transform the molecule so that it lies with its long axis horizontal."
  (multiple-value-bind (moment-of-inertia center-of-geometry)
      (chem:moment-of-inertia-tensor molecule)
    (chem:translate-all-atoms molecule (geom:v* center-of-geometry -1.0))
    (multiple-value-bind (eigen-vector-matrix eigen-values)
        (geom:eigen-system moment-of-inertia)
      (unless (< (first eigen-values) (second eigen-values) (third eigen-values))
        (multiple-value-setq (eigen-vector-matrix eigen-values)
          (organize-moment-of-inertia-eigen-system eigen-vector-matrix eigen-values)))
      (let ((transposed-transform (geom:transposed3x3 eigen-vector-matrix)))
        (chem:apply-transform-to-atoms molecule transposed-transform)))))



(defgeneric do-sketch2d (matter &key accumulate-coordinates)
  (:documentation "Return an edited molecule that looks like a chemdraw sketch of the molecule. 
The coordinates are all pressed into the X-Y plane and some hydrogens are added and lone-pairs removed."))

(defmethod do-sketch2d ((molecule chem:molecule) &key accumulate-coordinates transform)
  (let* ((sketch (sketch2d-molecule molecule :accumulate-coordinates accumulate-coordinates))
         (dynamics (dynamics sketch))
         (scoring-function (dynamics:scoring-function dynamics)))
    ;; Check for problems
    #+(or)
    (let* ((temp-coordinates (copy-seq (dynamics:coordinates dynamics)))
           (edited-molecule (chem:get-graph (dynamics:scoring-function dynamics))))
      (multiple-value-bind (number-of-problem-areas frozen)
          (identify-problem-areas dynamics)
        (if (> number-of-problem-areas 0)
            (let ((atom-table (chem:atom-table scoring-function)))
              ;; From here we use the coordinates in the atoms
              (format t "There were ~d problem areas~%" number-of-problem-areas)
              ;; Save the coordinates in the atoms
              (chem:load-coordinates-into-vector scoring-function temp-coordinates)
              (chem:reset-sketch-function scoring-function)
              ;; Jostle the non-frozen atoms
              (randomize-atoms atom-table :frozen frozen :width 5.0)
              (chem:load-coordinates-into-vector scoring-function (dynamics:coordinates dynamics))
              (let ((dynamics (sketch2d-dynamics dynamics
                                                 :accumulate-coordinates accumulate-coordinates
                                                 :frozen frozen
                                                 :bond-length +stage4-bond-length+)))
                (multiple-value-bind (new-number-of-problem-areas new-frozen)
                    (identify-problem-areas dynamics)
                  (cond
                    ((= new-number-of-problem-areas number-of-problem-areas)
                     (chem:save-coordinates-from-vector scoring-function temp-coordinates))
                    ((> new-number-of-problem-areas number-of-problem-areas)
                     (chem:save-coordinates-from-vector scoring-function temp-coordinates))
                    (t))
                  (values (chem:get-graph (dynamics:scoring-function dynamics)) dynamics)))))))
    ;; Look for bad angles
    (let* ((energy-function (dynamics:scoring-function dynamics))
           (molecule (chem:get-graph energy-function))
           (fixed (fix-bad-angles sketch)))
      (when fixed
        (chem:load-coordinates-into-vector energy-function (dynamics:coordinates dynamics))
        (loop for i from 0 below 5000
              do (advance-simulation dynamics))
        (dynamics:write-coordinates-back-to-matter dynamics)))
    (let ((result-molecule (chem:get-graph (dynamics:scoring-function dynamics))))
      (align-molecule-horizontal result-molecule)
      (when transform
        (check-type transform geom:m4)
        (chem:apply-transform-to-atoms (molecule sketch) transform))
      sketch)))

(defmethod do-sketch2d ((aggregate chem:aggregate) &key accumulate-coordinates transform)
  (if (= (chem:content-size aggregate) 1)
      (do-sketch2d (chem:content-at aggregate 0) :accumulate-coordinates accumulate-coordinates :transform transform)
      (error "sketch2d only accepts a molecule or an aggregate with a single molecule")))

(defun sketch2d (matter &key accumulate-coordinates transform)
  (do-sketch2d matter :accumulate-coordinates accumulate-coordinates :transform transform))

(defun similar-sketch2d (molecule other-sketch2d &key accumulate-coordinates)
  "Generates a sketch of MOLECULE that determines a common subgraph with the molecule in OTHER-SKETCH2D 
and aligns the new sketch the same way.  Matches atoms using element and name."
  (check-type molecule chem:molecule)
  (check-type other-sketch2d sketch2d)
  (let ((other-molecule (original-molecule other-sketch2d)))
    (multiple-value-bind (equiv diff1 diff2)
        (molecule-graph:compare-molecules molecule other-molecule :atom-match-callback #'molecule-graph:element-and-name-match)
      (let ((new-sketch (setup-simulation molecule :accumulate-coordinates accumulate-coordinates)))
        ;; We need hash-tables for the original atoms to the sketch atoms
        (let ((other-original-to-sketch (make-hash-table))
              (new-original-to-sketch (make-hash-table)))
          (maphash (lambda (key value)
                     (setf (gethash value new-original-to-sketch) key))
                   (sketch-atoms-to-original new-sketch))
          (maphash (lambda (key value)
                     (setf (gethash value other-original-to-sketch) key))
                   (sketch-atoms-to-original other-sketch2d))
          (let* ((dynamics (dynamics new-sketch))
                 (energy-function (dynamics:scoring-function dynamics))
                 (atom-table (chem:atom-table energy-function))
                 (number-of-atoms (chem:get-number-of-atoms atom-table))
                 (frozen (make-array number-of-atoms :element-type 'bit :initial-element 0)))
            (format t "There are ~a atoms in atom-table~%" number-of-atoms)
            (loop for equiv-pair in equiv
                  for new-atom = (car equiv-pair)
                  for other-atom = (cdr equiv-pair)
                  for new-sketch-atom = (gethash new-atom new-original-to-sketch)
                  for other-sketch-atom = (gethash other-atom other-original-to-sketch)
                  for other-pos = (chem:get-position other-sketch-atom)
                  do (unless (or (= 1 (chem:get-atomic-number other-sketch-atom))
                                 (eq :lp (chem:get-element other-sketch-atom)))
                       (let ((new-sketch-atom-index (chem:get-atom-index-or-nil atom-table new-sketch-atom)))
                         ;; Some atoms from the original molecule are not in the sketch molecule
                         ;; so if it's missing - then ignore it.
                         (when new-sketch-atom-index 
                           (chem:set-position new-sketch-atom other-pos)
                           (setf (aref frozen new-sketch-atom-index) 1)))))
            (format t "Using frozen: ~a~%" frozen)
            (chem:load-coordinates-into-vector energy-function (dynamics:coordinates dynamics))
            (sketch2d-dynamics dynamics :frozen frozen :unfreeze nil
                                        :accumulate-coordinates accumulate-coordinates)
            new-sketch
            ))))))
                   
