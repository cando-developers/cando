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

(defconstant +oozp-scale+ 1.0)
(defparameter *stage1-flatten-force-components* (list 0.85 0.85 0.85))
(defparameter *stage3-flatten-force-components* (list 0.85 0.85 0.65))
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
(defparameter *tertiary-carbon* (chem:make-chem-info-graph (chem:compile-smarts "[#6X4H1:0]")))
(defparameter *lp-oxygen* (chem:make-chem-info-graph (chem:compile-smarts "[#8X2:0]")))
(defparameter *lp-nitrogen* (chem:make-chem-info-graph (chem:compile-smarts "[#7X2:0]")))
(defparameter *lp-sulfur* (chem:make-chem-info-graph (chem:compile-smarts "[#16X2:0]")))

(defclass sketch-nonbond-force-field () ())

(defmethod chem:find-atom-type-position ((force-field sketch-nonbond-force-field) type) 0)


(defparameter *quaternary-center* (chem:make-chem-info-graph (chem:compile-smarts "[*X4:0](~[*:1])(~[*:2])(~[*:3])~[*:4]")))
(defparameter *three-bond-center* (chem:make-chem-info-graph (chem:compile-smarts "[*X3:0](~[*:1])(~[*:2])~[*:3]")))
(defparameter *linear-bond-center* (chem:make-chem-info-graph (chem:compile-smarts "[*X2:0](#[*:1])-[*:2]")))
(defparameter *cumulated-double-bond-center* (chem:make-chem-info-graph (chem:compile-smarts "[*X2:0](=[*:1])=[*:2]")))

(defmethod prepare-stage1-sketch-function (sketch-function &key (bond-length +stage1-bond-length+))
  "Generate bond and angle energy terms for a 2D sketch."
  (let* ((molecule (chem:get-matter sketch-function))
         (atom-table (chem:atom-table sketch-function))
         (nonbond-energy (chem:get-sketch-nonbond-component sketch-function))
         (bond-energy (chem:get-stretch-component sketch-function)))
    (chem:set-scale-sketch-nonbond nonbond-energy +stage1-scale-sketch-nonbond+)
    (loop for ia1 from 0 below (1- (chem:get-number-of-atoms atom-table))
          for atom1 = (chem:elt-atom atom-table ia1)
          do (loop for ia2 from (1+ ia1) below (chem:get-number-of-atoms atom-table)
                   for atom2 = (chem:elt-atom atom-table ia2)
                   do (chem:add-sketch-nonbond-term nonbond-energy atom-table atom1 atom2 +stage1-sketch-nonbond-force+)))
    (chem:map-bonds 'nil
                    (lambda (a1 a2 bond-order)
                      (chem:add-stretch-term bond-energy atom-table a1 a2 +first-bond-force+ bond-length))
                    molecule)
    sketch-function))


(defun add-angle-terms (sketch-function bond-force bond-length)
  (let* ((molecule (chem:get-matter sketch-function))
         (angle120-bond-length (* (sqrt 3.0) bond-length))
         (angle-quat-bond-length (* (sqrt 2.0) bond-length))
         (atom-table (chem:atom-table sketch-function))
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
          do (chem:add-stretch-term bond-energy atom-table atom1 atom2 bond-force angle120-bond-length)
             (chem:add-stretch-term bond-energy atom-table atom2 atom3 bond-force angle120-bond-length)
             (chem:add-stretch-term bond-energy atom-table atom3 atom1 bond-force angle120-bond-length))
    ;; Gather up all the unique quat-centers
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
          do (push (cons atom1 atom2) double-nonbond-pairs)
             (push (cons atom1 atom3) double-nonbond-pairs)
             (push (cons atom1 atom4) double-nonbond-pairs)
             (push (cons atom2 atom3) double-nonbond-pairs)
             (push (cons atom2 atom4) double-nonbond-pairs)
             (push (cons atom3 atom4) double-nonbond-pairs))
    ;; Double the nonbond constant for quat center angle outer atom pairs
    (chem:walk-sketch-nonbond-terms
     sketch-nonbond-component
     (lambda (index atom1 atom2 ia1 ia2 constant)
       (let ((key1 (cons atom1 atom2))
             (key2 (cons atom2 atom1)))
         (when (or (member key1 double-nonbond-pairs :test #'equal)
                   (member key2 double-nonbond-pairs :test #'equal))
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
  (let ((atom-table (chem:atom-table sketch-function))
        (oozp-component (chem:get-out-of-zplane-component sketch-function)))
    (cando:do-atoms (atom molecule)
      (chem:add-out-of-zplane-term oozp-component atom-table atom scale 0.0)
      )))

(defun system-preparation (molecule)
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
             idx))
    (let ((mol-copy (chem:matter-copy molecule)))
      (let ((mol-graph (chem:make-molecule-graph-from-molecule mol-copy)))
        (remove-connected-hydrogens mol-copy (find-groups *methyl* mol-graph))
        (remove-connected-hydrogens mol-copy (find-groups *hydroxyl* mol-graph))
        (remove-connected-hydrogens mol-copy (find-groups *thiol* mol-graph))
        (remove-connected-hydrogens mol-copy (find-groups *primary-ammonium* mol-graph))
        (remove-connected-hydrogens mol-copy (find-groups *primary-amine* mol-graph))
        (remove-connected-hydrogens mol-copy (find-groups *secondary-carbon* mol-graph) 1)
        (remove-connected-hydrogens mol-copy (find-groups *tertiary-carbon* mol-graph) 1)
        (let ((idx 0))
          (setf idx (add-lone-pair mol-copy (find-groups *lp-oxygen* mol-graph) idx))
          (setf idx (add-lone-pair mol-copy (find-groups *lp-nitrogen* mol-graph) idx))
          (setf idx (add-lone-pair mol-copy (find-groups *lp-sulfur* mol-graph) idx)))
        (chem:setf-force-field-name mol-copy :.hidden.sketch2d-1)
        (cando:jostle mol-copy 20.0 t)
        mol-copy))))


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

(defparameter *edited-mol* nil)
(defun setup-simulation (molecule &key accumulate-coordinates frozen)
  (let ((edited-mol (system-preparation molecule)))
    (setf *edited-mol* edited-mol)
    (cando:do-atoms (atom edited-mol)
      (chem:set-type atom :sketch))
    (let* ((dummy-sketch-nonbond-ff (make-instance 'sketch-nonbond-force-field))
           (sketch-function (chem:make-sketch-function edited-mol dummy-sketch-nonbond-ff))
           (atom-table (chem:atom-table sketch-function))
           (dynamics (dynamics:make-atomic-simulation sketch-function
                                                      :accumulate-coordinates accumulate-coordinates)))
      (randomize-atoms atom-table :frozen frozen)
      dynamics)))

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
         (molecule (chem:get-matter energy-function))
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

(defun sketch2d-dynamics (dynamics &key accumulate-coordinates frozen (bond-length +stage1-bond-length+))
  (let* ((sketch-function (dynamics:scoring-function dynamics))
         (atom-table (chem:atom-table sketch-function))
         (energy-stretch (chem:get-stretch-component sketch-function))
         (energy-sketch-nonbond (chem:get-sketch-nonbond-component sketch-function)))
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
    (add-flatten-function sketch-function (chem:get-matter sketch-function) *stage3-flatten-force*)
    (if frozen
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
    (chem:walk-stretch-terms 
     energy-stretch
     (lambda (index atom1 atom2 atom1-I1 atom2-I2 kb r0)
       (chem:modify-stretch-term-kb energy-stretch index +stage4-bond-force+)
       (chem:modify-stretch-term-r0 energy-stretch index +stage4-bond-length+)))
    (add-angle-terms sketch-function +stage4-bond-force+ +stage4-bond-length+)
    (loop for i from 0 below 5000
          for cutoff from 80.0 downto 1.0 by (/ 79.0 5000.0)
          do (chem:set-long-distance-cutoff energy-sketch-nonbond (* cutoff +stage4-bond-length+))
             (advance-simulation dynamics))
    ;; stage 5
    (chem:enable (chem:get-point-to-line-restraint-component sketch-function))
    (chem:set-ignore-hydrogens-and-lps energy-sketch-nonbond t)
    (loop for i from 0 below 5000
          do (advance-simulation dynamics))
    (dynamics:write-coordinates-back-to-matter dynamics)
    dynamics
    ))

(defun sketch2d-molecule (molecule &key accumulate-coordinates)
  (let* ((dynamics (setup-simulation molecule :accumulate-coordinates accumulate-coordinates))
         (scoring-function (dynamics:scoring-function dynamics)))
    (sketch2d-dynamics dynamics :accumulate-coordinates accumulate-coordinates)
    dynamics))

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


(defclass sketch2d ()
  ((molecule :initarg :molecule :accessor molecule)
   (dynamics :initarg :dynamics :accessor dynamics)))


(defgeneric do-sketch2d (matter &key accumulate-coordinates)
  (:documentation "Return an edited molecule that looks like a chemdraw sketch of the molecule. 
The coordinates are all pressed into the X-Y plane and some hydrogens are added and lone-pairs removed."))

(defmethod do-sketch2d ((molecule chem:molecule) &key accumulate-coordinates)
  (let* ((dynamics (sketch2d-molecule molecule :accumulate-coordinates accumulate-coordinates))
         (scoring-function (dynamics:scoring-function dynamics)))
    ;; Check for problems
    (let* ((temp-coordinates (copy-seq (dynamics:coordinates dynamics)))
           (edited-molecule (chem:get-matter (dynamics:scoring-function dynamics))))
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
                  (values (chem:get-matter (dynamics:scoring-function dynamics)) dynamics)))))))
    (let ((result-molecule (chem:get-matter (dynamics:scoring-function dynamics))))
      (align-molecule-horizontal result-molecule)
      (values result-molecule dynamics))))

(defmethod do-sketch2d ((aggregate chem:aggregate) &key accumulate-coordinates)
  (if (= (chem:content-size aggregate) 1)
      (do-sketch2d (chem:content-at aggregate 0) :accumulate-coordinates accumulate-coordinates)
      (error "sketch2d only accepts a molecule or an aggregate with a single molecule")))

(defun sketch2d (matter &key accumulate-coordinates)
  (multiple-value-bind (edited-mol dynamics)
      (do-sketch2d matter :accumulate-coordinates accumulate-coordinates)
    (make-instance 'sketch2d :molecule edited-mol
                             :dynamics dynamics)))



