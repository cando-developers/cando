(in-package :sketch2d)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (push :debug-sketch2d *features*))

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

(defclass annotation ()
  ((text :initarg :text :reader text)
   (atom-names :initarg :atom-names :reader atom-names)))

(defclass sketch2d ()
  ((molecule :initarg :molecule :accessor molecule)
   (double-bond-restraints :initarg :double-bond-restraints :accessor double-bond-restraints)
   #+debug-sketch2d(debug-info :initform (make-hash-table) :accessor debug-info)
   (verbose :initarg :verbose :accessor verbose)
   (rings :initarg :rings :accessor rings)
   (aromaticity-info :initarg :aromaticity-info :accessor aromaticity-info)
   (chiral-infos :initform nil :initarg :chiral-infos :accessor chiral-infos)
   (dynamics :initarg :dynamics :accessor dynamics)
   (sketch-atoms-to-original :initarg :sketch-atoms-to-original :accessor sketch-atoms-to-original)
   (original-molecule :initarg :original-molecule :accessor original-molecule)
   (annotations :initform nil :initarg :annotations :reader annotations)))

(defparameter *max-tries* 3)
(defconstant +oozp-scale+ 1.0)
(defconstant +150-degrees-rad+ (* 150.0 0.0174533))
(defparameter *stage1-flatten-force-components* (list 0.85 0.85 0.85))
(defparameter *stage3-flatten-force-components* (list 0.85 0.85 0.25))
(defconstant +stage1-nonbond-constant+ 0.04)
(defconstant +stage2-nonbond-constant+ 1.0)
(defconstant +first-bond-force+ 0.1)
(defconstant +stage1-bond-length+ 5.0)
(defconstant +stage1-bond-length-coordination+ 8.0)
(defconstant +stage1-double-bond-constraint-length+ (* 2.0 +stage1-bond-length+))
(defconstant +stage1-scale-sketch-nonbond+ 0.04)
(defconstant +stage1-sketch-nonbond-force+ 1.0)
(defparameter *stage3-flatten-force* 1.0)
(defconstant +stage4-bond-length+ 5.0) ; 1.5
(defconstant +stage4-stretch-bond-length-scale-factor+ (/ +stage4-bond-length+ +stage1-bond-length+))
(defconstant +stage4-bond-force+ (* 2.0 +first-bond-force+))
(defparameter *stage4-angle-stretch-force* 0.3)
(defconstant +angle-120+ (* 0.0174533 120.0))
(defconstant +angle-180+ (* 0.0174533 180.0))

(defparameter *methyl-smarts* #+nosmarts nil #-nosmarts(chem:compile-smarts "[#6H3:0]"))
(defparameter *methyl* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph *methyl-smarts*))
(defparameter *hydroxyl* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[#8H1:0]")))
(defparameter *thiol* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[#16H1:0]")))
(defparameter *primary-ammonium* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[#7H3:0]")))
(defparameter *primary-amine* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[#7H2:0]")))
(defparameter *secondary-carbon* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[#6X4H2:0]")))
(defparameter *secondary-carbon-missing-h* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[#6X2H0:0](-*)-*")))
(defparameter *tertiary-carbon* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[#6X4H1:0]")))
(defparameter *tertiary-amine-with-h* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[#7X4H1:0]")))
(defparameter *lp-oxygen* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[#8X2:0]")))
(defparameter *lp-nitrogen* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[#7X2:0]")))
(defparameter *lp-sulfur* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[#16X2:0]")))
(defparameter *bad-angle* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[*X3:0](~[*X1:1])(~[*:2])~[*:3]")))
(defparameter *ring6-across-length* 3.16)
(defparameter *ring3-long-length* 2.898)
(defparameter *ring3-short-length* 1.5)
(defparameter *exo-ring3* (chem:make-chem-info-graph (chem:compile-smarts "[*:0]1(-[*:1])[*:2][*:3]1")))
(defparameter *ring4-long-length* 2.77)
(defparameter *ring4-short-length* 2.123)
(defparameter *exo-ring4* (chem:make-chem-info-graph (chem:compile-smarts "[*:0]1(-[*:1])[*:2][*][*:3]1")))
(defparameter *exo-ring5* (chem:make-chem-info-graph (chem:compile-smarts "[*:0]1(-[*:1])[*:2][*][*:3][*:4]1")))
(defparameter *exo-ring6* (chem:make-chem-info-graph (chem:compile-smarts "[*:0]1(-[*:1])[*:2][*][*:3][*:4][*:5]1")))

(defparameter *bis-amino-acid-alpha* #+nosmarts nil #-nosmarts(chem:compile-smarts "[C:0]1(-[C:1]=[O])CCCN1"))

(defclass sketch-nonbond-force-field () ())

(defmethod chem:find-atom-type-position ((force-field sketch-nonbond-force-field) type)
  (declare (ignore force-field type))
  0)


(defparameter *quaternary-center* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[*X4:0](~[*:1])(~[*:2])(~[*:3])~[*:4]")))
(defparameter *three-bond-center* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[*X3:0](~[*:1])(~[*:2])~[*:3]")))
(defparameter *linear-bond-center* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[*X2:0](#[*:1])-[*:2]")))
(defparameter *cumulated-double-bond-center* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[*X2:0](=[*:1])=[*:2]")))

(defun add-stretch (bond-energy a1 a2 a1ci a2ci bond-force bond-length)
  #+(or)(format t "Adding sketch-stretch-term ~a ~a ~a ~a~%" (chem:get-name a1) (chem:get-name a2) bond-force bond-length)
  (chem:add-sketch-stretch-term bond-energy a1ci a2ci bond-force bond-length))

(defun prepare-stage1-sketch-function (sketch-function sketch)
  "Generate bond energy terms for a 2D sketch."
  (let* ((molecule (chem:get-graph sketch-function))
         (atom-table (chem:node-table sketch-function))
         (nonbond-energy (chem:get-sketch-nonbond-component sketch-function))
         (bond-energy (chem:get-stretch-component sketch-function)))
    (chem:set-scale-sketch-nonbond nonbond-energy +stage1-scale-sketch-nonbond+)
    (loop for ia1 from 0 below (1- (chem:get-number-of-atoms atom-table))
          for atom1 = (chem:elt-atom atom-table ia1)
          for atom1-coordinate-index = (chem:get-coordinate-index-times3 atom-table atom1)
          for freeze1 = (or (= (chem:get-atomic-number atom1) 1)
                            (eq (chem:get-element atom1) :lp))
          do (loop for ia2 from (1+ ia1) below (chem:get-number-of-atoms atom-table)
                   for atom2 = (chem:elt-atom atom-table ia2)
                   for atom2-coordinate-index = (chem:get-coordinate-index-times3 atom-table atom2)
                   for freeze2 = (or (= (chem:get-atomic-number atom2) 1)
                                     (eq (chem:get-element atom2) :lp))
                   for freeze-flag = (if (or freeze1 freeze2) 1 0)
                   do (chem:add-sketch-nonbond-term nonbond-energy
                                                    atom1-coordinate-index
                                                    atom2-coordinate-index
                                                    freeze-flag
                                                    +stage1-sketch-nonbond-force+)))
    (chem:map-bonds
     'nil
     (lambda (a1 a2 bond-order bond)
       (declare (ignore bond-order))
       (let ((atom1-coord-index (chem:get-coordinate-index-times3 atom-table a1))
             (atom2-coord-index (chem:get-coordinate-index-times3 atom-table a2)))
         (let* ((valence1 (chem:matter-get-property-or-default a1 :valence 1))
                (valence2 (chem:matter-get-property-or-default a2 :valence 1))
                (bond-length (if (or (> valence1 4)
                                     (> valence2 4))
                                 +stage1-bond-length-coordination+
                                 +stage1-bond-length+)))
           #+(or)(format t "a1 ~a a2 ~a   valence1 ~a valence2 ~a~%" (chem:get-name a1) (chem:get-name a2) valence1 valence2)
           (add-stretch bond-energy
                                    a1 a2
                                    atom1-coord-index atom2-coord-index
                                    +first-bond-force+ bond-length))))
     molecule)
    (loop for double-bond-restraint in (double-bond-restraints sketch)
          do (with-slots (top-left-atom bottom-left-atom left-atom right-atom top-right-atom bottom-right-atom) double-bond-restraint
               (let ((top-left-index (chem:get-coordinate-index-times3 atom-table top-left-atom))
                     (bottom-left-index (chem:get-coordinate-index-times3 atom-table bottom-left-atom))
                     (top-right-index (chem:get-coordinate-index-times3 atom-table top-right-atom))
                     (bottom-right-index (chem:get-coordinate-index-times3 atom-table bottom-right-atom)))
                 (add-stretch bond-energy
                              top-left-atom top-right-atom
                              top-left-index top-right-index
                              +first-bond-force+
                              (* 2.0 +stage1-bond-length+))
                 (add-stretch bond-energy
                              bottom-left-atom bottom-right-atom
                              bottom-left-index bottom-right-index
                              +first-bond-force+
                              (* 2.0 +stage1-bond-length+)))))
    sketch-function))

(defun debug-set-unique-atom-names (matter inc)
  "Assign each atom in the sketch a unique name"
  (let ((idx 0))
    (chem:do-residues (res matter)
      (chem:do-atoms (atm res)
        (let* ((element (chem:get-element atm))
               (name-str (format nil "~a~a" (string element) idx)))
          (when inc (incf idx))
          (chem:set-name atm (intern name-str :keyword)))))))

(defun ordered-key (ia1 ia2)
  (if (< ia1 ia2)
      (cons ia1 ia2)
      (cons ia2 ia1)))

(defun filter-using-first-atom (matches)
  (let ((unique-first-atom-matches (make-hash-table)))
    (loop for match in matches
          for first-atom = (elt match 0)
          do (setf (gethash first-atom unique-first-atom-matches) match))
    (let ((unique-matches nil))
      (maphash (lambda (key match)
                 (declare (ignore key))
                 (push match unique-matches))
               unique-first-atom-matches)
      unique-matches)))

(defun describe-terms (sketch-function name)
  (format t "------ ~a~%" name)
  (let* ((atom-table (chem:node-table sketch-function))
         (energy-nonbond (chem:get-sketch-nonbond-component sketch-function))
         (nonbond-scale (chem:get-scale-sketch-nonbond energy-nonbond))
         (energy-stretch (chem:get-stretch-component sketch-function)))
    (format t "Stretch terms~%")
    (chem:walk-sketch-stretch-terms
     energy-stretch
     (lambda (index ia1 ia2 kb r0)
       (let ((atom1 (chem:elt-atom atom-table (/ ia1 3)))
             (atom2 (chem:elt-atom atom-table (/ ia2 3))))
         (format t "stretch: ~a ~a ~a ~a~%" (chem:get-name atom1) (chem:get-name atom2) kb r0))))
    (format t "Nonbond terms scale: ~a~%" nonbond-scale)
    (chem:walk-sketch-nonbond-terms
     energy-nonbond
     (lambda (index freeze-flags ia1 ia2 constant)
       (declare (ignore index freeze-flags))
       (let ((atom1 (chem:elt-atom atom-table (/ ia1 3)))
             (atom2 (chem:elt-atom atom-table (/ ia2 3))))
         (format t "nonbond: ~a ~a ~a~%" (chem:get-name atom1) (chem:get-name atom2) constant))))))

(defun add-angle-terms (sketch sketch-function bond-force bond-length)
  (let* ((molecule (chem:get-graph sketch-function))
         (angle120-bond-length (* (sqrt 3.0) bond-length))
         (atom-table (chem:node-table sketch-function))
         (bond-energy (chem:get-stretch-component sketch-function))
         (sketch-nonbond-component (chem:get-sketch-nonbond-component sketch-function))
         (mol-graph (chem:make-molecule-graph-from-molecule molecule))
         #+(or)(linear-bond-centers (chem:boost-graph-vf2 *linear-bond-center* mol-graph))
         #+(or)(cumulated-double-bond-centers (chem:boost-graph-vf2 *cumulated-double-bond-center* mol-graph))
         (three-bond-centers (chem:boost-graph-vf2 *three-bond-center* mol-graph))
         (unique-three-bond-centers (make-hash-table :test #'equal))
         (quat-centers (chem:boost-graph-vf2 *quaternary-center* mol-graph))
         (ring3-centers (filter-using-first-atom (chem:boost-graph-vf2 *exo-ring3* mol-graph)))
         (ring4-centers (filter-using-first-atom (chem:boost-graph-vf2 *exo-ring4* mol-graph)))
         (ring5-centers (filter-using-first-atom (chem:boost-graph-vf2 *exo-ring5* mol-graph)))
         (ring6-centers (filter-using-first-atom (chem:boost-graph-vf2 *exo-ring6* mol-graph)))
         (ring-centers (let ((ring-centers (make-hash-table)))
                         (loop for center in ring3-centers
                               do (setf (gethash center ring-centers) t))
                         (loop for center in ring4-centers
                               do (setf (gethash center ring-centers) t))
                         (loop for center in ring5-centers
                               do (setf (gethash center ring-centers) t))
                         (loop for center in ring6-centers
                               do (setf (gethash center ring-centers) t))
                         ring-centers))
         (unique-quat-centers (make-hash-table :test #'equal))
         double-nonbond-pairs
         (atom-ids (make-hash-table)))
    (let ((idx 0))
      (chem:do-atoms (atm molecule)
        (setf (gethash atm atom-ids) (incf idx))))
    (progn
      (loop for tbc in three-bond-centers
            for center = (elt tbc 0)
            for atom1 = (elt tbc 1)
            for atom2 = (elt tbc 2)
            for atom3 = (elt tbc 3)
            for center-id = (gethash center atom-ids)
            for atom1-id = (gethash atom1 atom-ids)
            for atom2-id = (gethash atom2 atom-ids)
            for atom3-id = (gethash atom3 atom-ids)
            for sorted-ids = (sort (list atom1-id atom2-id atom3-id) #'<)
            when (null (gethash center ring-centers)) ; ignore centers that are in 3-, 4-, 5-, or 6-member rings
              do (progn
                   #+(or)(progn
                           (format t "angle bond ~a ~a~%" atom1 atom2)
                           (format t "angle bond ~a ~a~%" atom2 atom3)
                           (format t "angle bond ~a ~a~%" atom3 atom1))
                   (setf (gethash (cons center-id sorted-ids) unique-three-bond-centers) (list center atom1 atom2 atom3))))
      (loop for tbc being the hash-values of unique-three-bond-centers
            for center = (elt tbc 0)
            for atom1 = (elt tbc 1)
            for atom2 = (elt tbc 2)
            for atom3 = (elt tbc 3)
            for atom1-coord-index = (chem:get-coordinate-index-times3 atom-table atom1)
            for atom2-coord-index = (chem:get-coordinate-index-times3 atom-table atom2)
            for atom3-coord-index = (chem:get-coordinate-index-times3 atom-table atom3)
            do (add-stretch bond-energy atom1 atom2 atom1-coord-index atom2-coord-index bond-force angle120-bond-length)
               (add-stretch bond-energy atom2 atom3 atom2-coord-index atom3-coord-index bond-force angle120-bond-length)
               (add-stretch bond-energy atom3 atom1 atom3-coord-index atom1-coord-index bond-force angle120-bond-length)
            ))
    #+(or)
    (loop for ring3 in ring3-centers
          for center = (elt ring3 0)
          for exo = (elt ring3 1)
          for left = (elt ring3 2)
          for right = (elt ring3 3)
          for exo-coord-index = (chem:get-coordinate-index-times3 atom-table exo)
          for left-coord-index = (chem:get-coordinate-index-times3 atom-table left)
          for right-coord-index = (chem:get-coordinate-index-times3 atom-table right)
          do (format t "exocyclic ring3 ~a -> ~a ~a ~a~%" exo center left right)
          do (add-stretch bond-energy exo left exo-coord-index left-coord-index bond-force *ring3-long-length*)
             (add-stretch bond-energy exo right exo-coord-index right-coord-index bond-force *ring3-long-length*)
             (add-stretch bond-energy left right left-coord-index right-coord-index bond-force *ring3-short-length*))
    #+(or)
    (loop for ring4 in ring4-centers
          for center = (elt ring4 0)
          for exo = (elt ring4 1)
          for left = (elt ring4 2)
          for right = (elt ring4 3)
          for exo-coord-index = (chem:get-coordinate-index-times3 atom-table exo)
          for left-coord-index = (chem:get-coordinate-index-times3 atom-table left)
          for right-coord-index = (chem:get-coordinate-index-times3 atom-table right)
          do (format t "exocyclic ring4 ~a -> ~a ~a ~a~%" exo center left right)
          do (add-stretch bond-energy exo left exo-coord-index left-coord-index bond-force *ring4-long-length*)
             (add-stretch bond-energy exo right exo-coord-index right-coord-index bond-force *ring4-long-length*)
             (add-stretch bond-energy left right left-coord-index right-coord-index bond-force *ring4-short-length*))
    #+(or)
    (loop for ring in (rings sketch)
          when (= (length ring) 6)
            do (let (
                     (a0 (elt ring 0))
                     (a1 (elt ring 1))
                     (a2 (elt ring 2))
                     (a3 (elt ring 3))
                     (a4 (elt ring 4))
                     (a5 (elt ring 5))
                     (i0 (chem:get-coordinate-index-times3 atom-table a0))
                     (i1 (chem:get-coordinate-index-times3 atom-table a1))
                     (i2 (chem:get-coordinate-index-times3 atom-table a2))
                     (i3 (chem:get-coordinate-index-times3 atom-table a3))
                     (i4 (chem:get-coordinate-index-times3 atom-table a4))
                     (i5 (chem:get-coordinate-index-times3 atom-table a5)))
                 (add-stretch bond-energy a0 a3 i0 i3 bond-force *ring6-across-length*)
                 (add-stretch bond-energy a1 a4 i1 i4 bond-force *ring6-across-length*)
                 (add-stretch bond-energy a2 a5 i2 i5 bond-force *ring6-across-length*)))
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
            for atom1-index = (chem:get-coordinate-index-times3 atom-table atom1)
            for atom2-index = (chem:get-coordinate-index-times3 atom-table atom2)
            for atom3-index = (chem:get-coordinate-index-times3 atom-table atom3)
            for atom4-index = (chem:get-coordinate-index-times3 atom-table atom4)
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
    (chem:do-atoms (atom molecule)
      (let ((atom-coord-index (chem:get-coordinate-index-times3 atom-table atom)))
      (chem:add-out-of-zplane-term oozp-component atom-coord-index scale 0.0))
      )))

(defun add-atom-annotations (molecule)
  (chem:do-atoms (atom molecule)
    (let ((hydrogens 0))
      (cond
        ((/= (chem:get-atomic-number atom) 1)
         (let ((bonds (chem:bonds-as-list atom)))
           (loop for bond in bonds
                 for other-atom = (chem:bond/get-other-atom bond atom)
                 when (= (chem:get-atomic-number other-atom) 1)
                   do (incf hydrogens))
           (chem:set-property atom :hydrogens hydrogens)
           #+(or)(format t "atom: ~a  bonds: ~a  valence: ~a~%" atom bonds (chem:get-valence atom))
           (chem:set-property atom :valence (chem:get-valence atom))))
        (t (chem:set-property atom :valence 1))))))

(defun add-angle-bonds (molecule)
  ;; Add special bonds to constrain angles - across all atoms that share a common neighbor
  (let ((angle-restraint-bonds nil))
    (chem:map-angles
     'nil
     (lambda (a1 a2 a3)
       (declare (ignore a2))
       (push (cons a1 a3) angle-restraint-bonds))
     molecule)
    (loop for angle-restraint in angle-restraint-bonds
          for a1 = (car angle-restraint)
          for a2 = (cdr angle-restraint)
          do (let ((bond (chem:bond-to a1 a2 :virtual-bond t nil)))
               (chem:bond/set-property bond :special-bond :angle-bond-restraint)))))

(defclass double-bond-atoms ()
  ((top-left-atom :initarg :top-left-atom :accessor top-left-atom)
   (bottom-left-atom :initarg :bottom-left-atom :accessor bottom-left-atom)
   (left-atom :initarg :left-atom :accessor left-atom)
   (right-atom :initarg :right-atom :accessor right-atom)
   (top-right-atom :initarg :top-right-atom :accessor top-right-atom)
   (bottom-right-atom :initarg :bottom-right-atom :accessor bottom-right-atom)
   ))

(defmethod print-object ((dba double-bond-atoms) stream)
  (print-unreadable-object (dba stream :type t)
    (with-slots (top-left-atom bottom-left-atom left-atom right-atom top-right-atom bottom-right-atom) dba
      (format stream "~a/~a ~a=~a ~a/~a"
              (chem:get-name top-left-atom)
              (chem:get-name bottom-left-atom)
              (chem:get-name left-atom)
              (chem:get-name right-atom)
              (chem:get-name top-right-atom)
              (chem:get-name bottom-right-atom)
              )
      )))

(defun double-bond-atoms-from-names (a1 a2 atom-names)
  (let ((left-atom-name (third atom-names))
        (right-atom-name (fourth atom-names)))
    (multiple-value-bind (left-atom right-atom)
        (if (and (eq (chem:get-name a1) left-atom-name)
                 (eq (chem:get-name a2) right-atom-name))
            (values a1 a2)
            (if (and (eq (chem:get-name a2) left-atom-name)
                     (eq (chem:get-name a1) right-atom-name))
                (values a2 a1)
                (error "double-bond-geometry names ~a ~a don't match double bond atoms ~a ~a" left-atom-name right-atom-name a1 a2)))
      (flet ((find-vinyl-atom (double vinyl-name)
               (loop for atm in (chem:atom/bonded-atoms-as-list double)
                     when (eq (chem:get-name atm) vinyl-name)
                       do (return-from find-vinyl-atom atm))
               (error "Could not find vinyl atom with name ~a on double bond atom ~a" vinyl-name double)))
        (let ((top-left-atom (find-vinyl-atom left-atom (first atom-names)))
              (bottom-left-atom (find-vinyl-atom left-atom (second atom-names)))
              (top-right-atom (find-vinyl-atom right-atom (fifth atom-names)))
              (bottom-right-atom (find-vinyl-atom right-atom (sixth atom-names))))
          (make-instance 'double-bond-atoms
                         :top-left-atom top-left-atom
                         :bottom-left-atom bottom-left-atom
                         :left-atom left-atom
                         :right-atom right-atom
                         :top-right-atom top-right-atom
                         :bottom-right-atom bottom-right-atom))))))

(defun calculate-extra-bonds-for-angles-and-double-bonds (molecule)
  (let ((double-bond-restraint-bonds nil))
    (chem:map-bonds
     'nil
     (lambda (a1 a2 bond-order bond)
       (declare (ignore bond-order))
       (let* ((atom-names (chem:bond/get-property bond :double-bond-geometry nil)))
         (when atom-names
           (let ((dba (double-bond-atoms-from-names a1 a2 atom-names)))
             (push dba double-bond-restraint-bonds)
             #+(or)
             (with-slots (top-left-atom bottom-left-atom left-atom right-atom top-right-atom bottom-right-atom) dba
               (format t "double-bond ~s ~s~%" left-atom right-atom)
               (format t "double-bond cis ~a ~a~%" top-left-atom top-right-atom)
               (format t "double-bond cis bottom ~a ~a~%" bottom-left-atom bottom-right-atom))
             ))))
     molecule)
    (values double-bond-restraint-bonds)))

(defun system-preparation (molecule &key verbose)
  (labels ((find-groups (chem-info-graph molecule-graph)
             (loop for match in (chem:boost-graph-vf2 chem-info-graph molecule-graph)
                   collect (elt match 0)))
           (remove-connected-hydrogens (mol sketch-atoms-to-original atoms &optional (number-to-remove 999) verbose)
             (declare (ignore verbose))
             (let ((connected-hydrogens nil))
               (chem:do-residues (res mol)
                 (chem:do-atoms (atm res)
                   (when (member atm atoms)
                     (let ((neighbors (chem:atom/bonded-atoms-as-list atm))
                           (number-removed 0))
                       (loop for natm in neighbors
                             until (>= number-removed number-to-remove)
                             when (eq (chem:get-element natm) :H)
                               do (progn
                                    (incf number-removed)
                                    (push (cons natm res) connected-hydrogens)))))))
               (if connected-hydrogens
                   (progn
                     (loop for (hatom . res) in connected-hydrogens
                           do (remhash hatom sketch-atoms-to-original)
                           do (chem:remove-atom-delete-bonds res hatom))))))
           (add-lone-pair (mol atoms idx atom-to-residue-map)
             (declare (ignore mol))
             (loop for atom in atoms
                   for lp-idx from idx
                   for lp-name = (intern (format nil "lp~a" lp-idx) :keyword)
                   for lp = (chem:make-atom lp-name :lp)
                   for res = (gethash atom atom-to-residue-map)
                   do (incf idx)
                      (chem:add-matter res lp)
                      (chem:bond-to atom lp :single-bond t nil)
                   )
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
      ;; Set names here
      #+(or)(when *show-names* (debug-set-unique-atom-names mol-copy t))
      (add-atom-annotations mol-copy)
      (let ((mol-graph (chem:make-molecule-graph-from-molecule mol-copy)))
        (remove-connected-hydrogens mol-copy sketch-atoms-to-original (find-groups *methyl* mol-graph) 3)
        (remove-connected-hydrogens mol-copy sketch-atoms-to-original (find-groups *hydroxyl* mol-graph) 1)
        (remove-connected-hydrogens mol-copy sketch-atoms-to-original (find-groups *thiol* mol-graph) 1)
        (remove-connected-hydrogens mol-copy sketch-atoms-to-original (find-groups *primary-ammonium* mol-graph) 3)
        (remove-connected-hydrogens mol-copy sketch-atoms-to-original (find-groups *primary-amine* mol-graph) 2)
        (remove-connected-hydrogens mol-copy sketch-atoms-to-original (find-groups *secondary-carbon* mol-graph) 1)
        (remove-connected-hydrogens mol-copy sketch-atoms-to-original (find-groups *tertiary-carbon* mol-graph) 1 t)
        (remove-connected-hydrogens mol-copy sketch-atoms-to-original (find-groups *tertiary-amine-with-h* mol-graph) 1)
        (let ((idx 0)
              (atom-to-residue-map (make-hash-table))
              (atoms-need-lps nil))
          ;; Map atoms to residues
          (chem:do-residues (res mol-copy)
            (chem:do-atoms (atm res)
              (setf (gethash atm atom-to-residue-map) res)))
          ;; Collect all atoms that need lone-pairs.
          (chem:do-atoms (atm mol-copy)
            (when (and (member (chem:get-element atm) '(:O :N :S)) (= (chem:number-of-bonds atm) 2))
              (push atm atoms-need-lps)))
          ;; Add lone pairs to atoms that need them
          (let ((missing-h (unique-secondary-carbon-missing-h mol-graph)))
            (setf idx (add-lone-pair mol-copy missing-h idx atom-to-residue-map)))
          (setf idx (add-lone-pair mol-copy atoms-need-lps idx atom-to-residue-map)))
        ;; Add extra bonds for angles
        (let ((double-bond-restraints (calculate-extra-bonds-for-angles-and-double-bonds mol-copy)))
          (chem:setf-force-field-name mol-copy :.hidden.sketch2d-1)
;;;        (cando:jostle mol-copy 20.0 t)
          (make-instance 'sketch2d
                         :molecule mol-copy
                         :verbose verbose
                         :double-bond-restraints double-bond-restraints
                         :rings rings
                         :aromaticity-info aromaticity-info
                         :original-molecule molecule
                         :sketch-atoms-to-original sketch-atoms-to-original))))))


(defun randomize-atoms (atom-table &key unfrozen from-zero (width 40.0) &aux (half-width (/ width 2.0)))
  "Randomly jostle atoms from their current positions"
  (flet ((jostle-atom (atom)
           (when from-zero
             (chem:set-position atom (geom:vec 0.0 0.0 0.0)))
           (let* ((cp (chem:get-position atom))
                  (pos (geom:vec
                        (+ (- (random width) half-width) (geom:get-x cp))
                        (+ (- (random width) half-width) (geom:get-y cp))
                        (+ (- (random width) half-width) (geom:get-z cp)))))
             (chem:set-position atom pos))))
    (loop for index from 0 below (chem:get-number-of-atoms atom-table)
          for atom = (chem:elt-atom atom-table index)
          do (if unfrozen
                 (when (= (elt unfrozen index) 1)
                   (jostle-atom atom))
                 (jostle-atom atom)))))

(defun randomize-coordinates (coordinates &key unfrozen from-zero (width 40.0) &aux (half-width (/ width 2.0)))
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
          do (if unfrozen
                 (when (= (elt unfrozen index) 1)
                   (jostle-atom coord-index))
                 (jostle-atom coord-index)))))

(defparameter *edited-mol* nil)
(defun setup-simulation (molecule &key accumulate-coordinates verbose)
  (let* ((sketch (system-preparation molecule :verbose verbose))
         (edited-mol (molecule sketch))
         (atom-types (make-hash-table)))
    (setf *edited-mol* edited-mol)
    (chem:do-atoms (atom edited-mol)
      (setf (gethash atom atom-types) :sketch))
    (let* ((dummy-sketch-nonbond-ff (make-instance 'sketch-nonbond-force-field))
           (sketch-function (chem:make-sketch-function edited-mol dummy-sketch-nonbond-ff atom-types))
           (dynamics (dynamics:make-atomic-simulation sketch-function
                                                      :accumulate-coordinates accumulate-coordinates)))
      (setf (dynamics sketch) dynamics)
      sketch)))

(defun advance-simulation (dynamics &key unfrozen)
  (dynamics:velocity-verlet-step
   dynamics
   :velocity-verlet-function #'chem:sketch-function-velocity-verlet-step
   :unfrozen unfrozen))



#|
Use the approach described at:
https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
to check if two line segments (bonds) overlap/intersect
|#
(defun orientation (p1 p2 p3)
  (let ((orient (- (* (- (geom:get-y p2) (geom:get-y p1))
                      (- (geom:get-x p3) (geom:get-x p2)))
                   (* (- (geom:get-y p3) (geom:get-y p2))
                      (- (geom:get-x p2) (geom:get-x p1))))))
    (cond
      ((= 0.0 orient) 0)
      ((< 0.0 orient) 1)
      (t 2))))

;;; Given three colinear points p, q, r, the function checks if
;;; point q lies on line segment 'pr'
(defun on-segment (p q r)
  (and (<= (geom:get-x q) (max (geom:get-x p) (geom:get-x r)))
       (>= (geom:get-x q) (max (geom:get-x p) (geom:get-x r)))
       (<= (geom:get-y q) (max (geom:get-y p) (geom:get-y r)))
       (>= (geom:get-y q) (max (geom:get-y p) (geom:get-y r)))))

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
  (let ((x1 (geom:get-x p1))
        (y1 (geom:get-y p1))
        (x2 (geom:get-x q1))
        (y2 (geom:get-y q1))
        (x3 (geom:get-x p2))
        (y3 (geom:get-y p2))
        (x4 (geom:get-x q2))
        (y4 (geom:get-y q2)))
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


(defun good-double-bond-geometry (double-bond-atoms)
  (flet ((angle-degrees (a1 a2 a3)
           (let ((v1 (chem:get-position a1))
                 (v2 (chem:get-position a2))
                 (v3 (chem:get-position a3)))
             (/ (geom:calculate-angle v1 v2 v3) 0.0174533)))
         (good-angle (angle)
           (let ((good (< (abs (- angle 120.0)) 20.0)))
             good)))
    (with-slots (top-left-atom bottom-left-atom
                 left-atom right-atom
                 top-right-atom bottom-right-atom) double-bond-atoms
      (let ((top-left-angle (angle-degrees top-left-atom left-atom right-atom))
            (bottom-left-angle (angle-degrees bottom-left-atom left-atom right-atom))
            (top-right-angle (angle-degrees top-right-atom right-atom left-atom))
            (bottom-right-angle (angle-degrees bottom-right-atom right-atom left-atom)))
        (and (good-angle top-left-angle)
             (good-angle bottom-left-angle)
             (good-angle top-right-angle)
             (good-angle bottom-right-angle))))))

(defun detect-bad-double-bonds (sketch)
  (loop for double-bond-restraint in (double-bond-restraints sketch)
        unless (good-double-bond-geometry double-bond-restraint)
          collect double-bond-restraint))

(defstruct stretch-term atom1 atom2)

(defun detect-intersecting-bonds (matter)
  (let ((stretches (chem:map-bonds
                    'vector
                    (lambda (a1 a2 order bond)
                      (declare (ignore order bond))
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
          unless (or (eq (chem:get-element atomp1) :h)
                     (eq (chem:get-element atomq1) :h))
            do (loop named inner
                     for idx2 from (1+ idx1) below (length stretches)
                     for str2 = (elt stretches idx2)
                     for atomp2 = (stretch-term-atom1 str2)
                     for atomq2 = (stretch-term-atom2 str2)
                     for p2 = (chem:get-position atomp2)
                     for q2 = (chem:get-position atomq2)
                     unless (or (eq (chem:get-element atomp2) :h)
                                (eq (chem:get-element atomq2) :h))
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

(defclass overlap-problem-area ()
  ((atoms :initarg :atoms :accessor atoms)))

(defgeneric size-of (problem-area))
(defmethod size-of ((problem-area overlap-problem-area))
  (length (atoms problem-area)))

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
    (make-instance 'overlap-problem-area
                   :atoms close-atoms)))


(defun generate-double-bond-problem-area (double-bond-restraint cutoff atom-table)
  (let* ((number-of-atoms (chem:get-number-of-atoms atom-table))
         (point (geom:v* (geom:v+ (chem:get-position (left-atom double-bond-restraint))
                                  (chem:get-position (right-atom double-bond-restraint))) 0.5))
         (double-bond-atoms (with-slots (top-left-atom bottom-left-atom
                                         left-atom right-atom
                                         top-right-atom bottom-right-atom) double-bond-restraint
                              (list top-left-atom bottom-left-atom
                                    left-atom right-atom
                                    top-right-atom bottom-right-atom)))
         close-atoms)
    (loop for ia from 0 below number-of-atoms
          for atom = (chem:elt-atom atom-table ia)
          for atom-pos = (chem:get-position atom)
          for delta = (geom:v- atom-pos point)
          for distance = (geom:vlength delta)
          do (when (or (member atom double-bond-atoms)
                       (< distance cutoff))
               (push (make-distance-atom :distance distance
                                         :index ia
                                         :atom atom)
                     close-atoms)))
    (make-instance 'overlap-problem-area
                   :atoms close-atoms)))

(defun identify-problem-areas (sketch)
  (let* ((dynamics (dynamics sketch))
         (cutoff (* 2.4 +stage4-bond-length+))
         (energy-function (dynamics:scoring-function dynamics))
         (molecule (chem:get-graph energy-function))
         (atom-table (chem:atom-table energy-function))
         #+(or)(energy-stretch (chem:get-stretch-component energy-function))
         (intersections (detect-intersecting-bonds molecule))
         (bad-double-bonds (detect-bad-double-bonds sketch)))
    (if (or intersections bad-double-bonds)
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
                     ((> (size-of problem-area)
                         (size-of worst-problem-area))
                      (setf worst-problem-area problem-area))))
          (loop for double-bond in bad-double-bonds
                for problem-area = (generate-double-bond-problem-area double-bond cutoff atom-table)
                do (cond
                     ((null worst-problem-area)
                      (setf worst-problem-area problem-area))
                     ((> (size-of problem-area)
                         (size-of worst-problem-area))
                      (setf worst-problem-area problem-area))))
          (if worst-problem-area
              (let* ((number-of-atoms (chem:get-number-of-atoms atom-table))
                     (unfrozen (make-array number-of-atoms :element-type 'bit :initial-element 0)))
                (loop for isunfrozen in (atoms worst-problem-area)
                      for index = (distance-atom-index isunfrozen)
                      for atom = (distance-atom-atom isunfrozen)
                      do (setf (elt unfrozen index) 1))
                (values (+ (length bad-double-bonds) (length intersections)) unfrozen))
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

(defun unfreeze (atom-table coordinates unfrozen)
  (let ((new-unfrozen (copy-seq unfrozen)))
    (loop for index from 0 below (length unfrozen)
          for atom = (chem:elt-atom atom-table index)
          do (when (= (elt unfrozen index) 1)
               (loop named inner
                     for bondi from 0 below (chem:number-of-bonds atom)
                     for neighbor = (chem:bonded-neighbor atom bondi)
                     for neighbor-index*3 = (chem:get-coordinate-index-times3 atom-table neighbor)
                     for neighbor-index = (/ neighbor-index*3 3)
                     do (when (= (elt unfrozen neighbor-index) 0)
                          (progn
                            (setf (elt new-unfrozen neighbor-index) 1)
                            (randomize-neighbor-near-unfrozen-atom coordinates neighbor-index*3 (* index 3))
                            (return-from inner nil))))))
    new-unfrozen))

(defun geometric-center-of-unfrozen-atoms (atom-table coordinates unfrozen)
  (declare (ignore atom-table))
  (let ((xpos 0.0)
        (ypos 0.0)
        (zpos 0.0)
        (count 0))
    (loop for coord-index from 0 below (length unfrozen) by 3
          for bit-index from 0
          do (when (= 1 (elt unfrozen bit-index))
               (setf xpos (+ (elt coordinates coord-index) xpos)
                     ypos (+ (elt coordinates (+ coord-index 1)) ypos)
                     zpos (+ (elt coordinates (+ coord-index 2)) zpos))
               (incf count)))
    (geom:vec (/ xpos (float count 1.0d0))
              (/ ypos (float count 1.0d0))
              (/ zpos (float count 1.0d0)))))

(defun transform-unfrozen-atoms-center-to (to-center atom-table coordinates unfrozen)
  (let* ((from-center (geometric-center-of-unfrozen-atoms atom-table coordinates unfrozen))
         (offset (geom:v- to-center from-center)))
    (loop for coord-index from 0 below (length unfrozen) by 3
          for bit-index from 0
          do (when (= 1 (elt unfrozen bit-index))
               (let ((xpos (+ (geom:get-x offset) (elt coordinates coord-index)))
                     (ypos (+ (geom:get-y offset) (elt coordinates (+ 1 coord-index))))
                     (zpos (+ (geom:get-z offset) (elt coordinates (+ 2 coord-index)))))
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
  (declare (ignore mono))
  (let ((vcenter (chem:get-position center))
        #+(or)(vmono (chem:get-position mono))
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
    (let (#+(or)(temp-index nil)
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
             (declare (ignore why))
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

(defun sketch2d-dynamics (sketch dynamics &key accumulate-coordinates (unfreeze t) unfrozen)
  (declare (ignore accumulate-coordinates))
  (let* ((sketch-function (dynamics:scoring-function dynamics))
         (atom-table (chem:node-table sketch-function))
         (energy-stretch (chem:get-stretch-component sketch-function))
         (energy-sketch-nonbond (chem:get-sketch-nonbond-component sketch-function)))
    (randomize-coordinates (dynamics:coordinates dynamics) :from-zero t :unfrozen unfrozen)
    (chem:disable (chem:get-point-to-line-restraint-component sketch-function))
    (apply #'chem:setf-velocity-scale sketch-function *stage1-flatten-force-components*)
    ;; stage 1  step 0 - 4999
    (progn
      (prepare-stage1-sketch-function sketch-function sketch)
      (chem:set-scale-sketch-nonbond energy-sketch-nonbond +stage1-nonbond-constant+)
      (when (verbose sketch) (describe-terms sketch-function "stage1"))
      (dotimes (i 5000) (advance-simulation dynamics :unfrozen unfrozen)))
    ;; stage 2  step 5000 - 9999
    (progn
      (chem:set-scale-sketch-nonbond energy-sketch-nonbond +stage2-nonbond-constant+)
      (dotimes (i 5000) (advance-simulation dynamics :unfrozen unfrozen)))
    ;; stage 3  step 10000-14999
    ;; Add the flattening force
    (progn
      (dynamics:write-coordinates-back-to-matter dynamics)
      (let* ((coordinates (dynamics:coordinates dynamics)))
        (chem:save-coordinates-from-vector sketch-function coordinates)
        (leap.align-axes:tool-orient-principle-axis-along-coordinate-axis (molecule sketch)) ; try aligning the shortest axis around Z
        (chem:load-coordinates-into-vector sketch-function coordinates))
      (apply #'chem:setf-velocity-scale sketch-function *stage3-flatten-force-components*)
      (add-flatten-function sketch-function (chem:get-graph sketch-function) *stage3-flatten-force*)
      (if (and unfrozen unfreeze)
          (let ((center (geometric-center-of-unfrozen-atoms atom-table (dynamics:coordinates dynamics) unfrozen)))
            ;; Four epochs of unfreezing frozen atoms connected to unfrozen atoms
            (loop for stage from 0 below 4
                  do (dotimes (i 1000) (advance-simulation dynamics :unfrozen unfrozen))
                  do (setf unfrozen (unfreeze atom-table (dynamics:coordinates dynamics) unfrozen)))
            (dotimes (i 1000) (advance-simulation dynamics :unfrozen unfrozen))
            (transform-unfrozen-atoms-center-to center atom-table (dynamics:coordinates dynamics) unfrozen))
          (dotimes (i 5000) (advance-simulation dynamics :unfrozen unfrozen))))
    ;; stage 4  step 15000-19999
    ;; Double the bond spring force
    ;; Drop the bond lengths to 1.5
    ;; Add the angle terms
    ;; Progressively reduce the nonbond scale
    (progn
      (chem:walk-sketch-stretch-terms
       energy-stretch
       (lambda (index atom1-I1 atom2-I2 kb r0)
         (declare (ignorable atom1-I1 atom2-I2 kb r0))
         (let ((new-r0 (* r0 +stage4-stretch-bond-length-scale-factor+)))
           #+(or)(format t "r0 = ~a  new-r0 = ~a~%" r0 new-r0)
           (chem:modify-sketch-stretch-term-kb energy-stretch index +stage4-bond-force+)
           (chem:modify-sketch-stretch-term-r0 energy-stretch index new-r0))))
      (add-angle-terms sketch sketch-function *stage4-angle-stretch-force* +stage4-bond-length+)
      (when (verbose sketch) (describe-terms sketch-function "stage4"))
      (loop for i from 0 below 5000
            for cutoff from 80.0 downto 1.0 by (/ 79.0 5000.0)
            do (chem:set-long-distance-cutoff energy-sketch-nonbond (* cutoff +stage4-bond-length+))
               (advance-simulation dynamics :unfrozen (unless unfreeze unfrozen)))
      )
    ;; stage 5 step 20000-24999
    (progn
      (chem:enable (chem:get-point-to-line-restraint-component sketch-function))
      (chem:set-freeze-flags energy-sketch-nonbond 1)
      (when (verbose sketch) (format t "sketch2d-dynamics unfreeze: ~a unfrozen: ~a~%" unfreeze unfrozen))
      (loop for i from 20000 below 25000
            do (advance-simulation dynamics :unfrozen (unless unfreeze unfrozen))))
    (dynamics:write-coordinates-back-to-matter dynamics)
    dynamics
    ))

(defun sketch2d-molecule (molecule &key accumulate-coordinates verbose)
  (let* ((sketch (setup-simulation molecule :accumulate-coordinates accumulate-coordinates :verbose verbose))
         (dynamics (dynamics sketch)))
    (sketch2d-dynamics sketch dynamics :accumulate-coordinates accumulate-coordinates)
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
      (setf (geom:at new-eigen-vectors 0 col0) (geom:at eigen-vector-matrix 0 0)
            (geom:at new-eigen-vectors 1 col0) (geom:at eigen-vector-matrix 1 0)
            (geom:at new-eigen-vectors 2 col0) (geom:at eigen-vector-matrix 2 0)
            (geom:at new-eigen-vectors 0 col1) (geom:at eigen-vector-matrix 0 1)
            (geom:at new-eigen-vectors 1 col1) (geom:at eigen-vector-matrix 1 1)
            (geom:at new-eigen-vectors 2 col1) (geom:at eigen-vector-matrix 2 1)
            (geom:at new-eigen-vectors 0 col2) (geom:at eigen-vector-matrix 0 2)
            (geom:at new-eigen-vectors 1 col2) (geom:at eigen-vector-matrix 1 2)
            (geom:at new-eigen-vectors 2 col2) (geom:at eigen-vector-matrix 2 2))
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


(defgeneric do-sketch2d (thing &key accumulate-coordinates max-tries verbose fixup))

(defmethod do-sketch2d (molecule &key accumulate-coordinates (max-tries *max-tries*) verbose fixup )
  (let ((count 1)
        problems
        unfrozen
        sketch)
    (setf sketch (sketch2d-molecule molecule :accumulate-coordinates accumulate-coordinates :verbose verbose))
    (unless fixup (return-from do-sketch2d sketch))
    (loop named problem-loop
          for try-index below max-tries
          do (multiple-value-setq (problems unfrozen)
               (identify-problem-areas sketch))
          do (progn
               (when (= problems 0)
                 (return-from do-sketch2d sketch))
               (incf count)
               (progn
                 (format t "There were ~a sketch2d problems for try ~a~%" problems count)
                 (finish-output))
               (sketch2d-dynamics sketch (dynamics sketch)
                                  :accumulate-coordinates accumulate-coordinates
                                  :unfreeze nil
                                  :unfrozen unfrozen)))
    (values sketch count)))

(defun insert (item lst compare key)
  (let ((item-key (funcall key item)))
    (if (null lst)
        (list item)
        (if (funcall compare item-key (funcall key (car lst)))
            (cons item lst)
            (cons (car lst) (insert item (cdr lst) compare key))))))

(defun insertion-sort (lst compare key)
  (if (null lst)
    lst
    (insert (car lst) (insertion-sort (cdr lst) compare key) compare key)))


(defmethod do-sketch2d ((aggregate chem:aggregate) &key accumulate-coordinates (max-tries *max-tries*) verbose fixup)
  (if (= (chem:content-size aggregate) 1)
      (do-sketch2d (chem:content-at aggregate 0) :accumulate-coordinates accumulate-coordinates :max-tries max-tries :verbose verbose
        :fixup fixup)
      (error "sketch2d only accepts a molecule or an aggregate with a single molecule")))

#+debug-sketch2d
(defclass sketch-debug-info ()
  ((label :initarg :label :accessor label)))

(defclass chiral-bond-info ()
  ((bond-type :initarg :bond-type :accessor bond-type)
   (bond-atom :initarg :bond-atom :accessor bond-atom)))

(defclass chiral-info ()
  ((chiral-sketch-atom :initarg :chiral-sketch-atom :accessor chiral-sketch-atom)
   (neighbors :initarg :neighbors :accessor neighbors)
   (chiral-bonds :initarg :chiral-bonds :accessor chiral-bonds)))

(defmethod print-object ((obj chiral-info) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a -> ~a" (chiral-sketch-atom obj) (chiral-bonds obj))))

(defclass atom-order ()
  ((original-atom :initarg :original-atom :accessor original-atom)
   (sketch-atom :initarg :sketch-atom :accessor sketch-atom)
   (priority :initarg :priority :accessor priority)))

(defmethod print-object ((obj atom-order) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (original-atom obj))
    (when (slot-boundp obj 'priority)
      (format stream " :priority ~a" (priority obj)))))


(defun select-best-chiral-neighbor (sketch-chiral-atom sketch-neighbors chiral-atoms)
  (cond
    ((multiple-value-bind (matches match)
         (chem:chem-info-match *bis-amino-acid-alpha* sketch-chiral-atom)
       (if matches
           (chem:tag match 1)
           nil)))
    (t (loop for neighbor in (reverse sketch-neighbors)
             for atm = (sketch-atom neighbor)
             for orig-atm = (original-atom neighbor)
             when (not (member orig-atm chiral-atoms))
               do (return-from select-best-chiral-neighbor atm))
       (sketch-atom (car sketch-neighbors)))))

(defun augment-sketch-with-stereochemistry (use-structure-configurations sketch2d cip stereochemistry-types configurations)
  ;; Gather the chiral atoms
  (let ((chiral-atoms (let (chirals)
                        (chem:do-atoms (atm (original-molecule sketch2d))
                          (let ((stereo-type (gethash atm stereochemistry-types)))
                            (when (eq stereo-type :chiral) (push atm chirals)))
                          chirals)))
        (original-to-sketch (make-hash-table))
        )
    ;; calculate original atom to sketch atom map
    (chem:do-atoms (sketch-atom (molecule sketch2d))
      (let ((original-atom (gethash sketch-atom (sketch-atoms-to-original sketch2d))))
        (setf (gethash original-atom original-to-sketch) sketch-atom)))
    ;;
    (loop for chiral-atom in chiral-atoms
          for chiral-sketch-atom = (gethash chiral-atom original-to-sketch)
          for bonds = (chem:bonds-as-list chiral-atom)
          for neighbors = (chem:atom/bonded-atoms-as-list chiral-atom)
          do (let ((neighbors (loop for bond in bonds
                                    for other-atom = (chem:bond/get-other-atom bond chiral-atom)
                                    collect other-atom))
                   (config (gethash chiral-atom configurations)))
               (cond
                 ((member config '(:r :s))
                  (let* ((sketch-chiral (gethash chiral-atom original-to-sketch))
                         (sketch-neighbors (loop for neighbor in neighbors
                                                 for sketch-neighbor = (gethash neighbor original-to-sketch)
                                                 for priority = (chem:get-relative-priority neighbor cip)
                                                 when sketch-neighbor
                                                   collect (make-instance 'atom-order
                                                                          :original-atom neighbor
                                                                          :sketch-atom sketch-neighbor
                                                                          :priority priority))))
                    (setf sketch-neighbors (sort sketch-neighbors #'> :key #'priority))
                    (let ((chiral-info (make-instance 'chiral-info
                                                      :chiral-sketch-atom sketch-chiral
                                                      :neighbors sketch-neighbors)))
;;; svg is a left handed coordinate system
                      (case (length (neighbors chiral-info))
                        (3 (let* ((atm-1 (sketch-atom (car sketch-neighbors)))
                                  (atm-ch sketch-chiral)
                                  (atm-2 (sketch-atom (cadr sketch-neighbors)))
                                  (vec1 (geom:v- (chem:get-position atm-1)
                                                 (chem:get-position atm-ch)))
                                  (vec2 (geom:v- (chem:get-position atm-2)
                                                 (chem:get-position atm-ch)))
                                  (cross (geom:vcross vec1 vec2))
                                  (draw-config (if (> (geom:get-z cross) 0.0) :r :s))
                                  (bond-type (if (eq draw-config config)
                                                 :single-wedge-begin
                                                 :single-hash-begin)))
                             #+(or)(format t "atm-ch: ~a  atm-1: ~a  atm-2: ~a (vz cross): ~,1f draw-config: ~a config: ~a~%"
                                     (chem:get-name atm-ch)
                                     (chem:get-name atm-1)
                                     (chem:get-name atm-2)
                                     (geom:get-z cross)
                                     draw-config
                                     config)
                             (setf (chiral-bonds chiral-info)
                                   (list (make-instance 'chiral-bond-info
                                                        :bond-type bond-type
                                                        :bond-atom (select-best-chiral-neighbor sketch-chiral sketch-neighbors chiral-atoms))))))
                        (4 (let* ((atm-ch sketch-chiral)
                                  #+(or)(atm-1 (sketch-atom (first sketch-neighbors)))
                                  #+(or)(atm-2 (sketch-atom (second sketch-neighbors)))
                                  #+(or)(vec1 (geom:v- (chem:get-position atm-1)
                                                 (chem:get-position atm-ch)))
                                  #+(or)(vec2 (geom:v- (chem:get-position atm-2)
                                                 (chem:get-position atm-ch)))
                                  #+(or)(cross12 (geom:vcross vec1 vec2))
                                  #+(or)(cross12z (geom:get-z cross12))
                                  (atm-3 (sketch-atom (third sketch-neighbors)))
                                  (atm-4 (sketch-atom (fourth sketch-neighbors)))
                                  (vec3 (geom:v- (chem:get-position atm-3)
                                                 (chem:get-position atm-ch)))
                                  (vec4 (geom:v- (chem:get-position atm-4)
                                                 (chem:get-position atm-ch)))
                                  (cross34 (geom:vcross vec3 vec4))
                                  (cross34z (geom:get-z cross34))
                                  (draw-config (if (> cross34z 0.0) :r :s))) ; Change this to flip the rendered config of the center
                             #+(or)(format t "Q atm-ch: ~4a  atm-1: ~3a  atm-2: ~3a cross12z: ~,1f~%"
                                     (chem:get-name atm-ch)
                                     (chem:get-name atm-1)
                                     (chem:get-name atm-2)
                                     cross12z)
                             #+(or)(format t "Q atom-ch: ~4a       atm-3: ~3a  atm-4: ~3a cross34z: ~,1f~%"
                                     (chem:get-name atm-ch)
                                     (chem:get-name atm-3)
                                     (chem:get-name atm-4)
                                     cross34z)
                             (multiple-value-bind (bond-type1 bond-type2)
                                 (if (eq draw-config config)
                                     (values :single-wedge-begin :single-hash-begin)
                                     (values :single-hash-begin :single-wedge-begin))
                               (setf (chiral-bonds chiral-info)
                                     (list (make-instance 'chiral-bond-info
                                                          :bond-type bond-type1
                                                          :bond-atom (sketch-atom (first sketch-neighbors)))
                                           (make-instance 'chiral-bond-info
                                                          :bond-type bond-type2
                                                          :bond-atom (sketch-atom (second sketch-neighbors)))))))))
                      (push chiral-info (chiral-infos sketch2d)))))
                 ((member config '(:left-handed :right-handed))
                  (let* ((sketch-chiral-atom (gethash chiral-atom original-to-sketch))
                         (sketch-neighbors (loop for neighbor in neighbors
                                                 for sketch-neighbor = (gethash neighbor original-to-sketch)
                                                 when sketch-neighbor
                                                   collect (make-instance 'atom-order
                                                                          :original-atom neighbor
                                                                          :sketch-atom sketch-neighbor))))
                    (let ((chiral-info (make-instance 'chiral-info
                                                      :chiral-sketch-atom sketch-chiral-atom
                                                      :neighbors sketch-neighbors)))
;;; svg is a left handed coordinate system
                      ;;                      (case (length (neighbors chiral-info))
                      (let* ((atm-1 (sketch-atom (first sketch-neighbors)))
                             (atm-ch sketch-chiral-atom)
                             (atm-2 (sketch-atom (second sketch-neighbors)))
                             (vec1 (geom:v- (chem:get-position atm-1)
                                            (chem:get-position atm-ch)))
                             (vec2 (geom:v- (chem:get-position atm-2)
                                            (chem:get-position atm-ch)))
                             (cross (geom:vcross vec1 vec2))
                             (draw-config (if (> (geom:get-z cross) 0.0) :right-handed :left-handed)))
                        (multiple-value-bind (bond-type1 bond-type2)
                            (if (eq draw-config config)
                                (values :single-wedge-begin :single-hash-begin)
                                (values :single-hash-begin :single-wedge-begin))
                          (setf (chiral-bonds chiral-info)
                                (list* (make-instance 'chiral-bond-info
                                                      :bond-type bond-type1
                                                      :bond-atom (sketch-atom (third sketch-neighbors)))
                                       (if (> (length (neighbors chiral-info)) 3)
                                           (list (make-instance 'chiral-bond-info
                                                                :bond-type bond-type2
                                                                :bond-atom (sketch-atom (fourth sketch-neighbors))))
                                           nil)))))
                      (push chiral-info (chiral-infos sketch2d)))))
                 (t (warn "Handle config ~a for atom ~a" config chiral-atom)))))))


(defgeneric sketch2d (molecule &key accumulate-coordinates verbose fixup show-names&allow-other-keys))

(defmethod sketch2d ((obj topology:topology) &key accumulate-coordinates verbose fixup show-names &allow-other-keys)
  (let* ((*show-names* (or *show-names* show-names))
         (molecule (topology:build-one-molecule-for-topology obj))
         (sketch2d (do-sketch2d molecule :accumulate-coordinates accumulate-coordinates :verbose verbose :fixup fixup))
         )
    (multiple-value-bind (stereochemistry-types configurations cips)
        (chem:calculate-stereochemistry molecule :use-structure t)
      (augment-sketch-with-stereochemistry nil sketch2d cips stereochemistry-types configurations)
      (let* ((plugs (topology:plugs obj))
             (annotations (loop for plug in (alexandria:hash-table-values plugs)
                                for plug-name = (topology:name plug)
                                for plug-bonds = (topology:plug-bonds plug)
                                for plug-atom-names = (loop for plug-bond across plug-bonds
                                                            collect (topology:atom-name plug-bond))
                                collect (make-instance 'annotation
                                                       :text (format nil "~s" plug-name)
                                                       :atom-names plug-atom-names))))
        (reinitialize-instance sketch2d :annotations annotations)
    sketch2d))))


(defmethod sketch2d ((molecule chem:molecule) &key accumulate-coordinates verbose fixup use-structure show-names &allow-other-keys)
  "Calculate a sketch2d from the molecule.  If use-structure is set then calculate the stereochemical configuration of atoms
using the positions of the atoms in the molecule, otherwise get the configuration from the _Configuration slot."
  (declare (optimize (debug 3)))
  (let ((*show-names* (or *show-names* show-names))
        (sketch2d (do-sketch2d molecule :accumulate-coordinates accumulate-coordinates :verbose verbose :fixup fixup)))
    (multiple-value-bind (stereochemistry-types configurations cips)
        (chem:calculate-stereochemistry molecule :use-structure use-structure)
      (augment-sketch-with-stereochemistry use-structure sketch2d cips stereochemistry-types configurations)
      sketch2d)))

(defgeneric similar-sketch2d (thing anchor-sketch2d &key accumulate-coordinates atom-match-callback))

(defmethod similar-sketch2d ((molecule chem:molecule) anchor-sketch2d &key use-structure accumulate-coordinates (atom-match-callback #'molecule-graph:element-match))
  "Generates a sketch of molecule that shares a common substructure
with the molecule in anchor-sketch2d
and aligns the new sketch the same way.
The use-structure keyword argument is treated the same as the sketch2d function does.
If atom-match-callback is NIL then all atoms can match each other.
Otherwise pass a function that takes two atoms and returns T if they are matchable."
  (check-type molecule chem:molecule)
  (check-type anchor-sketch2d sketch2d)
  (let ((anchor-molecule (original-molecule anchor-sketch2d)))
    (multiple-value-bind (stereochemistry-types configurations cips)
        (chem:calculate-stereochemistry molecule :use-structure use-structure)
      #+(or)(format t "About to compare molecules~%")
      (multiple-value-bind (equiv diff1 diff2)
          (molecule-graph.max-clique:compare-molecules-return-largest-connected-cluster-of-equivalent-atoms
           molecule anchor-molecule
           :atom-match-callback atom-match-callback
           :topological-constraint-theta 2
           :remove-partial-rings t)
        (declare (ignore diff1 diff2))
        #+(or)(progn
                (format t "mol1: ~s~%" molecule)
                (format t "mol2: ~s~%" anchor-molecule)
                (format t "similar-sketch2d equiv atoms: ~a~%" equiv)
                (format t "diff1: ~s~%" diff1)
                (format t "diff2: ~s~%" diff2))
        (let ((new-sketch (setup-simulation molecule :accumulate-coordinates accumulate-coordinates)))
          ;; We need hash-tables for the original atoms to the sketch atoms
          (let ((anchor-original-to-sketch (make-hash-table))
                (new-original-to-sketch (make-hash-table)))
            (maphash (lambda (key value)
                       (setf (gethash value new-original-to-sketch) key))
                     (sketch-atoms-to-original new-sketch))
            (maphash (lambda (key value)
                       (setf (gethash value anchor-original-to-sketch) key))
                     (sketch-atoms-to-original anchor-sketch2d))
            (let* ((dynamics (dynamics new-sketch))
                   (energy-function (dynamics:scoring-function dynamics))
                   (atom-table (chem:atom-table energy-function))
                   (number-of-atoms (chem:get-number-of-atoms atom-table))
                   (debug-info (make-hash-table))
                   (unfrozen (make-array number-of-atoms :element-type 'bit :initial-element 1)))
              #+(or)(format t "There are ~a atoms in atom-table~%" number-of-atoms)
              (loop for equiv-pair in equiv
                    for new-atom = (car equiv-pair)
                    for anchor-atom = (cdr equiv-pair)
                    for new-sketch-atom = (gethash new-atom new-original-to-sketch)
                    for anchor-sketch-atom = (gethash anchor-atom anchor-original-to-sketch)
                    for anchor-pos = (chem:get-position anchor-sketch-atom)
                    do (unless (or (= 1 (chem:get-atomic-number anchor-sketch-atom))
                                   (eq :lp (chem:get-element anchor-sketch-atom)))
                         (let ((new-sketch-atom-index (chem:get-atom-index-or-nil atom-table new-sketch-atom)))
                           ;; Some atoms from the original molecule are not in the sketch molecule
                           ;; so if it's missing - then ignore it.
                           (when new-sketch-atom-index
                             #+debug-sketch2d
                             (setf (gethash new-sketch-atom debug-info)
                                   (make-instance 'sketch-debug-info
                                                  :label (format nil "~a/~a" (chem:get-name new-atom) (chem:get-name anchor-atom))))
                             (chem:set-position new-sketch-atom anchor-pos)
                             (setf (aref unfrozen new-sketch-atom-index) 0)))))
              #+(or)(format t "Using unfrozen: ~a  length: ~a~%" frozen (length unfrozen))
              (chem:load-coordinates-into-vector energy-function (dynamics:coordinates dynamics))
              (sketch2d-dynamics new-sketch dynamics :unfrozen unfrozen :unfreeze nil
                                                     :accumulate-coordinates accumulate-coordinates)
              (augment-sketch-with-stereochemistry use-structure new-sketch cips stereochemistry-types configurations)
              #+debug-sketch2d(setf (debug-info new-sketch) debug-info)
              new-sketch)))))))

(defmethod similar-sketch2d ((aggregate chem:aggregate) anchor-sketch2d
                             &key accumulate-coordinates (atom-match-callback #'molecule-graph:element-match))
  (unless (= (chem:content-size aggregate) 1)
    (error "There can only be one molecule passed to similar-sketch2d - you passed an aggregate containing ~a molecules" (chem:content-size aggregate)))
  (similar-sketch2d (chem:content-at aggregate 0) anchor-sketch2d
                    :accumulate-coordinates accumulate-coordinates
                    :atom-match-callback atom-match-callback))
