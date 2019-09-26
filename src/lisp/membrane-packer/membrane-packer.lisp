(in-package :membrane-packer)

#| Working with hex-grids is tricky.
This web page is excellent: https://www.redblobgames.com/grids/hexagons

This code uses Offset coordinates odd-r
|#
(defparameter *x-expand* 60.0)
(defparameter *y-expand* 60.0)
(defparameter *z-expand* 10.0)

(defparameter *z-space* 1.5)
(defparameter *mutate-lipid-fraction* 0.08)
(defparameter *mutate-position-fraction* 0.3)
(defparameter *lipid-density* 0.0098)
(defparameter *random-tip-angle-degrees* 10.0)
(defparameter *close-distance* 1.2)


(defparameter *oddr-directions* #( #( #( -1 -1)
                                     #( -1  0)
                                     #( -1 +1)
                                     #(  0 +1)
                                     #( +1  0)
                                     #(  0 -1))
                                  #( #(  0 -1)
                                    #( -1  0)
                                    #(  0 +1)
                                    #( +1 +1)
                                    #( +1  0)
                                    #( +1 -1)))
  "Directions for oddr hex-grid coordinate system")

(defun shallow-copy-object (original)
  (let* ((class (class-of original))
         (copy (allocate-instance class)))
    (dolist (slot (mapcar #'clos:slot-definition-name (clos:class-slots class)))
      (when (slot-boundp original slot)
        (setf (slot-value copy slot)
              (slot-value original slot))))
    copy))

(defclass order-info ()
  ((order-residue :initarg :order-residue :accessor order-residue)
   (order-atom :initarg :order-atom :accessor order-atom)))
(cando:make-class-save-load order-info)

(defclass lipid-info ()
  ((lipid-molecule :initarg :lipid-molecule :accessor lipid-molecule :type vector
                   :documentation "A single lipid-molecule")
   (order-info :initarg :order-info :accessor order-info :type order-info
          :documentation "The ordering of the atoms. This is a vector order-info objects.
It matches the order of the coordinates.")
   (coordinates :initarg :coordinates :accessor coordinates :type vector
                :documentation "A vector of the coordinates. A vector of vectors of doubles.
The inner vectors are the coordinates of the atoms defined by the order with the x,y,z coordinates
all layed out in a linear vector.")
   (center-height :initarg :center-height :accessor center-height
                  :documentation "The average height of the center of geometry of the lipid")
   (density :initarg :density :accessor density
            :documentation "Density of lipids per square angstrom")))

(cando:make-class-save-load lipid-info)

(defclass optimized-lipid-info (lipid-info)
  ((octrees :initarg :octrees :accessor octrees)))

(defclass ga-thing ()
  ((transform :initform (geom:make-m4-translate (geom:vec 0.0 0.0 0.0))
              :initarg :transform :accessor transform)
   (inverse-transform :initform (geom:make-m4-translate (geom:vec 0.0 0.0 0.0))
                      :initarg :inverse-transform :accessor inverse-transform)
   (adjacent-to-solute :initform nil :initarg :adjacent-to-solute :accessor adjacent-to-solute)))

(defclass ga-solute (ga-thing)
  ((solute :initarg :solute :accessor solute)
   (octree :initarg :octree :accessor octree)))

(defclass ga-lipid (ga-thing)
  ((array-indices :initarg :array-indices :accessor array-indices)
   (lipid-info :initarg :lipid-info :accessor lipid-info)
   (lipid-index :initarg :lipid-index :accessor lipid-index)
   (lipid-id :initarg :lipid-id :accessor lipid-id)
   ))

(defmethod print-object ((obj ga-lipid) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream ":id ~a :lipid-index ~a" (lipid-id obj) (lipid-index obj))))


(defgeneric shallow-copy-ga-object (original))

(defmethod shallow-copy-ga-object ((original ga-solute))
  ;; don't need to copy this
  original)

(defmethod shallow-copy-ga-object ((original ga-lipid))
  (make-instance 'ga-lipid
                 :array-indices (array-indices original)
                 :lipid-id (lipid-id original)
                 :lipid-info (lipid-info original)
                 :lipid-index (lipid-index original)
                 :transform (transform original)
                 :inverse-transform (inverse-transform original)
                 :adjacent-to-solute (adjacent-to-solute original)
                 ))



(defclass range ()
  ((source-ga-thing :initarg :source-ga-thing :accessor source-ga-thing)
   (target-ga-thing :initarg :target-ga-thing :accessor target-ga-thing)
   (start :initarg :start :accessor start)
   (end :initarg :end :accessor end)))

(defmethod print-object ((obj range) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a" (source-ga-thing obj) (target-ga-thing obj))))

(defclass collisions ()
  ((indices :initarg :indices :accessor indices)
   (ranges :initarg :ranges :accessor ranges)))

(defun make-optimized-lipid-selector (lipid-selector)
  (loop for prob-lipid-info in lipid-selector
        for probability = (car prob-lipid-info)
        for lipid-info = (shallow-copy-object (cdr prob-lipid-info))
        for octrees = (make-array (length (coordinates lipid-info)))
        do (loop for coordinates-index from 0 below (length (coordinates lipid-info))
                 for coordinates = (elt (coordinates lipid-info) coordinates-index)
                 for bounding-cuboid = (chem:nvector-bounding-cuboid coordinates)
                 for octree = (chem:make-generic-octree
                               (geom:get-center bounding-cuboid)
                               ;; Pad each octree by enough space so that points that match outermost points in octree will
                               ;; be within the octree boundaries
                               (geom:v+ (geom:v* (geom:vec 1.0 1.0 1.0) *close-distance*)
                                        (geom:get-half-width bounding-cuboid)))
                 do (loop for index below (length coordinates) by 3
                          with pos = (geom:vec 0.0 0.0 0.0)
                          do (geom:vec-extract pos coordinates index)
                          do (chem:generic-octree-insert octree pos index))
                 do (setf (elt octrees coordinates-index) octree))
        collect (cons probability (change-class lipid-info 'optimized-lipid-info
                                                :octrees octrees))))

(defun carve-membrane (membrane center-of-geometry-test)
  "Return an aggregate that has removed lipids from a membrane that do not pass the center-of-geometry-test"
  (let ((membrane-copy (chem:matter-copy membrane)))
    (let ((remove-molecule-list nil))
      (cando:do-molecules (mol membrane-copy)
        (let ((center-of-geometry (chem:geometric-center mol)))
          (unless (funcall center-of-geometry-test center-of-geometry)
            (push mol remove-molecule-list))))
      (loop for mol in remove-molecule-list
            do (chem:remove-molecule membrane-copy mol))
      membrane-copy)))

(defconstant +max-bins+ 30)
(defun analyze-membrane (membrane x-shrink y-shrink)
  (let ((lipid-vec (chem:map-molecules 'vector #'identity membrane)))
    (unless (> (length lipid-vec) 0)
      (error "You must provide a list of lipids from a bilayer"))
    ;; First analyze displacement from Z=0 plane
    (let ((average-height (/ (loop for ii from 0 below (length lipid-vec)
                                   for mol = (elt lipid-vec ii)
                                   for center = (chem:geometric-center mol)
                                   for delta = (if (> (float-sign (geom:vz center)) 0.0)
                                                   (geom:vz center)
                                                   (- (geom:vz center)))
                                   sum delta)
                             (length lipid-vec))))
      ;; Now calculate the shortest center-to-center distanc
      ;; in the x-y plane
      (let* ((bounding-box (chem:bounding-box membrane))
             (x-width (- (chem:get-x-width bounding-box) x-shrink))
             (y-width (- (chem:get-y-width bounding-box) y-shrink))
             (area (* 2.0 x-width y-width))
             (density (/ (length lipid-vec) area)))
        (values average-height density)))))

(defun order-atoms-in-molecule (molecule)
  (let ((order-info (make-array 256 :adjustable t :fill-pointer 0)))
    (cando:do-residues (residue molecule)
      (chem:set-name residue (intern (string (chem:get-name residue)) :keyword)) ; rename residues
      (cando:do-atoms (atom residue)
        (vector-push-extend (make-instance 'order-info
                                           :order-residue residue
                                           :order-atom atom) order-info)))
    order-info))

(defun orders-match-p (order1 order2)
  (unless (= (length order1) (length order2))
    (error "The two orders have different numbers of atoms"))
  (loop for index from 0 below (length order1)
        for entry1 = (elt order1 index)
        for entry2 = (elt order2 index)
        for residue1 = (order-residue entry1)
        for atom1 = (order-atom entry1)
        for residue2 = (order-residue entry2)
        for atom2 = (order-atom entry2)
        always (and (eq (chem:get-name residue1) (chem:get-name residue2))
                    (eq (chem:get-name atom1) (chem:get-name atom2)))))

(defun extract-coordinates (residue-atom-order)
  (let ((coordinates (make-array (* 3 (length residue-atom-order)) :element-type 'double-float)))
    (loop for index from 0 below (length residue-atom-order)
          for order-info = (elt residue-atom-order index)
          for residue = (order-residue order-info)
          for atom = (order-atom order-info)
          for pos = (chem:get-position atom)
          for coord-index = (* 3 index)
          do (setf (elt coordinates coord-index) (geom:vx pos)
                   (elt coordinates (+ 1 coord-index)) (geom:vy pos)
                   (elt coordinates (+ 2 coord-index)) (geom:vz pos)))
    coordinates))
  
(defun canonicalize-lipids (membrane top-atom-names bottom-atom-names)
  (let ((top (geom:vec 0.0 0.0 0.0))
        (bottom (geom:vec 0.0 0.0 0.0))
        (molecule-list (chem:map-molecules 'list #'identity membrane)))
    ;; First make sure all residues/atoms are in the same order
    (let* ((lipid-molecule (first molecule-list))
           (order-info (order-atoms-in-molecule lipid-molecule)))
      (unless (loop for molecule in (rest molecule-list)
                    always (orders-match-p order-info (order-atoms-in-molecule molecule)))
        (error "The orders of all residues/atoms in all molecules is not identical"))
      (let ((aligned-molecules nil))
        (loop for cmol in molecule-list
              for mol = (chem:matter-copy cmol)
              for cog = (chem:geometric-center mol)
              for m4 = (geom:make-m4-translate (geom:v* cog -1.0))
              do (chem:apply-transform-to-atoms mol m4)
              do (push mol aligned-molecules))
        (flet ((average-atom-pos (mol atom-names)
                 (loop named top-calc
                       with ttop = (geom:vec 0.0 0.0 0.0)
                       for name in atom-names
                       for atom = (chem:first-atom-with-name mol name)
                       do (setf ttop (geom:v+ ttop (chem:get-position atom)))
                       finally (return-from top-calc (geom:v* ttop (/ 1.0 (length atom-names)))))))
          (loop with zdir = (geom:vec 0.0 0.0 1.0)
                for mol in aligned-molecules
                for top = (average-atom-pos mol top-atom-names)
                for bot = (average-atom-pos mol bottom-atom-names)
                for delta = (geom:vnormalized (geom:v- top bot))
                for dot = (geom:vdot delta zdir)
                for angle = (acos dot)
                for cross = (geom:vcross delta zdir)
                for cross-length = (geom:vlength cross)
                for cross-dir = (if (or (< cross-length 0.035) ; 2 deg
                                        (> cross-length 3.107)) ; 178deg
                                    (geom:vec 1.0 0.0 0.0)
                                    (geom:vnormalized cross))
                for transform = (geom:make-m4-rotate-axis (- angle) cross)
                do (chem:apply-transform-to-atoms mol transform)))
        ;; extract the coordinates from the aligned molecules
        (let ((coordinates (make-array (length aligned-molecules))))
          (loop for index from 0 below (length aligned-molecules)
                for aligned-moleculei = (elt aligned-molecules index)
                for orderi = (order-atoms-in-molecule aligned-moleculei)
                for coordsi = (extract-coordinates orderi)
                do (setf (elt coordinates index) coordsi))
          ;; Analyze the lipid height and density of the membrane
          (multiple-value-bind (height density)
              (analyze-membrane membrane -2.0 -2.0)
            (make-instance 'lipid-info
                           :lipid-molecule lipid-molecule
                           :order-info order-info
                           :coordinates coordinates
                           :center-height height
                           :density density)))))))

(defparameter *canonicalized-popc* #P"cando:lisp;membrane-packer;data;popc.cando")

(defparameter *popc* nil)

(eval-when (:load-toplevel :execute)
  (setf *popc* (cando:load-cando *canonicalized-popc*)))

  
(defun prepare-popc-lipids ()
  "Prepare the popc lipids"
  (cando-user:setup-amber-paths)
  (cando-user:easy-gaff)
  (cando-user:source "leaprc.protein.ff14SB")
  (cando-user:source "leaprc.lipid14")
  (cando-user:source "leaprc.water.tip3p")
  (cando-user:source "leaprc.gaff")
  (let* ((membrane-protein (cando-user:load-pdb "cando:lisp;membrane-packer;data;4s0v_popc_membrane_hid.pdb"))
         (memb (carve-membrane membrane-protein (lambda (cog) (< (geom:vx cog) -20.0))))
         lipids)
    (leap.set-box:set-box memb :centers 0.0)
    (let ((lipid-info (canonicalize-lipids memb
                                          '(:p31 :n31)
                                          '(:c116 :c118))))
      (cando-user:save-cando lipid-info *canonicalized-popc*))))


(defun prepare-lipids ()
  (prepare-popc-lipids))

(defclass ga-membrane ()
  ((bounding-box :initarg :bounding-box :accessor bounding-box)
   (array :initarg :array :accessor array)
   (lipid-selector :initarg :lipid-selector :accessor lipid-selector)
   (ga-solute :initarg :ga-solute :accessor ga-solute)))

(defun copy-ga-membrane (original)
  (let ((copy (shallow-copy-object original)))
    (setf (array copy) (make-array (array-dimensions (array original))))
    (loop for index from 0 below (array-total-size (array copy))
          for lipid-copy = (shallow-copy-ga-object (row-major-aref (array original) index))
          do (setf (row-major-aref (array copy) index) lipid-copy))
    copy))

(defun wrap-in-bounding-box (xypos bounding-box)
  (let* ((xi (geom:vx xypos))
         (yi (geom:vy xypos))
         (center (chem:get-bounding-box-center bounding-box))
         (half-x-width (/ (chem:get-x-width bounding-box) 2.0))
         (half-y-width (/ (chem:get-y-width bounding-box) 2.0))
         (xmin (- (geom:vx center) half-x-width))
         (ymin (- (geom:vy center) half-y-width))
         (xmax (+ (geom:vx center) half-x-width))
         (ymax (+ (geom:vy center) half-y-width)))
    (setf xi (cond
               ((< xi xmin) (+ xi (chem:get-x-width bounding-box)))
               ((>= xi xmax) (- xi (chem:get-x-width bounding-box)))
               (t xi)))
    (setf yi (cond
               ((< yi ymin) (+ yi (chem:get-y-width bounding-box)))
               ((>= yi xmax) (- yi (chem:get-y-width bounding-box)))
               (t yi)))
    (geom:vec xi yi (geom:vz xypos))))

(defun select-lipid (lipid-selector &optional rnd index)
"A LIPID-SELECTOR is an alist of ( probability . lipid-info ) pairs.
The probabilities are sorted in ascending order.
Example ( (0.2 . foo) (0.5 . bar) (1.0 . baz))
will select lipid foo with 20% chance, lipid bar with 30% chance and lipid baz with 50% chance."
  (let ((rnd (or rnd (random 1.0))))
    (loop for pair in lipid-selector
          for prob = (car pair)
          for lipid-info = (cdr pair)
          when (< rnd prob)
            do (let ((index (or index (random (length (coordinates lipid-info))))))
                 (return-from select-lipid (cons lipid-info index))))))

(defun make-random-ga-lipid (array-indices xypos z-offset lipid-selector lipid-id top-bottom)
  (let* ((lipid-info-index (select-lipid lipid-selector))
         (lipid-info (car lipid-info-index))
         (lipid-index (cdr lipid-info-index))
         (translate (ecase top-bottom
                      (:top (geom:v+ xypos (geom:vec 0.0 0.0 (+ z-offset(center-height lipid-info)))))
                      (:bottom (geom:v+ xypos (geom:vec 0.0 0.0 (- (+ z-offset (center-height lipid-info))))))))
         (position (ecase top-bottom
                     (:top (geom:make-m4-translate translate))
                     (:bottom (geom:m*m (geom:make-m4-translate translate)
                                        (geom:make-m4-rotate-x PI)))))
         (tweakx (geom:make-m4-rotate-x (* (- (random *random-tip-angle-degrees*) (* *random-tip-angle-degrees* 0.5)) 0.0174533)))
         (tweaky (geom:make-m4-rotate-y (* (- (random *random-tip-angle-degrees*) (* *random-tip-angle-degrees* 0.5)) 0.0174533)))
         (tweakz (geom:make-m4-rotate-z (random (* 2.0 PI))))
         (tweak (let* ((m1 (geom:m*m tweaky tweakz))
                       (m2 (geom:m*m tweakx m1)))
                  m2))
         (transform (geom:m*m position tweak)))
    (make-instance 'ga-lipid
                   :array-indices array-indices
                   :lipid-info lipid-info
                   :lipid-index lipid-index
                   :lipid-id lipid-id
                   :transform transform
                   :inverse-transform (geom:invert-transform transform))))

(defun make-molecule-from-ga-lipid (ga-lipid &key debug)
  "Build a molecule from the ga-lipid.  If DEBUG is T then change the molecule name and residue names so we can pick out the lipid-index"
  (let* ((lipid-info (lipid-info ga-lipid))
         (lipid-index (lipid-index ga-lipid))
         (lipid-id (lipid-id ga-lipid))
         (transform (transform ga-lipid))
         (coordinates (elt (coordinates lipid-info) lipid-index))
         (new-to-old-map (make-hash-table))
         (lipid-mol (chem:matter-copy (lipid-molecule lipid-info) new-to-old-map))
         (order-info (order-info lipid-info))
         (old-to-new-map (make-hash-table)))
    (when debug
      (let ((lipid-name (format nil "L~a" lipid-id)))
        (chem:set-name lipid-mol (intern lipid-name :keyword))
        (cando:do-residues (res lipid-mol)
          (chem:set-name res (intern (format nil "~a-~a" lipid-name (chem:get-name res)) :keyword)))))
    ;; Flip the atom map so we can take the order of old atoms and find the new atom
    (maphash (lambda (new-atom old-atom)
               (setf (gethash old-atom old-to-new-map) new-atom))
             new-to-old-map)
    (loop for index from 0 below (length order-info)
          for coord-index = (* 3 index)
          for residue-atom = (elt order-info index)
          for residue = (order-residue residue-atom)
          for old-atom = (order-atom residue-atom)
          for new-atom = (gethash old-atom old-to-new-map)
          for xpos = (elt coordinates coord-index)
          for ypos = (elt coordinates (+ 1 coord-index))
          for zpos = (elt coordinates (+ 2 coord-index))
          for pos = (geom:vec xpos ypos zpos)
          do (chem:set-position new-atom pos))
    (chem:apply-transform-to-atoms lipid-mol transform)
    lipid-mol))

(defclass scoring-cell ()
  ((ga-thing :initform nil :initarg :ga-thing :accessor ga-thing)
   (inverse-transform :initarg :inverse-transform :accessor inverse-transform)
   (close-contact-indices :initform (make-array 16 :element-type 'cons :adjustable t :fill-pointer 0)
                          :accessor close-contact-indices)))

(defclass scoring-membrane ()
  ((solute :initarg :solute :accessor solute)
   (solute-scoring-cell :initarg :solute-scoring-cell :accessor solute-scoring-cell)
   (array :initarg :array :accessor array)))

#+(or)(defun coordinates-from-ga-lipid (ga-lipid scoring-cell)
  "Generate the x,y,z coordinates for the ga-lipid into the coordinates of the scoring-cell.
We don't want to do any consing, so we write the coordinates into the coordinate nvector
in the scoring cell."
  (let* ((lipid-info (lipid-info ga-lipid))
         (lipid-index (lipid-index ga-lipid))
         (transform (transform ga-lipid))
         (coordinates (elt (coordinates lipid-info) lipid-index)))
    (setf (ga-thing scoring-cell) ga-lipid)
    (chem:apply-transform-to-coordinates (coordinates scoring-cell) coordinates transform)
    (let* ((bounding-cuboid (chem:nvector-bounding-cuboid (coordinates scoring-cell) (length coordinates)))
           (octree (chem:make-generic-octree (geom:get-center bounding-cuboid) (geom:get-half-width bounding-cuboid))))
      ;; Fill the octree
      (loop for index below (length coordinates) by 3
            with pos = (geom:vec)
            do (progn
                 (geom:vec-extract pos (coordinates scoring-cell) index)
                 (chem:generic-octree-insert octree pos index)))
      (setf (octree scoring-cell) octree))))

(defun build-scoring-membrane (ga-membrane)
  "Build an object that will be used to quickly score a ga-membrane.
We don't want to do ANY consing when running the genetic algorithm. So we 
allocate enough space for the coordinates for the largest lipid in each cell of the membrane."
  (let* ((dimensions (array-dimensions (array ga-membrane)))
         (scoring-array (make-array dimensions))
         (lipid-selector (lipid-selector ga-membrane))
         (largest-coordinate-size (loop for prob-lipid in lipid-selector
                                              for lipid-info = (cdr prob-lipid)
                                              maximize (length (elt (coordinates lipid-info) 0)))))
    (loop for index from 0 below (array-total-size scoring-array)
          for scoring-cell = (make-instance 'scoring-cell
                                            :coordinates (chem:make-nvector largest-coordinate-size))
          do (setf (row-major-aref scoring-array index) scoring-cell))    (make-instance 'scoring-membrane
                   :array scoring-array
                   :solute-scoring-cell
                   )))

(defun calculate-center-and-half-width (bounding-cuboid)
  (let ((center (geom:vec (* (+ (geom:get-min-x bounding-cuboid) (geom:get-max-x bounding-cuboid)) 0.5)
                          (* (+ (geom:get-min-y bounding-cuboid) (geom:get-max-y bounding-cuboid)) 0.5)
                          (* (+ (geom:get-min-z bounding-cuboid) (geom:get-max-z bounding-cuboid)) 0.5)))
        (half-width (geom:vec (* 0.5 (- (geom:get-max-x bounding-cuboid) (geom:get-min-x bounding-cuboid)))
                              (* 0.5 (- (geom:get-max-y bounding-cuboid) (geom:get-min-y bounding-cuboid)))
                              (* 0.5 (- (geom:get-max-z bounding-cuboid) (geom:get-min-z bounding-cuboid))))))
    (values center half-width)))


(defun lipid-protein-overlaps-p (bounding-box ga-lipid ga-solute &optional (cutoff-squared 1.0))
  (let ((collisions (make-collisions)))
    (do-score-cell bounding-box ga-lipid ga-solute collisions)
    (> (length (indices collisions)) 0)))

(defun wrapped-array-indices (xi yi xdim ydim)
  (values (mod xi xdim) (mod yi ydim)))

(defun wrapped-xypos (array xi yi xstart ystart xdir ypdir ymdir)
  (let* ((dims (array-dimensions array))
         (xdim (elt dims 0))
         (ydim (elt dims 1))
         (zdim (elt dims 2)))
    (let* ((xdelta (geom:v* xdir (mod xi xdim)))
           (ydelta (geom:v+ (geom:v* (geom:v+ ypdir ymdir) (floor (/ yi 2)))
                            (geom:v* ypdir (if (oddp yi) 1 0))))
           (xypos (geom:v+ (geom:vec xstart ystart 0.0) (geom:v+ xdelta ydelta))))
      xypos)))


(defun neighbors (membrane xi yi)
  (let* ((array (array membrane))
         (dims (array-dimensions array))
         (xdim (elt dims 0))
         (ydim (elt dims 1))
         (zdim (elt dims 2)))
    (let* ((yparity (mod yi 2))
           (dirs (elt *oddr-directions* yparity)))
      (loop for dir across dirs
            for xdir = (elt dir 0)
            for ydir = (elt dir 1)
            do (multiple-value-bind (xp yp)
                   (wrapped-array-indices (+ xi xdir) (+ yi ydir) xdim ydim)
                 (format t "~a ~a~%" xp yp))))))
                                              

(defun prepare-solute (solute input-bounding-box)
  (let* ((solute-copy (let ((sol (chem:matter-copy solute)))
                        (cando:move-geometric-center-to-origin sol)
                        #+(or)(cando:move-geometric-center-to-point sol (geom:vec half-x-width half-y-width 0.0))
                        sol))
         (bounding-cuboid (chem:matter-calculate-bounding-cuboid solute-copy 1.0))
         (half-width (geom:get-half-width bounding-cuboid))
         (width (geom:v* half-width 2.0))
         (x-width (+ (geom:vx width) *x-expand*))
         (y-width (+ (geom:vx width) *y-expand*))
         (z-width (+ (geom:vz width) *z-expand*))
         (bounding-box (chem:make-bounding-box (list x-width y-width z-width))))
    (when input-bounding-box
      (setf bounding-box input-bounding-box))
    (values solute bounding-box)))

(defun build-ga-membrane (solute &key input-bounding-box
                                   (lipid-selector (list (cons 1.0 *popc*))))
  (multiple-value-bind (solute-copy suggested-bounding-box)
      (prepare-solute solute input-bounding-box)
    (unless input-bounding-box
      (setf input-bounding-box suggested-bounding-box))
    (let* ((lipid-selector (make-optimized-lipid-selector lipid-selector))
           (nonbond-db (chem:compute-merged-nonbond-force-field-for-aggregate solute))
           (x-width (chem:get-x-width input-bounding-box))
           (y-width (chem:get-y-width input-bounding-box))
           (z-width (chem:get-z-width input-bounding-box))
           (z-offset *z-space*)
           (half-x-width (/ x-width 2.0))
           (half-y-width (/ y-width 2.0))
           (lipid-id 0)
           (bounding-cuboid (chem:matter-calculate-bounding-cuboid solute 1.0)))
      (multiple-value-bind (center half-width)
          (calculate-center-and-half-width bounding-cuboid)
        (let* ((solute-octree (chem:make-generic-octree center half-width))
               (hex-angle (* 0.0174533 60.0)) ; 60 degrees in radians
               (xstep (sqrt (/ 1.0 *lipid-density*)))
               (xdir (geom:v* (geom:vec 1.0 0.0 0.0) xstep))
               (ypdir (geom:v* (geom:vec (cos hex-angle) (sin hex-angle) 0.0) xstep))
               (ymdir (geom:v* (geom:vec (- (cos hex-angle)) (sin hex-angle) 0.0) xstep))
               (ystep (* (sin hex-angle) xstep))
               (xnum (floor (/ x-width xstep)))
               (ynum (let ((simple-ynum (floor (/ y-width ystep))))
                       (if (evenp simple-ynum) simple-ynum (1+ simple-ynum)))) ; need even number of cells in y dir
               (x-width (* xnum xstep))
               (y-width (* ynum ystep))
               (xstart (- (* x-width 0.5)))
               (ystart (- (* y-width 0.5)))
               (lipid-id 0)
               (array (make-array (list xnum ynum 2)))
               (node-array (make-array (list xnum ynum 2)))
               (ga-solute (make-instance 'ga-solute :solute solute :octree solute-octree))
               (bounding-box (chem:make-bounding-box (list x-width y-width z-width)
                                                     :angles-degrees (chem:get-bounding-box-angles-degrees input-bounding-box)
                                                     :center (chem:get-bounding-box-center input-bounding-box))))
          ;; Adjust the xwidth and ywidth of the bounding box
          ;; Populate the solute-octree
          (cando:do-molecules (mol solute)
            (let ((mol-list (list mol)))
              (cando:do-residues (res mol)
                (let ((res-list (list* res mol-list)))
                  (cando:do-atoms (atm res)
                    (let ((atm-list (list* atm res-list)))
                      (chem:generic-octree-insert solute-octree (chem:get-position atm) atm-list)
                      ))))))
          (loop for yindex from 0 below ynum
                do (loop for xindex from 0 below xnum
                         for wrapped-xypos = (wrapped-xypos array xindex yindex xstart ystart xdir ypdir ymdir)
                         for top-ga-lipid = (make-random-ga-lipid (list xindex yindex 0) wrapped-xypos z-offset lipid-selector (incf lipid-id) :top)
                         for bottom-ga-lipid = (make-random-ga-lipid (list xindex yindex 1) wrapped-xypos z-offset lipid-selector (incf lipid-id) :bottom)
                         do (let* ((top-solute-collides-p (lipid-protein-overlaps-p bounding-box top-ga-lipid ga-solute))
                                   (top-thing (if top-solute-collides-p
                                                  ga-solute
                                                  top-ga-lipid)))
                              (setf (aref array xindex yindex 0) top-thing))
                         do (let* ((bottom-solute-collides-p (lipid-protein-overlaps-p bounding-box bottom-ga-lipid ga-solute))
                                   (bottom-thing (if bottom-solute-collides-p
                                                     ga-solute
                                                     bottom-ga-lipid)))
                              (setf (aref array xindex yindex 1) bottom-thing))))
          ;; Figure out what lipids are adjacent to the solute
          (flet ((adjacent-to-solute-p (membrane xi yi zi)
                   (let* ((yparity (mod yi 2))
                          (zother (mod (1+ zi) 2))
                          (dirs (elt *oddr-directions* yparity)))
                     (loop for dir across dirs
                           for xdir = (elt dir 0)
                           for ydir = (elt dir 1)
                           do (multiple-value-bind (xp yp)
                                  (wrapped-array-indices (+ xi xdir) (+ yi ydir) xnum ynum)
                                (let ((entry (aref array xp yp zi)))
                                  (when (typep entry 'ga-solute) (return-from adjacent-to-solute-p t)))
                                (let ((entry (aref array xp yp zother)))
                                  (when (typep entry 'ga-solute) (return-from adjacent-to-solute-p t)))))
                     (let ((entry (aref array xi yi zother)))
                       (when (typep entry 'ga-solute) (return-from adjacent-to-solute-p t))))))
            (loop for zindex from 0 below 2
                  do (loop for yindex from 0 below ynum
                           do (loop for xindex from 0 below xnum
                                    for entry = (aref array xindex yindex zindex)
                                    for adjacent-to-solute = (adjacent-to-solute-p array xindex yindex zindex)
                                    do (setf (adjacent-to-solute entry) adjacent-to-solute)))))
          (make-instance 'ga-membrane
                         :bounding-box bounding-box
                         :array array
                         :lipid-selector lipid-selector
                         :ga-solute ga-solute))))))

(defun do-score-cell (bounding-box source target collisions)
  (let ((range-start (length (indices collisions)))
        (from-source-transform (transform source))
        (from-target-transform (transform target))
        (target-octree (cond
                         ((typep target 'ga-lipid)
                          (let ((target-lipid-info (lipid-info target))
                                (target-lipid-index (lipid-index target)))
                            (elt (octrees target-lipid-info) target-lipid-index)))
                         ((typep target 'ga-solute)
                          (octree target))
                         (t (error "Illegal target ~a" target))))
        (lipid-info (lipid-info source))
        (lipid-index (lipid-index source)))
    (let ((transform from-source-transform) ; (geom:m*m to-target-transform from-source-transform))
          (source-coordinates (elt (coordinates lipid-info) lipid-index)))
      (loop for coord-index below (length source-coordinates) by 3
            with source-pos = (geom:vec 0.0 0.0 0.0)
            do (progn
                 (geom:inplace-transform-nvector-point source-pos transform source-coordinates coord-index)
                 (when chem:*debug-octree*
                   (let ((atom-name (chem:get-name (order-atom (elt (order-info lipid-info) (/ coord-index 3))))))
                     (when (eq atom-name :c117)
                       (let ((chem:*debug-octree* t))
                         (format t "transform ~%~a~%" transform)
                         (format t "lipid-id ~a  atom ~a  position: ~a~%" (lipid-id source) atom-name source-pos)
                         (let* ((bb-origin (chem:generic-octree-origin target-octree))
                                (bb-half-dim (chem:generic-octree-half-dimension target-octree))
                                (min-point (geom:v- bb-origin bb-half-dim))
                                (max-point (geom:v+ bb-origin bb-half-dim)))
                           (format t "target octree min-point: ~a     max-point: ~a~%" min-point max-point))))))
                 (chem:generic-octree-get-points-within-cutoff
                  target-octree
                  (+ *close-distance* 0.003)
                  coord-index
                  source-pos
                  bounding-box
                  from-target-transform
                  (indices collisions)))))
    (let ((range-end (length (indices collisions))))
      ;; (when (> range-end range-start)
      (vector-push-extend (make-instance 'range
                                         :source-ga-thing source
                                         :target-ga-thing target
                                         :start range-start
                                         :end range-end)
                          (ranges collisions)))))

(defgeneric score-cell (bounding-box cur other collisions))

(defmethod score-cell (bounding-box (cur ga-solute) (other ga-solute) collisions) #| Do nothing |#)

(defmethod score-cell (bounding-box (source ga-lipid) (target ga-lipid) collisions)
  (do-score-cell bounding-box source target collisions))

(defmethod score-cell (bounding-box (source ga-solute) (target ga-lipid) collisions) #| Do nothing |#)

(defmethod score-cell (bounding-box source (target null) collisions) #| Do nothing |# )

(defmethod score-cell (bounding-box (source ga-lipid) (target ga-solute) collisions) #| Do nothing |# )

(defun make-collisions (&optional (xdim 16) (ydim 16) (zdim 2))
  (let* ((ranges (make-array (+ 1024 (* xdim ydim zdim)) :adjustable t :fill-pointer 0))
         (indices (make-array 2048 :adjustable t :fill-pointer 0))
         (collisions (make-instance 'collisions :indices indices :ranges ranges)))
    collisions))

(defclass lipid-pair ()
  ((source-cell :initarg :source-cell :accessor source-cell)
   (target-cell :initarg :target-cell :accessor target-cell)))

(defmethod print-object ((obj lipid-pair) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a -> ~a" (source-cell obj) (target-cell obj))))

(defclass work-list ()
  ((lipid-pairs :initarg :lipid-pairs :accessor lipid-pairs)
   (solute-tests :initarg :solute-tests :accessor solute-tests)))
  
(defun maybe-score-cell-work-list (sc-cur entry work-list &key debug)
  "Only score the cell if lipid-ids are in ascending order"
  (let ((id1 (lipid-id sc-cur))
        (id2 (lipid-id entry)))
    (when (< id1 id2)
      (when debug (format t "~a/" (array-indices entry)))
      (vector-push-extend (make-instance 'lipid-pair
                                         :source-cell (array-indices sc-cur)
                                         :target-cell (array-indices entry))
                          (lipid-pairs work-list)))))

(defun build-ga-membrane-work-list (ga-membrane &key debug)
  (let* ((array (array ga-membrane))
         (dims (array-dimensions array))
         (xnum (elt dims 0))
         (ynum (elt dims 1))
         (znum (elt dims 2))
         (work-list (make-instance 'work-list
                                   :lipid-pairs (make-array 1024 :adjustable t :fill-pointer 0)
                                   :solute-tests (make-array 256 :adjustable t :fill-pointer 0))))
    (flet ((aref-wrapped (xi yi zi xd yd)
             (multiple-value-bind (xp yp)
                 (wrapped-array-indices (+ xi xd) (+ yi yd) xnum ynum)
               (when debug (format t " ~a ~a ~a, " xp yp zi))
               (aref (array ga-membrane) xp yp zi))))
      (loop for zi below znum
            for zother = (mod (1+ zi) 2)
            do (loop for yi below ynum
                     for yparity = (mod yi 2)
                     for dirs = (elt *oddr-directions* yparity)
                     do (loop for xi below xnum
                              for sc-cur = (aref (array ga-membrane) xi yi zi)
                              when (typep sc-cur 'ga-lipid)
                                do (progn
                                     (when debug (format t "~a ~a ~a | dirs: ~a  xnum: ~a  ynum: ~a~%" xi yi zi dirs xnum ynum))
                                     (loop for dir across dirs
                                           for xdir = (elt dir 0)
                                           for ydir = (elt dir 1)
                                           for xd = (+ xi xdir)
                                           for yd = (+ yi ydir)
                                           do (multiple-value-bind (xp yp)
                                                  (wrapped-array-indices (+ xi xdir) (+ yi ydir) xnum ynum)
                                                (when debug (format t "~a,~a[~a,~a]{~a,~a}" xdir ydir xd yd xp yp))
                                                (let ((entry (aref array xp yp zi)))
                                                  (when (typep entry 'ga-lipid) (maybe-score-cell-work-list sc-cur entry work-list)))
                                                (let ((entry (aref array xp yp zother)))
                                                  (when (typep entry 'ga-lipid) (maybe-score-cell-work-list sc-cur entry work-list)))))
                                     (let ((entry (aref array xi yi zother)))
                                       (when (typep entry 'ga-lipid) (maybe-score-cell-work-list sc-cur entry work-list))))
                                   (when debug (format t "~%") (finish-output))
                                   (when (adjacent-to-solute sc-cur)
                                     (vector-push-extend (array-indices sc-cur) (solute-tests work-list)))))))
    work-list))

(defun score-work-list (work-list membrane &key debug)
  (let* ((array (array membrane))
         (dims (array-dimensions array))
         (xnum (elt dims 0))
         (ynum (elt dims 1))
         (znum (elt dims 2))
         (bounding-box (bounding-box membrane))
         (collisions (make-collisions xnum ynum znum)))
    (loop for lipid-pair across (lipid-pairs work-list)
          for source-cell = (source-cell lipid-pair)
          for target-cell = (target-cell lipid-pair)
          for source-xi = (elt source-cell 0)
          for source-yi = (elt source-cell 1)
          for source-zi = (elt source-cell 2)
          for target-xi = (elt target-cell 0)
          for target-yi = (elt target-cell 1)
          for target-zi = (elt target-cell 2)
          for sc-cur = (aref array source-xi source-yi source-zi)
          for sc-other = (aref array target-xi target-yi target-zi)
          do (score-cell bounding-box sc-cur sc-other collisions))
    (loop for lipid-cell across (solute-tests work-list)
          for source-xi = (elt lipid-cell 0)
          for source-yi = (elt lipid-cell 1)
          for source-zi = (elt lipid-cell 2)
          for sc-cur = (aref array source-xi source-yi source-zi)
          do (score-cell bounding-box sc-cur (ga-solute membrane) collisions))
    (values (/ (length (indices collisions)) 2) collisions)))

(defun score-work-list-atom-vectors (work-list atom-vectors &key debug no-result)
  (let ((close-contact-vec (make-array 1024 :adjustable t :fill-pointer 0))
        (cutoff-squared (* *close-distance* *close-distance*)))
    (let* ((array (array atom-vectors))
           (dims (array-dimensions array))
           (xnum (elt dims 0))
           (ynum (elt dims 1))
           (znum (elt dims 2))
           (bounding-box (bounding-box atom-vectors))
           (num 0)
           (close-contacts 0))
      (loop for lipid-pair across (lipid-pairs work-list)
            for source-cell = (source-cell lipid-pair)
            for target-cell = (target-cell lipid-pair)
            for source-xi = (elt source-cell 0)
            for source-yi = (elt source-cell 1)
            for source-zi = (elt source-cell 2)
            for target-xi = (elt target-cell 0)
            for target-yi = (elt target-cell 1)
            for target-zi = (elt target-cell 2)
            for sc-cur = (aref array source-xi source-yi source-zi)
            for sc-other = (aref array target-xi target-yi target-zi)
            for m1-atoms = (elt (atom-vectors atom-vectors) sc-cur)
            for m2-atoms = (elt (atom-vectors atom-vectors) sc-other)
            do (multiple-value-bind (one-num one-close-contacts)
                   (close-contacts-between-atom-vectors m1-atoms m2-atoms bounding-box close-contact-vec cutoff-squared num)
                 (setf num one-num)
                 (incf close-contacts one-close-contacts)))
      (loop for lipid-cell across (solute-tests work-list)
            for source-xi = (elt lipid-cell 0)
            for source-yi = (elt lipid-cell 1)
            for source-zi = (elt lipid-cell 2)
            for sc-cur = (aref array source-xi source-yi source-zi)
            for m1-atoms = (elt (atom-vectors atom-vectors) sc-cur)
            with solute-atoms = (elt (atom-vectors atom-vectors) (solute-atom-vector-index atom-vectors))
            do (multiple-value-bind (one-num one-close-contacts)
                   (close-contacts-between-atom-vectors m1-atoms solute-atoms bounding-box close-contact-vec cutoff-squared num)
                 (setf num one-num)
                 (incf close-contacts one-close-contacts)))
      (format t "There are ~d close-contacts of a possible ~a~%" close-contacts num)
      (if no-result
          nil close-contact-vec))))

(defun maybe-score-cell (bounding-box sc-cur entry collisions &key debug)
  "Only score the cell if lipid-ids are in ascending order"
  (let ((id1 (lipid-id sc-cur))
        (id2 (lipid-id entry)))
    (when (< id1 id2)
      (when debug (format t "~a/" (array-indices entry)))
      (score-cell bounding-box sc-cur entry collisions))))

(defun score-ga-membrane (ga-membrane &key debug)
  (let* ((array (array ga-membrane))
         (dims (array-dimensions array))
         (xnum (elt dims 0))
         (ynum (elt dims 1))
         (znum (elt dims 2))
         (bounding-box (bounding-box ga-membrane))
         (collisions (make-collisions xnum ynum znum)))
    (flet ((aref-wrapped (xi yi zi xd yd)
             (multiple-value-bind (xp yp)
                 (wrapped-array-indices (+ xi xd) (+ yi yd) xnum ynum)
               (when debug (format t " ~a ~a ~a, " xp yp zi))
               (aref (array ga-membrane) xp yp zi))))
      (loop for zi below znum
            for zother = (mod (1+ zi) 2)
            do (loop for yi below ynum
                     for yparity = (mod yi 2)
                     for dirs = (elt *oddr-directions* yparity)
                     do (loop for xi below xnum
                              for sc-cur = (aref (array ga-membrane) xi yi zi)
                              when (typep sc-cur 'ga-lipid)
                                do (progn
                                     (when debug (format t "~a ~a ~a | dirs: ~a  xnum: ~a  ynum: ~a~%" xi yi zi dirs xnum ynum))
                                     (loop for dir across dirs
                                           for xdir = (elt dir 0)
                                           for ydir = (elt dir 1)
                                           for xd = (+ xi xdir)
                                           for yd = (+ yi ydir)
                                           do (multiple-value-bind (xp yp)
                                                  (wrapped-array-indices (+ xi xdir) (+ yi ydir) xnum ynum)
                                                (when debug (format t "~a,~a[~a,~a]{~a,~a}" xdir ydir xd yd xp yp))
                                                (let ((entry (aref array xp yp zi)))
                                                  (when (typep entry 'ga-lipid) (maybe-score-cell bounding-box sc-cur entry collisions)))
                                                (let ((entry (aref array xp yp zother)))
                                                  (when (typep entry 'ga-lipid) (maybe-score-cell bounding-box sc-cur entry collisions)))))
                                     (let ((entry (aref array xi yi zother)))
                                       (when (typep entry 'ga-lipid) (maybe-score-cell bounding-box sc-cur entry collisions))))
                                   (when debug (format t "~%") (finish-output))
                                   (when (adjacent-to-solute sc-cur)
                                     (do-score-cell bounding-box sc-cur (ga-solute ga-membrane) collisions))))))
    (values (/ (length (indices collisions)) 2) collisions)))

(defun mutate-all-lipids (ga-membrane)
  "Mutate the lipid in every position in the membrane"
  (loop for row-major-index from 0 below (array-total-size (array ga-membrane))
        for new-lipid-info-index = (select-lipid (lipid-selector ga-membrane))
        for new-lipid-info = (car new-lipid-info-index)
        for new-lipid-index = (cdr new-lipid-info-index)
        for lipid-or-solute = (row-major-aref (array ga-membrane) row-major-index)
        do (unless (typep lipid-or-solute 'ga-solute)
             (setf (lipid-info lipid-or-solute) new-lipid-info
                   (lipid-index lipid-or-solute) new-lipid-index))))


(defun collisions-as-point-pairs (collisions)
  (let ((point-pairs (make-array 1024 :adjustable t :fill-pointer 0)))
    (loop for rangei from 0 below (length (ranges collisions))
          for range = (elt (ranges collisions) rangei)
          for target-ga-thing = (target-ga-thing range)
          for source-ga-thing = (source-ga-thing range)
          for start = (start range)
          for end = (end range)
          if (typep target-ga-thing 'ga-lipid)
            do (let* ((source-lipid-info (lipid-info source-ga-thing))
                      (source-lipid-index (lipid-index source-ga-thing))
                      (target-lipid-info (lipid-info target-ga-thing))
                      (target-lipid-index (lipid-index target-ga-thing))
                      (source-coordinates (elt (coordinates source-lipid-info) source-lipid-index))
                      (target-coordinates (elt (coordinates target-lipid-info) target-lipid-index)))
                 (loop for index from start below end by 2
                       for source-coord-index = (elt (indices collisions) index)
                       for target-coord-index = (elt (indices collisions) (1+ index))
                       for source-pos = (geom:vec 0.0 0.0 0.0)
                       for target-pos = (geom:vec 0.0 0.0 0.0)
                       do (progn
                            (geom:vec-extract-transformed source-pos source-coordinates source-coord-index (transform source-ga-thing))
                            (geom:vec-extract-transformed target-pos target-coordinates target-coord-index (transform target-ga-thing))
                            (let ((points (list source-pos target-pos)))
                              (vector-push-extend points point-pairs)))))
          else if (typep target-ga-thing 'ga-solute)
                 do (loop for index from start below end by 2
                          for source-coord-index = (elt (indices collisions) index)
                          for atom-res-mol = (elt (indices collisions) (1+ index))
                          for atom = (first atom-res-mol)
                          for solute-pos = (chem:get-position atom)
                          for source-pos = (geom:vec 0.0 0.0 0.0)
                          do (let ((point-pair (list source-pos solute-pos)))
                               (vector-push-extend point-pair point-pairs))))
    point-pairs))

(defun add-point-pairs-as-shape (widget point-pairs &key (color #(1 0 1)) (radius 4.0) (name "collisions"))
  (let ((shapes (make-array (* 3 (length point-pairs)))))
    (loop for index from 0 below (length point-pairs)
          for shape-index from 0 by 3
          for point-pair = (elt point-pairs index)
          for point1 = (first point-pair)
          for point2 = (second point-pair)
          for cyl = (vector "cylinder"
                               (vector (geom:vx point1) (geom:vy point1) (geom:vz point1))
                               (vector (geom:vx point2) (geom:vy point2) (geom:vz point2))
                               color
                               radius)
          for sphere1 = (vector "sphere"
                                (vector (geom:vx point1) (geom:vy point1) (geom:vz point1))
                                color
                                radius)
          for sphere2 = (vector "sphere"
                                (vector (geom:vx point2) (geom:vy point2) (geom:vz point2))
                                color
                                radius)
          do (setf (elt shapes shape-index) cyl)
          do (setf (elt shapes (+ shape-index 1)) sphere1)
          do (setf (elt shapes (+ shape-index 2)) sphere2))
    (let ((add-shape (find-symbol "ADD-SHAPE" :nglv)))
      (funcall add-shape widget shapes :name name))))

(defun add-collisions-as-shape (widget collisions &key (color #(1 0 1)) (radius 2.1))
  (let ((point-pairs (collisions-as-point-pairs collisions)))
    (add-point-pairs-as-shape widget point-pairs :color color :radius radius)))

(defun find-collision-range (collisions ga-thing1 ga-thing2)
  (loop for range across (ranges collisions)
        when (or (and (eq ga-thing1 (source-ga-thing range))
                      (eq ga-thing2 (target-ga-thing range)))
                 (and (eq ga-thing1 (target-ga-thing range))
                      (eq ga-thing2 (source-ga-thing range))))
          do (return-from find-collision-range range)))

(defclass atom-vectors ()
  ((array :initarg :array :accessor array)
   (solute-atom-vector-index :initarg :solute-atom-vector-index :accessor solute-atom-vector-index)
   (atom-vectors :initarg :atom-vectors :accessor atom-vectors)
   (bounding-box :initarg :bounding-box :accessor bounding-box)))

(defun build-aggregate-from-ga-membrane (ga-membrane &key debug)
  "Build an aggregate for the membrane and return it as well as a vector of vectors of atoms for
testing the scoring."
  (let ((solute (chem:matter-copy (solute (ga-solute ga-membrane)))))
    (chem:apply-transform-to-atoms solute (transform (ga-solute ga-membrane)))
    (let* ((aggregate (cando-user:combine solute))
           (dimensions (array-dimensions (array ga-membrane)))
           (xnum (elt dimensions 0))
           (ynum (elt dimensions 1))
           (znum (elt dimensions 2))
           (bounding-box (bounding-box ga-membrane))
           (lipid-molecules (make-array 1024 :adjustable t :fill-pointer 0))
           (atom-vectors (make-instance 'atom-vectors
                                        :array (make-array (list xnum ynum znum) :element-type t)
                                        :bounding-box bounding-box)))
      (chem:set-bounding-box aggregate (bounding-box ga-membrane))
      (loop for xindex from 0 below (elt dimensions 0)
            do (loop for yindex from 0 below (elt dimensions 1)
                     do (loop for zindex from 0 to 1
                              for ga-thing = (aref (array ga-membrane) xindex yindex zindex)
                              unless (typep ga-thing 'ga-solute)
                                do (let* ((ga-lipid ga-thing)
                                          (transform (transform ga-lipid))
                                          (lipid-info (lipid-info ga-lipid))
                                          (lipid-index (lipid-index ga-lipid)))
                                     (let ((lipid-mol (make-molecule-from-ga-lipid ga-lipid :debug debug)))
                                       (chem:add-molecule aggregate lipid-mol)
                                       (let ((lipid-mol-index (vector-push-extend lipid-mol lipid-molecules)))
                                         (setf (aref (array atom-vectors) xindex yindex zindex) lipid-mol-index)))))))
      ;; Build a vector of atom vectors
      (let ((molecule-atoms (make-array (1+ (length lipid-molecules))))
            (total-pairs 0))
        (loop for index below (length lipid-molecules)
              for lipid-molecule = (elt lipid-molecules index)
              for mol-list = (list lipid-molecule)
              for atoms = (make-array (chem:number-of-atoms lipid-molecule))
              for ai = -1
              do (cando:do-residues (res lipid-molecule)
                   (let ((res-list (list* res mol-list)))
                     (cando:do-atoms (atm res)
                       (let ((atm-list (list* atm res-list)))
                         (setf (elt atoms (incf ai)) atm-list)))))
              do (setf (elt molecule-atoms index) atoms))
        ;; add in the solute atoms at the end
        (let ((solute-atoms (make-array (chem:number-of-atoms solute)))
              (si -1))
          (cando:do-molecules (mol solute)
            (let ((mol-list (list mol)))
              (cando:do-residues (res mol)
                (let ((res-list (list* res mol-list)))
                  (cando:do-atoms (atm res)
                    (let ((atm-list (list* atm res-list)))
                      (setf (elt solute-atoms (incf si)) atm-list)))))))
          (setf (elt molecule-atoms (1- (length molecule-atoms))) solute-atoms))
        (setf (solute-atom-vector-index atom-vectors) (1- (length molecule-atoms)))
        (setf (atom-vectors atom-vectors) molecule-atoms)
        (values aggregate atom-vectors)))))

(defun simulate-membrane (aggregate)
  (let* ((energy-function (chem:make-energy-function aggregate))
         (dyn (dynamics:make-atomic-simulation energy-function
                                               :delta-t 1.0e-3
                                               :accumulate-coordinates t)))
    (chem:setf-velocity-scale energy-function 0.1 0.1 0.1)
    dyn))
         

(defun load-test-molecule ()
  (cando-user:setup-amber-paths)
  (cando-user:easy-gaff)
  (cando-user:source "leaprc.protein.ff14SB")
  (cando-user:source "leaprc.lipid14")
  (cando-user:source "leaprc.water.tip3p")
  (cando-user:source "leaprc.gaff")
  (let* ((prot0 (chem:moe-read-aggregate "/Users/meister/Development/cando-dev/extensions/cando/src/lisp/membrane-packer/data/merck_suvorexant_derivatives.moe"))
         (prot (chem:make-aggregate)))
    (chem:add-molecule prot (cando:mol prot0 0))
    (chem:add-molecule prot (cando:mol prot0 1))
    (cando:move-geometric-center-to-origin prot)
    (leap:assign-atom-types prot)
    prot
    ))

(defun test-everything1 ()
  (format t "Building an empty membrane~%")
  (let* ((prot0 (chem:moe-read-aggregate "/Users/meister/Development/cando-dev/extensions/cando/src/lisp/membrane-packer/data/merck_suvorexant_derivatives.moe"))
         (prot (chem:make-aggregate)))
    (chem:add-molecule prot (cando:mol prot0 0))
    (chem:add-molecule prot (cando:mol prot0 1))
    (cando:move-geometric-center-to-origin prot)
    (leap:assign-atom-types prot)
    (let ((memb (build-ga-membrane prot (chem:make-bounding-box (list 100.0 100.0 40.0)))))
      memb)))

(defun test-everything2 (aggregate)
  (simulate-membrane aggregate))

(defun test-everything (&optional (num 10))
  (test-everything0)
  (let ((a (test-everything1)))
    (let ((dyn (test-everything2 a)))
      (dynamics:with-dynamics (dyn)
        (loop for x from 0 below num
              do (format t "Dynamics step ~a~%" x)
              do (dynamics:velocity-verlet-step-limit-displacement dyn))))))

(defun add-close-contacts-as-shape (widget close-contacts &key (radius 2.0))
  "Render the close-contacts returned by close-contacts as shapes"
  (let ((point-pairs (make-array (length close-contacts) :fill-pointer 0 :adjustable t)))
    (loop for entry across close-contacts
          for dist = (first entry)
          for atom1-info = (second entry)
          for atom2-info = (third entry)
          for atom1 = (first atom1-info)
          for atom2 = (first atom2-info)
          for pos1 = (chem:get-position atom1)
          for pos2 = (chem:get-position atom2)
          for pair = (list pos1 pos2)
          do (vector-push-extend pair point-pairs))
    (add-point-pairs-as-shape widget point-pairs :color #(1 1 0) :radius radius :name "close-contacts")))


(defun close-contacts-between-atom-vectors (m1-atoms m2-atoms bounding-box close-contact-vec cutoff-squared num)
  (let ((close-contacts 0))
    (loop for i1 below (length m1-atoms)
          for atom1-res-mol = (elt m1-atoms i1)
          for atom1 = (first atom1-res-mol)
          do (loop for i2 below (length m2-atoms)
                   for atom2-res-mol = (elt m2-atoms i2)
                   for atom2 = (first atom2-res-mol)
                   for dist-squared = (chem:bounding-box-distance-squared-between-two-atoms bounding-box atom1 atom2)
                   do (incf num)
                   when (< dist-squared cutoff-squared)
                     do (progn
                          (vector-push-extend (list (sqrt dist-squared) atom1-res-mol atom2-res-mol) close-contact-vec)
                          (incf close-contacts))))
    (values num close-contacts)))

(defun close-contacts (atom-vectors &key (cutoff *close-distance*))
  "Build the atoms for the atom-vectors and look for close contacts between atoms generating the same
pair-list that would be tested using the octree approach. Return the results in a vector of
elements that look like (list dist-between-atoms atom1-res-mol atom2-res-mol).
The atom-res-mol is a list of the atom,residue and molecule."
  (let ((close-contact-vec (make-array 1024 :adjustable t :fill-pointer 0))
        (cutoff-squared (* cutoff cutoff))
        (molecule-atoms (atom-vectors atom-vectors))
        (total-pairs 0)
        (close-contacts 0))
    ;; Now count the number of interactions -> total-pairs
    (loop for index1 from 0 below (1- (length molecule-atoms))
          for m1-atoms = (elt molecule-atoms index1)
          for num-m1-atoms = (length m1-atoms)
          do (incf total-pairs (* num-m1-atoms
                                  (loop for index2 from (1+ index1) below (length molecule-atoms)
                                        sum (length (elt molecule-atoms index2))))))
    ;; Check for collisions between atoms in lipids
    (let ((bounding-box (bounding-box atom-vectors))
          (progress (cando-user:make-progress-bar :total total-pairs))
          (num 0))
      (loop for index1 from 0 below (1- (length molecule-atoms))
            for m1-atoms = (elt molecule-atoms index1)
            do (loop for index2 from (1+ index1) below (length molecule-atoms)
                     for m2-atoms = (elt molecule-atoms index2)
                     do (multiple-value-bind (one-num one-close-contacts)
                            (close-contacts-between-atom-vectors m1-atoms m2-atoms bounding-box close-contact-vec cutoff-squared num)
                          (incf close-contacts one-close-contacts)
                          (setf num one-num)
                          (cando-user:progress-advance progress num))))
      (format t "There are ~d close-contacts of a possible ~a~%" close-contacts num)
      close-contact-vec)))

(defun find-ga-lipid (ga-membrane lipid-id)
  "Return the ga-lipid in the GA-MEMBRANE that has the LIPID-ID"
  (loop for row-major-index below (array-total-size (array ga-membrane))
        for ga-thing = (row-major-aref (array ga-membrane) row-major-index)
        when (and (typep ga-thing 'ga-lipid)
                  (= (lipid-id ga-thing) lipid-id))
          do (return-from find-ga-lipid ga-thing)))
  
(defun focused-aggregate (ga-membrane lipids)
  (let ((aggregate (chem:make-aggregate))
        (atom-vectors (make-array (length lipids))))
    (chem:set-bounding-box aggregate (bounding-box ga-membrane))
    (loop for ga-thing in lipids
          for index from 0
          do (let* ((ga-lipid ga-thing)
                    (lipid-mol (make-molecule-from-ga-lipid ga-lipid :debug t)))
               (chem:add-molecule aggregate lipid-mol)
               (let ((mol-list (list lipid-mol))
                     (atoms (make-array (chem:number-of-atoms lipid-mol)))
                     (ai -1))
                 (cando:do-residues (res lipid-mol)
                   (let ((res-list (list* res mol-list)))
                     (cando:do-atoms (atm res)
                       (let ((atm-list (list* atm res-list)))
                         (setf (elt atoms (incf ai)) atm-list)))))
                 (setf (elt atom-vectors index) atoms))))
    (values aggregate (make-instance 'atom-vectors
                                     :atom-vectors atom-vectors
                                     :bounding-box (bounding-box ga-membrane)))))

(defun lots-of-ga-membranes (original num)
  (let ((membranes (make-array num)))
    (loop for index below num
          for new-membrane = (copy-ga-membrane original)
          do (mutate-all-lipids new-membrane)
          do (setf (aref membranes index) new-membrane))
    membranes))



(defun recombine-slice (memb1 memb2)
  (let* ((dims (array-dimensions (array memb1)))
         (xdim (elt dims 0))
         (ydim (elt dims 1))
         (zdim (elt dims 2)))
    (let ((dir (random 3))
          (from-to (random 2)))
      (multiple-value-bind (from to)
          (if (= from-to 0)
              (values memb1 (copy-ga-membrane memb2))
              (values memb2 (copy-ga-membrane memb1)))
        (multiple-value-bind (xlim ylim zlim)
            (case dir
              (0 (values (random xdim) ydim zdim))
              (1 (values xdim (random ydim) zdim))
              (2 (values xdim ydim (random zdim))))
          (loop for xi from 0 below xlim
                do (loop for yi from 0 below ylim
                         do (loop for zi from 0 below zlim
                                  when (typep (aref (array from) xi yi zi) 'ga-lipid)
                                    do (setf (aref (array to) xi yi zi) (aref (array from) xi yi zi))))))
        to))))


(defun mutate-lipids (orig-membrane &optional (fraction-mutations 0.03))
  (let ((membrane orig-membrane)
        (num-mutations (floor (* fraction-mutations (array-total-size (array orig-membrane))))))
    (loop for num from 0 below num-mutations
          for row-major-index = (random (array-total-size (array orig-membrane)))
          for lipid-or-solute = (row-major-aref (array membrane) row-major-index)
          when (typep lipid-or-solute 'ga-lipid)
            do (let* ((new-lipid-info-index (select-lipid (lipid-selector orig-membrane)))
                      (new-lipid-info (car new-lipid-info-index))
                      (new-lipid-index (cdr new-lipid-info-index)))
                 (setf (lipid-info lipid-or-solute) new-lipid-info
                       (lipid-index lipid-or-solute) new-lipid-index)))
    membrane))

(defun mutate-position (orig-membrane &key (delta-x 0.1) (delta-y 0.1) (fraction-mutations 0.1))
  (let ((membrane orig-membrane)
        (num-mutations (floor (* fraction-mutations (array-total-size (array orig-membrane))))))
    (loop for num from 0 below num-mutations
          for row-major-index = (random (array-total-size (array orig-membrane)))
          for ga-thing = (row-major-aref (array membrane) row-major-index)
          when (typep ga-thing 'ga-lipid)
            do (let* ((lipid-shift (geom:vec (- (random (* 2.0 delta-x)) delta-x)
                                              (- (random (* 2.0 delta-y)) delta-y)
                                              0.0))
                      (extra-transform (geom:make-m4-translate lipid-shift))
                      (transform (geom:m*m extra-transform (transform ga-thing))))
                 (setf (transform ga-thing) transform)))
    membrane))

(defun score-all (membs scoring-membrane)
  (format t "In score-all~%")
  (let ((scores (loop for index below (length membs)
                      for memb = (elt membs index)
                      collect (cons (score-ga-membrane memb scoring-membrane) memb))))
    (sort scores #'> :key #'car)))
       


(defun try-times (memb num)
  (let ((min-collisions 9999))
    (loop for index below num
          for score = (progn
                        (mutate-all-lipids memb)
                        (score-ga-membrane memb))
          do (setf min-collisions (min min-collisions score))
          do (format t "num = ~a   best-score ~a~%" index min-collisions))))


(defun generic-octree-as-shape (octree transform shape)
  (let ((origin (chem:generic-octree-origin octree))
        (half-dim (chem:generic-octree-half-dimension octree)))
    (flet ((calc-point (xd yd zd)
             (let ((xd (if (= xd 0) -1 1))
                   (yd (if (= yd 0) -1 1))
                   (zd (if (= zd 0) -1 1)))
               (geom:m*v transform
                         (geom:v+ origin (geom:vec (* (geom:vx half-dim) xd)
                                                   (* (geom:vy half-dim) yd)
                                                   (* (geom:vz half-dim) zd))))))
           (draw-line (p1 p2 &optional (color #(0.5 0.5 0.5)))
             (vector "cylinder"
                     (vector (float (geom:vx p1) 1.0s0) (float (geom:vy p1) 1.0s0) (float (geom:vz p1) 1.0s0))
                     (vector (float (geom:vx p2) 1.0s0) (float (geom:vy p2) 1.0s0) (float (geom:vz p2) 1.0s0))
                     color
                     0.02)))
      (let ((p000 (calc-point 0 0 0))
            (p001 (calc-point 0 0 1))
            (p010 (calc-point 0 1 0))
            (p011 (calc-point 0 1 1))
            (p100 (calc-point 1 0 0))
            (p101 (calc-point 1 0 1))
            (p110 (calc-point 1 1 0))
            (p111 (calc-point 1 1 1)))
        (vector-push-extend (draw-line p000 p001) shape)
        (vector-push-extend (draw-line p000 p010) shape)
        (vector-push-extend (draw-line p000 p100) shape)
        (vector-push-extend (draw-line p001 p011) shape)
        (vector-push-extend (draw-line p001 p101) shape)
        (vector-push-extend (draw-line p010 p011) shape)
        (vector-push-extend (draw-line p010 p110) shape)
        (vector-push-extend (draw-line p100 p101) shape)
        (vector-push-extend (draw-line p100 p110) shape)
        (vector-push-extend (draw-line p011 p111) shape)
        (vector-push-extend (draw-line p101 p111) shape)
        (vector-push-extend (draw-line p110 p111) shape)))))

(defun find-ga-lipid-position (ga-lipid residue-name atom-name)
  (let* ((lipid-index (lipid-index ga-lipid))
         (lipid-info (lipid-info ga-lipid))
         (transform (transform ga-lipid))
         (coordinates (elt (coordinates lipid-info) lipid-index))
         (order-infos (order-info lipid-info)))
    (loop for index below (length order-infos)
          for coord-index = (* 3 index)
          for order-info = (elt order-infos index)
          for res = (order-residue order-info)
          for atom = (order-atom order-info)
          with pos = (geom:vec 0.0 0.0 0.0)
          when (and (eq (chem:get-name res) residue-name)
                    (eq (chem:get-name atom) atom-name))
            do (progn
                 (geom:vec-extract-transformed pos coordinates coord-index transform)
                 (return-from find-ga-lipid-position pos)))))

(defun find-ga-lipid-octree (ga-lipid residue-name atom-name)
  (let* ((lipid-index (lipid-index ga-lipid))
         (lipid-info (lipid-info ga-lipid))
         (transform (transform ga-lipid))
         (octree (elt (octrees lipid-info) lipid-index))
         (find-octree nil)
         (stack (make-array 32 :fill-pointer 0)))
    (chem:walk-generic-octree octree
                              (lambda (ot)
                                (let ((depth (chem:generic-octree-depth ot)))
                                  (setf (fill-pointer stack) (1+ depth))
                                  (setf (elt stack depth) ot))
                                (if (and (chem:is-leaf-node ot)
                                         (chem:generic-octree-data-boundp ot))
                                    (let* ((coord-index (chem:generic-octree-data ot))
                                           (index (/ coord-index 3))
                                           (order-info (elt (order-info lipid-info) index))
                                           (residue (order-residue order-info))
                                           (atom (order-atom order-info)))
                                      (when (and (eq (chem:get-name residue) residue-name)
                                                 (eq (chem:get-name atom) atom-name))
                                        (format t "generic-octree-data -> residue: ~a   atom: ~a~%" residue atom)
                                        (return-from find-ga-lipid-octree (values (chem:generic-octree-data ot) stack))))
                                    (format t "Skipping octree ~a  data: ~a origin: ~a~%" ot (and (chem:is-leaf-node ot)
                                                                                                  (chem:generic-octree-data-boundp ot)
                                                                                                  (chem:generic-octree-data ot))
                                            (chem:generic-octree-origin ot)))))))

(defun ga-lipid-octree-as-shape (ga-lipid)
  (let* ((shape (make-array 1024 :adjustable t :fill-pointer 0))
         (lipid-index (lipid-index ga-lipid))
         (lipid-info (lipid-info ga-lipid))
         (transform (transform ga-lipid))
         (octree (elt (octrees lipid-info) lipid-index)))
    (chem:walk-generic-octree octree
                              (lambda (ot)
                                (generic-octree-as-shape ot transform shape)))
    shape))

(defun evolve (solute &rest args &key input-bounding-box
                                   (lipid-selector (list (cons 1.0 *popc*)))
                                   (population-size 200)
                                   (select-fraction 0.1)
                                   (number-of-generations 100)
                                   (parallel nil))
  (let* ((template-membrane (build-ga-membrane solute :input-bounding-box input-bounding-box :lipid-selector lipid-selector))
         (generation 0))
    ;; Fill the population at start
    (let ((membranes (loop for index below population-size
                           for membrane = (let ((copy-membrane (copy-ga-membrane template-membrane)))
                                            (mutate-all-lipids copy-membrane)
                                            copy-membrane)
                           collect membrane)))
      (let (scored-membranes
            (work-list (build-ga-membrane-work-list (first membranes))))
        (loop named evolution
              for unsorted-scored-membranes = (if parallel
                                                  (lparallel:pmap 'vector
                                                                  (lambda (memb)
                                                                    (multiple-value-bind (num-collisions collisions)
                                                                        (score-work-list work-list memb)
                                                                      (list num-collisions memb collisions)))
                                                                  membranes)
                                                  (map 'vector
                                                       (lambda (memb)
                                                         (multiple-value-bind (num-collisions collisions)
                                                             (score-work-list work-list memb)
                                                           (list num-collisions memb collisions)))
                                                       membranes))
              for sorted-scored-membranes = (sort unsorted-scored-membranes #'< :key #'first)
              for best-score = (car (elt sorted-scored-membranes 0))
              do (setf scored-membranes sorted-scored-membranes)
              do (format t "Generation ~3a ... least collisions -> ~3a~%" generation (loop for index from 0 below 20
                                                                                           collect (car (elt scored-membranes index))))
              do (finish-output)
              do (when (= (first (elt scored-membranes 0)) 0) (return-from evolution t))
              do (incf generation)
              do (when (> generation number-of-generations) (return-from evolution nil))
                 ;; Scramble the bad lipids
#|              do (loop for mutant-index below (length scored-membranes)
                       for score-memb-collisions = (elt scored-membranes mutant-index)
                       for score = (first score-memb-collisions) 
                       for membrane = (second score-memb-collisions) 
                       for collisions = (third score-memb-collisions) 
                       do (loop for range across (ranges collisions) 
                                for source = (source-ga-thing range)   
                                for target = (target-ga-thing range)   
                                when (typep source 'ga-lipid) 
                                  do (let ((lipid-info (lipid-info source))) 
                                       (setf (lipid-index source) (random (length (coordinates lipid-info)))))
                                when (typep target 'ga-lipid) 
                                  do (let ((lipid-info (lipid-info target))) 
                                       (setf (lipid-index target) (random (length (coordinates lipid-info)))))) 
                       )
|#
                 ;; generate recombined membranes
              do (loop for mutant-index below (length membranes) 
                             for parent1-index = (random (floor (* population-size select-fraction))) 
                             for parent2-index = (random (floor (* population-size select-fraction))) 
                             for parent1 = (second (elt scored-membranes parent1-index)) 
                             for parent2 = (second (elt scored-membranes parent2-index)) 
                             for mutant = (recombine-slice parent1 parent2) 
                             do (setf (elt membranes mutant-index) mutant)) 
                 ;; Do some random mutations of *mutate-lipid-fraction* 
              do (loop for mutant-index below (length membranes) 
                       for mutate-p = (< (random 1.0) *mutate-lipid-fraction*) 
                       for membrane = (elt membranes mutant-index) 
                       when mutate-p                               
                         do (mutate-lipids membrane))              
                 ;; Do some random shifting of *mutate-position-fraction* 
              do (loop for mutant-index below (length membranes) 
                       for mutate-p = (< (random 1.0) *mutate-position-fraction*) 
                       for membrane = (elt membranes mutant-index) 
                       when mutate-p                               
                         do (mutate-position membrane)))
        scored-membranes))))


(defun pack (solute &rest args &key input-bounding-box
                                 (lipid-selector (list (cons 1.0 *popc*)))
                                 (number-of-generations 100)
                                 (parallel nil))
  (apply 'evolve solute :input-bounding-box input-bounding-box args))
