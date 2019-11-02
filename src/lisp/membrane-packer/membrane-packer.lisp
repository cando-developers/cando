(in-package :membrane-packer)

#| Working with hex-grids is tricky.
This web page is excellent: https://www.redblobgames.com/grids/hexagons

This code uses Offset coordinates odd-r
|#
(defparameter *x-expand* 25.0) ; 40.0)
(defparameter *y-expand* 25.0) ; 40.0)
(defparameter *z-expand* 10.0)

(defparameter *1fs-time-units* (/ 1.0 48.89)
  "One time unit is 48.89 fs - if we want 1 fs then we want 1/48.89 of that")
(defparameter *select-fraction* 0.1)
(defparameter *fraction-collisions-randomize* 0.4)
(defparameter *z-space* 1.5)
(defparameter *z-height* 10.0)
(defparameter *mutate-lipid-fraction* 0.08)
(defparameter *mutate-position-fraction* 0.3)
(defparameter *population-size* 200)
(defparameter *lipid-density* 0.0098) ; 0.0098)
(defparameter *random-tip-angle-degrees* 10.0)
(defparameter *lipid-radius* 5.2)
(defparameter *lipid-vdw-radius-nanometers* 0.62)
(defparameter *lipid-vdw-edep* 1.0)
(defparameter *close-distance* 1.8) ; 1.8)
(defparameter *max-close-distance* 3.0)
(defparameter *close-packed-circle-density* 0.9
  "The maximum density that equal circles can be packed in a plane")
(defparameter *membrane-force-field-name* :system-membrane)


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

(defun setup-system-membrane-force-field ()
  (let ((force-field (chem:make-force-field))
        (nonbond-force-field (chem:make-ffnonbond-db))
        (nonbond (chem:make-ffnonbond :lipid0
                                      :radius-nanometers *lipid-vdw-radius-nanometers*
                                      :epsilon-kj *lipid-vdw-edep*
                                      :mass 10.0
                                      )))
    (chem:ffnonbond-db-add nonbond-force-field nonbond)
    (chem:set-ffnonbond-db force-field nonbond-force-field)
    (leap.core:add-force-field-or-modification force-field :force-field-name *membrane-force-field-name*
                                                           :combined-force-field-class-name 'chem:combined-force-field)))

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
   (adjacent-to-solute :initform nil :initarg :adjacent-to-solute :accessor adjacent-to-solute)))


(defclass ga-lipid (ga-thing)
  ((array-index :initarg :array-index :accessor array-index)
   (lipid-info :initarg :lipid-info :accessor lipid-info)
   (lipid-index :initarg :lipid-index :accessor lipid-index)
   (lipid-id :initarg :lipid-id :accessor lipid-id)
   (top-bottom :initarg :top-bottom :accessor top-bottom)
   ))

(defun lipid-center (lipid)
  (let* ((lipid-info-ai (lipid-info lipid))
         (lipid-index-ai (lipid-index lipid))
         (octree-ai (elt (octrees lipid-info-ai) lipid-index-ai))
         (transform-ai (transform lipid))
         (origin-ai (chem:generic-octree-origin octree-ai))
         (coord-ai (geom:m*v transform-ai origin-ai)))
    coord-ai))

;;; Dummy values for cells that contain solute
(defclass ga-solute (ga-thing)
  ((octree :initarg :octree :accessor octree)
   (solute :initarg :solute :accessor solute)
   (array-index :initform nil :initarg :array-index :accessor array-index)
   (lipid-info :initform nil :initarg :lipid-info :accessor lipid-info)
   (lipid-index :initform nil :initarg :lipid-index :accessor lipid-index)
   (lipid-id :initform nil :initarg :lipid-id :accessor lipid-id)
   ))

(defclass scored-membrane ()
  ((score :initarg :score :accessor score)
   (membrane :initarg :membrane :accessor membrane)
   (collisions :initarg :collisions :accessor collisions)))


(defclass atm-res-mol ()
  ((atm :initarg :atm :accessor atm)
   (res :initarg :res :accessor res)
   (mol :initarg :mol :accessor mol)))

(defmethod print-object ((obj atm-res-mol) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a ~a" (chem:get-name (atm obj)) (chem:get-name (res obj)) (chem:get-name (mol obj)))))

(defclass atom-vectors ()
  ((thing-array :initarg :thing-array :accessor thing-array)
   (solute-atom-vector-index :initarg :solute-atom-vector-index :accessor solute-atom-vector-index)
   (atom-vectors :initarg :atom-vectors :accessor atom-vectors)
   (bounding-box :initarg :bounding-box :accessor bounding-box)))


(defmethod print-object ((obj ga-lipid) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream ":id ~a :lipid-index ~a" (lipid-id obj) (lipid-index obj))))


(defgeneric shallow-copy-ga-object (original))

(defmethod shallow-copy-ga-object ((original ga-solute))
  ;; don't need to copy this
  original)

(defmethod shallow-copy-ga-object ((original ga-lipid))
  (make-instance 'ga-lipid
                 :array-index (array-index original)
                 :lipid-id (lipid-id original)
                 :lipid-info (lipid-info original)
                 :lipid-index (lipid-index original)
                 :transform (transform original)
                 :adjacent-to-solute (adjacent-to-solute original)))

(defclass range ()
  ((source-ga-thing :initarg :source-ga-thing :accessor source-ga-thing)
   (target-ga-thing :initarg :target-ga-thing :accessor target-ga-thing)
   (start :initarg :start :accessor start)
   (end :initarg :end :accessor end)))

(defmethod print-object ((obj range) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream ":source ~a :target ~a :start ~a :end ~a" (source-ga-thing obj) (target-ga-thing obj)
            (start obj) (end obj))))

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
                               (geom:v+ (geom:v* (geom:vec 1.0 1.0 1.0) *max-close-distance*)
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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Calculate the number of lipid molecules that can
;;; fit in the bounding box and avoiding the protein.
;;;
;;; This creates a bitmap and paints it with the solute.
;;; Then it treats lipids as a model sphere with about 5.2A radius.
;;; Then it drags a lipid circle through the box and if the
;;; lipid circle doesn't contain any solute points - then paint
;;; the lipid circle with 1.
;;; Then add up the area of the bits in the lipid space and
;;; multiply it by the fraction of space filled by close-packed
;;; spheres on a hexagonal grid. The number of lipids is then
;;; estimated by dividing this area by the area of a lipid.

(defclass bit-field ()
  ((bit-array :initarg :bit-array :accessor bit-array)
   (bit-spacing :initarg :bit-spacing :accessor bit-spacing)
   (min-corner :initarg :min-corner :accessor min-corner)
   (bounding-box :initarg :bounding-box :accessor bounding-box)
   ))

(defun z-height (bit-field)
  (geom:vz (min-corner bit-field)))

(defun max-corner (bit-field)
  (let* ((min-corner (min-corner bit-field))
         (dims (array-dimensions (bit-array bit-field)))
         (xdim (first dims))
         (ydim (second dims))
         (diagonal (geom:v* (geom:vec xdim ydim 0.0) (bit-spacing bit-field))))
    (geom:v+ min-corner diagonal)))
    
(defun make-bit-field (bounding-box z-height bit-spacing)
  (let* ((center (chem:get-bounding-box-center bounding-box))
         (widths (chem:get-bounding-box-widths bounding-box))
         (half-width (geom:v* widths 0.5))
         (min-corner-bb (geom:v- center half-width))
         (min-corner (geom:vec (geom:vx min-corner-bb) (geom:vy min-corner-bb) z-height))
         (bit-array (make-array (list (round (/ (geom:vx widths) bit-spacing))
                                      (round (/ (geom:vy widths) bit-spacing)))
                                :element-type 'bit :initial-element 0)))
    (make-instance 'bit-field
                   :bit-array bit-array
                   :bit-spacing bit-spacing
                   :min-corner min-corner)))

(defun count-circle-bits (bit-field position radius &optional (value 1))
  "Return the count of the number of bits that are on in the circle"
  (let* ((min-corner (min-corner bit-field))
         (bit-spacing (bit-spacing bit-field))
         (iradius (round (/ radius bit-spacing)))
         (iradius-sq (* iradius iradius))
         (x (geom:vx position))
         (y (geom:vy position))
         (x-start-index (round (/ (- x (geom:vx min-corner) radius) bit-spacing)))
         (y-start-index (round (/ (- y (geom:vy min-corner) radius) bit-spacing)))
         (bit-array (bit-array bit-field))
         (dims (array-dimensions bit-array))
         (xdim (elt dims 0))
         (ydim (elt dims 1))
         (count 0))
    (loop for xi from x-start-index
          for xir from (- iradius) to iradius
          do (loop for yi from y-start-index
                   for yir from (- iradius) to iradius
                   when (< (+ (* xir xir) (* yir yir)) iradius-sq)
                     do (let ((xi (mod xi xdim))
                              (yi (mod yi ydim)))
                          (incf count (if (= value (aref bit-array xi yi)) 1 0)))))
    count))

(defun paint-circle-in-bit-field (bit-field position radius &optional (value 1))
  "Fill in the circle"
  (let* ((min-corner (min-corner bit-field))
         (bit-spacing (bit-spacing bit-field))
         (iradius (round (/ radius bit-spacing)))
         (iradius-sq (* iradius iradius))
         (x (geom:vx position))
         (y (geom:vy position))
         (x-start-index (round (/ (- x (geom:vx min-corner) radius) bit-spacing)))
         (y-start-index (round (/ (- y (geom:vy min-corner) radius) bit-spacing)))
         (bit-array (bit-array bit-field))
         (dims (array-dimensions bit-array))
         (xdim (elt dims 0))
         (ydim (elt dims 1)))
    (loop for xi from x-start-index
          for xir from (- iradius) to iradius
          do (loop for yi from y-start-index
                   for yir from (- iradius) to iradius
                   when (< (+ (* xir xir) (* yir yir)) iradius-sq)
                     do (let ((xi (mod xi xdim))
                              (yi (mod yi ydim)))
                          (setf (aref bit-array xi yi) value))))))

(defun dump-bit-field (bit-field)
  (loop for y below (second (array-dimensions (bit-array bit-field)))
        do (loop for x below (first (array-dimensions (bit-array bit-field)))
                 if (= (aref (bit-array bit-field) x y) 1)
                   do (format t "*")
                 else
                   do (format t "."))
        do (terpri)))

(defun paint-sphere-in-bit-field (bit-field position radius)
  (let* ((z-height (geom:vz (min-corner bit-field)))
         (z-delta (abs (- (geom:vz position) z-height))))
    (when (< z-delta radius)
      (let* ((angle (acos (/ z-delta radius)))
             (circle-radius (* radius (sin angle))))
        (paint-circle-in-bit-field bit-field position circle-radius)))))
   
(defun paint-protein-at-height (bounding-box solute z-height &key (resolution 0.5))
  (let* ((widths (chem:get-bounding-box-widths bounding-box))
         (xbits (round (/ (geom:vx widths) resolution)))
         (ybits (round (/ (geom:vy widths) resolution)))
         (min-corner (chem:min-corner bounding-box))
         (zcorner (geom:vec (geom:vx min-corner) (geom:vy min-corner) z-height))
         (bit-field (make-instance 'bit-field
                                   :bit-spacing resolution
                                   :min-corner zcorner
                                   :bounding-box bounding-box
                                   :bit-array (make-array (list xbits ybits)
                                                          :element-type 'bit
                                                          :initial-element 0)))
         (total 0)
         (intersect 0))
    (cando:do-atoms (atm solute)
      (incf total)
      (let* ((pos (chem:get-position atm))
             (element (chem:get-element atm))
             (radius (+ (chem:vdw-radius-for-element element) (chem:vdw-radius-for-element :H)))
             (from-z (abs (- (geom:vz pos) z-height))))
        (when (< (abs (- (geom:vz pos) z-height)) radius)
          ;; sphere intercepts bitmap
          (incf intersect)
          (paint-sphere-in-bit-field bit-field pos radius))))
    bit-field))

(defun excluded-lipid (bit-field)
  (let* ((new-bit-field (make-bit-field (bounding-box bit-field) (z-height bit-field) (bit-spacing bit-field)))
         (dims (array-dimensions (bit-array bit-field)))
         (xdim (elt dims 0))
         (ydim (elt dims 1))
         (min-corner (min-corner bit-field))
         (bit-spacing (bit-spacing bit-field))
         (xdir (geom:vec bit-spacing 0.0 0.0))
         (ydir (geom:vec 0.0 bit-spacing 0.0)))
    (loop for xi below xdim
          do (loop for yi below ydim
                   for posx = (geom:v+ min-corner (geom:v* xdir xi))
                   for pos = (geom:v+ posx (geom:v* ydir yi))
                   for count = (count-circle-bits bit-field pos *lipid-radius*)
                   when (= count 0)
                     do (paint-circle-in-bit-field new-bit-field pos *lipid-radius*)))
    new-bit-field))

(defun calculate-how-many-lipids-fit (bounding-box solute z-height &key (resolution 0.5))
  (let* ((protein-bit-field (paint-protein-at-height bounding-box solute z-height))
         (lipid-bit-field (excluded-lipid protein-bit-field))
         (bit-array (bit-array lipid-bit-field))
         (lipid-bits (loop for index below (array-total-size bit-array)
                           when (= (row-major-aref bit-array index) 1)
                             count 1))
         (lipid-area (* lipid-bits resolution resolution))
         (one-lipid-area (* PI *lipid-radius* *lipid-radius*))
         (num-lipids (round (/ (* lipid-area *close-packed-circle-density*) one-lipid-area))))
    (format t "one-lipid-area -> ~a  lipid-area -> ~a~%" one-lipid-area lipid-area)
    num-lipids))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Carry out a molecular dynamics simulation of lipids modeled as
;;; spheres avoiding the protein.
;;;

(defun random-aggregate (num)
  (let ((agg (chem:make-aggregate))
        (bbwidth 20.0)
        (width 10.0))
    (chem:set-bounding-box agg (chem:make-bounding-box (list bbwidth bbwidth bbwidth)))
    (loop for index from 0 below num
          for name = (intern (format nil "L~a" index) :keyword)
          for mol = (chem:make-molecule name)
          for res = (chem:make-residue name)
          for atm = (chem:make-atom name :|Ar|)
          for pos = (geom:vec (- (random width) (* width 0.5))
                              (- (random width) (* width 0.5))
                              (- (random width) (* width 0.5)))
          do (chem:setf-force-field-name mol :system-membrane)
          do (chem:set-position atm pos)
          do (chem:set-type atm :lipid0)
          do (chem:add-matter res atm)
          do (chem:add-matter mol res)
          do (chem:add-matter agg mol))
    agg))


(defun look-for-close-contacts (dynamics frame)
  (let* ((distances nil)
         (coordinates (elt (dynamics:accumulated-coordinates dynamics)
                           (if frame
                               frame
                               (1- (length (dynamics:accumulated-coordinates dynamics)))))))
    (loop for i1 below (/ (length coordinates) 3)
          for ic1 = (* 3 i1)
          for a1 = (chem:elt-atom (chem:atom-table (dynamics:scoring-function dynamics)) i1)
          for v1 = (geom:vec (elt coordinates ic1)
                             (elt coordinates (+ ic1 1))
                             (elt coordinates (+ ic1 2)))
          do (loop for i2 from (1+ i1) below (/ (length coordinates) 3)
                   for ic2 = (* 3 i2)
                   for a2 = (chem:elt-atom (chem:atom-table (dynamics:scoring-function dynamics)) i2)
                   for v2 = (geom:vec (elt coordinates ic2)
                                      (elt coordinates (+ ic2 1))
                                      (elt coordinates (+ ic2 2)))
                   do (let* ((vdelta (geom:v- v1 v2))
                             (vlen (geom:vlength vdelta)))
                        (push (list vlen a1 a2) distances))))
    (sort distances #'< :key #'car)))
    
                          
   
(defun show-well (agg energy &optional index)
  (let* ((mol0 (cando:mol agg 1))
         (res0 (cando:res mol0 0))
         (atm0 (cando:atm res0 0))
         (mol1 (cando:mol agg 0))
         (res1 (cando:res mol1 0))
         (atm1 (cando:atm res1 0)))
    (when index (chem:set-position atm0 (geom:vec index 0 0)))
    (format t "~a ~a~%" (abs (- (geom:vx (chem:get-position atm0)) (geom:vx (chem:get-position atm1)))) (chem:calculate-energy-and-force energy))))

(defun build-simple-membrane-aggregate (membrane z-height top-bottom &optional (lipid-radius *lipid-radius*))
  (let* ((solute (solute (ga-solute membrane)))
         (agg (chem:make-aggregate)))
    (cando:do-atoms (atm solute)
      (let* ((zpos (geom:vz (chem:get-position atm)))
             (zdelta (abs (- zpos z-height))))
        (when (< zdelta lipid-radius)
          (let ((natm (chem:make-atom (chem:get-name atm) (chem:get-element atm)))
                (res (chem:make-residue :solute))
                (mol (chem:make-molecule :solute)))
            (chem:set-position natm (chem:get-position atm))
            (chem:set-type natm (chem:get-type atm))
            (chem:add-matter res natm)
            (chem:add-matter mol res)
            (chem:add-matter agg mol)))))
    (loop for lipid across (lipids membrane)
          for index from 0
          when (eq (top-bottom lipid) top-bottom)
            do (let* ((center (lipid-center lipid))
                      (name (intern (format nil "LIPID~a" index) :keyword))
                      (atm (chem:make-atom name :|Ar|))
                      (res (chem:make-residue name))
                      (mol (chem:make-molecule name)))
                 ;;(chem:set-position atm center)
                 (chem:set-position atm (geom:v+ (geom:vec 0.0 0.0 0.0) (geom:vec (random 40.0) (random 40.0) z-height)))
                 (chem:set-type atm :lipid0)
                 (chem:setf-force-field-name mol *membrane-force-field-name*)
                 (chem:add-matter res atm)
                 (chem:add-matter mol res)
                 (chem:add-matter agg mol)))
    (chem:aggregate-set-bounding-box agg (bounding-box membrane))
    agg))


(defun build-simple-aggregate-energy-function (simple-aggregate &optional (lipid-radius *lipid-radius*))
  (cando:do-atoms (atm simple-aggregate)
    (unless (chem:get-type atm)
      (format t "Could not get type for ~a~%" atm)))
  (let* ((energy-function (chem:make-energy-function :matter simple-aggregate
                                                     :assign-types nil
                                                     :use-excluded-atoms t
                                                     :bounding-box (chem:bounding-box simple-aggregate)))
         (atom-table (chem:atom-table energy-function)))
    (chem:set-nonbond-cutoff (chem:get-nonbond-component energy-function) 16.0)
    (chem:set-bounding-box energy-function (chem:bounding-box simple-aggregate))
    (format t "atom-table -> ~a~%" atom-table)
    (format t "First atom: ~a~%" (chem:elt-atom atom-table 0))
    (let* ((num-atoms (chem:get-number-of-atoms atom-table))
           (frozen (make-array (* num-atoms 3) :element-type 'bit :initial-element 0)))
      (loop for atom-index below num-atoms
            for coord-index from 0 by 3
            for atom = (chem:elt-atom atom-table atom-index)
            if (eq (chem:get-type atom) :lipid0)
              do (setf (elt frozen (+ 2 coord-index)) 1)
            else
              do (setf (elt frozen coord-index) 1
                       (elt frozen (+ coord-index 1)) 1
                       (elt frozen (+ coord-index 2)) 1))
      (values energy-function frozen))))

(defun build-simple-membrane-energy-function (membrane z-height top-bottom &optional (lipid-radius *lipid-radius*))
  (leap:assign-atom-types (solute (ga-solute membrane)))
  (let ((simple-aggregate (build-simple-membrane-aggregate membrane z-height top-bottom lipid-radius)))
    (build-simple-aggregate-energy-function simple-aggregate lipid-radius)))

(defun setup-simple-membrane-dynamics (energy-function bounding-box)
  (let ((dyn (dynamics:make-atomic-simulation energy-function
                                              :delta-t *1fs-time-units*
                                              :accumulate-coordinates t)))
    (chem:setf-velocity-scale energy-function 0.9 0.9 0.9)
    dyn))


(defun advance-simple-membrane-dynamics (dyn bounding-box &optional frozen)
  ;; Bring stuff back into the periodic box.
  (dynamics:with-dynamics (dyn)
    (let ((pos (geom:vec 0 0 0)))
      (dynamics:velocity-verlet-step dyn :frozen frozen)
      (loop for index below (length (dynamics:coordinates dyn)) by 3
            do (geom:vec-extract pos (dynamics:coordinates dyn) index)
            do (let ((new-pos (chem:normalize-position bounding-box pos)))
                 (geom:vec-put (dynamics:coordinates dyn) new-pos index))))))

(defun simulate-simple-membrane-energy-function (energy-function bounding-box frozen num)
  (let ((dyn (setup-simple-membrane-dynamics energy-function bounding-box)))
    ;; Bring stuff back into the periodic box.
    (dynamics:with-dynamics (dyn)
      (loop for x from 0 below num
            do (advance-simple-membrane-dynamics dyn bounding-box frozen)))
    dyn))

(defun simulate-simple-membrane (membrane num)
  (multiple-value-bind (energy-function frozen)
      (membrane-packer::build-simple-membrane-energy-function membrane 10 :top)
    (simulate-simple-membrane-energy-function energy-function (bounding-box membrane) frozen num)))
        

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
   (lipids :initarg :lipids :accessor lipids)
   (lipid-selector :initarg :lipid-selector :accessor lipid-selector)
   (ga-solute :initarg :ga-solute :accessor ga-solute)))

(defun copy-ga-membrane (original)
  (let ((copy (shallow-copy-object original)))
    (setf (lipids copy) (make-array (array-dimensions (lipids original))))
    (loop for index from 0 below (array-total-size (lipids copy))
          for ga-thing = (row-major-aref (lipids original) index)
          for lipid-copy = (progn
                             (unless ga-thing
                               (error "ga-thing is nil"))
                             (shallow-copy-ga-object ga-thing))
          do (setf (row-major-aref (lipids copy) index) lipid-copy))
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

(defun make-random-ga-lipid (array-index xypos z-offset lipid-selector lipid-id top-bottom)
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
                   :array-index array-index
                   :lipid-info lipid-info
                   :lipid-index lipid-index
                   :lipid-id lipid-id
                   :transform transform
                   :top-bottom top-bottom)))

(defun make-molecule-from-ga-lipid (ga-lipid &key debug)
  "Build a molecule from the ga-lipid.  If DEBUG is T then change the molecule name and residue names so we can pick out the lipid-index.
   Return the molecule and the vector of atoms that correspond to the molecule order."
  (let* ((lipid-info (lipid-info ga-lipid))
         (lipid-index (lipid-index ga-lipid))
         (lipid-id (lipid-id ga-lipid))
         (transform (transform ga-lipid))
         (coordinates (elt (coordinates lipid-info) lipid-index))
         (new-to-old-map (make-hash-table))
         (lipid-mol (chem:matter-copy (lipid-molecule lipid-info) new-to-old-map))
         (order-info (order-info lipid-info))
         (old-atom-to-order-index (make-hash-table)))
    (chem:set-name lipid-mol (intern (format nil "mol-~a" lipid-id) :keyword))
    ;; Build a map from old-atoms to index in order-info
    (loop for index from 0 below (length order-info)
          for residue-atom = (elt order-info index)
          for old-atom = (order-atom residue-atom)
          do (setf (gethash old-atom old-atom-to-order-index) index))
    ;; If we are debugging - then rename the lipids
    (when debug
      (let ((lipid-name (format nil "L~a" lipid-id)))
        (chem:set-name lipid-mol (intern lipid-name :keyword))
        (cando:do-residues (res lipid-mol)
          (chem:set-name res (intern (format nil "~a-~a" lipid-name (chem:get-name res)) :keyword)))))
    ;; Build a new atom order
    (let ((index 0)
          (atom-vec (make-array (chem:number-of-atoms lipid-mol))))
      (cando:do-residues (new-res lipid-mol)
        (cando:do-atoms (new-atom new-res)
          (let* ((old-atom (gethash new-atom new-to-old-map))
                 (old-index (gethash old-atom old-atom-to-order-index))
                 (old-coord-index (* 3 old-index))
                 (pos (geom:vec 0.0 0.0 0.0)))
            (geom:vec-extract-transformed pos coordinates old-coord-index transform)
            (chem:set-position new-atom pos)
            (setf (elt atom-vec (prog1 index (incf index))) (make-instance 'atm-res-mol :atm new-atom :res new-res :mol lipid-mol)))))
      (values lipid-mol atom-vec))))

(defun calculate-center-and-half-width (bounding-cuboid)
  (let ((center (geom:vec (* (+ (geom:get-min-x bounding-cuboid) (geom:get-max-x bounding-cuboid)) 0.5)
                          (* (+ (geom:get-min-y bounding-cuboid) (geom:get-max-y bounding-cuboid)) 0.5)
                          (* (+ (geom:get-min-z bounding-cuboid) (geom:get-max-z bounding-cuboid)) 0.5)))
        (half-width (geom:vec (* 0.5 (- (geom:get-max-x bounding-cuboid) (geom:get-min-x bounding-cuboid)))
                              (* 0.5 (- (geom:get-max-y bounding-cuboid) (geom:get-min-y bounding-cuboid)))
                              (* 0.5 (- (geom:get-max-z bounding-cuboid) (geom:get-min-z bounding-cuboid))))))
    (values center half-width)))


(defun wrapped-array-indices (xi yi xdim ydim)
  (values (mod xi xdim) (mod yi ydim)))

(defun wrapped-xypos (xi yi xstart ystart xdir ypdir ymdir xdim)
  (let* ((xdelta (geom:v* xdir (mod xi xdim)))
         (ydelta (geom:v+ (geom:v* (geom:v+ ypdir ymdir) (floor (/ yi 2)))
                          (geom:v* ypdir (if (oddp yi) 1 0))))
         (xypos (geom:v+ (geom:vec xstart ystart 0.0) (geom:v+ xdelta ydelta))))
    xypos))


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

(defun lipid-protein-overlaps-p (bounding-box ga-lipid ga-solute &optional (close-distance *close-distance*))
  (let ((collisions (make-collisions)))
    (do-score-cell bounding-box ga-lipid ga-solute collisions close-distance)
    (let ((solute-collisions (if (> (length (indices collisions)) 0)
                                 (length (indices collisions))
                                 nil)))
      solute-collisions)))


(defun simulate-crude-lipids (bounding-box lipids solute &key top-bottom z-height)
  (let ((agg (chem:make-aggregate)))
    (loop for index below (length lipids)
          for lipid = (elt lipids index)
          for transform = (transform lipid)
          for position = (geom:vec 0.0 0.0 0.0)
          for transformed-position = (geom:m*v transform position)
          for mol = (chem:make-molecule)
          for name = (intern (format nil "DU~a" (lipid-index lipid)) :keyword)
          for res = (chem:make-residue name)
          for atom = (chem:make-atom name :du)
          do (chem:add-matter agg mol)
          do (chem:add-matter mol res)
          do (chem:add-matter res atom)
          do (chem:set-position atom transformed-position)
          do (chem:set-type atom :lipid))
    agg))


(defun fill-in-lipids (lipids solute-colliders solute top-bottom &key z-height expected-lipids)
  (let ((sorted-solute-colliders (sort solute-colliders #'< :key #'car)))
    (loop for solute-collider in sorted-solute-colliders
          for lipid-index from (length lipids) below expected-lipids
          for collisions = (car solute-collider)
          for ga-lipid = (cdr solute-collider)
          do (progn
               (format t "Adding ")
               (push ga-lipid lipids))
          do (format t "#collisions ~a  for solute ~a~%" collisions (lipid-id ga-lipid))))
  lipids)

(defparameter *top-agg* nil)

(defun build-ga-membrane (solute &key input-bounding-box
                                   (lipid-selector (list (cons 1.0 *popc*)))
                                   (z-height *z-height*))
  (let (top-lipids bottom-lipids)
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
                 (ga-solute (make-instance 'ga-solute :solute solute :octree solute-octree)) ;;  :solute solute :octree solute-octree))
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
            (let (top-solute-colliders
                  bottom-solute-colliders)
              (loop for yindex from 0 below ynum
                    do (loop for xindex from 0 below xnum
                             for wrapped-xypos = (wrapped-xypos xindex yindex xstart ystart xdir ypdir ymdir xnum)
                             for top-ga-lipid = (make-random-ga-lipid (list xindex yindex 0) wrapped-xypos z-offset lipid-selector (incf lipid-id) :top)
                             for bottom-ga-lipid = (make-random-ga-lipid (list xindex yindex 1) wrapped-xypos z-offset lipid-selector (incf lipid-id) :bottom)
                             do (let ((top-solute-collides-p (lipid-protein-overlaps-p bounding-box top-ga-lipid ga-solute)))
                                  (if top-solute-collides-p
                                      (progn
                                        (push (cons top-solute-collides-p top-ga-lipid) top-solute-colliders))
                                      (push top-ga-lipid top-lipids)))
                             do (let ((bottom-solute-collides-p (lipid-protein-overlaps-p bounding-box bottom-ga-lipid ga-solute)))
                                  (if bottom-solute-collides-p
                                      (progn
                                        (push (cons bottom-solute-collides-p bottom-ga-lipid) bottom-solute-colliders)
                                        )
                                      (push bottom-ga-lipid bottom-lipids)))))
              (format t "There are ~a top-lipids~%" (length top-lipids))
              (let ((total-top-lipids (calculate-how-many-lipids-fit bounding-box (solute ga-solute) z-height)))
                (setf top-lipids (fill-in-lipids top-lipids top-solute-colliders solute :top
                                             :z-height z-height
                                             :expected-lipids total-top-lipids))
                (format t "    there should be ~a~%" total-top-lipids))
              (format t "There are ~a bottom-lipids~%" (length bottom-lipids))
              (let ((total-bottom-lipids (calculate-how-many-lipids-fit bounding-box (solute ga-solute) (- z-height))))
                (setf bottom-lipids (fill-in-lipids bottom-lipids bottom-solute-colliders solute :bottom
                                             :z-height z-height
                                             :expected-lipids total-bottom-lipids))
                (format t "    there should be ~a~%" total-bottom-lipids))
              (setf *top-agg* (simulate-crude-lipids bounding-box top-lipids solute :top-bottom :top :z-height z-height))
              (simulate-crude-lipids bounding-box bottom-lipids solute :top-bottom :bottom :z-height z-height))
            (let ((lipids (make-array (+ (length top-lipids) (length bottom-lipids))
                                      :initial-contents (append top-lipids bottom-lipids))))
              (make-instance 'ga-membrane
                             :bounding-box bounding-box
                             :lipids lipids
                             :lipid-selector lipid-selector
                             :ga-solute ga-solute))))))))

(defun do-score-cell (bounding-box source target collisions close-distance)
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
                  close-distance
                  coord-index
                  source-pos
                  (indices collisions)
                  :bounding-box bounding-box
                  :octree-transform from-target-transform))))
    (let ((range-end (length (indices collisions))))
      ;; (when (> range-end range-start)
      (vector-push-extend (make-instance 'range
                                         :source-ga-thing source
                                         :target-ga-thing target
                                         :start range-start
                                         :end range-end)
                          (ranges collisions))
      (/ (- range-end range-start) 2))))

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

(defun octree-info-for-lipid-id (ga-membrane lipid-id)
  (let* ((lipids (lipids ga-membrane))
         (bounding-box (bounding-box ga-membrane))
         (lipid-pos (position-if (lambda (lipid) (= (lipid-id lipid) lipid-id)) (lipids ga-membrane)))
         (lipid (elt (lipids ga-membrane) lipid-pos))
         (lipid-info-ai (lipid-info lipid))
         (lipid-index-ai (lipid-index lipid))
         (octree-ai (elt (octrees lipid-info-ai) lipid-index-ai))
         (transform-ai (transform lipid))
         (origin-ai (chem:generic-octree-origin octree-ai))
         (coord-ai (geom:m*v transform-ai origin-ai)))
    (format t "Octree center: ~a     half-dimensions: ~a~%" coord-ai (chem:generic-octree-half-dimension octree-ai))))

(defun overlapping-octrees-p (ga-membrane lipid-id-ai lipid-id-bi)
  (let* ((lipids (lipids ga-membrane))
         (bounding-box (bounding-box ga-membrane))
         (lipid-pos (position-if (lambda (lipid) (= (lipid-id lipid) lipid-id-ai)) (lipids ga-membrane)))
         (lipid (elt (lipids ga-membrane) lipid-pos))
         (lipid-info-ai (lipid-info lipid))
         (lipid-index-ai (lipid-index lipid))
         (octree-ai (elt (octrees lipid-info-ai) lipid-index-ai))
         (transform-ai (transform lipid))
         (origin-ai (chem:generic-octree-origin octree-ai))
         (coord-ai (geom:m*v transform-ai origin-ai))
         (shapes (make-array 128 :adjustable t :fill-pointer 0))
         (color (vector 0.0 1.0 1.0)))
    (let* ((bi (position-if (lambda (lipid) (= (lipid-id lipid) lipid-id-bi)) (lipids ga-membrane)))
           (sc-bi (elt lipids bi))
           (lipid-info-bi (lipid-info sc-bi))
           (lipid-index-bi (lipid-index sc-bi))
           (octree-bi (elt (octrees lipid-info-bi) lipid-index-bi))
           (transform-bi (transform sc-bi))
           (origin-bi (chem:generic-octree-origin octree-bi))
           (coord-bi (geom:m*v transform-bi origin-bi)))
      (chem:octree-transformed-intersects-in-bounding-box octree-ai transform-ai octree-bi transform-bi bounding-box))))
  
(defun overlapping-octrees-as-shape (ga-membrane lipid-id &key (radius 0.5))
  (let* ((lipids (lipids ga-membrane))
         (bounding-box (bounding-box ga-membrane))
         (lipid-pos (position-if (lambda (lipid) (= (lipid-id lipid) lipid-id)) (lipids ga-membrane)))
         (lipid (elt (lipids ga-membrane) lipid-pos))
         (lipid-info-ai (lipid-info lipid))
         (lipid-index-ai (lipid-index lipid))
         (octree-ai (elt (octrees lipid-info-ai) lipid-index-ai))
         (transform-ai (transform lipid))
         (origin-ai (chem:generic-octree-origin octree-ai))
         (coord-ai (geom:m*v transform-ai origin-ai))
         (shapes (make-array 128 :adjustable t :fill-pointer 0))
         (color (vector 0.0 1.0 1.0)))
    (loop for bi from 0 below (length lipids)
          for sc-bi = (elt lipids bi)
          for lipid-info-bi = (lipid-info sc-bi)
          for lipid-index-bi = (lipid-index sc-bi)
          for octree-bi = (elt (octrees lipid-info-bi) lipid-index-bi)
          for transform-bi = (transform sc-bi)
          for origin-bi = (chem:generic-octree-origin octree-bi)
          for coord-bi = (geom:m*v transform-bi origin-bi)
          Unless (= (lipid-id sc-bi) lipid-id)
            do (when (chem:octree-transformed-intersects-in-bounding-box octree-ai transform-ai octree-bi transform-bi bounding-box)
                 (vector-push-extend
                  (vector "cylinder"
                          (vector (geom:vx coord-ai) (geom:vy coord-ai) (geom:vz coord-ai))
                          (vector (geom:vx coord-bi) (geom:vy coord-bi) (geom:vz coord-bi))
                          color
                          radius)
                  shapes)))
    shapes))



(defun score-ga-membrane (ga-membrane &key debug (close-distance *close-distance*))
  (let* ((lipids (lipids ga-membrane))
         (bounding-box (bounding-box ga-membrane))
         (collisions (make-collisions))
         (scores 0)
         (lipid-scores 0)
         (solute-scores 0))
    (loop for ai below (1- (length lipids))
          for sc-ai = (elt lipids ai)
          for lipid-info-ai = (lipid-info sc-ai)
          for lipid-index-ai = (lipid-index sc-ai)
          for octree-ai = (elt (octrees lipid-info-ai) lipid-index-ai)
          for transform-ai = (transform sc-ai)
          do (loop for bi from (1+ ai) below (length lipids)
                   for sc-bi = (elt lipids bi)
                   for lipid-info-bi = (lipid-info sc-bi)
                   for lipid-index-bi = (lipid-index sc-bi)
                   for octree-bi = (elt (octrees lipid-info-bi) lipid-index-bi)
                   for transform-bi = (transform sc-bi)
                   do (incf scores)
                   do (when (chem:octree-transformed-intersects-in-bounding-box octree-ai transform-ai octree-bi transform-bi bounding-box)
                        (incf lipid-scores)
                        (do-score-cell bounding-box sc-ai sc-bi collisions close-distance)))
          do (when (chem:octree-transformed-intersects-in-bounding-box octree-ai transform-ai (octree (ga-solute ga-membrane)) (transform (ga-solute ga-membrane)) bounding-box)
               (incf solute-scores)
               (let ((num (do-score-cell bounding-box sc-ai (ga-solute ga-membrane) collisions close-distance))))))
    (values (/ (length (indices collisions)) 2) collisions)))

(defun mutate-all-lipids (ga-membrane)
  "Mutate the lipid in every position in the membrane"
  (loop for row-major-index from 0 below (array-total-size (lipids ga-membrane))
        for new-lipid-info-index = (select-lipid (lipid-selector ga-membrane))
        for new-lipid-info = (car new-lipid-info-index)
        for new-lipid-index = (cdr new-lipid-info-index)
        for lipid-or-solute = (row-major-aref (lipids ga-membrane) row-major-index)
        do (unless (typep lipid-or-solute 'ga-solute)
             (setf (lipid-info lipid-or-solute) new-lipid-info
                   (lipid-index lipid-or-solute) new-lipid-index))))


(defun mutate-all-lipid-orientations (ga-membrane)
  "Mutate the lipid in every position in the membrane"
  (loop for row-major-index from 0 below (array-total-size (lipids ga-membrane))
        for new-lipid-info-index = (select-lipid (lipid-selector ga-membrane))
        for new-lipid-info = (car new-lipid-info-index)
        for new-lipid-index = (cdr new-lipid-info-index)
        for lipid-or-solute = (row-major-aref (lipids ga-membrane) row-major-index)
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
                          for atom = (atm atom-res-mol)
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

(defun dump-scores (scored-membranes)
  (let ((step 0))
    (mapc (lambda (m) (format t "~4d " (score m))
            (when (= (mod (incf step) 30) 0)
              (terpri) (finish-output)))
          scored-membranes))
  nil)

(defun compare-two-membranes (membrane1 membrane2)
  (loop for index from 1 below (length (lipids membrane1))
        for ga-thing1 = (elt (lipids membrane1) index)
        for ga-thing2 = (elt (lipids membrane2) index)
        for lipid-index1 = (if (typep ga-thing1 'ga-lipid) (lipid-index ga-thing1) 0)
        for lipid-index2 = (if (typep ga-thing2 'ga-lipid) (lipid-index ga-thing2) 0)
        do (if (/= lipid-index1 lipid-index2)
               (format t "*")
               (format t "."))
        do (when (= (mod index 40) 0)
             (terpri)
             (finish-output))))

(defun check-aggregate-atom-vector-atom-order (aggregate atom-vectors)
  (let ((aggregate-atoms (make-array 1024 :fill-pointer 0 :adjustable t)))
    (cando:do-molecules (mol aggregate)
      (cando:do-residues (res mol)
        (cando:do-atoms (atm res)
          (vector-push-extend atm aggregate-atoms))))
    (let ((atom-vector-atoms (make-array 1024 :fill-pointer 0 :adjustable t)))
      (loop for atom-vector across (atom-vectors atom-vectors)
            do (loop for atm-res-mol across atom-vector
                     do (vector-push-extend (atm atm-res-mol) atom-vector-atoms)))
      (format t "Checking atom order~%")
      (let ((checked (loop named check
                           for agg-index below (length aggregate-atoms)
                           for av-index below (length atom-vector-atoms)
                           for agg-atom = (elt aggregate-atoms agg-index)
                           for av-atom = (elt atom-vector-atoms av-index)
                           unless (eq agg-atom av-atom)
                             do (progn
                                  (format t "Mismatch in atoms at index ~a~%" agg-index)
                                  (return-from check nil))
                           finally (return-from check t))))
        (format t "aggregate-atoms length -> ~a~%" (length aggregate-atoms))
        (format t "atom-vector-atoms length -> ~a~%" (length atom-vector-atoms))
        (when (/= (length aggregate-atoms) (length atom-vector-atoms))
          (format t "There is a mismatch in the number of atoms~%"))
        (format t "Done~%")))))


(defun build-aggregate-from-ga-membrane (ga-membrane &key debug)
  "Build an aggregate for the membrane and return it as well as a vector of vectors of atoms for
testing the scoring."
  (let ((solute (chem:matter-copy (solute (ga-solute ga-membrane)))))
    (chem:apply-transform-to-atoms solute (transform (ga-solute ga-membrane)))
    (let* ((aggregate (chem:make-aggregate)) ; make an empty aggregate
           (lipids (lipids ga-membrane))
           (index -1)
           (bounding-box (bounding-box ga-membrane))
           (lipid-molecules (make-array 1024 :adjustable t :fill-pointer 0))
           (atom-vectors (make-instance 'atom-vectors
                                        :thing-array (make-array (length lipids) :element-type t)
                                        :bounding-box bounding-box)))
      (chem:set-bounding-box aggregate (bounding-box ga-membrane))
      (let ((molecule-atoms (make-array 256 :fill-pointer 0 :adjustable t)))
        (loop for index below (length lipids)
              for ga-thing = (elt lipids index)
              unless (typep ga-thing 'ga-solute)
                do (let* ((ga-lipid ga-thing)
                          (transform (transform ga-lipid))
                          (lipid-info (lipid-info ga-lipid))
                          (lipid-index (lipid-index ga-lipid)))
                     (multiple-value-bind (lipid-mol lipid-atom-vector)
                         (make-molecule-from-ga-lipid ga-lipid :debug debug)
                       (chem:add-molecule aggregate lipid-mol)
                       (let ((lipid-mol-index (vector-push-extend lipid-mol lipid-molecules)))
                         (setf (elt (thing-array atom-vectors) index) lipid-mol-index))
                       (vector-push-extend lipid-atom-vector molecule-atoms))))
        ;; add in the solute atoms at the end
        (let ((solute-atoms (make-array (chem:number-of-atoms solute)))
              (si -1))
          (cando:do-molecules (mol solute)
            (let ((mol-list (list mol)))
              (chem:add-molecule aggregate mol)
              (cando:do-residues (res mol)
                (let ((res-list (list* res mol-list)))
                  (cando:do-atoms (atm res)
                    (let ((atm-list (list* atm res-list)))
                      (setf (elt solute-atoms (incf si)) (make-instance 'atm-res-mol :atm atm :res res :mol mol))))))))
          (setf (solute-atom-vector-index atom-vectors) (length molecule-atoms))
          (vector-push-extend solute-atoms molecule-atoms))
        (setf (atom-vectors atom-vectors) molecule-atoms)
        ;; (check-aggregate-atom-vector-atom-order aggregate atom-vectors)
        (values aggregate atom-vectors)))))

(defun simulate-membrane (aggregate num)
  (let* ((energy-function (chem:make-energy-function :matter aggregate))
         (dyn (dynamics:make-atomic-simulation energy-function
                                               :delta-t 0.01
                                               :accumulate-coordinates t)))
    (chem:setf-velocity-scale energy-function 0.1 0.1 0.1)
    (dynamics:with-dynamics (dyn)
      (loop for x from 0 below num
            do (format t "Dynamics step ~a~%" x)
            do (dynamics:velocity-verlet-step-limit-displacement dyn)))
    dyn))

(defun geometric-center-atom-vector (atom-vector)
  (unless (> (length atom-vector) 0)
    (error "Empty atom vector"))
  (let ((sum (geom:vec 0.0 0.0 0.0)))
    (loop for atom-res-mol across atom-vector
          for atom = (atm atom-res-mol)
          for pos = (chem:get-position atom)
          do (setf sum (geom:v+ sum pos)))
    (geom:v* sum (/ 1.0 (length atom-vector)))))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;;
;;;




(defclass rigid-body-dynamics ()
  ((membrane :initarg :membrane :accessor membrane)
   (trajectory :initarg :trajectory :accessor trajectory)
   (atom-vectors :initarg :atom-vectors :accessor atom-vectors)
   (energy-function :initarg :energy-function :accessor energy-function)
   (velocity-scale :initarg :velocity-scale :accessor velocity-scale)
   (limit-displacement :initarg :limit-displacement :accessor limit-displacement)
   (coordinates :initarg :coordinates :accessor coordinates)
   (velocity :initarg :velocity :accessor velocity)
   (force :initarg :force :accessor force)
   (force-dt :initarg :force-dt :accessor force-dt)
   (delta-t-over-mass :initarg :delta-t-over-mass :reader delta-t-over-mass)
   (frozen :initarg :frozen :accessor frozen)))


(defun build-rigid-body-energy-function (membrane &key debug)
  (multiple-value-bind (aggregate atom-vectors)
      (build-aggregate-from-ga-membrane membrane :debug debug)
    (let ((end-atom-vec (make-array (length (atom-vectors atom-vectors)) :element-type 'ext:byte32))
          (nonbond-db (chem:compute-merged-nonbond-force-field-for-aggregate aggregate))
          (centers (make-array (length (atom-vectors atom-vectors)))))
      (loop for index from 0 below (length (atom-vectors atom-vectors))
            for atom-vector = (elt (atom-vectors atom-vectors) index)
            for end-atom-index = (length atom-vector) then (+ end-atom-index (length atom-vector))
            do (setf (elt end-atom-vec index) end-atom-index)
            do (setf (elt centers index) (geometric-center-atom-vector atom-vector)))
      (let ((nonbond-component (chem:make-energy-rigid-body-nonbond end-atom-vec))
            (nb-index 0))
        (loop for index below (length (atom-vectors atom-vectors))
              for atom-vector = (elt (atom-vectors atom-vectors) index)
              for center = (elt centers index)
              for sum = (geom:vec 0.0 0.0 0.0)
              do (loop for ai from 0 below (length atom-vector)
                       for atom-res-mol = (elt atom-vector ai)
                       for atom = (atm atom-res-mol)
                       for pos = (geom:v- (chem:get-position atom) center)
                       for type = (chem:get-type atom)
                       for type-index = (chem:find-type-index nonbond-db type)
                       for ffnonbond = (chem:get-ffnonbond-using-type-index nonbond-db type-index)
                       for atom-radius = (chem:get-radius-angstroms ffnonbond)
                       for epsilon = (chem:get-epsilon-k-cal ffnonbond)
                       do (setf sum (geom:v+ sum pos))
                       do (chem:energy-rigid-body-nonbond-set-term
                           nonbond-component
                           (prog1 nb-index (incf nb-index))
                           atom
                           atom-radius  ; atom-radius
                           epsilon      ; epsilon
                           0.0          ; charge
                           pos))
              )
        (let ((energy (chem:make-rigid-body-energy-function (length (atom-vectors atom-vectors)) (bounding-box membrane))))
          (chem:rigid-body-energy-function-add-term energy nonbond-component)
          (let ((start-pos (make-array (chem:get-nvector-size energy) :element-type 'double-float)))
            (loop for index from 0 below (chem:get-nvector-size energy) by 7
                  for center-index from 0
                  for center = (elt centers center-index)
                  do (progn
                       (setf (elt start-pos (+ index 0)) 1.0)
                       (setf (elt start-pos (+ index 1)) 0.0)
                       (setf (elt start-pos (+ index 2)) 0.0)
                       (setf (elt start-pos (+ index 3)) 0.0)
                       (setf (elt start-pos (+ index 4)) (geom:vx center))
                       (setf (elt start-pos (+ index 5)) (geom:vy center))
                       (setf (elt start-pos (+ index 6)) (geom:vz center))))
            (chem:save-coordinates-from-vector energy start-pos)
            (values energy aggregate atom-vectors start-pos)))))))

(defun build-rigid-body-minimizer (membrane &key debug)
  (multiple-value-bind (energy aggregate atom-vectors start-pos)
      (build-rigid-body-energy-function membrane :debug debug)
    (let* ((frozen (make-array (* 7 (length (atom-vectors atom-vectors))) :element-type 'bit :initial-element 0))
           (minimizer (chem:make-minimizer energy))
           (solute-atom-vector-index (solute-atom-vector-index atom-vectors)))
      (chem:enable-print-intermediate-results minimizer 1 2)
      (chem:set-maximum-number-of-steepest-descent-steps minimizer 100)
      (chem:set-steepest-descent-tolerance minimizer 1.0)
      (chem:set-maximum-number-of-conjugate-gradient-steps minimizer 0)
      (chem:set-maximum-number-of-truncated-newton-steps minimizer 0)
      #+(or)
      (progn
        (loop for index from (* 7 solute-atom-vector-index)
                below (* 7 (1+ solute-atom-vector-index))
              do (setf (elt frozen index) 1)) ; freeze the solute
        (loop for index from 0 below (length frozen) by 7
              do (setf (elt frozen (+ 0 index)) 1 ; quaternion a
                       (elt frozen (+ 1 index)) 1 ; quaternion b
                       (elt frozen (+ 2 index)) 1 ; quaternion c
                       (elt frozen (+ 3 index)) 1 ; quaternion d
                       (elt frozen (+ 6 index)) 1 ; dir z
                       ))
        (chem:minimizer-set-frozen minimizer frozen))
      minimizer)))

(defun build-rigid-body-dynamics (membrane &key debug)
  (multiple-value-bind (energy aggregate atom-vectors start-pos)
      (build-rigid-body-energy-function membrane :debug debug)
    (let ((frozen (make-array (length (atom-vectors atom-vectors)) :element-type 'bit :initial-element 0))
          (trajectory (dynamics:make-trajectory aggregate)))
      (setf (elt frozen (solute-atom-vector-index atom-vectors)) 1)
      (make-instance 'rigid-body-dynamics
                     :membrane membrane
                     :atom-vectors atom-vectors
                     :trajectory trajectory
                     :energy-function energy
                     :velocity-scale (make-array 7 :element-type 'double-float :initial-contents '(0.00 0.00 0.00 0.00 0.1 0.1 0.0))
                     :limit-displacement (make-array 7 :element-type 'double-float :initial-contents '(0.05 0.05 0.05 0.05 0.1 0.1 0.0))
                     :coordinates start-pos
                     :velocity (make-array (chem:get-nvector-size energy) :element-type 'double-float :initial-element 0.0d0)
                     :force (make-array (chem:get-nvector-size energy) :element-type 'double-float :initial-element 0.0d0)
                     :force-dt (make-array (chem:get-nvector-size energy) :element-type 'double-float :initial-element 0.0d0)
                     :delta-t-over-mass (make-array (chem:get-nvector-size energy) :element-type 'double-float :initial-element 1.0d0)
                     :frozen frozen
                     ))))

(defun advance-rigid-body-dynamics (dynamics)
  (let* ((energy-function (energy-function dynamics))
         (rbcoords (coordinates dynamics)))
    (chem:rigid-body-velocity-verlet-step-limit-displacement
     energy-function                      ; scoring-func
     (coordinates dynamics)               ; position
     (velocity dynamics)                  ; velocity
     (force dynamics)                     ; force
     (force-dt dynamics)                  ; force-dt
     (delta-t-over-mass dynamics)         ; delta-t-over-mass
     1.0                                  ; delta-t
     (frozen dynamics)                    ; frozen
     (velocity-scale dynamics)            ; velocity-scale
     (limit-displacement dynamics))       ; limit-displacement
    (chem:rigid-body-energy-function-normalize-position energy-function rbcoords)
    (let* ((terms (chem:rigid-body-energy-function-terms energy-function))
           (nonbond-term (find-if (lambda (term)
                                    (let ((result (typep term 'chem:energy-rigid-body-nonbond)))
                                      result)) terms))
           (coords (make-array (* 3 (chem:energy-rigid-body-nonbond-number-of-nonbond-atoms nonbond-term))
                               :element-type 'single-float :adjustable t :fill-pointer 0)))
      (chem:write-nonbond-atom-coordinates-to-complex-vector-float nonbond-term rbcoords coords)
      (let ((copy-coords (copy-seq coords)))
        (vector-push-extend copy-coords (dynamics:coordinates (trajectory dynamics)))))))


(defun rigid-body-origins (num)
  (let ((coordinates (make-array (* 3 num) :initial-element 0.0 :element-type 'double-float))
        (indices (make-array num)))
    (loop for index from 0 below num
          for coord-index = (* 3 index)
          for coord-end = (* 3 (1+ index))
          do (setf (elt indices index) coord-end))
    (values coordinates indices)))

(defun rigid-body-origins-as-shape (coordinates)
  (let ((shapes (make-array 256 :fill-pointer 0 :adjustable t)))
    (loop for index3 from 0 below (length coordinates) by 3
          for xpos = (elt coordinates index3)
          for ypos = (elt coordinates (+ 1 index3))
          for zpos = (elt coordinates (+ 2 index3))
          for part = (vector "sphere"
                             (vector xpos ypos zpos)
                             (vector 1 0 0)
                             4.0)
          do (vector-push-extend part shapes))
    shapes))

(defun rigid-body-axes (num)
  (let ((coordinates (make-array (* 4 3 num) :initial-element 0.0 :element-type 'double-float))
        (indices (make-array num))
        (line-len 10.0))
    (loop for index from 0 below num
          for coord-index = (* 4 3 index)
          for coord-end = (* 4 3 (1+ index))
          do (setf (elt indices index) coord-end)
          do (setf (elt coordinates (+ 0 0 coord-index)) 0.0
                   (elt coordinates (+ 1 0 coord-index)) 0.0
                   (elt coordinates (+ 2 0 coord-index)) 0.0
                   (elt coordinates (+ 0 3 coord-index)) line-len
                   (elt coordinates (+ 1 3 coord-index)) 0.0
                   (elt coordinates (+ 2 3 coord-index)) 0.0
                   (elt coordinates (+ 0 6 coord-index)) 0.0
                   (elt coordinates (+ 1 6 coord-index)) line-len
                   (elt coordinates (+ 2 6 coord-index)) 0.0
                   (elt coordinates (+ 0 9 coord-index)) 0.0
                   (elt coordinates (+ 1 9 coord-index)) 0.0
                   (elt coordinates (+ 2 9 coord-index)) line-len))
    (values coordinates indices)))

(defun rigid-body-axes-as-shape (coordinates)
  (let ((shapes (make-array 256 :fill-pointer 0 :adjustable t))
        (arrow-radius 0.3))
    (loop for index3 from 0 below (length coordinates) by 12
          for xpos-origin = (elt coordinates (+ 0 0 index3))
          for ypos-origin = (elt coordinates (+ 1 0 index3))
          for zpos-origin = (elt coordinates (+ 2 0 index3))
          for xpos-xaxis = (elt coordinates (+ 0 3 index3))
          for ypos-xaxis = (elt coordinates (+ 1 3 index3))
          for zpos-xaxis = (elt coordinates (+ 2 3 index3))
          for xpos-yaxis = (elt coordinates (+ 0 6 index3))
          for ypos-yaxis = (elt coordinates (+ 1 6 index3))
          for zpos-yaxis = (elt coordinates (+ 2 6 index3))
          for xpos-zaxis = (elt coordinates (+ 0 9 index3))
          for ypos-zaxis = (elt coordinates (+ 1 9 index3))
          for zpos-zaxis = (elt coordinates (+ 2 9 index3))
          for part1 = (vector "arrow"
                              (vector xpos-origin ypos-origin zpos-origin)
                              (vector xpos-xaxis ypos-xaxis zpos-xaxis)
                              (vector 1 0 0)
                              arrow-radius)
          for part2 = (vector "arrow"
                              (vector xpos-origin ypos-origin zpos-origin)
                              (vector xpos-yaxis ypos-yaxis zpos-yaxis)
                              (vector 0 1 0)
                              arrow-radius)
          for part3 = (vector "arrow"
                              (vector xpos-origin ypos-origin zpos-origin)
                              (vector xpos-zaxis ypos-zaxis zpos-zaxis)
                              (vector 0 0 1)
                              arrow-radius)
          do (vector-push-extend part1 shapes)
          do (vector-push-extend part2 shapes)
          do (vector-push-extend part3 shapes))
    shapes))

(defun rigid-body-coordinates-as-shape (dynamics)
  (let* ((energy-function (energy-function dynamics))
         (rbcoords (coordinates dynamics))
         (shapes (make-array 256 :adjustable t :fill-pointer 0)))
    (loop for index from 0 below (length rbcoords) by 7
          for xpos = (elt rbcoords (+ index 4))
          for ypos = (elt rbcoords (+ index 5))
          for zpos = (elt rbcoords (+ index 6))
          for pos = (geom:vec xpos ypos zpos)
          for shape = (vector "sphere"
                              (vector xpos ypos zpos)
                              (vector 1 0 0)
                              3.0)
          do (vector-push-extend shape shapes))
    shapes))

(defun rigid-body-dynamics-origins (dynamics)
  (let* ((energy-function (energy-function dynamics))
         (rbcoords (coordinates dynamics))
         (origins (rigid-body-origins (/ (length rbcoords) 7)))
         (shapes (make-array 256 :adjustable t :fill-pointer 0)))
    (loop for index from 0 below (length rbcoords) by 7
          for xpos = (elt rbcoords (+ index 4))
          for ypos = (elt rbcoords (+ index 5))
          for zpos = (elt rbcoords (+ index 6))
          for pos = (geom:vec xpos ypos zpos)
          for shape = (vector "sphere"
                              (vector xpos ypos zpos)
                              (vector 1 0 0)
                              3.0)
          do (vector-push-extend shape shapes))
    shapes))

(defun rigid-body-dynamics-origins-as-shape-using-nonbond-code (dynamics)
  (let* ((energy-function (energy-function dynamics))
         (terms (chem:rigid-body-energy-function-terms energy-function))
         (nonbond-term (find-if (lambda (term)
                                  (let ((result (typep term 'chem:energy-rigid-body-nonbond)))
                                    result)) terms))
         (rbcoords (coordinates dynamics)))
    (multiple-value-bind (coordinates end-indices)
        (rigid-body-origins (/ (length rbcoords) 7))
      (let ((transformed-origins (make-array (length coordinates) :element-type 'single-float :fill-pointer 0 :adjustable t)))
        (chem:write-rigid-body-coordinates-to-complex-vector-float nonbond-term rbcoords end-indices coordinates transformed-origins)
        (rigid-body-origins-as-shape transformed-origins)))))

(defun rigid-body-dynamics-axes-as-shape-using-nonbond-code (dynamics)
  (let* ((energy-function (energy-function dynamics))
         (terms (chem:rigid-body-energy-function-terms energy-function))
         (nonbond-term (find-if (lambda (term)
                                  (let ((result (typep term 'chem:energy-rigid-body-nonbond)))
                                    result)) terms))
         (rbcoords (coordinates dynamics)))
    (multiple-value-bind (coordinates end-indices)
        (rigid-body-axes (/ (length rbcoords) 7))
      (let ((transformed-axes (make-array (length coordinates) :element-type 'single-float :fill-pointer 0 :adjustable t)))
        (chem:write-rigid-body-coordinates-to-complex-vector-float nonbond-term rbcoords end-indices coordinates transformed-axes)
        (rigid-body-axes-as-shape transformed-axes)))))

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

#+(or)
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
          for atom1 = (atm atom1-info)
          for atom2 = (atm atom2-info)
          for pos1 = (chem:get-position atom1)
          for pos2 = (chem:get-position atom2)
          for pair = (list pos1 pos2)
          do (vector-push-extend pair point-pairs))
    (add-point-pairs-as-shape widget point-pairs :color #(1 1 0) :radius radius :name "close-contacts")))


(defun close-contacts-between-atom-vectors (m1-atoms m2-atoms bounding-box close-contact-vec cutoff-squared num)
  (let ((close-contacts 0))
    (loop for i1 below (length m1-atoms)
          for atom1-res-mol = (elt m1-atoms i1)
          for atom1 = (atm atom1-res-mol)
          do (loop for i2 below (length m2-atoms)
                   for atom2-res-mol = (elt m2-atoms i2)
                   for atom2 = (atm atom2-res-mol)
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
  (loop for index below (length (lipids ga-membrane))
        for ga-thing = (elt (lipids ga-membrane) index)
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


(defclass membrane-strip ()
  ((dir-fn :initarg :dir-fn :accessor dir-fn)
   (min-val :initarg :min-val :accessor min-val)
   (max-val :initarg :max-val :accessor max-val)
   (inside :initarg :inside :accessor inside)
   (inside-counts :initform (make-hash-table) :initarg :inside-counts :accessor inside-counts)
   (outside :initarg :outside :accessor outside)
   (outside-counts :initform (make-hash-table) :initarg :outside-counts :accessor outside-counts)))

(defmethod print-object ((object membrane-strip) stream)
  (print-unreadable-object (object stream :type t)
    (format stream ":dir-fn ~a :min-val ~a :max-val ~a" (dir-fn object) (min-val object) (max-val object))))

(defun classify-strip (membrane selection-function dir-fn min-val max-val)
  (let (inside
        outside
        (inside-counts (make-array (length (lipid-selector membrane)) :initial-element 0))
        (outside-counts (make-array (length (lipid-selector membrane)) :initial-element 0)))
    (loop for one across (lipids membrane)
          for lipid-info = (lipid-info one)
          for count-index = (position lipid-info (lipid-selector membrane) :key #'cdr)
          if (funcall selection-function one)
            do (progn
                 (push one inside)
                 (incf (aref inside-counts count-index)))
          else
            do (progn
                 (push one outside)
                 (incf (aref outside-counts count-index))))
    (make-instance 'membrane-strip
                   :dir-fn dir-fn
                   :min-val min-val
                   :max-val max-val
                   :inside inside
                   :outside outside
                   :inside-counts inside-counts
                   :outside-counts outside-counts)))

(defun compatible-strips-p (strip1 strip2)
  (and (= (length (inside strip1)) (length (inside strip2)))
       (= (length (outside strip1)) (length (outside strip2)))
       (equalp (inside-counts strip1) (inside-counts strip2))
       (equalp (outside-counts strip1) (outside-counts strip2))))

(defun identify-strip (memb1 memb2)
  (let ((dir (random 2))
        (from-to (random 2))
        (bounding-box (bounding-box memb1))
        )
    (let ((min-corner (chem:min-corner bounding-box))
          (max-corner (chem:max-corner bounding-box))
          (bound0 (random 1.0))
          (bound1 (random 1.0)))
      (multiple-value-bind (bound0 bound1)
          (if (< bound0 bound1)
              (values bound0 bound1)
              (values bound1 bound0))
        (unless (or (> (- bound1 bound0) 0.3)
                    (< (- bound1 bound0) 0.7))
          (return-from identify-strip nil))
        (let* ((dir-fn (ecase dir
                         (0 'geom:vx)
                         (1 'geom:vy)))
               (val0 (+ (funcall dir-fn min-corner) (* bound0 (funcall dir-fn (chem:get-bounding-box-widths bounding-box)))))
               (val1 (+ (funcall dir-fn min-corner) (* bound1 (funcall dir-fn (chem:get-bounding-box-widths bounding-box))))))
          (let* ((selector
                   (flet ((lipid-position (lipid)
                            (let* ((lipid-info (lipid-info lipid))
                                   (lipid-index (lipid-index lipid))
                                   (transform (transform lipid))
                                   (octree (elt (octrees lipid-info) lipid-index))
                                   (origin (chem:generic-octree-origin octree))
                                   (coord (geom:m*v transform origin)))
                              coord)))
                     (lambda (lipid)
                       (let* ((pos (lipid-position lipid))
                              (val (funcall dir-fn pos)))
                         (and (>= val val0) (< val val1))))))
                 (strip1 (classify-strip memb1 selector dir-fn val0 val1))
                 (strip2 (classify-strip memb2 selector dir-fn val0 val1)))
            (when (compatible-strips-p strip1 strip2)
              (values strip1 strip2))))))))

(defun cross-two-membranes (memb1 memb2)
  (let ((cmemb1 (copy-ga-membrane memb1))
        (cmemb2 (copy-ga-membrane memb2)))
    (loop for index below 100
          do (multiple-value-bind (strip1 strip2)
                 (identify-strip cmemb1 cmemb2)
               (when strip1
                 (let ((lipids (if (= (random 2) 0)
                                   (coerce (append (inside strip1) (outside strip2)) 'vector)
                                   (coerce (append (outside strip1) (inside strip2)) 'vector))))
                   (setf (lipids cmemb1) lipids))
                 (return-from cross-two-membranes cmemb1)))))
  memb1)

(defun mutate-lipids (orig-membrane &optional (fraction-mutations 0.03))
  (let ((membrane (copy-ga-membrane orig-membrane))
        (num-mutations (floor (* fraction-mutations (length (lipids orig-membrane))))))
    (loop for num from 0 below num-mutations
          for index = (random (length (lipids orig-membrane)))
          for lipid-or-solute = (aref (lipids membrane) index)
          when (typep lipid-or-solute 'ga-lipid)
            do (let* ((new-lipid-info-index (select-lipid (lipid-selector orig-membrane)))
                      (new-lipid-info (car new-lipid-info-index))
                      (new-lipid-index (cdr new-lipid-info-index)))
                 (setf (lipid-info lipid-or-solute) new-lipid-info
                       (lipid-index lipid-or-solute) new-lipid-index)))
    membrane))

(defun mutate-position (orig-membrane &key (delta-x 0.04) (delta-y 0.04) (fraction-mutations 0.1))
  (let ((membrane (copy-ga-membrane orig-membrane))
        (num-mutations (floor (* fraction-mutations (length (lipids orig-membrane))))))
    (loop for num from 0 below num-mutations
          for index = (random (length (lipids orig-membrane)))
          for ga-thing = (aref (lipids membrane) index)
          when (typep ga-thing 'ga-lipid)
            do (let* ((lipid-shift (geom:vec (- (random (* 2.0 delta-x)) delta-x)
                                              (- (random (* 2.0 delta-y)) delta-y)
                                              0.0))
                      (extra-transform (geom:make-m4-translate lipid-shift))
                      (transform (geom:m*m extra-transform (transform ga-thing))))
                 (setf (transform ga-thing) transform)))
    membrane))


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

(defun score-membranes (membranes &key (parallel nil) sort (close-distance *close-distance*))
  (check-type membranes list)
  (let ((unsorted-scored-membranes
          (if parallel
              (progn
                (cando:lazy-setup-lparallel)
                (lparallel:pmap 'list
                                (lambda (memb)
                                  (multiple-value-bind (num-collisions collisions)
                                      (score-ga-membrane memb :close-distance close-distance)
                                    (make-instance 'scored-membrane
                                                   :score num-collisions
                                                   :membrane memb
                                                   :collisions collisions)))
                                membranes))
              (map 'list
                   (lambda (memb)
                     (multiple-value-bind (num-collisions collisions)
                         (score-ga-membrane memb :close-distance close-distance)
                       (make-instance 'scored-membrane
                                      :score num-collisions
                                      :membrane memb
                                      :collisions collisions)))
                   membranes))))
    (if sort
        (sort unsorted-scored-membranes #'< :key #'score)
        unsorted-scored-membranes)))

(defun mutate-collisions-one-membrane (original-membrane
                                       membrane-collisions
                                       &optional
                                         (fraction-collisions-randomize
                                          *fraction-collisions-randomize*))
  (let ((new-membrane (copy-ga-membrane original-membrane)))
    (loop for range across (ranges membrane-collisions) 
          for source = (source-ga-thing range)
          for target = (target-ga-thing range)
          for start = (start range)
          for end = (end range)
          when (and (/= start end) (< (random 1.0) *fraction-collisions-randomize*))
            do (progn
                 (when (typep source 'ga-lipid)
                   (let* ((source-lipid source)
                          (address (position source-lipid (lipids original-membrane)))
                          (lipid-info (lipid-info source-lipid))
                          (new-lipid (aref (lipids new-membrane) address)))
                     (setf (lipid-index new-lipid) (random (length (coordinates lipid-info))))))
                 (when (typep target 'ga-lipid)
                   (let* ((target-lipid source)
                          (address (position target-lipid (lipids original-membrane)))
                          (lipid-info (lipid-info target-lipid))
                          (new-lipid (aref (lipids new-membrane) address)))
                     (setf (lipid-index new-lipid) (random (length (coordinates lipid-info))))))))
    new-membrane))
  
(defun mutate-collisions (membranes collisions &key (fraction-collisions-randomize *fraction-collisions-randomize*))
  (loop for index below (length membranes)
        for membrane in membranes
        for membrane-collisions in collisions
        collect (mutate-collisions-one-membrane membrane membrane-collisions fraction-collisions-randomize)))

(defun cross-membranes (membranes &optional (select-fraction *select-fraction*))
  (loop with population-size = (length membranes)
        for index below population-size
        for parent1-index = (random (floor (* population-size select-fraction))) 
        for parent2-index = (random (floor (* population-size select-fraction))) 
        for parent1 = (elt membranes parent1-index)
        for parent2 = (elt membranes parent2-index)
        for mutant = (cross-two-membranes parent1 parent2) 
        collect mutant))


(defun randomly-mutate-membranes (membranes mutate-lipid-fraction)
  (loop for membrane in membranes
        for mutate-p = (< (random 1.0) mutate-lipid-fraction)
        collect (if mutate-p
                    (mutate-lipids membrane)
                    membrane)))

(defun randomly-shift-lipids-in-membranes (membranes mutate-position-fraction)
  (loop for membrane in membranes
        for mutate-p = (< (random 1.0) mutate-position-fraction)
        collect (if mutate-p
                    (mutate-position membrane)
                    membrane)))

(defun build-population (template-membrane &key (population-size *population-size*))
  (loop for index below population-size
        for membrane = (let ((copy-membrane (copy-ga-membrane template-membrane)))
                         (mutate-all-lipids copy-membrane)
                         copy-membrane)
        collect membrane))

(defun build-population-orientation-only (template-membrane &key (population-size *population-size*))
  (loop for index below population-size
        for membrane = (let ((copy-membrane (copy-ga-membrane template-membrane)))
                         (mutate-position copy-membrane)
                         copy-membrane)
        collect membrane))

(defun evolve (solute &rest args &key input-bounding-box
                                   (lipid-selector (list (cons 1.0 *popc*)))
                                   (population-size 200)
                                   (select-fraction 0.1)
                                   (number-of-generations 100)
                                   (parallel nil))
  (let* ((template-membrane (build-ga-membrane solute :input-bounding-box input-bounding-box :lipid-selector lipid-selector))
         (generation 0))
    ;; Fill the population at start
    (let* ((membranes (build-population template-membrane :population-size population-size))
           scored-membranes)
      (loop named evolution
            for sorted-scored-membranes = (score-membranes membranes :sort t :parallel parallel)
            for best-score = (score (first sorted-scored-membranes))
            do (setf scored-membranes sorted-scored-membranes)
            do (progn
                 (format t "Generation ~3a size ~a~%" generation (length scored-membranes))
                 (format t "  ... least collisions -> ~3a~%" (loop for index from 0 below 10
                                                                   collect (score (elt scored-membranes index))))
                 (finish-output))
            do (when (= (score (first scored-membranes)) 0) (return-from evolution t))
            do (incf generation)
            do (when (> generation number-of-generations) (return-from evolution scored-membranes))
               ;; Copy the scored membranes into membranes
            do (setf membranes (mapcar #'membrane scored-membranes))
               ;; Mutate lipid pairs that have bad contacts.
            do (format t "About to mutate-collisions~%")
            do (let ((collisions (mapcar #'collisions scored-membranes)))
                 (setf membranes (mutate-collisions membranes collisions)))
               ;; generate recombined membranes by mixing parts of one top membrane with another
            do  (let ((new-membranes (cross-membranes membranes select-fraction)))
                  (setf membranes new-membranes))
               ;; Do some random mutations of *mutate-lipid-fraction* 
            do (format t "About to randomly-mutate-membranes~%")
            do (let ((new-membranes (randomly-mutate-membranes membranes *mutate-lipid-fraction*)))
                 (setf membranes new-membranes))
               ;; Do some random shifting of *mutate-position-fraction* 
            do (format t "About to randomly-shift-lipids-in-membranes~%")
            do (let ((new-membranes (randomly-shift-lipids-in-membranes membranes *mutate-position-fraction*)))
                 (setf membranes new-membranes)))
      scored-membranes)))


#+(or)
(defun evolve2 (template-membrane &rest args &key input-bounding-box
                                               (lipid-selector (list (cons 1.0 *popc*)))
                                               (population-size 200)
                                               (select-fraction 0.1)
                                               (number-of-generations 100)
                                               (parallel nil)
                                               (close-distance *close-distance*))
  (let* ((generation 0))
    ;; Fill the population at start
    (let* ((membranes (build-population-orientation-only template-membrane :population-size population-size))
           scored-membranes
           (work-list (build-ga-membrane-work-list (first membranes))))
      (loop named evolution
            for sorted-scored-membranes = (score-membranes work-list membranes :sort t :close-distance close-distance)
            for best-score = (score (first sorted-scored-membranes))
            do (setf scored-membranes sorted-scored-membranes)
            do (progn
                 (format t "Generation ~3a size ~a~%" generation (length scored-membranes))
                 (format t "  ... least collisions -> ~3a~%" (loop for index from 0 below 10
                                                                   collect (score (elt scored-membranes index))))
                 (finish-output))
            do (when (= (score (first scored-membranes)) 0) (return-from evolution t))
            do (incf generation)
            do (when (> generation number-of-generations) (return-from evolution scored-membranes))
               ;; Copy the scored membranes into membranes
            do (setf membranes (mapcar #'membrane scored-membranes))
               ;; Mutate lipid pairs that have bad contacts.
               ;; Do some random shifting of *mutate-position-fraction* 
            do (let ((new-membranes (randomly-shift-lipids-in-membranes membranes *mutate-position-fraction*)))
                 (setf membranes new-membranes)))
      scored-membranes)))


(defun pack (solute &rest args &key input-bounding-box
                                 (lipid-selector (list (cons 1.0 *popc*)))
                                 (number-of-generations 100)
                                 (parallel nil))
  (apply 'evolve solute :input-bounding-box input-bounding-box args))
