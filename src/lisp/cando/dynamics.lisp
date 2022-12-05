(in-package :dynamics)


(defclass trajectory ()
  ((coordinates :initform (make-array 16 :adjustable t :fill-pointer 0) :initarg :coordinates :accessor coordinates)
   (matter :initarg :matter :accessor matter)
   (number-of-atoms :initarg :number-of-atoms :accessor number-of-atoms))
  (:documentation "Maintain a vector of coordinates for a matter. 
This is to store a molecular dynamics trajectory or a trajectory from energy minimization.
Methods are specialized on this class in cando-nglview.lisp."))

(defclass simulation ()
  ((scoring-function :initarg :scoring-function :accessor scoring-function)
   (coordinates :initarg :coordinates :accessor coordinates)
   (velocity :initarg :velocity :accessor velocity)
   (forces :initarg :forces :accessor forces :type 'vector)
   (temp-forces :initarg :temp-forces :accessor temp-forces :type 'vector)
   (delta-t :initarg :delta-t :accessor delta-t :type 'double-float)
   (delta-t-over-mass :initarg :delta-t-over-mass :accessor delta-t-over-mass :type (vector))
   (current-time :initform 0.0 :accessor current-time)
   (accumulate-coordinates :initform nil :initarg :accumulate-coordinates :accessor accumulate-coordinates)
   (accumulated-coordinates :initform (make-array 16 :adjustable t :fill-pointer 0) :accessor accumulated-coordinates)))


(defgeneric make-trajectory (structure))

(defmethod make-trajectory ((structure chem:matter))
  (make-instance 'trajectory :matter structure
                                   :number-of-atoms (chem:number-of-atoms structure)))

(defmethod make-trajectory ((dynamics dynamics:simulation))
  (let* ((scoring-function (dynamics:scoring-function dynamics))
         (matter (chem:get-matter scoring-function))
         (trajectory (dynamics:accumulated-coordinates dynamics)))
    (unless (> (length trajectory) 1)
      (error "The trajectory must have more than one coordinate set"))
    (cond
      ((typep matter 'chem:aggregate) nil)
      ((typep matter 'chem:molecule)
       (let ((agg (chem:make-aggregate)))
         (chem:add-matter agg matter)
         (setf matter agg)))
      (t (error "matter must be aggregate or molecule")))
    (make-instance 'trajectory :matter matter
                                     :number-of-atoms (chem:number-of-atoms matter)
                                     :coordinates trajectory)))

(defmethod make-trajectory ((trajectory trajectory))
  (let ((matter (dynamics:matter trajectory))
        (coordinates (dynamics:coordinates trajectory)))
    (unless (> (length trajectory) 1)
      (error "The trajectory must have more than one coordinate set"))
    (make-instance 'trajectory :matter matter
                                     :number-of-atoms (chem:number-of-atoms matter)
                                     :trajectory coordinates)))

(defun do-accumulate-coordinates (simulation coordinates)
  "Convert coordinates into a vector of single-floats and append it to accumulated-coordinates"
  (let ((single-float-coordinates (make-array (length coordinates) :element-type 'single-float :adjustable nil)))
    (loop for index from 0 below (length coordinates)
          for sf-val = (float (elt coordinates index) 1.0s0)
          do (setf (elt single-float-coordinates index) sf-val))
    (vector-push-extend single-float-coordinates (accumulated-coordinates simulation))))

(defun make-atomic-simulation (scoring-function &key (delta-t 1.0) accumulate-coordinates)
  (let* ((vector-size (chem:get-nvector-size scoring-function))
         (coordinates (make-array vector-size :element-type 'double-float))
         (velocity (make-array vector-size :initial-element 0.0 :element-type 'double-float))
         (forces (make-array vector-size :element-type 'double-float))
         (temp-forces (make-array vector-size :element-type 'double-float))
         (delta-t-over-mass (make-array (/ vector-size 3) :element-type 'double-float))
         (atom-table (chem:atom-table scoring-function))
         (number-of-atoms (chem:get-number-of-atoms atom-table)))
    (loop for index from 0 below number-of-atoms
          for atom = (chem:elt-atom atom-table index)
          for element = (chem:get-element atom)
          for mass = (let ((mass (chem:atomic-weight-for-element element)))
                       (if (< mass 0.1)
                           1.0
                           mass))
          for dt-over-m = (/ delta-t mass)
          do (setf (aref delta-t-over-mass index) (float dt-over-m 1.0d0)))
    (chem:load-coordinates-into-vector scoring-function coordinates)
    (chem:evaluate-energy-force scoring-function coordinates t forces)
    (let ((sim (make-instance 'simulation
                              :scoring-function scoring-function
                              :coordinates coordinates
                              :velocity velocity
                              :forces forces
                              :temp-forces temp-forces
                              :delta-t delta-t
                              :delta-t-over-mass delta-t-over-mass
                              :accumulate-coordinates accumulate-coordinates)))
      (when accumulate-coordinates
        (do-accumulate-coordinates sim coordinates))
      sim)))


(defun make-sketch-simulation (scoring-function &key (delta-t 1.0) accumulate-coordinates)
  "Build a simulation object that carries out dynamics on a non-physical system - like the layout of a graph"
  (let* ((vector-size (chem:get-nvector-size scoring-function))
         (coordinates (make-array vector-size :element-type 'double-float))
         (velocity (make-array vector-size :initial-element 0.0 :element-type 'double-float))
         (forces (make-array vector-size :element-type 'double-float))
         (temp-forces (make-array vector-size :element-type 'double-float))
         (delta-t-over-mass (make-array (/ vector-size 3) :element-type 'double-float))
         (node-table (chem:node-table scoring-function))
         (node-table-size (chem:node-table-size node-table)))
    (loop for index from 0 below node-table-size
          for atom = (chem:node-table-node-at-index node-table index)
          for dt-over-m = delta-t
          do (setf (aref delta-t-over-mass index) dt-over-m))
    (chem:load-coordinates-into-vector scoring-function coordinates)
    (chem:evaluate-energy-force scoring-function coordinates t forces)
    (let ((sim (make-instance 'simulation
                              :scoring-function scoring-function
                              :coordinates coordinates
                              :velocity velocity
                              :forces forces
                              :temp-forces temp-forces
                              :delta-t delta-t
                              :delta-t-over-mass delta-t-over-mass
                              :accumulate-coordinates accumulate-coordinates)))
      (when accumulate-coordinates
        (do-accumulate-coordinates sim coordinates))
      sim)))


(defun dynamics-start (dynamics)
  (chem:load-coordinates-into-vector (scoring-function dynamics)
                                     (coordinates dynamics)))

(defun dynamics-end (dynamics)
  (chem:save-coordinates-from-vector (scoring-function dynamics)
                                     (coordinates dynamics)))


(defmacro with-dynamics ((dynamics) &body body)
  `(progn
     (dynamics-start ,dynamics)
     ,@body
     (dynamics-end ,dynamics)
     ,dynamics))

(defun velocity-verlet-step (dynamics &key (velocity-verlet-function #'chem:velocity-verlet-step) frozen)
  (funcall velocity-verlet-function
           (scoring-function dynamics)
           (coordinates dynamics)
           (velocity dynamics)
           (forces dynamics)
           (temp-forces dynamics)
           (delta-t-over-mass dynamics)
           (delta-t dynamics)
           frozen)
  (when (accumulate-coordinates dynamics)
    (do-accumulate-coordinates dynamics (coordinates dynamics)))
  (incf (current-time dynamics) (delta-t dynamics)))


(defun velocity-verlet-step-limit-displacement (dynamics
                                                &key (velocity-verlet-function #'chem:velocity-verlet-step-limit-displacement)
                                                  (limit-displacement (geom:vec 0.5 0.5 0.5))
                                                  frozen)
  (funcall velocity-verlet-function
           (scoring-function dynamics)
           (coordinates dynamics)
           (velocity dynamics)
           (forces dynamics)
           (temp-forces dynamics)
           (delta-t-over-mass dynamics)
           (delta-t dynamics)
           frozen
           limit-displacement)
  (when (accumulate-coordinates dynamics)
    (do-accumulate-coordinates dynamics (coordinates dynamics)))
  (incf (current-time dynamics) (delta-t dynamics)))

(defun write-coordinates-back-to-matter (simulation)
  (let ((scoring-function (scoring-function simulation))
        (coordinates (coordinates simulation)))
    (chem:save-coordinates-from-vector scoring-function coordinates)))
