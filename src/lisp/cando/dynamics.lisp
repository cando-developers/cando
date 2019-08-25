(in-package :dynamics)

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

(defun velocity-verlet-step (simulation &key velocity-verlet-function frozen)
  (funcall velocity-verlet-function
           (scoring-function simulation)
           (coordinates simulation)
           (velocity simulation)
           (forces simulation)
           (temp-forces simulation)
           (delta-t-over-mass simulation)
           (delta-t simulation)
           frozen)
  (when (accumulate-coordinates simulation)
    (do-accumulate-coordinates simulation (coordinates simulation)))
  (incf (current-time simulation) (delta-t simulation)))

(defun write-coordinates-back-to-matter (simulation)
  (let ((scoring-function (scoring-function simulation))
        (coordinates (coordinates simulation)))
    (chem:save-coordinates-from-vector scoring-function coordinates)))
