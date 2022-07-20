(in-package :leap.set-box)

(defun set-box (aggregate enclosure &optional buffer)
"    setBox solute enclosure [ buffer ]
      UNIT                         _solute_
      \"vdw\" OR \"centers\"           _enclosure_
      object                       _buffer_

The setBox command creates a periodic box around the _solute_ UNIT, turning
it into a periodic system for the simulation programs.  It does not add
any solvent to the system. The choice of \"vdw\" or \"centers\" determines
whether the box encloses all entire atoms or just all atom centers -
use \"centers\" if the system has been previously equilibrated as a
periodic box. See the solvateBox command for a description of the
buffer object, which extends either type of box by an arbitrary amount.
"
  (check-type aggregate chem:aggregate)

  (let ((buffer-vec (make-array 3)))
    (if (and (not (eq :vdw enclosure))
             (not (eq :centers enclosure)))
        (error "~s - Expected :vdw or :centers for second argument~%" enclosure))
    (calculate-bounding-box aggregate (eq :vdw enclosure))
    (when buffer
      (if (listp buffer)
          (progn
            (if (/= (list-length buffer) 3)
                (error "Expected 3 non-negative floating point numbers {x y z} for third argument~%"))
            (loop for i from 0 below 3
                  for number = (nth i buffer)
                  do (if (not (numberp number))
                         (error "Expected 3 non-negative floating point numbers {x y z} for third argument~%"))
                  do (setf (aref buffer-vec i) (* number 2.0))
                  do (if (< (aref buffer-vec i) 0.0)
                         (error "Expected 3 non-negative floating point numbers {x y z} for third argument~%"))))
          (progn
            (if (not (numberp buffer))
                (error "Expected 3 non-negative floating point numbers {x y z} for third argument~%"))
            (setf (aref buffer-vec 0) (* buffer 2.0))
            (if (< (aref buffer-vec 0) 0.0)
                (error "Expected 3 non-negative floating point numbers {x y z} for third argument~%"))
            (setf (aref buffer-vec 1) (aref buffer-vec 0))
            (setf (aref buffer-vec 2) (aref buffer-vec 0))))
      (let ((bounding-box (chem:bounding-box aggregate)))
        (setf (chem:bounding-box aggregate)
              (chem:make-bounding-box (list (+ (chem:get-x-width bounding-box) (aref buffer-vec 0))
                                            (+ (chem:get-y-width bounding-box) (aref buffer-vec 1))
                                            (+ (chem:get-z-width bounding-box) (aref buffer-vec 2)))))))
    (let ((bounding-box (chem:bounding-box aggregate)))
      (format t "Box dimensions ~a ~a ~a~%"
              (chem:get-x-width bounding-box)
              (chem:get-y-width bounding-box)
              (chem:get-z-width bounding-box)))))

(defun calculate-bounding-box (aggregate vdw)
  (let ((cuboid (geom:make-bounding-cuboid)))
    (chem:map-atoms nil
                    (lambda (atom)
                      (let* ((radius (cond ((not vdw)
                                            0)
                                           ((plusp (chem:get-vdw-radius atom))
                                            (chem:get-vdw-radius atom))
                                           (t
                                            (chem:vdw-radius-for-element (chem:get-element atom)))))
                             (corner (geom:vec radius radius radius)))
                        (geom:expand-to-encompass-point cuboid
                                                        (geom:v+ (chem:get-position atom) corner))
                        (geom:expand-to-encompass-point cuboid
                                                        (geom:v- (chem:get-position atom) corner))))
                    aggregate)
    (setf (chem:bounding-box aggregate)
          (chem:make-bounding-box (list (geom:get-extent-x cuboid)
                                        (geom:get-extent-y cuboid)
                                        (geom:get-extent-z cuboid))
                                  :center (list (geom:get-x (geom:get-center cuboid))
                                                (geom:get-y (geom:get-center cuboid))
                                                (geom:get-z (geom:get-center cuboid)))))))
