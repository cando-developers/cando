(in-package :leap.set-box)

(defun set-box (aggregate enclosure &optional buffer)
  "The set-box command creates a bounding-box around the solute."
  (check-type aggregate chem:aggregate)
  (let ((buffer-vec (make-array 3))
        (dx 0)
        (dy 0)
        (dz 0))
    (if (and (not (eq :vdw enclosure))
             (not (eq :centers enclosure)))
        (error "~s - Expected :vdw or :centersfor second argument~%" enclosure))
    (if (eq :vdw enclosure)
        (center-unit-by-radii aggregate)
        (set-unit-box-by-centers aggregate))
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
        (setf dx (+ (chem:get-x-width bounding-box) (aref buffer-vec 0)))
        (setf dy (+ (chem:get-y-width bounding-box) (aref buffer-vec 1)))
        (setf dz (+ (chem:get-z-width bounding-box) (aref buffer-vec 2))))
      (chem:set-bounding-box aggregate (chem:make-bounding-box (list dx dy dz))))
    (let ((bounding-box (chem:bounding-box aggregate)))
      (format t "Box dimensions ~a ~a ~a~%"
              (chem:get-x-width bounding-box)
              (chem:get-y-width bounding-box)
              (chem:get-z-width bounding-box)))))

(defun center-unit-by-radii (aggregate)
  (check-type aggregate chem:aggregate)
  (let ((atom-default-radius 1.5)
        (ifirst 1)
        (x 0.0)
        (y 0.0)
        (z 0.0)
        (x-max 0.0)
        (y-max 0.0)
        (z-max 0.0)
        (x-min 0.0)
        (y-min 0.0)
        (z-min 0.0)
        (radius 0.0))
    (chem:map-atoms
     nil
     (lambda (a)
       (setf radius (chem:get-vdw-radius a))
       (if (< radius 0.1)
           (if (eq (chem:get-element-as-symbol a) :H)
               (setf radius 1.0)
               (progn
                 #+(or)(format t "(using default radius ~a for ~a)~%"
                               atom-default-radius (chem:get-name a))
                 (setf radius atom-default-radius))))
       (chem:set-vdw-radius a radius)
       (setf x (+ (geom:vx (chem:get-position a)) radius))
       (setf y (+ (geom:vy (chem:get-position a)) radius))
       (setf z (+ (geom:vz (chem:get-position a)) radius))
       (if (=  ifirst 1)
           (setf x-max x
                 y-max y
                 z-max z
                 ifirst 0)
           (progn
             (if (> x x-max)
                 (setf x-max x))
             (if (> y y-max)
                 (setf y-max y))
             (if (> z z-max)
                 (setf z-max z))))
       (setf x (- (geom:vx (chem:get-position a)) radius))
       (setf y (- (geom:vy (chem:get-position a)) radius))
       (setf z (- (geom:vz (chem:get-position a)) radius))
       (if (=  ifirst 1)
           (setf x-min x
                 y-min y
                 z-min z
                 ifirst 0)
           (progn
             (if (< x x-min)
                 (setf x-min x))
             (if (< y y-min)
                 (setf y-min y))
             (if (< z z-min)
                 (setf z-min z)))))
     aggregate)

    ;;Define center of bounding box
    (setf x (+ x-min (* (- x-max x-min) 0.5)))
    (setf y (+ x-min (* (- x-max x-min) 0.5)))
    (setf z (+ x-min (* (- x-max x-min) 0.5)))

    ;;Translate center to origin
    (let* ((translate-mol (geom:vec x y z))
           (mol-transform (geom:make-m4-translate translate-mol)))
      (chem:apply-transform-to-atoms aggregate mol-transform)
      (chem:set-bounding-box aggregate (chem:make-bounding-box (list (- x-max x-min) (- y-max y-min) (- z-max z-min)))))))

(defun set-unit-box-by-centers (aggregate)
  (let ((atom-default-radius 1.5)
        (ifirst 1)
        (x 0.0)
        (y 0.0)
        (z 0.0)
        (x-max 0.0)
        (y-max 0.0)
        (z-max 0.0)
        (x-min 0.0)
        (y-min 0.0)
        (z-min 0.0)
        (radius 0.0))
    (chem:map-atoms
     nil
     (lambda (a)
       (setf x (geom:vx (chem:get-position a)))
       (setf y (geom:vy (chem:get-position a)))
       (setf z (geom:vz (chem:get-position a)))
       (if (=  ifirst 1)
           (setf x-max x
                 x-min x
                 y-max y
                 y-min y
                 z-max z
                 z-min z
                 ifirst 0))
       (if (> x x-max)
           (setf x-max x)
           (if (< x x-min)
               (setf x-min x)))
       (if (> y y-max)
           (setf y-max y)
           (if (< y y-min)
               (setf y-min y)))
       (if (> z z-max)
           (setf z-max z)
           (if (< z z-min)
               (setf z-min z))))
     aggregate)
    (chem:set-bounding-box aggregate
                           (chem:make-bounding-box (list (- x-max x-min) (- y-max y-min) (- z-max z-min))))))
