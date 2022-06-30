(in-package :leap.align-axes)

(defun tool-orient-principle-axis-along-coordinate-axis (unit)
  (let* ((unit-center (chem:geometric-center unit))
         (transform-to-origin (geom:make-m4-translate (geom:v* unit-center -1.0)))
         (eigen-values (make-array 3 :element-type 'double-float :initial-contents '(1.0d0 2.0d0 3.0d0)))
         (diagonalize (geom:make-matrix-identity))
         (moment (geom:make-matrix-identity)))
     ;;Translate the GEOMETRIC center of the UNIT to the origin
    (chem:apply-transform-to-atoms unit transform-to-origin)
    ;;Now calculate the moment of inertia in terms of GEOMETRY
    ;;meaning assign a mass of 1 to ALL ATOMS
    (setf moment (math-op-moment-of-inertia unit))
    ;;Diagonalize the moment of inertia matrix
    (math-op-diagonalize moment eigen-values diagonalize)
    (when (chem:verbose 2)
      (format t "eigen-values ~s~%" eigen-values)
      (format t "moment: ~s~%" moment)
      (format t "diagonalize: ~%~s~%" diagonalize))
    (let* ((diagonalize-0 (geom:vec (geom:at-row-col-get diagonalize 0 0)
                                    (geom:at-row-col-get diagonalize 0 1)
                                    (geom:at-row-col-get diagonalize 0 2)))
           (diagonalize-1 (geom:vec (geom:at-row-col-get diagonalize 1 0)
                                    (geom:at-row-col-get diagonalize 1 1)
                                    (geom:at-row-col-get diagonalize 1 2)))
           (diagonalize-2 (geom:vec (geom:at-row-col-get diagonalize 2 0)
                                    (geom:at-row-col-get diagonalize 2 1)
                                    (geom:at-row-col-get diagonalize 2 2)))
           (pos (geom:vcross diagonalize-0 diagonalize-1))
           (dot (geom:vdot pos diagonalize-2)))
      ;;If the handedness of the matrix is wrong then change it
      (if (< dot 0.0)
          (progn
            (geom:at-row-col-put diagonalize 2 0 (* (geom:at-row-col-get diagonalize 2 0) -1.0))
            (geom:at-row-col-put diagonalize 2 1 (* (geom:at-row-col-get diagonalize 1 0) -1.0))
            (geom:at-row-col-put diagonalize 2 2 (* (geom:at-row-col-get diagonalize 2 0) -1.0))))
      ;; Get the largest eigenvalue along z
      (let* ((xeigen (elt eigen-values 0))
             (yeigen (elt eigen-values 1))
             (zeigen (elt eigen-values 2))
             (rotate-largest-to-zaxis (cond
                                       ((and (> xeigen yeigen)
                                             (> xeigen zeigen))
                                        (geom:make-m4-rotate-y (* 90.0 0.0174533)))
                                       ((and (> yeigen xeigen)
                                             (> yeigen zeigen))
                                        (geom:make-m4-rotate-x (* 90.0 0.0174533)))
                                       (t
                                        (geom:make-m4-rotate-z 0.0)) ; identity
                                       )))
        (when (chem:verbose 2)
          (format t "rotate-largest-to-zaxis:~%~s~%" rotate-largest-to-zaxis))
        (setf diagonalize (geom:m*m rotate-largest-to-zaxis diagonalize)
              eigen-values (geom:m*v rotate-largest-to-zaxis eigen-values)))
      (when (chem:verbose 2)
        (format t "after eigen-values ~s~%" eigen-values)
        (format t "after diagonalize: ~%~s~%" diagonalize))
      ;; Apply the matrix to the UNIT
      (chem:map-atoms
       'nil
       (lambda (atom)
         (let ((atom-position (chem:get-position atom)))
           (geom:set-all3 atom-position
                          (+ (* (geom:at-row-col-get diagonalize 0 0) (geom:get-x atom-position))
                             (* (geom:at-row-col-get diagonalize 1 0) (geom:get-y atom-position))
                             (* (geom:at-row-col-get diagonalize 2 0) (geom:get-z atom-position)))
                          (+ (* (geom:at-row-col-get diagonalize 0 1) (geom:get-x atom-position))
                             (* (geom:at-row-col-get diagonalize 1 1) (geom:get-y atom-position))
                             (* (geom:at-row-col-get diagonalize 2 1) (geom:get-z atom-position)))
                          (+ (* (geom:at-row-col-get diagonalize 0 2) (geom:get-x atom-position))
                             (* (geom:at-row-col-get diagonalize 1 2) (geom:get-y atom-position))
                             (* (geom:at-row-col-get diagonalize 2 2) (geom:get-z atom-position))))
           #+(or)(geom:set-all3 atom-position
                          (geom:get-x (geom:m*v diagonalize atom-position))
                          (geom:get-y (geom:m*v diagonalize atom-position))
                          (geom:get-z (geom:m*v diagonalize atom-position)))
           (chem:set-position atom atom-position)))
       unit))))


    
(defun math-op-moment-of-inertia (unit)
  (let ((moment (geom:make-matrix-identity))
        (moment-array (make-array 16 :element-type 'double-float))
        (sxx 0.0)
        (syy 0.0)
        (szz 0.0)
        (sxy 0.0)
        (sxz 0.0)
        (syz 0.0))
    (chem:map-atoms
     'nil
     (lambda (atom)
       (let* ((position (chem:get-position atom))
              (x2 (* (geom:get-x position) (geom:get-x position)))
              (y2 (* (geom:get-y position) (geom:get-y position)))
              (z2 (* (geom:get-z position) (geom:get-z position))))
         (setf sxx (+ sxx y2 z2)
               syy (+ syy x2 z2)
               szz (+ szz x2 y2)
               sxy (+ sxy (* (geom:get-x position) (geom:get-y position)))
               sxz (+ sxz (* (geom:get-x position) (geom:get-z position)))
               syz (+ syz (* (geom:get-y position) (geom:get-z position))))))
     unit)
    (setf (aref moment-array 0) sxx
          (aref moment-array 1) (- sxy)
          (aref moment-array 2) (- sxz)
          (aref moment-array 4) (- sxy)
          (aref moment-array 5) syy
          (aref moment-array 6) (- syz)
          (aref moment-array 8) (- sxz)
          (aref moment-array 9) (- syz)
          (aref moment-array 10) szz)
    (geom:set-from-vector moment moment-array)
    moment))

(defun math-op-diagonalize (ma eigen-values diagonalize)
  (let* ((vector-size 3)
        (b-vector (make-array vector-size :element-type 'double-float))
        (d-vector (make-array vector-size :element-type 'double-float))
        (z-vector (make-array vector-size :element-type 'double-float))
        (thresh 0.0)
        (h 0.0)
        (tt 0.0)
        (theta 0.0)
        (c 0.0)
        (s 0.0)
        (tau 0.0)
        (steps 0)
         )
    (setf d-vector eigen-values)
     (loop for i from 0 below vector-size
                          ;;(format t "i ~a j ~a value ~a~%" i j (geom:at-row-col-get ma i j))))
          do (setf (aref b-vector i) (geom:at-row-col-get ma i i)
                   (aref d-vector i) (geom:at-row-col-get ma i i)
                   (aref z-vector i) 0.0d0))
    (loop for i from 1 below 51
          for sm = 0.0
          do (loop for j from 0 below (- vector-size 1)
                   do (loop for k from (+ j 1) below vector-size
                            do (setf sm (+ sm (abs (geom:at-row-col-get ma j k))))))       
          do (if (= sm 0.0) (return-from math-op-diagonalize))
          do (if (< i 4)
                 (setf thresh (/ (* 0.2 sm) (* vector-size vector-size)))
                 (setf thresh 0.0))
          do (loop for j from 0 below (- vector-size 1)
                   do (loop for k from (+ j 1) below vector-size
                            for g = (* 100.0 (abs  (geom:at-row-col-get ma j k)))
                            do (if (and (> i 4)
                                        (= (+ g (abs (aref d-vector j))) (abs (aref d-vector j)))
                                        (= (+ g (abs (aref d-vector k))) (abs (aref d-vector k))))
                                    (geom:at-row-col-put ma j k 0.0)
                                    (if (> (abs (geom:at-row-col-get ma j k)) thresh)
                                       (progn
                                         (setf h (- (aref d-vector k) (aref d-vector j)))
                                         (if (= (abs h) (+ (abs h) g))
                                             (setf tt (/ (geom:at-row-col-get ma j k) h))
                                             (progn
                                               (setf theta (* 0.5 (/ h  (geom:at-row-col-get ma j k))))
                                               (setf tt (/ 1.0 (+ (abs theta) (sqrt (+ 1.0 (* theta theta))))))
                                               (if (< theta 0.0) (setf tt (- tt)))))
                                         (setf c (/ 1.0 (sqrt (+ 1 (* tt tt))))
                                               s (* tt c)
                                               tau (/ s (+ 1.0 c))
                                               h (* tt  (geom:at-row-col-get ma j k))
                                               (aref z-vector j) (- (aref z-vector j) h)
                                               (aref z-vector k) (+ (aref z-vector k) h)
                                               (aref d-vector j) (- (aref d-vector j) h)
                                               (aref d-vector k) (+ (aref d-vector k) h))
                                         (geom:at-row-col-put ma j k 0.0)
                                         (loop for l from 0 below j
                                               do (progn
                                                    (setf g  (geom:at-row-col-get ma l j)
                                                          h  (geom:at-row-col-get ma l k))
                                                    (geom:at-row-col-put ma l j (- g (* s (+ h (* g tau)))))
                                                    (geom:at-row-col-put ma l k (+ h (* s (- g (* h tau)))))))
                                         (loop for l from (+ j 1) below k
                                               do (progn
                                                    (setf g (geom:at-row-col-get ma j l)
                                                          h (geom:at-row-col-get ma l k))
                                                    (geom:at-row-col-put ma j l (- g (* s (+ h (* g tau)))))
                                                    (geom:at-row-col-put ma l k (+ h (* s (- g (* h tau)))))))
                                          (loop for l from (+ k 1) below vector-size
                                               do (progn
                                                    (setf g (geom:at-row-col-get ma j l)
                                                          h (geom:at-row-col-get ma k l))
                                                    (geom:at-row-col-put ma j l (- g (* s (+ h (* g tau)))))
                                                    (geom:at-row-col-put ma k l (+ h (* s (- g (* h tau)))))))
                                          (loop for l from 0 below vector-size
                                               do (progn
                                                    (setf g (geom:at-row-col-get diagonalize l j)
                                                          h (geom:at-row-col-get diagonalize l k))
                                                    (geom:at-row-col-put diagonalize l j (- g (* s (+ h (* g tau)))))
                                                    (geom:at-row-col-put diagonalize l k (+ h (* s (- g (* h tau)))))))
                                         (incf steps))))))
           (loop for j from 0 below vector-size
                   do (setf (aref b-vector j) (+ (aref b-vector j) (aref z-vector j))
                            (aref d-vector j) (aref b-vector j)
                            (aref z-vector j) 0.0d0)))))
