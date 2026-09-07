#+tests(in-package #:clasp-tests)

;;; Shared finite-difference support for two-atom pair-potential regressions.
;;;
;;; The Rosetta pair components currently expose their generated second
;;; derivatives through H*d, but do not populate a supplied Hessian matrix.
;;; Applying H to every Cartesian basis vector therefore gives every Hessian
;;; column and lets us check all 36 entries of a two-atom potential.

(defun pairwise-derivatives--nvector-ref (vector index)
  "Return scalar INDEX from the packed XYZ NVECTOR VECTOR."
  (let* ((base (* 3 (floor index 3)))
         (component (mod index 3))
         (value (geom:vec-array vector base)))
    (ecase component
      (0 (geom:get-x value))
      (1 (geom:get-y value))
      (2 (geom:get-z value)))))

(defun pairwise-derivatives--nvector-set (vector index new-value)
  "Set scalar INDEX in the packed XYZ NVECTOR VECTOR to NEW-VALUE."
  (let* ((base (* 3 (floor index 3)))
         (component (mod index 3))
         (value (geom:vec-array vector base))
         (x (geom:get-x value))
         (y (geom:get-y value))
         (z (geom:get-z value)))
    (ecase component
      (0 (setf x new-value))
      (1 (setf y new-value))
      (2 (setf z new-value)))
    (geom:vec-set (geom:vec x y z) vector base)))

(defun pairwise-derivatives--error-ratio
    (actual expected absolute-tolerance relative-tolerance)
  "Return a normalized error; values at or below one are within tolerance."
  (/ (abs (- actual expected))
     (+ absolute-tolerance
        (* relative-tolerance
           (max (abs actual) (abs expected))))))

(defun pairwise-derivatives--set-oblique-separation (atom1 atom2 distance)
  "Place ATOM1 and ATOM2 DISTANCE apart along an oblique unit vector."
  (chem:set-position atom1 (geom:vec 0.0d0 0.0d0 0.0d0))
  ;; (2/3,-1/3,2/3) has unit length and makes every mixed Cartesian
  ;; Hessian entry nonzero away from stationary points.
  (chem:set-position atom2
                     (geom:vec (* distance (/ 2.0d0 3.0d0))
                               (* distance (/ -1.0d0 3.0d0))
                               (* distance (/ 2.0d0 3.0d0)))))

(defun pairwise-derivatives--save-position (position size)
  "Copy the SIZE scalar coordinates in POSITION to a double-float array."
  (let ((saved (make-array size :element-type 'double-float)))
    (dotimes (coordinate size saved)
      (setf (aref saved coordinate)
            (pairwise-derivatives--nvector-ref position coordinate)))))

(defun pairwise-derivatives--restore-position (position saved)
  "Restore POSITION from the scalar coordinate array SAVED."
  (dotimes (coordinate (length saved))
    (pairwise-derivatives--nvector-set
     position coordinate (aref saved coordinate))))

(defun pairwise-derivatives--set-vector (vector values)
  "Copy scalar VALUES into packed XYZ NVECTOR VECTOR."
  (dotimes (coordinate (length values))
    (pairwise-derivatives--nvector-set
     vector coordinate (aref values coordinate))))

(defun pairwise-derivatives--set-offset-position
    (position origin direction scale)
  "Set POSITION to ORIGIN + SCALE*DIRECTION."
  (dotimes (coordinate (length origin))
    (pairwise-derivatives--nvector-set
     position coordinate
     (+ (aref origin coordinate)
        (* scale (aref direction coordinate))))))

(defun pairwise-derivatives--central-hd-error
    (energy-function position origin direction analytic-hdvec step)
  "Compare ANALYTIC-HDVEC with -d(force)/dt along DIRECTION."
  (let* ((size (length origin))
         (force-plus (chem:make-nvector size))
         (force-minus (chem:make-nvector size))
         (maximum-error 0.0d0))
    (unwind-protect
         (progn
           (pairwise-derivatives--set-offset-position
            position origin direction step)
           (chem:evaluate-energy-force
            energy-function position :calc-force t :force force-plus)
           (pairwise-derivatives--set-offset-position
            position origin direction (- step))
           (chem:evaluate-energy-force
            energy-function position :calc-force t :force force-minus)
           ;; Cando force is -gradient, so H*d is -d(force)/dt.
           (dotimes (coordinate size)
             (let ((numerical-hdvec
                     (- (/ (- (pairwise-derivatives--nvector-ref
                                force-plus coordinate)
                              (pairwise-derivatives--nvector-ref
                                force-minus coordinate))
                           (* 2.0d0 step)))))
               (setf maximum-error
                     (max maximum-error
                          (pairwise-derivatives--error-ratio
                           (pairwise-derivatives--nvector-ref
                            analytic-hdvec coordinate)
                           numerical-hdvec
                           2.0d-5 2.0d-6))))))
      (pairwise-derivatives--restore-position position origin))
    maximum-error))

(defun pairwise-derivatives--finite-difference-errors
    (energy-function position &key (step 1.0d-5))
  "Check force, every Hessian entry, and a dense H*d by finite differences.

The component under test must be the only enabled component of ENERGY-FUNCTION.
Return four values: energy, maximum normalized force error, maximum normalized
Hessian-column error, and maximum normalized dense H*d error."
  (let* ((size (chem:get-nvector-size energy-function))
         (dense-direction
           #(0.31d0 -0.47d0 0.23d0 -0.19d0 0.41d0 -0.37d0))
         (origin (pairwise-derivatives--save-position position size))
         (dvec (chem:make-nvector size))
         (analytic-force (chem:make-nvector size))
         (numerical-force (chem:make-nvector size))
         (analytic-hdvec (chem:make-nvector size)))
    (unless (= size (length dense-direction))
      (error "The pairwise derivative regression expected ~d coordinates, not ~d"
             (length dense-direction) size))
    (unwind-protect
         (progn
           (pairwise-derivatives--set-vector dvec dense-direction)
           (let ((energy
                   (chem:scoring-function/evaluate-all
                    energy-function position
                    :calc-force t
                    :force analytic-force
                    :calc-diagonal-hessian t
                    :calc-off-diagonal-hessian t
                    :hdvec analytic-hdvec
                    :dvec dvec)))
             (chem:evaluate-finite-difference-force
              energy-function position :force numerical-force :delta step)
             (let ((maximum-force-error 0.0d0)
                   (maximum-hessian-error 0.0d0))
               (dotimes (coordinate size)
                 (setf maximum-force-error
                       (max maximum-force-error
                            (pairwise-derivatives--error-ratio
                             (pairwise-derivatives--nvector-ref
                              analytic-force coordinate)
                             (pairwise-derivatives--nvector-ref
                              numerical-force coordinate)
                             2.0d-6 2.0d-7))))
               (let ((dense-hdvec-error
                       (pairwise-derivatives--central-hd-error
                        energy-function position origin dense-direction
                        analytic-hdvec step)))
                 ;; Each basis-vector product is one complete Hessian column.
                 (dotimes (column size)
                   (let ((basis (make-array size
                                            :element-type 'double-float
                                            :initial-element 0.0d0))
                         (basis-force (chem:make-nvector size))
                         (basis-hdvec (chem:make-nvector size)))
                     (setf (aref basis column) 1.0d0)
                     (pairwise-derivatives--set-vector dvec basis)
                     (chem:scoring-function/evaluate-all
                      energy-function position
                      :calc-force t
                      :force basis-force
                      :calc-diagonal-hessian t
                      :calc-off-diagonal-hessian t
                      :hdvec basis-hdvec
                      :dvec dvec)
                     (setf maximum-hessian-error
                           (max maximum-hessian-error
                                (pairwise-derivatives--central-hd-error
                                 energy-function position origin basis
                                 basis-hdvec step)))))
                 (values energy maximum-force-error
                         maximum-hessian-error dense-hdvec-error)))))
      (pairwise-derivatives--restore-position position origin))))
