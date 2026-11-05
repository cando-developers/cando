(in-package #:topology)

;;; Keep the related widths together. Foldamer methods select their policy
;;; explicitly by calling the helpers below with these or other widths.
(defparameter +shape-dihedral-bin-degrees+ 10)
(defparameter +shape-cp-phase-bin-degrees+ 10)
(defparameter +shape-cp-polar-bin-degrees+ 10)
(defparameter +shape-cp-amplitude-bin-width+ 0.1d0)

(defun check-shape-component-values (kind values &key bins)
  (let ((size (ecase kind (:dihedral 1) (:ring5 2) (:ring6 3))))
    (unless (and (binned-shape-key-proper-list-p values)
                 (= size (length values))
                 (every (if bins #'integerp
                            (lambda (value)
                              (and (realp value)
                                   (<= (- most-positive-double-float)
                                       value most-positive-double-float))))
                        values))
      (error "Invalid ~s ~a: ~s" kind (if bins "bins" "parameters") values))))

(defun shape-angular-bin-count (width period)
  (unless (and (realp width) (plusp width) (zerop (mod period width)))
    (error "Bin width ~s must divide angular period ~s" width period))
  (round period width))

(defun bin-shape-component-with-widths
    (kind parameters &key (dihedral-width +shape-dihedral-bin-degrees+)
                           (phase-width +shape-cp-phase-bin-degrees+)
                           (polar-width +shape-cp-polar-bin-degrees+)
                           (amplitude-width +shape-cp-amplitude-bin-width+))
  "Shared numerical implementation; physical angles are radians."
  (check-shape-component-values kind parameters)
  (shape-angular-bin-count polar-width 180)
  (unless (and (realp amplitude-width) (plusp amplitude-width))
    (error "Amplitude bin width must be positive"))
  (labels ((circular (angle width)
             (mod (round (* angle (/ 180d0 pi)) width)
                  (shape-angular-bin-count width 360)))
           (amplitude (value)
             (when (minusp value) (error "Negative CP amplitude ~s" value))
             (round value amplitude-width)))
    (ecase kind
      (:dihedral (list (circular (first parameters) dihedral-width)))
      (:ring5 (list (amplitude (first parameters))
                    (circular (second parameters) phase-width)))
      (:ring6
       (unless (<= 0 (second parameters) pi)
         (error "CP polar angle must be in [0, pi]"))
       (list (amplitude (first parameters))
             (round (* (second parameters) (/ 180d0 pi)) polar-width)
             (circular (third parameters) phase-width))))))

(defun unbin-shape-component-with-widths
    (kind bins &key (dihedral-width +shape-dihedral-bin-degrees+)
                     (phase-width +shape-cp-phase-bin-degrees+)
                     (polar-width +shape-cp-polar-bin-degrees+)
                     (amplitude-width +shape-cp-amplitude-bin-width+))
  "Decode canonical bins to representative geometry, with angles in radians."
  (check-shape-component-values kind bins :bins t)
  (shape-angular-bin-count polar-width 180)
  (unless (and (realp amplitude-width) (plusp amplitude-width))
    (error "Amplitude bin width must be positive"))
  (labels ((angle (bin width)
             (unless (<= 0 bin (1- (shape-angular-bin-count width 360)))
               (error "Noncanonical circular bin ~s for width ~s" bin width))
             (* bin width (/ pi 180d0)))
           (amplitude (bin)
             (when (minusp bin) (error "Negative amplitude bin ~s" bin))
             (* bin amplitude-width)))
    (ecase kind
      (:dihedral (list (angle (first bins) dihedral-width)))
      (:ring5 (list (amplitude (first bins)) (angle (second bins) phase-width)))
      (:ring6
       (unless (<= 0 (second bins) (shape-angular-bin-count polar-width 180))
         (error "Noncanonical polar-angle bin ~s" (second bins)))
       (list (amplitude (first bins)) (* (second bins) polar-width (/ pi 180d0))
             (angle (third bins) phase-width))))))
