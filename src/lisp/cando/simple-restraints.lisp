(in-package :simple-restraints)

;;;
;;; Restrain atoms to points in space
;;;

;;; Different patterns of points in space
(defun circle-points (radius points &key (z 0.0))
  "Generate points on a circle"
  (loop for aindex from 0 below points
     for angrad = (* 360.0 (/ (float aindex) points) 0.0174533)
     collect (geom:make-v3 (* (cos angrad) radius)
                                 (* (sin angrad) radius) z)))


;;;
;;; Anchor atoms to points
;;; First always clear the anchor restraints
