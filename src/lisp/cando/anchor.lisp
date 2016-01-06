(in-package :anchor)

;;;
;;; Restrain atoms to points in space
;;;

;;; Different patterns of points in space
(defun circle-points (radius points &key (z 0.0))
  "Generate points on a circle"
  (loop for aindex from 0 below points
     for angrad = (* 360.0 (/ (float aindex) points) 0.0174533)
     collect (geom:make-ovector3 (* (cos angrad) radius)
                                 (* (sin angrad) radius) z)))


(defun zig-zag-cylinder-points (radius num height)
  "* Arguments
- radius :: The radius of the cylinder.
- num :: The number of points on each end of the cylinder.
- height :: The height of the cylinder.
* Description
Generate points that zig-zag from one end of a cylinder to the other in the pattern top-(bottom-bottom-top-top)n-bottom Where n = (floor (/ 1- num) 2)."
  (let (cy
	(top-circle (anchor:circle-points radius num :z (/ height 2)))
	(bottom-circle (anchor:circle-points radius num :z (- (/ height 2)))))
    (push (pop top-circle) cy)
    (dotimes (i (floor (/ (1- num) 2)))
      (push (pop bottom-circle) cy)
      (push (pop bottom-circle) cy) 
      (push (pop top-circle) cy)
      (push (pop top-circle) cy))
    (push (pop bottom-circle) cy)
    (nreverse cy)))

;;;
;;; Anchor atoms to points
;;; First always clear the anchor restraints
(defun anchor-atom (atom pos &optional (weight 1.0))
  (let ((restraint-anchor (core:make-cxx-object 'chem:restraint-anchor :atom atom :position pos :weight weight)))
    (chem:clear-restraints atom)
    (chem:add-restraint atom restraint-anchor)
    restraint-anchor))

(defun apply-restraints (atoms points &optional (weight 1.0))
  (loop for a in atoms
     for p in points
     do (anchor-atom a p weight)))

(defun on-circle (atoms &key (radius 15.0) (weight 1.0))
  "Anchor atoms on a circle"
  (let ((points (circle-points radius (length atoms))))
    (apply-restraints atoms points weight)))

(defun on-points (atoms points &key (weight 1.0))
  (apply-restraints atoms points weight))
