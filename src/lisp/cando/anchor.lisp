;;;
;;;    File: anchor.lisp
;;;

;; Open Source License
;; Copyright (c) 2016, Christian E. Schafmeister
;; Permission is hereby granted, free of charge, to any person obtaining a copy
;; of this software and associated documentation files (the "Software"), to deal
;; in the Software without restriction, including without limitation the rights
;; to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
;; copies of the Software, and to permit persons to whom the Software is
;; furnished to do so, subject to the following conditions:
;; The above copyright notice and this permission notice shall be included in
;; all copies or substantial portions of the Software.
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
;; IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
;; OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
;; THE SOFTWARE.
;;  
;; This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.

;; -^-
(in-package :anchor)

;;;
;;; Restrain atoms to points in space
;;;

;;; Different patterns of points in space
(defun circle-points (radius points &key (z 0.0))
  "* Arguments
- radius :: Number in angstroms
- points :: Number of points
- z :: Height along z-axis.
* Description
Return POINTS points on a circle of RADIUS where the normal to the circle points along the Z-axis and the points are all offset from the origin by Z."
  (loop for aindex from 0 below points
     for angrad = (* 360.0 (/ (float aindex) points) 0.0174533)
     collect (geom:vec (* (cos angrad) radius)
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

(defun top-bottom-cylinder-points (radius num height)
  "* Arguments
- radius :: The radius of the cylinder.
- num :: The number of points on each end of the cylinder.
- height :: The height of the cylinder.
* Description
Generate (* 2 num) points that first go around the circle on the top of the cylinder and then around the bottom  Where n = (floor (/ 1- num) 2)."
  (let (cy
	(top-circle (anchor:circle-points radius (/ num 2) :z (/ height 2)))
	(bottom-circle (anchor:circle-points radius (/ num 2) :z (- (/ height 2)))))
    (dotimes (i (/ num 2))
      (push (pop top-circle) cy))
    (dotimes (i (/ num 2))
      (push (pop bottom-circle) cy))
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

(defun clear-all-atom-restraints (matter)
  (chem:map-atoms
   nil
   (lambda (a) (chem:clear-restraints a))
   matter))
