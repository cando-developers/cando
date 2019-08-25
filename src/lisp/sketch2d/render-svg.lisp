(in-package :sketch2d)

(defparameter *xscale* 600)
(defparameter *yscale* 600)
(defparameter *perpendicular-fraction* 0.15)
(defparameter *parallel-fraction* 0.1)

(defun sketch-bonds (molecule)
  (let* ((scene (cl-svg:make-svg-toplevel 'cl-svg:svg-1.2-toplevel :height 100 :width 100)))
    (chem:map-bonds
     nil
     (lambda (a1 a2 order)
       (unless (or (= (chem:get-atomic-number a1) 1)
                   (= (chem:get-atomic-number a2) 1)
                   (eq (chem:get-element a1) :lp)
                   (eq (chem:get-element a2) :lp))
         (let* ((p1 (chem:get-position a1))
                (p2 (chem:get-position a2))
                (df (geom:v* (geom:vnormalized (geom:v- p1 p2)) 40.0)))
           (cl-svg:draw scene (:path :d (cl-svg:path
                                         (cl-svg:move-to 50 50)
                                         (cl-svg:line-to-r (geom:vx df) (geom:vy df)))
                                     :fill "none" :stroke "blue" :stroke-width 1)))))
     molecule)
    (with-output-to-string (sout) (cl-svg:stream-out sout scene))))

(defstruct draw-bond a1 a2 order)

(defun gather-bonds-and-atoms (molecule)
  "Return the bonds and atoms that we need to render"
  (let (bonds
        (atoms (make-hash-table)))
    (chem:map-bonds
     nil
     (lambda (a1 a2 bo)
       (cond
         ((eq (chem:get-element a1) :lp)
          (setf (gethash a2 atoms) t))
         ((eq (chem:get-element a2) :lp)
          (setf (gethash a1 atoms) t))
         ((= (chem:get-atomic-number a1) 1)
          (setf (gethash a2 atoms) t))
         ((= (chem:get-atomic-number a2) 1)
          (setf (gethash a1 atoms) t))
         (t (push (make-draw-bond :a1 a1 :a2 a2 :order bo) bonds)
            (setf (gethash a1 atoms) t
                  (gethash a2 atoms) t))))
     molecule)
    (values bonds atoms)))


(defun calculate-side-bond (x1 y1 x2 y2 &key (towards 1.0))
  "Calculate the start and end of a side-bond.  If towards is a vector it must be to the side of 
the main bond and the calculated bond will be on the same side.
If towards is a number then it should be 1.0 or -1.0.
This will place the calculated bond on one or the other side of the x1,y1-x2,y2 bond."
  (let* ((p1 (geom:vec x1 y1 0.0))
         (p2 (geom:vec x2 y2 0.0))
         (center (geom:v* (geom:v+ p1 p2) 0.5)) ; center of bond
         (d1-2 (geom:v- p2 p1)) ; bond direction vector
         (dz (geom:vec 0.0 0.0 1.0)) ; z-axis
         (dcross (geom:vnormalized (geom:vcross dz d1-2))) ; perpendicular to bond
         (dtowards (if (numberp towards)
                       (geom:v* dcross towards)
                       (geom:v- towards center))) ; if no towards then use pcross
         (dot (float-sign (geom:vdot dcross dtowards) 1.0))
         (len1-2 (geom:vlength d1-2))
         (offset-perpendicular (geom:v* dcross (* *perpendicular-fraction* len1-2 dot)))
         (offset-parallel (geom:v* d1-2 *parallel-fraction*))
         (p1o (geom:v+ (geom:v+ p1 offset-parallel) offset-perpendicular))
         (p2o (geom:v+ (geom:v- p2 offset-parallel) offset-perpendicular)))
    (values (geom:vx p1o) (geom:vy p1o) (geom:vx p2o) (geom:vy p2o))))
  

(defun draw-bond (scene x1 y1 x2 y2)
  (cl-svg:draw scene (:line :x1 x1 :y1 y1 :x2 x2 :y2 y2
                            :fill "none" :stroke "white" :stroke-width 4))
  (cl-svg:draw scene (:line :x1 x1 :y1 y1 :x2 x2 :y2 y2
                            :fill "none" :stroke "black" :stroke-width 2)))

(defun calculate-towards (a1 a2 aromatic-rings rings)
  1.0)

(defun generate-svg (molecule &key (width 600))
  (let* ((rings (chem:identify-rings molecule))
         (chem:*current-rings* rings)
         (aromatic-rings (chem:identify-aromatic-rings molecule :mdl))
         (bbox (chem:matter-bounding-box molecule 0.0))
         (xwidth (* (- (geom:get-max-x bbox) (geom:get-min-x bbox)) 1.1))
         (xscale (/ width xwidth))
         (ywidth (* (- (geom:get-max-y bbox) (geom:get-min-y bbox)) 1.1))
         (yscale (* (/ ywidth xwidth) xscale))
         (x-mol-center (/ (+ (geom:get-max-x bbox) (geom:get-min-x bbox)) 2.0))
         (y-mol-center (/ (+ (geom:get-max-y bbox) (geom:get-min-y bbox)) 2.0))
         (x-scene-center (/ xwidth 2.0))
         (y-scene-center (/ ywidth 2.0))
         (scene (cl-svg:make-svg-toplevel 'cl-svg:svg-1.2-toplevel :width (* xwidth xscale) :height (* ywidth yscale))))
    (flet ((transform-point (x y)
             (let ((x1 (* (+ (- x x-mol-center) x-scene-center) xscale))
                   (y1 (* (+ (- y y-mol-center) y-scene-center) yscale)))
               (values x1 y1))))
      (multiple-value-bind (bonds atoms)
          (gather-bonds-and-atoms molecule)
        (loop for bond in bonds
              for a1 = (draw-bond-a1 bond)
              for a2 = (draw-bond-a2 bond)
              for order = (draw-bond-order bond)
              for p1 = (chem:get-position a1)
              for p2 = (chem:get-position a2)
              do (multiple-value-bind (x1 y1)
                     (transform-point (geom:vx p1) (geom:vy p1))
                   (multiple-value-bind (x2 y2)
                       (transform-point (geom:vx p2) (geom:vy p2))
                     (draw-bond scene x1 y1 x2 y2)
                     (cond
                       ((eq order :double-bond)
                        (let ((towards (calculate-towards a1 a2 aromatic-rings rings)))
                          (multiple-value-bind (xs1 ys1 xs2 ys2)
                              (calculate-side-bond x1 y1 x2 y2 :towards towards)
                            (draw-bond scene xs1 ys1 xs2 ys2))))
                       ((eq order :triple-bond)
                        (multiple-value-bind (xs1 ys1 xs2 ys2)
                            (calculate-side-bond x1 y1 x2 y2 :towards 1.0)
                          (draw-bond scene xs1 ys1 xs2 ys2))
                        (multiple-value-bind (xs1 ys1 xs2 ys2)
                            (calculate-side-bond x1 y1 x2 y2 :towards -1.0)
                          (draw-bond scene xs1 ys1 xs2 ys2)))
                       (t nil)))))
        (maphash (lambda (atom value)
                   (declare (ignore value))
                   (unless (eq (chem:get-element atom) :c)
                     (multiple-value-bind (xs1 ys1)
                         (let* ((pos (chem:get-position atom))
                                (x (geom:vx pos))
                                (y (geom:vy pos)))
                           (transform-point x y))
                       (cl-svg:draw scene (:circle :cx xs1 :cy ys1 :r 11 :fill "white"))
                       (cl-svg:text scene (:x xs1 :y (+ ys1 00.5) :text-anchor "middle" :alignment-baseline "middle")
                                    (string (chem:get-element atom))))))
                 atoms)))
    (with-output-to-string (sout) (cl-svg:stream-out sout scene))))
    
