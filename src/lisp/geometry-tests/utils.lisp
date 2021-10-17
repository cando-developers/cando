
(in-package :geometry-tests)


(defclass plane ()
  ((centroid :initarg :centroid :accessor centroid)
   (normal :initarg :normal :accessor normal)))

(defmethod print-object ((object plane) stream)
  (print-unreadable-object (object stream :type t)
    (format stream ":centroid ~a :normal ~a" (centroid object) (normal object))))



;; Constructs a plane from a collection of points
;; so that the summed squared distance to all points is minimized
;; From: https://www.ilikebigbits.com/2015_03_04_plane_from_points.html
(defun plane-from-points (points)
  "Constructs a plane from a collection of points
so that the summed squared distance to all points is minimized.
Returns (values t centroid normalized-dir) or (values nil)"
  (if (< (length points) 3)
      nil
      (let ((sum (geom:vec 0.0 0.0 0.0))
            (xx 0.0d0) (xy 0.0d0) (xz 0.0d0)
            (yy 0.0d0) (yz 0.0d0) (zz 0.0d0))
        (declare (double-float xx xy xz yy yz zz))
        (loop for vec in points
              do (setf sum (geom:v+ sum vec)))
        (let ((centroid (geom:v* sum (/ 1.0 (length points)))))
          ;; Calculate full 3x3 covariance matrix, excluding symmetries
          (loop for p in points
                for r = (geom:v- p centroid)
                for rx double-float = (geom:vx r)
                for ry double-float = (geom:vy r)
                for rz double-float = (geom:vz r)
                do (incf xx (infix:infix rx * rx))
                do (incf xy (infix:infix rx * ry))
                do (incf xz (infix:infix rx * rz))
                do (incf yy (infix:infix ry * ry))
                do (incf yz (infix:infix ry * rz))
                do (incf zz (infix:infix rz * rz)))
          (let ((det-x (infix:infix yy * zz - yz * yz))
                (det-y (infix:infix xx * zz - xz * xz))
                (det-z (infix:infix xx * yy - xy * xy)))
            (declare (double-float det-x det-y det-z))
            (let ((det-max (max det-x det-y det-z)))
              (declare (double-float det-max))
              (if (<= det-max 0.0d0)
                  nil
                  ;; Pick path with best conditioning
                  (let ((dir (cond
                               ((= det-max det-x)
                                (geom:vec det-x (infix:infix xz * yz - xy * zz) (infix:infix xy * yz - xz * yy)))
                               ((= det-max det-y)
                                (geom:vec (infix:infix xz * yz - xy * zz) det-y (infix:infix xy * xz - yz * xx)))
                               ((= det-max det-z)
                                (geom:vec (infix:infix xy * yz - xz * yy) (infix:infix xy * xz - yz * xx) det-z)))))
                    (make-instance 'plane :centroid centroid :normal (geom:vnormalized dir))))))))))

(defun perpendicular (direction)
  (let ((vx (geom:vx direction))
        (vy (geom:vy direction))
        (vz (geom:vz direction)))
    (let* ((perp-try1 (geom:vec vz vz (- (- vx) vy)))
           (perp-try2 (geom:vec (- (- vy) vz) vx vx))
           (select-try (and (/= vz 0.0) (/= (- vx) vy)))
           (perp (if select-try perp-try1 perp-try2)))
      perp)))


(defun plane-shapes (plane)
  (let* ((centroid (centroid plane))
         (normal (normal plane))
         (normal-pos (geom:v+ centroid (geom:v* normal 2.0))))
    `((:type "sphere"
       :position ,(list (geom:vx centroid) (geom:vy centroid) (geom:vz centroid))
       :color (0 1 1) :radius 0.4)
      (:type "arrow"
       :position1 ,(list (geom:vx centroid) (geom:vy centroid) (geom:vz centroid))
       :position2 ,(list (geom:vx normal-pos) (geom:vy normal-pos) (geom:vz normal-pos))
       :color (0 0 1) :radius 0.2))))


(defun ring-planes (rings)
  (loop for ring in rings
        for points = (mapcar (lambda (atm) (chem:get-position atm)) ring)
        for plane = (plane-from-points points)))


(defun ring-plane-shapes (rings)
  (loop for ring in rings
        for points = (mapcar (lambda (atm) (chem:get-position atm)) ring)
        for plane = (plane-from-points points)
        for plane-shapes = (plane-shapes plane)
        append plane-shapes into shapes
        finally (return shapes)))


(defun log-picks (view matter)
  (let ((stage (loop named find-stage
                     for child in (jw:widget-children view)
                     when (typep child 'ngl:stage)
                       do (return-from find-stage child))))
    (flet ((do-pick (inst pick)
             (let ((px (gethash "x" (gethash "atom" pick)))
                   (py (gethash "y" (gethash "atom" pick)))
                   (pz (gethash "z" (gethash "atom" pick)))
                   )
               (format t "inst = ~a~%" inst)
               (format t "(hash-table-plist pick) = ~s~%" (alexandria:hash-table-plist pick))
               (format t "(hash-table-plist (gethash \"atom\" pick)) -> ~s~%" (alexandria:hash-table-plist (gethash "atom" pick)))
               (format t "(gethash \"x\" (gethash \"atom\" pick)) -> ~s~%" px)
               (format t "(gethash \"y\" (gethash \"atom\" pick)) -> ~s~%" py)
               (format t "(gethash \"z\" (gethash \"atom\" pick)) -> ~s~%" pz)
               (format t "(hash-table-plist (gethash \"position\" pick) -> ~s~%" (alexandria:hash-table-plist (gethash "position" pick)))
               (let ((compare nil))
                 (cando:do-atoms (atm matter)
                   (let* ((dx (abs (- px (geom:vx (chem:get-position atm)))))
                          (dy (abs (- px (geom:vy (chem:get-position atm)))))
                          (dz (abs (- px (geom:vz (chem:get-position atm)))))
                          (dist (sqrt (+ (* dx dx) (* dy dy) (* dz dz)))))
                     (push (cons dist atm) compare)))
                 (let ((sorted-compare (sort compare #'< :key #'car)))
                   (loop for one in sorted-compare
                         do (format t "~s~%" one)))))))
      (setf (ngl:on-pick stage) (list (lambda (inst pick) (do-pick inst pick)))))))
