
(in-package :geometry-tests)


(defclass ring-plane ()
  ((centroid :initarg :centroid :accessor centroid)
   (normal :initarg :normal :accessor normal)
   (distance-to-farthest-point :initarg :distance-to-farthest-point :accessor distance-to-farthest-point)
   (atoms :initarg :atoms :accessor atoms)
  ))

(defmethod print-object ((object ring-plane) stream)
  (print-unreadable-object (object stream :type t)
    (format stream ":centroid ~a :normal ~a" (centroid object) (normal object))))



;; Constructs a plane from a collection of atoms
;; so that the summed squared distance to all atoms is minimized
;; From: https://www.ilikebigbits.com/2015_03_04_plane_from_points.html
(defun plane-from-atoms (atoms)
  "Constructs a plane from a collection of atoms
so that the summed squared distance to all atoms is minimized.
Returns the plane or nil"
  (if (< (length atoms) 3)
      nil
      (let ((points (mapcar #'chem:get-position atoms))
            (sum (geom:vec 0.0 0.0 0.0))
            (xx 0.0d0) (xy 0.0d0) (xz 0.0d0)
            (yy 0.0d0) (yz 0.0d0) (zz 0.0d0))
        (declare (double-float xx xy xz yy yz zz))
        (loop for pnt in points
              do (setf sum (geom:v+ sum pnt)))
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
                    ;; Find the ring point that is farthest from the centroid.
                    ;;  Bonds that intersect the ring closer than the farthest point are intersecting the ring.
                    (let ((longest-distance (loop for pnt in points
                                                  for distance = (geom:vdistance pnt centroid)
                                                  maximize distance)))
                      (make-instance 'ring-plane :centroid centroid
                                                 :normal (geom:vnormalized dir)
                                                 :atoms atoms
                                                 :distance-to-farthest-point longest-distance))))))))))

(defun bonds-crossing-ring-plane (plane matter)
  (let (collisions
        (l0 (geom:vec 0 0 0))
        (l1 (geom:vec 0 0 0))
        (l (geom:vec 0 0 0))
        (p (geom:vec 0 0 0))
        (delta (geom:vec 0 0 0))
        (p0 (centroid plane))
        (n (normal plane)))
    (declare (dynamic-extent l0 l1 l p delta))
    ;; Loop over all bonds 
    (chem:map-bonds
     'nil
     (lambda (a0 a1 order)
       (declare (ignoreable order))
       (chem:get-position! a0 l0)
       (chem:get-position! a1 l1)
       (geom:v-! l1 l0 l)
       (unless (or (member a0 (atoms plane))
                   (member a1 (atoms plane)))
         (let* ((p0-l0 (geom:v- p0 l0))
                (l.n (geom:vdot l n)))
           (declare (double-float l.n))
           ;; if l.n == 0.0 then the bond is perpendicular to the normal and
           ;;  cannot be intersecting the ring
           (unless (= l.n 0.0)
             (let* ((top (geom:vdot p0-l0 n))
                    (d (/ top l.n)))
               (declare (double-float top d))
               ;; If d < 0.0 or d > 1.0 then the bond cannot intersect the ring.
               (when (<= 0.0 d 1.0)
                 (geom:v*! l d p)
                 (geom:v+! p l0 p)
                 (geom:v-! p p0 delta)
                 (let ((distance (geom:vlength delta)))
                   (declare (double-float distance))
                   ;; If distance to intersection is shorter than distance to farthest atom
                   ;; then the intersection is in the ring.
                   (when (< distance (distance-to-farthest-point plane))
                     (push (list plane a0 a1) collisions)))))))))
     matter)
    collisions))


(defun perpendicular (direction)
  (let ((vx (geom:vx direction))
        (vy (geom:vy direction))
        (vz (geom:vz direction)))
    (let* ((perp-try1 (geom:vec vz vz (- (- vx) vy)))
           (perp-try2 (geom:vec (- (- vy) vz) vx vx))
           (select-try (and (/= vz 0.0) (/= (- vx) vy)))
           (perp (if select-try perp-try1 perp-try2)))
      perp)))

(defun segment-crosses-plane (plane l0 l1)
  "Return T if the segment defined by l0 and l1 crosses the plane"
  (let* ((l (geom:v- l1 l0)) ; allocate l on the stack
         (p0-l0 (geom:v- (centroid plane) l0)) ; allocate p0-l0 on the stack
         (l.normal (geom:vdot l (normal plane))))
    (declare (dynamic-extent l p0-l0))
    (if (= l.normal 0.0)
        nil
        (let ((d (/ (geom:vdot (geom:vdot p0-l0 (normal plane))) l.normal)))
          (<= 0.0 d 1.0)))))


(defun segment-crosses-plane (plane l0 l1)
  "Return T if the segment defined by l0 and l1 crosses the plane"
  (let* ((l (geom:v- l1 l0)) ; allocate l on the stack
         (p0-l0 (geom:v- p0 l0)) ; allocate p0-l0 on the stack
         (l.normal (geom:vdot l (normal plane))))
    (declare (dynamic-extent l p0-l0))
    (if (= l.normal 0.0)
        nil
        (let ((d (/ (geom:vdot (geom:vdot p0-l0 (normal plane))) l.normal)))
          (<= 0.0 d 1.0)))))


(defun segment-crosses-plane-temps (plane l0 l1 l p0-l0)
  "Return T if the segment defined by l0 and l1 crosses the plane"
  (geom:v-! l l1 l0)
  (geom:v-! p0-l0 p0 l0)
  (let ((l.normal (geom:vdot l (normal plane))))
    (if (= l.normal 0.0)
        nil
        (let ((d (/ (geom:vdot (geom:vdot p0-l0 (normal plane))) l.normal)))
          (<= 0.0 d 1.0)))))

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
        for plane = (plane-from-atoms ring)))


(defun ring-plane-shapes (rings)
  (loop for ring in rings
        for plane = (plane-from-atoms ring)
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
