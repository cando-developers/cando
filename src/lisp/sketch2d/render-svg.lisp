(in-package :sketch2d)

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
         (magnitude-towards (if (numberp towards)
                                (abs towards)
                                1.0))
         (dtowards (if (numberp towards)
                       (geom:v* dcross towards)
                       (geom:v- towards center))) ; if no towards then use pcross
         (dot (float-sign (geom:vdot dcross dtowards) 1.0))
         (len1-2 (geom:vlength d1-2))
         (offset-perpendicular (geom:v* (geom:v* dcross (* *perpendicular-fraction* len1-2 dot)) magnitude-towards))
         (offset-parallel (geom:v* d1-2 *parallel-fraction*))
         (p1o (geom:v+ (geom:v+ p1 offset-parallel) offset-perpendicular))
         (p2o (geom:v+ (geom:v- p2 offset-parallel) offset-perpendicular)))
    (values (geom:vx p1o) (geom:vy p1o) (geom:vx p2o) (geom:vy p2o))))

(defparameter *font-pixels* 16.0)

(defun draw-bond (scene l1 x1 y1 l2 x2 y2)
  (let* ((p1e (geom:vec x1 y1 0.0))
         (p2e (geom:vec x2 y2 0.0))
         (phalf (geom:v* (geom:v+ p1e p2e) 0.5))
         (delta (geom:vnormalized (geom:v- p2e p1e)))
         (p1bg (geom:v+ p1e (geom:v* delta *font-pixels*)))
         (p2bg (geom:v+ p2e (geom:v* delta *font-pixels*)))
         )
    (format t "delta * 16 -> ~a~%" (geom:v* delta 16.0))
    (cl-svg:draw scene (:line :x1 (geom:vx p1bg) :y1 (geom:vy p1bg)
                              :x2 (geom:vx p2bg) :y2 (geom:vy p2bg)
                              :fill "none" :stroke "white" :stroke-width 4))
    (cl-svg:draw scene (:line :x1 x1 :y1 y1 :x2 x2 :y2 y2
                              :fill "none" :stroke "black" :stroke-width 2
                              :stroke-linecap "round"))))

(defun calculate-towards (a1 a2 aromatic-rings rings)
  1.0)

(defun draw-atom-text (scene xs1 ys1 label)
  (format t "text at ~a ~a~%" xs1 ys1)
  (cl-svg:text scene (:x xs1 :y ys1 :text-anchor "middle" :alignment-baseline "middle") label))

(defun draw-hydrogens (scene xs1 ys1 num-hydrogens)
  (when (> num-hydrogens 0)
    (cl-svg:text scene (:x xs1 :y ys1 :dx 12 :text-anchor "middle" :alignment-baseline "middle") "H")
    (when (> num-hydrogens 1)
      (cl-svg:text scene (:x xs1 :y ys1 :dx 24 :text-anchor "middle" :alignment-baseline "middle") 
                   (format nil "~a" num-hydrogens)))))


(defclass atom-node ()
  ((the-atom :initarg :the-atom :accessor the-atom)
   (render :initform nil :accessor render)
   (label :initform nil :accessor label)
   (bond-nodes :initform nil :initarg :bond-nodes :accessor bond-nodes)
   (element-label :initarg :element-label :accessor element-label)
   (pos :initarg :pos :accessor pos)
   (hydrogens :initarg :hydrogens :accessor hydrogens)
   (hydrogens-side :initarg :hydrogens-position :accessor hydrogens-position)
   (charge :initarg :charge :accessor charge)))

(defclass bond-node ()
  ((bond-order :initarg :bond-order :accessor bond-order)
   (atom-node1 :initarg :atom-node1 :accessor atom-node1)
   (atom-node2 :initarg :atom-node2 :accessor atom-node2)
   (pos1 :initarg :pos1 :accessor :pos1)
   (pos2 :initarg :pos2 :accessor :pos2)))

(defclass sketch ()
  ((molecule :initarg :molecule :accessor molecule)
   (rings :initarg :rings :accessor rings)
   (atom-nodes :initarg :atom-nodes :accessor atom-nodes)
   (bond-nodes :initarg :bond-nodes :accessor bond-nodes)))

(defun get-atom-node (atom atom-ht)
  (let ((found (gethash atom atom-ht)))
    (if found
        found
        (let ((new-atom-node (make-instance 'atom-node :the-atom atom)))
          (setf (gethash atom atom-ht) new-atom-node)
          new-atom-node))))

(defun generate-sketch (molecule)
  "Return the bonds and atoms that we need to render"
  (let (bonds
        (atoms (make-hash-table)))
    (chem:map-bonds
     nil
     (lambda (a1 a2 bo)
       (push (make-instance 'bond-node
                            :bond-order bo
                            :atom1 (get-atom-node a1 atoms)
                            :atom2 (get-atom-node a2 atoms))
             bonds))
     molecule)
    (let (atom-list)
      (maphash (lambda (atom dummy)
                 (declare (ignore dummy))
                 (push atom atom-list))
               atoms)
      ;; Tell the atom-nodes their bond-nodes
      (loop for bond-node in bonds
            do (push bond-node (bond-nodes (atom-node1 bond-node)))
               (push bond-node (bond-nodes (atom-node2 bond-node))))
      (let* ((rings (chem:identify-rings molecule)))
        (error "Handle rings differently - they need to be discovered on the original molecule and mapped to edited-mol")
        (loop for ring in rings
              do (loop for atom in ring
                       for atom-node = (gethash atom atoms)
                       do (push ring (rings atom-node))))
        (make-instance 'sketch
                       :molecule molecule
                       :rings rings
                       :aromatic-rings aromatic-rings
                       :atom-nodes atom-list
                       :bond-nodes bonds)))))

(defun optimize-sketch (sketch)
  ;; Set up which atoms we want to render and which ones we want to label
  (loop for atom in (atom-nodes sketch)
        unless (or (eq :lp (chem:get-element atom))
                   (= 1 (chem:get-atomic-number atom)))
          do (setf (renderp atom) t)
        when (not (eq (chem:get-element atom) :c))
          do (setf (labelp atom) t))
  ;; Calculate a direction vector for each atom that
  ;;  points in the average direction of the remaining
  ;;  that will be rendered.  This is to hint what side
  ;;  double bonds should render.
  (loop for atom-node in (atom-nodes sketch)
        for atom = (the-atom atom-node)
        for bond-nodes = (bond-nodes atom-node)
        for dir = (let ((dir (geom:vec 0.0 0.0 0.0))
                        (count 0))
                    (loop for bond-node in bond-nodes
                          for other-atom-node = (other-atom-node bond-node atom-node)
                          for other-atom = (the-atom other-atom-node)
                          when (renderp other-atom-node)
                            do (setf dir (geom:v+ dir (chem:get-position other-atom))
                                     count (1+ count)))
                    (if (> count 0)
                        (geom:v* dir (/ 1.0 count))
                        dir))
        do (setf (bond-weight-dir atom-node) dir))
  ;; Calculate the side for double bonds
  ;;  1) If double bond is in one aromatic ring put it in the ring, if in two rings then 
  ;;  2)
  ;; For each bond draw the bond lines
  (loop for bond-node in (bond-nodes sketch)
        for atom-node1 = (atom-node1 bond-node)
        for atom-node2 = (atom-node2 bond-node)
        for atom1 = (the-atom atom-node1)
        for atom2 = (the-atom atom-node2)
        for pos1 = (chem:get-position atom1)
        for pos2 = (chem:get-position atom2)
        ;; temporary - just set pos of atoms
        do (setf (pos1 bond-node) pos1)
           (setf (pos2 bond-node) pos2))
  )


(defun render-sketch (scene sketch)
  (loop for bond-node in (bond-nodes sketch)
        do (render-node scene bond-node))
  (loop for atom-node in (atom-nodes sketch)
        do (render-node scene atom-node)))


(defun svg (sketch2d &key (width 1000) (xbuffer 0.1) (ybuffer 0.1) )
  (let* ((molecule (molecule sketch2d))
         (rings (chem:identify-rings molecule))
         (chem:*current-rings* rings)
         (aromatic-rings (chem:identify-aromatic-rings molecule :mdl))
         (bbox (chem:matter-bounding-box molecule 0.0))
         (xscale 20.0)
         (yscale 20.0)
         (xviewport (* (- (geom:get-max-x bbox) (geom:get-min-x bbox)) (+ 1.0 (* 2.0 xbuffer)) xscale))
         (yviewport (* (- (geom:get-max-y bbox) (geom:get-min-y bbox)) (+ 1.0 (* 2.0 ybuffer)) yscale))
         (x-mol-center (/ (+ (geom:get-max-x bbox) (geom:get-min-x bbox)) 2.0))
         (y-mol-center (/ (+ (geom:get-max-y bbox) (geom:get-min-y bbox)) 2.0))
         (x-viewport-center (/ xviewport 2.0))
         (y-viewport-center (/ yviewport 2.0))
         (height (* width (/ yviewport xviewport)))
         ;; It looks like for jupyter we need to specify width and height
         ;; Here it says to ignore them?  https://css-tricks.com/scale-svg/
         (scene (cl-svg:make-svg-toplevel 'cl-svg:svg-1.2-toplevel :width (round xviewport) :height (round yviewport) ; :width "1000" :height "auto" ; :width width :height height
                                                                   :viewport (format nil "0 0 ~d ~d" (round xviewport) (round yviewport)))))

    (flet ((transform-point (x y)
             (let ((x1 (+ (* xscale (- x x-mol-center)) x-viewport-center))
                   (y1 (+ (* yscale (- y y-mol-center)) y-viewport-center)))
               (values x1 y1))))
      #+(or)(cl-svg:draw scene (:rect :x 0 :y 0 :width xviewport :height yviewport :stroke "green" :fill "green"))
      (let ((sketch (generate-sketch molecule)))
        (optimize-sketch sketch)
        (render-sketch scene sketch))
      #+(or)(loop for bond in (bonds sketch)
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
      #+(or)
      (maphash (lambda (atom value)
                 (declare (ignore value))
                 (let ((properties (chem:properties atom)))
                   (unless (eq (chem:get-element atom) :c)
                     (multiple-value-bind (xs1 ys1)
                         (let* ((pos (chem:get-position atom))
                                (x (geom:vx pos))
                                (y (geom:vy pos)))
                           (transform-point x y))
                       (draw-atom-text scene xs1 ys1 (string (chem:get-element atom)))
                       (let ((hydrogens (getf properties :hydrogens)))
                         (when (integerp hydrogens)
                           (draw-hydrogens scene xs1 ys1 hydrogens)))))))
               atoms)))
  (with-output-to-string (sout) (cl-svg:stream-out sout scene)))
