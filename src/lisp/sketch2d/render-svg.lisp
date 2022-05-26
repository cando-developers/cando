(in-package :sketch2d)

(defparameter *number-of-hashes-in-hashed-bond* 5)
(defparameter *show-all* nil)
(defparameter *show-names* nil)
(defparameter *perpendicular-fraction* 0.15)
(defparameter *parallel-fraction* 0.1)
(defparameter *dirz* (geom:vec 0.0 0.0 1.0))
(defparameter *pi-bond-width* 5)
(defparameter *wedge-bond-scale* 0.7
  "Scale *pi-bond-width* for wedge bonds")
(defparameter *character-pts* 12)
(defparameter *label-offset* 9)
(defparameter *lower-text* 2)
(defparameter *shorten-pi-bond* 4)
(defparameter *font-pixels* 16.0)
(defparameter *right-h* 11)
(defparameter *right-subh* 10)
(defparameter *left-h* -11)
(defparameter *left-hn* -9)
(defparameter *left-subh* -9)
(defparameter *sub-hn-dy* 6)

(defmacro render-dbg (fmt &rest args)
  nil)

(defclass line ()
  ((p1 :initarg :p1 :accessor p1)
   (p2 :initarg :p2 :accessor p2)
   (class :initform "" :initarg :class :accessor %class)
   (style :initform nil :initarg :style :accessor style :type (member nil :dash :single-wedge-begin :single-hash-begin))
   (color :initform "black" :initarg :color :accessor color)))

(defclass polygon ()
  ((points :initarg :points :accessor points)
   (class :initform "" :initarg :class :accessor %class)
   (color :initform "black" :initarg :color :accessor color)))

(defclass atom-node ()
  ((atom :initarg :atom :accessor atom)
   (debug-info :initarg :debug-info :accessor debug-info)
   (renderp :initform nil :accessor renderp)
   (heavyp :initform nil :accessor heavyp)
   (labelp :initform nil :accessor labelp)
   (bond-nodes :initform nil :initarg :bond-nodes :accessor bond-nodes)
   (rings :initform nil :accessor rings)
   (aromaticity-info :initform nil :accessor aromaticity-info)
   (bond-weight-dir :accessor bond-weight-dir)
   (terminal :initform nil :accessor terminal)
   (element-label :initarg :element-label :accessor element-label)
   (pos :initarg :pos :accessor pos)
   (hydrogens :initform 0 :initarg :hydrogens :accessor hydrogens)
   (hydrogens-dir :initarg :hydrogens-dir :accessor hydrogens-dir)
   (charge :initarg :charge :accessor charge)))

(defclass bond-node ()
  ((bond-order :initarg :bond-order :accessor bond-order)
   (atom-node1 :initarg :atom-node1 :accessor atom-node1)
   (atom-node2 :initarg :atom-node2 :accessor atom-node2)
   (lines :initarg :lines :accessor lines)))

(defclass sketch-svg ()
  ((molecule :initarg :molecule :accessor molecule)
   (sketch2d* :initarg :sketch2d* :accessor sketch2d*)
   (height :initarg :height :accessor height)
   (width :initarg :width :accessor width)
   (sketch-atoms-to-original :initarg :sketch-atoms-to-original :accessor sketch-atoms-to-original)
   (rings :initarg :rings :accessor rings)
   (aromaticity-info :initarg :aromaticity-info :accessor aromaticity-info)
   (atoms-to-nodes :initarg :atoms-to-nodes :accessor atoms-to-nodes)
   (atom-nodes :initarg :atom-nodes :accessor atom-nodes)
   (bond-nodes :initarg :bond-nodes :accessor bond-nodes)
   (before-render :initarg :before-render :accessor before-render)
   (show-names :initarg :show-names :accessor show-names)
   (after-render :initarg :after-render :accessor after-render)
   (scene :initarg :scene :accessor scene)
   (character-pts :initform *character-pts* :reader character-pts)))


(defun calculate-bond-geometry (atom-node1 atom-node2)
  (let* ((pos1 (pos atom-node1))
         (pos2 (pos atom-node2))
         (delta12 (geom:v- pos2 pos1))
         (len12 (geom:vlength delta12))
         (dir12 (geom:v* delta12 (/ 1.0 len12)))
         (left (geom:v* (geom:vnormalized (geom:vcross dir12 *dirz*)) *pi-bond-width*))
         (start (if (labelp atom-node1)
                    *label-offset*
                    0.0))
         (stop (if (labelp atom-node2)
                   (- len12 *label-offset*)
                   len12)))
    (values pos1 pos2 delta12 len12 dir12 left start stop)))


(defun calculate-side-bond (p1 p2 offset1 offset2 &key (towards 1.0))
  "Calculate the start and end of a side-bond.  If towards is a vector it must be to the side of 
the main bond and the calculated bond will be on the same side.
If towards is a number then it should be 1.0 or -1.0.
This will place the calculated bond on one or the other side of the x1,y1-x2,y2 bond."
  (let* ((center (geom:v* (geom:v+ p1 p2) 0.5)) ; center of bond
         (d1-2 (geom:v- p2 p1))         ; bond direction vector
         (dz (geom:vec 0.0 0.0 1.0))    ; z-axis
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
    (values p1o p2o)))

(defgeneric draw-bond (scene shape))

(defmethod draw-bond (scene (line line))
  (let ((p1 (p1 line))
        (p2 (p2 line)))
    (cl-svg:draw scene (:line :x1 (geom:get-x p1) :y1 (geom:get-y p1)
                              :x2 (geom:get-x p2) :y2 (geom:get-y p2)
                              :class (format nil "~A~@[ ~(~A~)~]" (%class line) (style line))))))


(defmethod draw-bond (scene (obj polygon))
  (cl-svg:draw scene (:polygon :points (format nil "~{~5f~^ ~}" (mapcan (lambda (vec)
                                                                          (list (geom:get-x vec) (geom:get-y vec)))
                                                                        (points obj)))
                               :class (%class obj))))

(defun draw-atom-name (scene atom-node show-names)
  (when show-names
    (if (and (<= (chem:get-atomic-number (atom atom-node)) 1)
             (not (eq show-names :light)))
        (return-from draw-atom-name)
        (let* ((pos (pos atom-node))
               (xs1 (geom:get-x pos))
               (ys1 (geom:get-y pos))
               (label #+debug-sketch2d (string (if (debug-info atom-node)
                                                   (label (debug-info atom-node))
                                                   (chem:get-name (atom atom-node))))
                      #-debug-sketch2d (string (chem:get-name (atom atom-node)))))
          (cl-svg:text scene (:x xs1 :y (- ys1 *lower-text*)
                              :class (format nil "atom name ~A" (chem:get-name (atom atom-node))))
                       label)))))

(defun draw-atom-text (scene atom-node)
  (let* ((pos (pos atom-node))
         (xs1 (geom:get-x pos))
         (ys1 (geom:get-y pos))
         (label (string (chem:get-element (atom atom-node)))))
    (cl-svg:text scene (:x xs1 :y (+ ys1 *lower-text*)
                        :class (format nil "atom element ~A" (chem:get-name (atom atom-node))))
                 label)
    (when (and (hydrogens atom-node) (> (hydrogens atom-node) 0))
      (multiple-value-bind (hdx hdy hsubdx hsubdy)
          (cond
            ((= 1 (hydrogens atom-node))
             (case (hydrogens-dir atom-node)
               (:right (values (* *right-h* (length label)) 0 0 0))
               (:left (values *left-h* 0 0 0))
               (otherwise (error "Illegal hydrogens-dir ~a" (hydrogens-dir atom-node)))))
            ((> (hydrogens atom-node) 1)
             (case (hydrogens-dir atom-node)
               (:right (values (* *right-h* (length label)) 0 (+ (* *right-h* (length label)) *right-subh*) *sub-hn-dy*))
               (:left (values (+ *left-hn* *left-subh*) 0 *left-subh* *sub-hn-dy*))
               (otherwise (error "Illegal hydrogens-dir ~a" (hydrogens-dir atom-node))))))
        (cl-svg:text scene (:x xs1 :y (+ ys1 *lower-text*) :dx hdx :dy hdy :class "atom element hydrogen") "H")
        (when (> (hydrogens atom-node) 1)
          (cl-svg:text scene (:x xs1 :y (+ ys1 *lower-text*) :dx hsubdx :dy hsubdy :class "atom element hydrogen") (format nil "~a" (hydrogens atom-node))))))))

(defun draw-hydrogens (scene xs1 ys1 num-hydrogens)
  (case num-hydrogens
    (0)
    (1 (cl-svg:text scene (:x xs1 :y ys1 :dx 12 :class "atom element hydrogen") "H"))
    (otherwise
      (cl-svg:text scene (:x xs1 :y ys1 :dx 24 :class "atom element hydrogen")
                   (format nil "~a" num-hydrogens)))))


(defun get-atom-node (atom atom-ht sketch2d)
  (let ((found (gethash atom atom-ht)))
    (if found
        found
        (let* ((atom-debug-info (gethash atom (debug-info sketch2d)))
               (new-atom-node (make-instance 'atom-node :atom atom :debug-info atom-debug-info)))
          (setf (gethash atom atom-ht) new-atom-node)
          new-atom-node))))

(defun generate-sketch (sketch2d width height)
  "Return the bonds and atoms that we need to render"
  (let ((molecule (molecule sketch2d))
        bonds
        (atoms (make-hash-table)))
    (chem:map-bonds
     nil
     (lambda (a1 a2 bo)
       (push (make-instance 'bond-node
                            :bond-order bo
                            :atom-node1 (get-atom-node a1 atoms sketch2d)
                            :atom-node2 (get-atom-node a2 atoms sketch2d))
             bonds))
     molecule)
    (let (atom-node-list)
      (maphash (lambda (atom dummy)
                 (declare (ignore dummy))
                 (push (get-atom-node atom atoms sketch2d) atom-node-list))
               atoms)
      ;; Tell the atom-nodes their bond-nodes
      (loop for bond-node in bonds
            do (push bond-node (bond-nodes (atom-node1 bond-node)))
               (push bond-node (bond-nodes (atom-node2 bond-node))))
      (let* ((rings (rings sketch2d))
             (aromaticity-info (aromaticity-info sketch2d)))
        (loop for ring in rings
              do (loop for atom in ring
                       for atom-node = (gethash atom atoms)
                       do (push ring (rings atom-node))))
        (make-instance 'sketch-svg
                       :molecule molecule
                       :sketch2d* sketch2d
                       :width width
                       :height height
                       :sketch-atoms-to-original (sketch-atoms-to-original sketch2d)
                       :rings rings
                       :atoms-to-nodes atoms
                       :aromaticity-info aromaticity-info
                       :atom-nodes atom-node-list
                       :bond-nodes bonds)))))

(defun other-atom-node (bond-node atom-node)
  (cond
    ((eq (atom-node1 bond-node) atom-node)
     (atom-node2 bond-node))
    ((eq (atom-node2 bond-node) atom-node)
     (atom-node1 bond-node))
    (t (error "other-atom-node did not match ~a in ~a" atom-node bond-node))))

(defun layout-sketch (sketch transform)
  (loop for atom-node in (atom-nodes sketch)
        for atom = (atom atom-node)
        do (setf (pos atom-node) (funcall transform (chem:get-position atom)))))

(defun center-of-ring (atoms-to-nodes ring)
  (let ((pos (geom:vec 0.0 0.0 0.0))
        (count 0))
    (loop for atom in ring
          for atom-node = (gethash atom atoms-to-nodes)
          do (setf pos (geom:v+ pos (pos atom-node))
                   count (1+ count)))
    (if (> count 0)
        (geom:v* pos (/ 1.0 count))
        (error "There was no ring"))))

(defun draw-double-bond (atoms-to-nodes pos dir left start stop atom-node1 atom-node2 rings &optional aromaticp)
  (let ((in-rings (loop for ring in rings
                        when (and (member (atom atom-node1) ring)
                                  (member (atom atom-node2) ring))
                          collect ring))
        (pi-side :centered)
        (atom1-side (geom:vdot left (bond-weight-dir atom-node1)))
        (atom2-side (geom:vdot left (bond-weight-dir atom-node2)))
        (bond-class (format nil "bond ~A ~A" (chem:get-name (atom atom-node1)) (chem:get-name (atom atom-node2)))))
    (render-dbg "atom1-side atom2-side ~a ~a~%" atom1-side atom2-side)
    (if (< (abs atom1-side) 20.0) (setf atom1-side 0.0))
    (if (< (abs atom2-side) 20.0) (setf atom2-side 0.0))
    (render-dbg "after atom1-side atom2-side ~a ~a~%" atom1-side atom2-side)
    (let* ((pi-side
             (cond
               ((= 1 (length in-rings))
                (let* ((ring-center (center-of-ring atoms-to-nodes (first in-rings)))
                       (towards-ring-center (geom:v- ring-center (geom:v+ pos (geom:v* dir (/ (- stop start) 2.0)))))
                       (side (float-sign (geom:vdot towards-ring-center left))))
                  (render-dbg "pi-side 1~%")
                  (if (> side 0)
                      :left
                      :right))
                )
               ((< (* atom1-side atom2-side) 0.0)
                (render-dbg "pi-side 2~%")
                :left                ; bond-weight-dir are conflicting
                )
               ((= 0.0 atom1-side atom2-side)
                (render-dbg "pi-side 3~%")
                (if (or (terminal atom-node1)
                        (terminal atom-node2))
                    :centered
                    :left               ; pick one - it looks better
                    )
                )
               ((> (abs atom1-side) 20.0)
                (render-dbg "pi-side 4~%")
                (if (> (float-sign atom1-side) 0.0)
                    :left
                    :right)
                )
               ((> (abs atom2-side) 0.0)
                (render-dbg "pi-side 5~%")
                (if (> (float-sign atom2-side) 0.0)
                    :left
                    :right)
                )
               )
             ))
      (multiple-value-bind (pi1 pi2 shorten1 shorten2)
          (case pi-side
            (:centered (values 0.5 -0.5 0.0 0.0))
            (:left (values 0.0 1.0
                           (if (labelp atom-node1) 0.0 *shorten-pi-bond*)
                           (if (labelp atom-node2) 0.0 *shorten-pi-bond*)))
            (:right (values 0.0 -1.0
                            (if (labelp atom-node1) 0.0 *shorten-pi-bond*)
                            (if (labelp atom-node2) 0.0 *shorten-pi-bond*))))
        (list
         #+(or)(make-instance 'line
                              :p1 pos
                              :p2 (geom:v+ pos (geom:v* left 10.0))
                              :color "blue")
         (make-instance 'line
                        :p1 (geom:v+ (geom:v* left pi1) (geom:v+ pos (geom:v* dir start)))
                        :p2 (geom:v+ (geom:v* left pi1) (geom:v+ pos (geom:v* dir stop)))
                        :class bond-class)
         (if aromaticp
             (make-instance 'line
                            :style :dash
                            :p1 (geom:v+ (geom:v* left pi2) (geom:v+ pos (geom:v* dir (+ start shorten1))))
                            :p2 (geom:v+ (geom:v* left pi2) (geom:v+ pos (geom:v* dir (- stop shorten2))))
                            :class bond-class)
             (make-instance 'line
                            :p1 (geom:v+ (geom:v* left pi2) (geom:v+ pos (geom:v* dir (+ start shorten1))))
                            :p2 (geom:v+ (geom:v* left pi2) (geom:v+ pos (geom:v* dir (- stop shorten2))))
                            :class bond-class)))))))

(defun lookup-chiral-info (sketch atom-node1 atom-node2)
  (let ((atom1 (atom atom-node1))
        (atom2 (atom atom-node2)))
    (loop for ci in (chiral-infos (sketch2d* sketch))
          for chiral-sketch-atom = (chiral-sketch-atom ci)
          for bond-atoms = (chiral-bonds ci)
          do (loop for ba in bond-atoms
                   for bond-atom = (bond-atom ba)
                   when (and (eq chiral-sketch-atom atom1)
                             (eq bond-atom atom2))
                     do (return-from lookup-chiral-info
                          (values ci atom-node1 atom-node2 (bond-type ba)))
                   when (and (eq chiral-sketch-atom atom2)
                             (eq bond-atom atom1))
                     do (return-from lookup-chiral-info
                          (values ci atom-node2 atom-node1 (bond-type ba)))))
    (values nil)))

(defun calculate-bond (bond-node sketch)
  (Let* ((atom-node1 (atom-node1 bond-node))
         (atom-node2 (atom-node2 bond-node))
         (bond-class (format nil "bond ~A ~A" (chem:get-name (atom atom-node1)) (chem:get-name (atom atom-node2)))))
    (when (or *show-all*
              (and (renderp atom-node1)
                   (renderp atom-node2)))
      (let ((bond-order (bond-order bond-node)))
        (labels ((calc-pos (pos dir dist)
                   (geom:v+ pos (geom:v* dir dist))))
          (if (member bond-order (list :single-bond :single-wedge-begin :single-hash-begin))
              (multiple-value-bind (chiral-info from-atom-node to-atom-node bond-type)
                  (lookup-chiral-info sketch atom-node1 atom-node2)
                (if chiral-info
                    (multiple-value-bind (pos1 pos2 delta12 len12 dir12 left start stop)
                        (calculate-bond-geometry from-atom-node to-atom-node)
                      (let* ((p1 pos1)
                             (p2 (calc-pos pos1 dir12 stop))
                             (c2 (geom:v+ p2 (geom:v* left *wedge-bond-scale* )))
                             (c3 (geom:v+ p2 (geom:v* left (- *wedge-bond-scale*)))))
                        (case bond-type
                          (:single-wedge-begin
                           #+(or)(when (/= 0.0 start)
                             (warn "The start of a wedge bond is not zero - that means we need a more complex polygon than a triangle start: ~a" start))
                           (list (make-instance 'polygon :points (list p1 c2 c3) :class (concatenate 'string bond-class " wedge"))))
                          (:single-hash-begin
                           #+(or)(when (/= 0.0 start)
                             (warn "The start of a wedge bond is not zero - that means we need a more complex polygon than a triangle start: ~a" start))
                           (let ((hashes (loop for index from 1 below *number-of-hashes-in-hashed-bond*
                                               for frac-offset = (/ (float index) *number-of-hashes-in-hashed-bond*)
                                               for abs-offset = (* len12 frac-offset)
                                               for plus-left = (geom:v* left frac-offset)
                                               for minus-left = (geom:v* left (- frac-offset))
                                               for offset = (calc-pos pos1 dir12 abs-offset)
                                               for cc2 = (geom:v+ plus-left offset)
                                               for cc3 = (geom:v+ minus-left offset)
                                               collect (make-instance 'line :p1 cc2 :p2 cc3 :class bond-class))))
                             hashes))
                          (otherwise (warn "single-bond Handle bond type: ~a~%" (bond-type chiral-info))
                           (list (make-instance 'polygon :points (list p1 c2 c3) :class bond-class))))))
                    (multiple-value-bind (pos1 pos2 delta12 len12 dir12 left start stop)
                        (calculate-bond-geometry atom-node1 atom-node2)
                      (list (make-instance 'line :p1 (calc-pos pos1 dir12 start) ; A simple single-line bond
                                                 :p2 (calc-pos pos1 dir12 stop)
                                                 :class bond-class)))))
              (multiple-value-bind (pos1 pos2 delta12 len12 dir12 left start stop)
                  (calculate-bond-geometry atom-node1 atom-node2)
                (cond
                  ((eq :unknown-order-bond bond-order)
                   (list (make-instance 'line :p1 (calc-pos pos1 dir12 start)
                                              :p2 (calc-pos pos1 dir12 stop)
                                              :class bond-class)))
                  ((eq :double-bond bond-order)
                   (draw-double-bond (atoms-to-nodes sketch)
                                     pos1 dir12
                                     left
                                     start stop
                                     atom-node1 atom-node2
                                     (rings sketch)))
                  ((eq :aromatic-bond bond-order)
                   (draw-double-bond (atoms-to-nodes sketch)
                                     pos1 dir12
                                     left
                                     start stop
                                     atom-node1 atom-node2
                                     (rings sketch)
                                     t))
                  ((eq :triple-bond bond-order)
                   (let ((p1 (calc-pos pos1 dir12 start))
                         (p2 (calc-pos pos1 dir12 stop)))
                     (list
                      (make-instance 'line :p1 p1 :p2 p2 :class bond-class)
                      (make-instance 'line :p1 (geom:v+ p1 left) :p2 (geom:v+ p2 left) :class bond-class)
                      (make-instance 'line :p1 (geom:v- p1 left) :p2 (geom:v- p2 left) :class bond-class))))
                  (t (warn "In render-svg calculate-bond handle bond ~a" bond-order))))))))))

(defun calculate-bonds (sketch)
  (let ((lines nil))
    (loop for bond-node in (bond-nodes sketch)
          for lines = (calculate-bond bond-node sketch)
          do (setf (lines bond-node) lines))))
                 
(defun optimize-sketch (sketch)
  ;; Set up which atoms we want to render and which ones we want to label
  (loop for atom-node in (atom-nodes sketch)
        for atom = (atom atom-node)
        unless (or (eq :lp (chem:get-element atom))
                   (= 1 (chem:get-atomic-number atom)))
          do (setf (renderp atom-node) t
                   (heavyp atom-node) t)
        when (not (eq (chem:get-element atom) :c))
          do (setf (labelp atom-node) t))
  ;; Calculate a direction vector for each atom that
  ;;  points in the average direction of the remaining
  ;;  that will be rendered.  This is to hint what side
  ;;  double bonds should render.
  (loop for atom-node in (atom-nodes sketch)
        for atom = (atom atom-node)
        for bond-nodes = (bond-nodes atom-node)
        for dir = (let ((dir (geom:vec 0.0 0.0 0.0))
                        (count 0))
                    (loop for bond-node in bond-nodes
                          for other-atom-node = (other-atom-node bond-node atom-node)
                          when (heavyp other-atom-node)
                            do (setf dir (geom:v+ dir (geom:v- (pos other-atom-node) (pos atom-node)))
                                     count (1+ count)))
                    (if (> count 0)
                        (geom:v* dir (/ 1.0 count))
                        dir))
        do (setf (bond-weight-dir atom-node) dir))
  ;; identify terminal atoms
  (loop for atom-node in (atom-nodes sketch)
        for atom = (atom atom-node)
        for bond-nodes = (bond-nodes atom-node)
        for other-bonds-count = 0
        do (loop for bond-node in bond-nodes
                 for other-atom-node = (other-atom-node bond-node atom-node)
                 when (heavyp other-atom-node)
                   do (incf other-bonds-count))
        do (when (= 1 other-bonds-count)
             (setf (terminal atom-node) t)))
  ;; figure out attached hydrogens
  (loop for atom-node in (atom-nodes sketch)
        for properties = (chem:properties (atom atom-node))
        for hydrogens = (getf properties :hydrogens)
        for bond-weight-dir = (bond-weight-dir atom-node)
        do (setf (hydrogens atom-node) hydrogens)
        do (if (> (geom:get-x bond-weight-dir) 0.0)
               (setf (hydrogens-dir atom-node) :left)
               (setf (hydrogens-dir atom-node) :right)))
  (calculate-bonds sketch)
  )


(defgeneric render-node (scene node))

(defmethod render-node (scene (node atom-node))
  (draw-atom-name scene node *show-names*)
  (when (and (or *show-all* (renderp node)) (labelp node))
    (draw-atom-text scene node)))

(defmethod render-node (scene (node bond-node))
  (loop for line in (lines node)
        do (draw-bond scene line)))


(defvar *svg-stylesheet*
  ".name {
     font: italic 9px sans-serif;
     fill: red;
     text-anchor: left;
     alignment-baseline: left;
   }
   .element {
     text-anchor: middle;
     alignment-baseline: middle;
   }
   .bond {
     fill: none;
     stroke: black;
     stroke-width: 2;
     stroke-linecap: round;
   }
   .dash {
     stroke-dasharray: 5, 5;
   }
   .bond.wedge {
     fill: black;
     stroke: black;
     stroke-width: 1;
   }")


(cl-svg::define-element-maker :filter "filter" '(:id))
(cl-svg::define-defs-group-maker make-filter :filter)

(cl-svg::define-element-maker :feMorphology "feMorphology" '())
(cl-svg::define-element-maker :feFlood "feFlood" '())
(cl-svg::define-element-maker :feComposite "feComposite" '())
(cl-svg::define-element-maker :feMerge "feMerge" '())
(cl-svg::define-element-maker :feMergeNode "feMergeNode" '())


(defmacro feMerge ((&rest opts) &body inputs)
  (let ((feMerge (gensym "feMerge")))
    `(let ((,feMerge (cl-svg:draw* (:feMerge))))
       (dolist (input (list ,@inputs))
         (if input
           (cl-svg:draw ,feMerge (:feMergeNode :in input))
           (cl-svg:draw ,feMerge (:feMergeNode))))
       ,feMerge)))


(defun render-sketch (scene sketch)
  (make-filter scene (:id "highlight" :x "-1" :y "-1" :width "3" :height "3")
    (cl-svg:draw* (:feMorphology :operator "dilate" :radius "3" :result "dilate"))
    (cl-svg:draw* (:feFlood :flood-color "rgb(237,212,0)"))
    (cl-svg:draw* (:feComposite :in2 "dilate" :operator "in"))
    (feMerge () nil "SourceGraphic"))
  (cl-svg:style scene *svg-stylesheet*)
  (loop for bond-node in (bond-nodes sketch)
        for count from 0
        do (render-dbg "bond-node #~a~%" count)
        do (render-node scene bond-node))
  (loop for atom-node in (atom-nodes sketch)
        do (render-node scene atom-node)))


(defun svg (sketch2d &key (toplevel t) (width 1000) (xbuffer 0.1) (ybuffer 0.1) before-render after-render (id "") show-names (scale 20) (margin 40))
  "Generate SVG to render the molecule.  Pass a BEFORE-RENDER function or AFTER-RENDER function to add info to the structure.
Each of these functions take two arguments, the svg-scene and the sketch-svg. 
The caller provided functions should use cl-svg to render additional graphics."
  (let* ((molecule (molecule sketch2d))
         (xmin most-positive-single-float)
         (ymin most-positive-single-float)
         (xmax most-negative-single-float)
         (ymax most-negative-single-float))
    (flet ((transform-point (atomic-pos)
             (let ((x (* scale (geom:get-x atomic-pos)))
                   (y (* scale (geom:get-y atomic-pos))))
               (setf xmin (min xmin x))
               (setf xmax (max xmax x))
               (setf ymin (min ymin y))
               (setf ymax (max ymax y))
               (geom:vec x y 0.0))))
      (let ((sketch (generate-sketch sketch2d 0 0)))
        (layout-sketch sketch #'transform-point)
        (setf xmin (floor (- xmin margin))
              ymin (floor (- ymin margin))
              xmax (ceiling (+ xmax margin))
              ymax (ceiling (+ ymax margin))
              (width sketch) (- xmax xmin)
              (height sketch) (- ymax ymin))
        (optimize-sketch sketch)
        (render-dbg "About to render-sketch  bond-nodes ~a~%" (length (bond-nodes sketch)))
        ;; It looks like for jupyter we need to specify width and height
        ;; Here it says to ignore them?  https://css-tricks.com/scale-svg/
        (let* ((view-box (format nil "~a ~a ~a ~a" xmin ymin (width sketch) (height sketch)))
               (width (format nil "~apx" (width sketch)))
               (height (format nil "~apx" (height sketch)))
               (scene (if toplevel
                        (cl-svg:make-svg-toplevel 'cl-svg:svg-1.2-toplevel
                                                  :id id
                                                  :width width
                                                  :height height
                                                  "viewBox" view-box)
                        (make-instance 'cl-svg::svg-element
                                       :name "svg"
                                       :attributes (list :id id
                                                         :width width
                                                         :height height
                                                         "viewBox" view-box)))))
          (setf (scene sketch) scene
                (before-render sketch) before-render
                (show-names sketch) show-names
                (after-render sketch) after-render)
          (values sketch))))))

(defun render-svg-scene (sketch-svg)
  (let ((scene (scene sketch-svg)))
    (when (before-render sketch-svg)
      (funcall (before-render sketch-svg) scene sketch-svg))
    (let ((*show-names* (or *show-names* (show-names sketch-svg))))
      (render-sketch scene sketch-svg))
    (when (after-render sketch-svg)
      (funcall (after-render sketch-svg) scene sketch-svg))
    scene))

(defun render-svg-to-string (sketch-svg)
  (let ((scene (render-svg-scene sketch-svg)))
    (with-output-to-string (sout) (cl-svg:stream-out sout scene))))



