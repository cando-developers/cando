
;; Variables/Constants

(defconstant spin-speed 5)
(defvar slices 40)    ; number of hori/verti slices on spheres
(defvar model-type 'solid)
(defvar curr-mol)
(defvar view-rotx 20)
(defvar view-roty 30)
(defvar view-rotz 0)
(defvar view-scale 1.0)
(defvar view-style :line)
(defvar motion nil)
(defvar show-light-source nil)
(defvar light-theta 0)
(defvar light-phi 0)
(defvar light-dist 5)
(defconstant light-spin-speed .0872664625) ; 5 degrees in radians
(defvar light-r .8)
(defvar light-g .8)
(defvar light-b .8)
(defconstant light-colour-vel .1)
(defvar walk-mode nil)
(defvar *view-aggregate*)


(defclass lines (glut:window)
  ((buffers :accessor buffers :initform nil)
   (draw-calls :accessor draw-calls :initform nil)
   (a :accessor a :initform 0.0)
   (count :initform 1)
   (t0 :initform 0))
  (:default-initargs :width 1024 :height 1024 :title "lines.lisp"
                     :mode '(:double :rgb :multisample :depth)))

(defconstant +index-bytes+ 4) ;; 32bit indices


(defun gen-vertex-array ()
  (cffi:with-foreign-object (array '%gl:uint 1)
    (%gl::gen-vertex-arrays-apple 1 array)
    (cffi:mem-aref array '%gl:uint 0)))


(defun make-vaos ()
  (let ((vao (gen-vertex-array))
        (vbo (gl:gen-buffer))
        (ivbo (gl:gen-buffer)))
    (list vao vbo ivbo)))

(defun fill-vao (vertices indices vao vbo ivbo)
  (%gl:bind-vertex-array-apple vao)
  (gl:bind-buffer :element-array-buffer ivbo)
  (%gl:buffer-data :element-array-buffer (* (length indices) +index-bytes+)
                   (cffi:null-pointer)
                   :static-draw)
  (gl:with-mapped-buffer (p :element-array-buffer :write-only)
    (loop for i from 0
          for v across indices
          do (setf (cffi:mem-aref p :unsigned-int i) v)))
  (gl:bind-buffer :array-buffer vbo)
  (%gl:buffer-data :array-buffer (* (length vertices) 4)
                   (cffi:null-pointer)
                   :static-draw)
  (gl:with-mapped-buffer (p :array-buffer :write-only)
    (loop for i from 0
          for v across vertices
          do (setf (cffi:mem-aref p :float i) v)))
  ;; if using shaders:
  ;; (gl:vertex-attrib-pointer 0 3 :float nil 0 (cffi:null-pointer))
  ;; (gl:enable-vertex-attrib-array 0)
  ;; no shaders
  (%gl:vertex-pointer 3 :float 0 (cffi:null-pointer))
  (gl:enable-client-state :vertex-array)

  (%gl:bind-vertex-array-apple 0)
  (gl:bind-buffer :array-buffer 0)
  (gl:bind-buffer :element-array-buffer 0)
  )


(defun pair (l)
  (let* ((x1 (random 100.0))
         (y1 (random 100.0))
         (z1 (random 100.0))
         (x2 (+ x1 (- (random l) (/ l 2))))
         (y2 (+ y1 (- (random l) (/ l 2))))
         (z2 (+ z1 (- (random l) (/ l 2)))))
    (list (list x1 y1 z1) (random 12) (list x2 y2 z2) (random 12))))

(defun translate-data-lines (w)
  (let* ((nlines (expt 2 14) ;;(expt 2 23)
                 )
         (nverts (* nlines 3))
         (verts (make-array (* 3 nverts) ;; 3 floats per vertex
                            :element-type 'single-float))
         (indices nil)
         (interned-verts (make-hash-table :test 'equalp))
         (l 5.0)
         (i 0)
         (colors (make-array 12)))
    (loop for i below 12
          do (setf (aref colors i)
                   (make-array 1024 :adjustable t
                                    :fill-pointer 0)))
    (flet ((v (p)
             (or (gethash p interned-verts )
                 (progn
                   (setf (aref verts (+ (* i 3) 0)) (first p))
                   (setf (aref verts (+ (* i 3) 1)) (second p))
                   (setf (aref verts (+ (* i 3) 2)) (third p))
                   (prog1 i
                     (setf (gethash p interned-verts) i)
                     (incf i))))))
      (loop repeat nlines
            for (p1 c1 p2 c2) = (pair l)
            for pm = (mapcar (lambda (x y) (/ (+ x y) 2))
                             p1 p2)
            for s = (v p1)
            for m = (v pm)
            for e = (v p2)
            do (vector-push-extend s (aref colors c1))
               (vector-push-extend m (aref colors c1))
               (vector-push-extend m (aref colors c2))
               (vector-push-extend e (aref colors c2))))
    (setf indices (apply #'concatenate
                         '(simple-array (unsigned-byte 32) (*))
                         (coerce colors 'list)))
    (setf (draw-calls w)
          (loop with start = 0
                for i below 12
                for count = (length (aref colors i))
                collect (list (list (random 1.0) (random 1.0) (random 1.0))
                              start count)
                do (incf start count)))
    (apply #'fill-vao verts indices (buffers w)))
  )

(defun count-bonds (agg)
  (let ((count 0))
    (chem:map-bonds
     nil
     (lambda (x y o b)
       (declare (ignore b))
       (incf count))
     agg)
    count))

(defmacro lim-between (sym bot top)
  `(setq ,sym (max ,bot (min ,top ,sym))))

(defvar *element-color-index*)
(defvar *indexed-color*)
(eval-when (:compile-toplevel :load-toplevel :execute)
  (setf *element-color-index* (make-hash-table :test #'eql))
  (setf *indexed-color* (make-hash-table :test #'eql))
  (setf (gethash 0 *indexed-color*) '(1.0 1.0 1.0))
  (setf (gethash 1 *indexed-color*) '(0.0 1.0 1.0))
  (setf (gethash 2 *indexed-color*) '(0.0 0.0 1.0))
  (setf (gethash 3 *indexed-color*) '(1.0 0.0 0.0))
  (setf (gethash 4 *indexed-color*) '(1.0 1.0 0.0))
                                               
  (setf (gethash :H *element-color-index*) (cons 0 nil))
  (setf (gethash :C *element-color-index*) (cons 1 nil))
  (setf (gethash :N *element-color-index*) (cons 2 nil))
  (setf (gethash :O *element-color-index*) (cons 3 nil))
  (setf (gethash :P *element-color-index*) (cons 4 nil)))
(defparameter *max-colors* 5)


(defun vlist (p) (list (geom:get-x p) (geom:get-y p) (geom:get-z p)))


(defun translate-data-aggregate (w)
  (format t "In translate-data-aggregate~%")
  (let* ((nlines (* (count-bonds *view-aggregate*) 2))
         (nverts (* nlines 3))
         (verts (make-array (* 3 nverts) ;; 3 floats per vertex
                            :element-type 'single-float))
         (indices nil)
         (interned-verts (make-hash-table :test 'equalp))
         (l 5.0)
         (i 0)
         (colors (make-array *max-colors*)))
    (format t "In translate-data-aggregate step 2~%")
    (loop for i below *max-colors*
       do (setf (aref colors i)
                (make-array 1024 :adjustable t
                            :fill-pointer 0)))
    (flet ((v (p)
             (or (gethash p interned-verts )
                 (progn
                   (setf (aref verts (+ (* i 3) 0)) (first p))
                   (setf (aref verts (+ (* i 3) 1)) (second p))
                   (setf (aref verts (+ (* i 3) 2)) (third p))
                   (prog1 i
                     (setf (gethash p interned-verts) i)
                     (incf i))))))
      (chem:map-bonds
       nil
       (lambda (x y o b)
         (declare (ignore b))
         (let* ((p1 (vlist (chem:get-position x)))
                (c1 (car (gethash (chem:get-element x) *element-color-index*)))
                (p2 (vlist (chem:get-position y)))
                (c2 (car (gethash (chem:get-element y) *element-color-index*)))
                (s (v p1))
                (pm (mapcar (lambda (x y) (/ (+ x y) 2))
                            p1 p2))
                (m (v pm))
                (e (v p2)))
           (vector-push-extend s (aref colors c1))
           (vector-push-extend m (aref colors c1))
           (vector-push-extend m (aref colors c2))
           (vector-push-extend e (aref colors c2))))
       *view-aggregate*))
    (setf indices (apply #'concatenate
                         '(simple-array (unsigned-byte 32) (*))
                         (coerce colors 'list)))
    (setf (draw-calls w)
          (loop with start = 0
             for i below *max-colors*
             for count = (length (aref colors i))
             for color = (gethash i *indexed-color*)
             collect (list color
                           start count)
             do (incf start count)))
    (apply #'fill-vao verts indices (buffers w)))
  )

(defun translate-data (w)
  (if *view-aggregate*
      (translate-data-aggregate w)
      (translate-data-lines w)))

(defmethod glut:display-window :before ((w lines))
  (setf (buffers w) (make-vaos))
  (translate-data w)
  (gl:enable :depth-test :multisample
             ;:line-smooth :sample-alpha-to-coverage
             )
  (gl:clear-color 0 0 0 0))



(defun print-frame-rate (window)
  "Prints the frame rate every ~5 seconds."
  (with-slots (count t0) window
    (incf count)
    (let ((time (get-internal-real-time)))
      (when (= t0 0)
        (setq t0 time))
      (when (>= (- time t0) (* 5 internal-time-units-per-second))
        (let* ((seconds (/ (- time t0) internal-time-units-per-second))
               (fps (/ count seconds)))
          (format *terminal-io* "~D frames in ~3,1F seconds = ~6,3F FPS~%"
                  count seconds fps))
        (with-simple-restart (continue "continue")
          (format t "~s samples~%" (gl:get* :samples)))
        (setq t0 time)
        (setq count 0)))))


(defmethod glut:display ((w lines))
  (gl:enable :multisample :line-smooth :sample-alpha-to-coverage)
  #++(gl:disable ;:sample-alpha-to-coverage
   :line-smooth
   )
  (glut:enable-event w :idle)
  (gl:clear-color 0.2 0.0 0.0 0.0)
  (gl:clear :color-buffer :depth-buffer)
  (gl:load-identity)
  (gl:translate 0 0 -2)
  (gl:with-pushed-matrix
      (gl:rotate view-rotx 1 0 0)
    (gl:rotate view-roty 0 1 0)
    (gl:rotate view-rotz 0 0 1)
    ;;(gl:rotate (mod (/ (get-internal-real-time) 100.0) 360) 0 1 0)
    (gl:scale (* view-scale 0.02) (* view-scale 0.02) (* view-scale 0.02))
    (gl:translate -50 -50 -50)

    (when (buffers w)
      (gl:line-width 2)
      (%gl:bind-vertex-array-apple (first (buffers w)))
      (loop for (color start count) in (draw-calls w)
         do (when *once*
              (format t "~s ~s ~s~%" color start count))
                                        ;  (setf color (append color (list 0.6)))
           (apply #'gl:color color)
           (%gl:draw-elements :lines count :unsigned-int
                              (* start +index-bytes+)
                              ))
      (%gl:bind-vertex-array-apple 0)))
  (setf *once* nil)
  (glut:swap-buffers)
  (print-frame-rate w))

(defparameter *once* t)


#+(or)(defmethod glut:idle ((w lines))
  (glut:post-redisplay))

(defmethod glut:reshape ((w lines) width height)
  (gl:viewport 0 0 width height)
  (gl:matrix-mode :projection)
  (gl:load-identity)
  (glu:perspective 50 (/ width height) 0.5 200)
  (gl:matrix-mode :modelview)
  (gl:load-identity))

(defmethod glut:keyboard ((w lines) key x y)
  (declare (ignore x y))
  (case key
    (#\q
     (setf (a w) (float (mod (+ (a w) (* 2 pi ) -0.02) (* 2 pi)) 1.0)))
    (#\e
     (setf (a w) (float (mod (+ (a w) 0.02) (* 2 pi)) 1.0)))
    (#\space
     (translate-data w))
    (#\Esc
     (glut:destroy-current-window))))

(defun lines ()
;;  (glut:set-option :multisample 4)
  (glut:display-window (make-instance 'lines)))




;;; ------------------------------------------------------------

(defvar origclick)
(defvar origrot)

(defmethod glut:mouse ((window glut:window) button state x y)
  (cond
    ((eq button :left-button)
     (if (eq state :down)
	 (progn (setf motion :rotate)
		(setf origrot (list view-rotx view-roty))
		(setf origclick (list x y)))
	 (setf origclick ()
	       motion nil)))
    ((eq button :right-button)
     (if (eq state :down)
	 (setf motion :scale
	       origscale view-scale
	       origclick (list x y))
	 (setf origclick ()
	       motion nil)))))

(defmethod glut:motion ((window glut:window) x y)
  (case motion
    (:rotate
     (setf view-rotx (+ (car origrot) (- y (cadr origclick))))
     (setf view-roty (+ (cadr origrot) (- x (car origclick))))
     (glut:post-redisplay))
    (:scale
     (let* ((delta (- y (second origclick)))
            (scale-factor (expt 2.0 (* delta 0.01))))
       (setf view-scale (* origscale scale-factor))
       (lim-between view-scale 0.01 10.0)))))


(defun random-interval (bot top)
  (+ (* (- top bot) (/ (random 100000) 100000.0)) bot))

(defvar view-rotx-vel 0)
(defvar view-roty-vel 0)
(defvar view-rotz-vel 0)

(defvar last-update 0)
(defvar counter 0)

(defmethod glut:idle ((window glut:window))
  (swank::process-requests t)
  (if walk-mode
    (progn
      (incf counter)

      (if (< (+ last-update internal-time-units-per-second) (get-internal-real-time))
        (progn
          (format t "~a frames per second with ~a slices.~%" counter slices)
          (setq counter 0)
          (setq last-update (get-internal-real-time))))

      (incf view-rotx-vel (random-interval -.1 .1))
      (incf view-roty-vel (random-interval -.1 .1))
      (incf view-rotz-vel (random-interval -.1 .1))
      (lim-between view-rotx-vel -2 2)
      (lim-between view-roty-vel -2 2)
      (lim-between view-rotz-vel -2 2)
      (incf view-rotx view-rotx-vel)
      (incf view-roty view-roty-vel)
      (incf view-rotz view-rotz-vel)

      (incf light-r (random-interval -.02 .02))
      (incf light-g (random-interval -.02 .02))
      (incf light-b (random-interval -.02 .02))
      (lim-between light-r 0 1)
      (lim-between light-g 0 1)
      (lim-between light-b 0 1)))
  (glut:post-redisplay))





#++
(lines)
