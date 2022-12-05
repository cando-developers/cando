(in-package :modeller)

(defparameter unitx (geom:vec 1.0 0.0 0.0))
(defparameter unity (geom:vec 0.0 1.0 0.0))
(defparameter unitz (geom:vec 0.0 0.0 1.0))



(defclass scene ()
  ((shapes :initform nil :initarg shapes :accessor shapes)))

(defun add-shape (scene shape)
  (push shape (shapes scene)))

(defclass shape ()
  ((vertices :initform (make-array 9 :element-type 'single-float :adjustable t :fill-pointer 0)
             :initarg :vertices :accessor vertices)
   (faces :initform (make-array 9 :element-type 'ext:byte64 :adjustable t :fill-pointer 0)
          :initarg :faces :accessor faces)))

(defun make-cylinder (center axis radius length &optional (divisions 36))
  (let* ((axis (geom:vnormalized axis))
         (half-length (/ length 2.0))
         (top (geom:v* axis half-length))
         (bottom (geom:v* axis (- half-length)))
         (perp (geom:vcross axis unitx))
         (perplen (geom:vlength perp)))
    (when (< perplen 0.02)
      (setf perp (geom:vcross axis unity)))
    (setf perp (geom:vnormalized perp))
    (let* ((cross (geom:vcross axis perp))
           (vertices (make-array (+ 2 (* 2 divisions)) :element-type 'geom:vector3))
           (faces (make-array (* 3 4 divisions) :element-type 'ext:byte32 :fill-pointer 0))
           (top-center-idx (* 2 divisions))
           (bottom-center-idx (1+ top-center-idx)))
      (flet ((add-face (faces i0 i1 i2)
               (vector-push i0 faces)
               (vector-push i1 faces)
               (vector-push i2 faces)))
        (loop for angle from 0.0 below 360.0 by (/ 360.0 divisions)
              for index from 0 below (* 2 divisions) by 2
              for angle-rad = (* 0.0174533 angle)
              for pos = (geom:v+ (geom:v* perp (* (cos angle-rad) radius))
                                 (geom:v* cross (* (sin angle-rad) radius)))
              for top-pos = (geom:v+ (geom:v+ pos top) center)
              for bottom-pos = (geom:v+ (geom:v+ pos bottom) center)
              do (setf (elt vertices index) top-pos)
              do (setf (elt vertices (1+ index)) bottom-pos))
        (setf (elt vertices top-center-idx) (geom:v+ top center)
              (elt vertices bottom-center-idx) (geom:v+ bottom center))
        (loop with divisions*2 = (* 2 divisions)
              for idx below divisions*2 by 2
              for top-cur-idx = idx
              for bottom-cur-idx = (1+ idx)
              for top-next-idx = (multiple-value-bind (div mod)
                                     (floor (+ 2 top-cur-idx) divisions*2)
                                   (declare (ignore div))
                                   mod)
              for bottom-next-idx = (1+ top-next-idx)
              do (add-face faces top-center-idx top-cur-idx top-next-idx)
                 (add-face faces top-cur-idx bottom-cur-idx top-next-idx)
                 (add-face faces bottom-cur-idx bottom-next-idx top-next-idx)
                 (add-face faces bottom-center-idx bottom-next-idx bottom-cur-idx)))
      (make-instance 'shape
                     :vertices vertices
                     :faces faces))))


(defun make-icosphere (center radius &optional (order 3))
  (let* ((f (/ (+ 1.0 (expt 5.0 0.5)) 2.0))
         (-f (- f))
         (tt (expt 4.0 order))
         (vertices-length (floor (* (+ (* 10 tt) 2) 3)))
         (vertices (make-array vertices-length :element-type 'single-float :initial-element 0.0s0))
         (triangle-indices (list
                            0 11 5 0 5 1 0 1 7 0 7 10 0 10 11
                            11 10 2 5 11 4 1 5 9 7 1 8 10 7 6
                            3 9 4 3 4 2 3 2 6 3 6 8 3 8 9
                            9 8 1 4 9 5 2 4 11 6 2 10 8 6 7))
         (triangles (make-array (length triangle-indices) :element-type 'ext:byte32
                                                          :initial-contents triangle-indices))
         (triangles-prev nil)
         (vv 12)
         (mid-cache (if (= order 0)
                        nil
                        (make-hash-table :test 'eql)))
         )
    (loop for val in (list
                      -1.0s0 f 0.0s0 1.0s0 f 0.0s0 -1.0s0 -f 0.0s0 1.0s0 -f 0.0s0 
                      0.0s0 -1.0s0 f 0.0s0 1.0s0 f 0.0s0 -1.0s0 -f 0.0s0 1.0s0 -f
                      f 0.0s0 -1.0s0 f 0.0s0 1.0s0 -f 0.0s0 -1.0s0 -f 0.0s0 1.0s0 )
          for index from 0
          do (setf (elt vertices index) val))
    (flet ((add-mid-point (aa bb)
             (let* ((key (+ (floor (* (+ aa bb) (+ aa bb 1) 0.5)) (min aa bb)))
                    (ikey (gethash key mid-cache)))
               (when ikey
                 (remhash key mid-cache)
                 (return-from add-mid-point ikey))
               (setf (gethash key mid-cache) vv)
               (loop for kk below 3
                     for aa-val = (elt vertices (+ (* 3 aa) kk))
                     for bb-val = (elt vertices (+ (* 3 bb) kk))
                     for new-val = (+ aa-val bb-val)
                     for index = (+ (* 3 vv) kk)
                     do (if (>= index (length vertices))
                            (warn "index ~a is >= (length vertices) ~a" index (length vertices))
                            (setf (elt vertices index) new-val)))
               (setf ikey vv)
               (incf vv)
               ikey))
           (hypot (&rest vals)
             (let ((sum-squares (loop for val in vals
                                      sum (* val val))))
               (sqrt sum-squares))))
      (setf triangles-prev (copy-seq triangles))
      (loop for ii below order
            do (setf triangles (make-array (* 4 (length triangles-prev))
                                           :element-type 'ext:byte32))
            do (loop for kk below (length triangles-prev) by 3
                     for v1 = (elt triangles-prev kk)
                     for v2 = (elt triangles-prev (+ kk 1))
                     for v3 = (elt triangles-prev (+ kk 2))
                     for aa = (add-mid-point v1 v2)
                     for bb = (add-mid-point v2 v3)
                     for cc = (add-mid-point v3 v1)
                     for tt = (* 4 kk)
                     do (setf (elt triangles (+ tt 0)) v1
                              (elt triangles (+ tt 1)) aa
                              (elt triangles (+ tt 2)) cc
                              (elt triangles (+ tt 3)) v2
                              (elt triangles (+ tt 4)) bb
                              (elt triangles (+ tt 5)) aa
                              (elt triangles (+ tt 6)) v3
                              (elt triangles (+ tt 7)) cc
                              (elt triangles (+ tt 8)) bb
                              (elt triangles (+ tt 9)) aa
                              (elt triangles (+ tt 10)) bb
                              (elt triangles (+ tt 11)) cc))
            do (setf triangles-prev (copy-seq triangles)))
      (let ((vertices-vecs (make-array (/ (length vertices) 3))))
        (loop for ii below (length vertices) by 3
              for vi from 0
              for mm = (/ 1.0 (hypot (elt vertices ii) (elt vertices (+ ii 1)) (elt vertices (+ ii 2))))
              for xx = (* (elt vertices ii) mm radius)
              for yy = (* (elt vertices (+ 1 ii)) mm radius)
              for zz = (* (elt vertices (+ 2 ii)) mm radius)
              do (setf (elt vertices-vecs vi) (geom:v+ center (geom:vec xx yy zz))))
        (make-instance 'shape
                       :vertices vertices-vecs
                       :faces triangles)))))

(defun make-cuboid (center xsize ysize zsize)
  (let* ((half-xsize (/ xsize 2.0))
         (half-ysize (/ ysize 2.0))
         (half-zsize (/ zsize 2.0))
         (xn (+ (- half-xsize) (geom:vx center)))
         (xx (+ half-xsize (geom:vx center)))
         (yn (+ (- half-ysize) (geom:vy center)))
         (yx (+ half-ysize (geom:vy center)))
         (zn (+ (- half-zsize) (geom:vz center)))
         (zx (+ half-zsize (geom:vz center))))
    (let ((vertices (list (geom:vec xn yn zn) ;  0
                          (geom:vec xn yx zn) ;  1
                          (geom:vec xx yn zn) ;  2
                          (geom:vec xx yx zn) ;  3
                          (geom:vec xn yn zx) ;  4
                          (geom:vec xn yx zx) ;  5
                          (geom:vec xx yn zx) ;  6
                          (geom:vec xx yx zx) ;  7
                          ))
          (faces (list 0 1 2
                       1 3 2
                       0 2 4
                       2 6 4
                       2 7 6
                       2 3 7
                       0 4 1
                       1 4 5
                       1 5 3
                       3 5 7
                       4 6 5
                       5 6 7
                       )))
      (make-instance 'shape
                     :vertices (make-array (length vertices) :element-type 'geom:v3 :initial-contents vertices)
                     :faces (make-array (length faces) :element-type 'ext:byte32 :initial-contents faces)))))


(defun generate-obj (scene stream)
  (let ((shape-vertex-start (make-hash-table))
        (abs-vertex-index 1))
    (loop for shape in (shapes scene)
          do (setf (gethash shape shape-vertex-start) abs-vertex-index)
          do (loop for vertex-index below (length (vertices shape))
                   for vertex = (elt (vertices shape) vertex-index)
                   do (format stream "v ~8,3f ~8,3f ~8,3f~%" (geom:vx vertex) (geom:vy vertex) (geom:vz vertex))
                   do (incf abs-vertex-index)))
    (loop for shape in (shapes scene)
          for vertex-start-index = (gethash shape shape-vertex-start)
          do (loop for face-index from 0 below (length (faces shape)) by 3
                   for face0 = (+ (elt (faces shape) face-index) vertex-start-index)
                   for face1 = (+ (elt (faces shape) (+ 1 face-index)) vertex-start-index)
                   for face2 = (+ (elt (faces shape) (+ 2 face-index)) vertex-start-index)
                   do (format stream "f ~a ~a ~a~%" face0 face1 face2)))))


(defun save-obj (scene filename)
  (with-open-file (fout filename :direction :output)
    (generate-obj scene fout)))



(defstruct facet
  ni nj nk
  v1x v1y v1z
  v2x v2y v2z
  v3x v3y v3z
  )

(defun safe-normal (vec1 vec2 vec3)
  (let ((vec21 (geom:v- vec2 vec1))
        (vec31 (geom:v- vec3 vec1)))
    (let ((vcross (geom:vcross vec21 vec31)))
      (when (> (geom:vlength vcross) 0.0000001)
        (return-from safe-normal (geom:vnormalized vcross)))))
  (let ((vec32 (geom:v- vec3 vec2))
        (vec12 (geom:v- vec1 vec2)))
    (let ((vcross (geom:vcross vec32 vec12)))
      (when (> (geom:vlength vcross) 0.0000001)
        (return-from safe-normal (geom:vnormalized vcross)))))
  (error "triangle is too acute to calculate normal ~a ~a ~a" vec1 vec2 vec3))


(defun generate-stl (scene stream)
  (let ((number-of-triangles (loop for shape in (shapes scene)
                                   sum (/ (length (faces shape)) 3)))
        (header (make-array 80 :element-type 'ext:byte8 :initial-element (char-int #\space))))
    (flet ((write-header (header start msg)
             (loop for idx below (length msg)
                   do (setf (elt header (+ idx start)) (char-int (elt msg idx))))))
      (write-header header 0 "UNITS=mm"))
    (write-sequence header stream)
    (sys:write-c-uint32 stream number-of-triangles)
    (loop for shape in (shapes scene)
          do (loop with vertices = (vertices shape)
                   with faces = (faces shape)
                   for index from 0 below (length faces) by 3
                   for f1 = (elt faces index)
                   for f2 = (elt faces (+ 1 index))
                   for f3 = (elt faces (+ 2 index))
                   for v1 = (elt vertices f1)
                   for v2 = (elt vertices f2)
                   for v3 = (elt vertices f3)
                   for vecnormal = (safe-normal v1 v2 v3)
                   do (flet ((safe-write-single-float (stream val)
                               (sys:write-ieee-single-float stream (float val 1.0s0))))
                        (safe-write-single-float stream (geom:vx vecnormal))
                        (safe-write-single-float stream (geom:vy vecnormal))
                        (safe-write-single-float stream (geom:vz vecnormal))
                        (safe-write-single-float stream (geom:vx v1))
                        (safe-write-single-float stream (geom:vy v1))
                        (safe-write-single-float stream (geom:vz v1))
                        (safe-write-single-float stream (geom:vx v2))
                        (safe-write-single-float stream (geom:vy v2))
                        (safe-write-single-float stream (geom:vz v2))
                        (safe-write-single-float stream (geom:vx v3))
                        (safe-write-single-float stream (geom:vy v3))
                        (safe-write-single-float stream (geom:vz v3))
                        (sys:write-c-uint16 stream 0))))))

(defun save-stl (scene filename)
  (with-open-file (fout filename :direction :output :element-type '(unsigned-byte 8))
    (generate-stl scene fout)))


(defun plate-384 (&key (post-radius 1.0) (post-height 2.5)
                    (xbegin 0) (xend 24)
                    (ybegin 0) (yend 16)
                    (plate-height 2.21))
  (let* ((scene (make-instance 'scene))
         (plate-width-x 150.0)
         (plate-width-y 100.0)
         (plate (make-cuboid (geom:vec 0 0 (/ plate-height 2.0)) plate-width-x plate-width-y plate-height))
         (upper-outer-x-width 107.24) ; distance from upper first last hole x outside
         (upper-inner-x-width 99.89)
         (upper-outer-y-width 71.26) ; distance from upper first last hole y outside
         (upper-inner-y-width 63.84)
         (post-x-start-center-offset (+ 0.2 (/ (+ upper-outer-x-width upper-inner-x-width) 4.0)))
         (post-y-start-center-offset (+ 0.2 (/ (+ upper-outer-y-width upper-inner-y-width) 4.0)))
         (post-start-x (- post-x-start-center-offset))
         (post-stop-x post-x-start-center-offset)
         (post-start-y (- post-y-start-center-offset))
         (post-stop-y post-y-start-center-offset)
         (xnum 24)
         (ynum 16)
         (post-inc-x (/ (- post-stop-x post-start-x) (1- xnum)))
         (post-inc-y (/ (- post-stop-y post-start-y) (1- ynum))))
    (loop for xi from xbegin below xend
          for xpos = (+ post-start-x (* post-inc-x xi))
          do (loop for yi from ybegin below yend
                   for ypos = (+ post-start-y (* post-inc-y yi))
                   for zpos = (+ plate-height (/ post-height 2.0))
                   for post = (make-cylinder (geom:vec xpos ypos zpos)
                                             unitz
                                             post-radius
                                             post-height)
                   do (push post (shapes scene))))
    (push plate (shapes scene))
    scene))

(defun bottom-plate-384 () 
  (plate-384 :plate-height 3.0
             :post-radius 1.1
             :post-height 2.5))

(defun top-plate-384 ()
  (plate-384 :plate-height 3.0
             :post-radius 1.2
             :post-height 5.0))

(defun bottom-plate-2-col ()
    (plate-384 :plate-height 3.0
               :post-radius 1.1
               :post-height 2.5
               :xbegin 11
               :xend 13))

(defun top-plate-2-col ()
    (plate-384 :plate-height 3.0
               :post-radius 1.2
               :post-height 5.0
               :xbegin 11
               :xend 13))

(defun spacer ()
  (let* ((width-x 135.0)
         (width-y 110.0)
         (width-z 15.0)
         (posts-x 13)
         (posts-y 11)
         (thickness 3.0)
         (post-radius 2.0)
         (top-z-center (- width-z (/ thickness 2.0)))
         (top-z-bottom (- width-z thickness))
         (bottom-z-top thickness)
         (limit-x (- (/ width-x 2.0) post-radius))
         (limit-y (- (/ width-y 2.0) post-radius))
         (inc-x (/ (* 2.0 limit-x) (1- posts-x)))
         (inc-y (/ (* 2.0 limit-y) (1- posts-y)))
         (bottom-z-center (/ thickness 2.0))
         (scene (make-instance 'scene))
         (top-plate (make-cuboid (geom:vec 0 0 top-z-center) width-x width-y thickness))
         (bottom-plate (make-cuboid (geom:vec 0 0 bottom-z-center) width-x width-y thickness)))
    (loop for xidx from 0 below posts-x
          for xpos from (- limit-x) by inc-x
          do (loop for yidx from 0 below posts-y
                   for ypos from (- limit-y) by inc-y
                   for post = (make-cylinder (geom:vec xpos ypos (/ width-z 2.0)) unitz post-radius (- top-z-bottom bottom-z-top))
                   do (push post (shapes scene))))
    (push top-plate (shapes scene))
    (push bottom-plate (shapes scene))
    scene))

(defun slide ()
  (let* ((scene (make-instance 'scene))
         (plate-height 1.1)
         (plate-width-x 75.7)
         (plate-width-y 25.8)
         (plate (make-cuboid (geom:vec 0 0 (/ plate-height 2.0)) plate-width-x plate-width-y plate-height))
         )
    (push plate (shapes scene))
    scene))

(defun shim-96-well ()
  (let* ((scene (make-instance 'scene))
         (height 10)
         (outer-x-width 127.0 ); 127.8) ; 127.0) ; 126.7) ; 126.4)
         (outer-y-width 85.3) ; 85.5) ; 85.2) ; 84.9)
         (inner-x-width 122.9) ; 122.7) ; 122.4) ; 122.1)
         (inner-y-width 81.3) ; 80.9) ; 80.7) ; 80.4) ; 80.1)
         (x-thickness (/ (- outer-x-width inner-x-width) 2.0))
         (y-thickness (/ (- outer-y-width inner-y-width) 2.0))
         (x-center (+ (/ inner-x-width 2.0) (/ x-thickness 2.0)))
         (y-center (+ (/ inner-y-width 2.0) (/ y-thickness 2.0))))
    (push (make-cuboid (geom:vec x-center 0 0) x-thickness outer-y-width height) (shapes scene))
    (push (make-cuboid (geom:vec (- x-center) 0 0) x-thickness outer-y-width height) (shapes scene))
    (push (make-cuboid (geom:vec 0 y-center 0) outer-x-width y-thickness height) (shapes scene))
    (push (make-cuboid (geom:vec 0 (- y-center) 0) outer-x-width y-thickness height) (shapes scene))
    scene))

(defun generate-all (&optional (path "/tmp/shapes/"))
  (ensure-directories-exist path)
  (save-stl (shim-96-well) (merge-pathnames #P"shim.stl" path))
  (save-stl (slide) (merge-pathnames #P"slide.stl" path))
  (save-stl (spacer) (merge-pathnames #P"spacer.stl" path))
  (save-stl (bottom-plate-384) (merge-pathnames #P"bottom-plate-384.stl" path))
  (save-stl (top-plate-384) (merge-pathnames #P"top-plate-384.stl" path))
  (save-stl (bottom-plate-2-col) (merge-pathnames #P"bottom-plate-2-col.stl" path))
  (save-stl (top-plate-2-col) (merge-pathnames #P"top-plate-2-col.stl" path)))
