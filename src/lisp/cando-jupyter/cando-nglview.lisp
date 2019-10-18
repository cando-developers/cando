
(in-package :cando-jupyter)

(defun pick-history (widget)
  (funcall (find-symbol "PICK-HISTORY" :NGLV) widget))

(defun safe-add-shape (widget parts &key name)
  (funcall (find-symbol "ADD-SHAPE" :NGLV) widget parts :name name))


(defconstant +width+ 1.0)
(defconstant +half-width+ (/ +width+ 2.0))
(defun build-hydrogens (matter)
  (chem:map-atoms nil
                  (lambda (a)
                    (when (eq (chem:get-element a) :H)
                      (let* ((neighbors (chem:bonds-as-list a))
                             (_         (unless (eql (length neighbors) 1)
                                          (error "There are too many neighbors: ~a" neighbors)))
                             (neighbor (chem:get-other-atom (first neighbors) a))
                             (neighbor-pos (chem:get-position neighbor))
                             (pos (geom:vec
                                   (+ (- (random +width+) +half-width+) (geom:vx neighbor-pos))
                                   (+ (- (random +width+) +half-width+) (geom:vy neighbor-pos))
                                   (+ (- (random +width+) +half-width+) (geom:vz neighbor-pos)))))
                        (chem:set-position a pos))))
                  matter))

(defun npicked (widget n)
    (with-output-to-string (sout)
  (mapcar (lambda (a) 
                  (format sout ".~a "(cdr (assoc "atomname" (cdr (car a)) :test #'equal)))) 
          (subseq (pick-history widget) 0 n))))


(defun pick-property (pick property)
  (cdr (assoc property (cdr (car pick)) :test #'equal)))

(defun atom-with-name (matter name)
  (let (atom)
    (chem:map-atoms 'nil (lambda (a)
                           (when (eq (chem:get-name a) name)
                             (setf atom a)))
                    matter)
    atom))

(defun atom-map (widget n matter-from matter-to)
  (unless (>= (length (pick-history widget)) (* 2 n))
    (error "There aren't enough atoms selected to create a map for ~a pairs." n))
  (let ((map (loop for (pick-to pick-from) on (pick-history widget) by #'cddr
                   for atom-from = (atom-with-name matter-from (intern (pick-property pick-from "atomname") :keyword))
                   for atom-to = (atom-with-name matter-to (intern (pick-property pick-to "atomname") :keyword))
                   collect (cons atom-from atom-to))))
    (subseq map 0 n)))

(defun named-atom-map (widget n matter-from matter-to)
  (let ((atom-map (atom-map widget n matter-from matter-to)))
    (loop for (from . to) in atom-map
          collect (cons (chem:get-name from) (chem:get-name to)))))


(defun anchor-named-atom-map (name-atom-map matter-from matter-to)
  (loop for (from-name . to-name) in name-atom-map
        for from-atom = (chem:first-atom-with-name matter-from from-name)
        for to-atom = (chem:first-atom-with-name matter-to to-name)
        for to-position = (chem:get-position to-atom)
        do (cando:anchor-atom from-atom to-position)))


(defun coord-to-vector (c)
  (vector (float (geom:vx c) 1.0s0)
          (float (geom:vy c) 1.0s0)
          (float (geom:vz c) 1.0s0)))

(defun arrow (from-vec to-vec &optional (radius 0.2))
  (let ((from-coord (coord-to-vector from-vec))
        (to-coord (coord-to-vector to-vec)))
    (vector "arrow" from-coord to-coord #(1 0 1) radius)))
(defun cartoon-atom-map-arrows (atom-map from-matter to-matter)
  (let (arrows)
    (loop for pair in atom-map
          for from-name = (car pair)
          for to-name = (cdr pair)
          for from-atom = (atom-with-name from-matter from-name)
          for to-atom = (atom-with-name to-matter to-name)
          for from-pos = (chem:get-position from-atom)
          for from-coord = (coord-to-vector from-pos)
          for to-pos = (chem:get-position to-atom)
          for to-coord = (coord-to-vector to-pos)
          for arrow = (vector "arrow" from-coord to-coord #(1 0 1) 0.2)
;;;          do (format t "Names:  ~a -> ~a   atoms:  ~a -> ~a~%" from-name to-name from-atom to-atom)
;;;          do (format t "   coords: ~a -> ~a~%" from-coord to-coord)
;;;          do (format t "    arrow: ~a~%" arrow)
             collect arrow
          )))

(defun cartoon-atom-map (widget atom-map from-matter to-matter)
  (let ((shape (cartoon-atom-map-arrows atom-map from-matter to-matter)))
    (safe-add-shape widget (coerce shape 'vector) :name "arrows")
    nil))



#++(defun carboxylic-acid-atoms (matter unique-name)
  (let ((c (chem:first-atom-with-name matter unique-name))
	(carb (core:make-cxx-object 'chem:chem-info)))
    (chem:compile-smarts carb "C1(~O2)~O3")
    (or (chem:matches carb c)
	(error "The atom ~a is not a carboxylic acid carbon" c))
    (let* ((m (chem:get-match carb))
	   (o2 (chem:get-atom-with-tag m :2))
	   (o3 (chem:get-atom-with-tag m :3)))
      (list c o2 o3))))


(defun center-on (matter &optional (dest (geom:vec 0.0 0.0 0.0)))
  (let ((dest (cond
                ((typep dest 'geom:v3) dest)
                ((consp dest) (geom:vec (first dest) (second dest) (third dest)))
                (t (error "Convert ~a to vec" dest)))))
    (let ((transform (geom:make-m4-translate 
                      (geom:v- dest (chem:geometric-center matter) ))))
      (chem:apply-transform-to-atoms matter transform))))

(defun rotate-x (matter angle-degrees)
  (let ((transform (geom:make-m4-rotate-x (* 0.0174533 angle-degrees))))
    (chem:apply-transform-to-atoms matter transform)))

(defun rotate-y (matter angle-degrees)
  (let ((transform (geom:make-m4-rotate-y (* 0.0174533 angle-degrees))))
    (chem:apply-transform-to-atoms matter transform)))

(defun rotate-z (matter angle-degrees)
  (let ((transform (geom:make-m4-rotate-z (* 0.0174533 angle-degrees))))
    (chem:apply-transform-to-atoms matter transform)))

(defun distance-two-positions (p1 p2)
  (geom:vlength (geom:v- p1 p2)))

(defun picked-atoms (history agg &optional (num 1))
  "Return a list of the num most recently picked atoms"
    (loop for picked in (reverse (subseq history 0 num))
           collect (let* ((atom (cdr (assoc "atom" picked :test #'string=)))
                           (x (cdr (assoc "x" atom :test #'string=)))
                           (y (cdr (assoc "y" atom :test #'string=)))
                           (z (cdr (assoc "z" atom :test #'string=)))
                           (pos (geom:vec x y z))
                            resulta)
                         (chem:map-atoms nil (lambda (a) 
                                 (when (< (distance-two-positions pos (chem:get-position a)) 0.01) 
                                     (setf resulta a))) agg)
                        resulta)))

(defclass cando-structure (nglv:structure)
  ((%matter :initarg :matter :accessor matter)))

(defmethod nglv:ext ((self cando-structure))
  "mol2")

(defmethod nglv:get-structure-string ((self cando-structure))
  (check-type self cando-structure)
  (progn
    (cl-jupyter:logg 2 "Generating mol2 as string~%")
    (chem:aggregate-as-mol2-string (matter self) t))
  #++(progn
       (cl-jupyter:logg 2 "Saving structure to /tmp/structure.mol2~%")
       (cando:save-mol2 (matter self) "/tmp/structure.mol2" :use-sybyl-types t)
       (with-open-file (stream "/tmp/structure.mol2" :direction :input)
	 (let* ((entire-file (make-string (+ (file-length stream) 2)
					  :initial-element #\newline)))
	   (read-sequence entire-file stream)
	   (close stream)
	   entire-file))))

(defclass cando-trajectory (nglv:trajectory nglv:structure dynamics:trajectory)
  ())

(defmethod initialize-instance :after ((self cando-trajectory) &key)
  (setf (gethash "cando" nglv:*BACKENDS*) 'cando-trajectory)
  (values))

(defmethod nglv:ext ((self cando-trajectory))
  "mol2")

(defmethod nglv:get-structure-name ((self cando-trajectory))
  (chem:get-name (dynamics:matter self)))

(defmethod nglv:get-structure-string ((self cando-trajectory))
  (chem:aggregate-as-mol2-string (dynamics:matter self) t))


(defmethod nglv:append-coordinates ((self cando-trajectory) coordinates)
  (vector-push-extend coordinates (dynamics:coordinates self)))

(defmethod nglv:get-coordinates ((self cando-trajectory) index)
  (aref (dynamics:coordinates self) index))

(defmethod nglv:n-frames ((self cando-trajectory))
  (length (dynamics:coordinates self)))

(defmethod show ((trajectory cando-trajectory) &rest kwargs &key &allow-other-keys)
  (apply #'nglv:make-nglwidget :structure trajectory kwargs))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; amber-netcdf-trajectory
;;;

(defclass amber-netcdf-trajectory (nglv:trajectory nglv:structure)
  ((netcdf :initarg :netcdf :accessor netcdf)
   (number-of-atoms :initarg :number-of-atoms :accessor number-of-atoms)
   (matter :initarg :matter :accessor matter)
   (number-of-frames :initarg :number-of-frames :accessor number-of-frames)
   (current-coordinates :initform nil :initarg :current-coordinates :accessor current-coordinates)))

(defun make-amber-netcdf-trajectory (&key netcdf matter)
  (let* ((number-of-atoms (chem:number-of-atoms matter))
         (coords (static-vectors:make-static-vector (+ 10 (* number-of-atoms 3)) :element-type 'single-float))
         (number-of-frames (netcdf:get-dimension netcdf "frame")))
    (make-instance 'amber-netcdf-trajectory
                   :netcdf netcdf
                   :number-of-atoms number-of-atoms
                   :matter matter
                   :number-of-frames number-of-frames
                   :current-coordinates coords)))

(defmethod initialize-instance :after ((self cando-trajectory) &key)
  (setf (gethash "amber-netcdf" nglv:*BACKENDS*) 'amber-netcdf-trajectory)
  (values))

(defmethod nglv:ext ((self amber-netcdf-trajectory))
  "mol2")

(defmethod nglv:get-structure-name ((self amber-netcdf-trajectory))
  (chem:get-name (matter self)))

(defmethod nglv:get-structure-string ((self amber-netcdf-trajectory))
  (chem:aggregate-as-mol2-string (matter self) t))

(defmethod nglv:get-coordinates ((self amber-netcdf-trajectory) index)
  (netcdf:get-vara-float (netcdf self) "coordinates"
                                 (vector index 0 0)
                                 (vector 1 (number-of-atoms self) 3)
                                 (current-coordinates self))
  (current-coordinates self))

(defmethod nglv:n-frames ((self amber-netcdf-trajectory))
  (number-of-frames self))

(defmethod show ((trajectory amber-netcdf-trajectory) &rest kwargs &key &allow-other-keys)
  (apply #'nglv:make-nglwidget :structure trajectory kwargs))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Change the way that code cells are evaluated to allow
;;; leap syntax to be used.
;;;
;;;

(defun parse-lisp-or-leap (code)
  (let ((code (string-trim (list #\space #\tab) (copy-seq code))))
    (if (or (char= (char code 0) #\() (char= (char code 0) #\*))
        (let ((sexp (read-from-string (format nil "(progn ~A~%)" code))))
          `(core:call-with-stack-top-hint
            (lambda ()
              ,sexp)))
        (let ((ast (architecture.builder-protocol:with-builder ('list)
                     (handler-bind ((esrap:esrap-parse-error
                                      (lambda (c)
                                        (format t "cando-nglview: Encountered error ~s while parsing ~s~%" c code))))
                       (esrap:parse 'leap.parser::leap code)))))
          `(core:call-with-stack-top-hint
            (lambda ()
              (leap.core:evaluate 'list ',ast leap.core:*leap-env*)))))))

(eval-when (:load-toplevel :execute)
  (when (find-symbol "*READ-CODE-HOOK*" :cl-jupyter)
    (setf cl-jupyter:*read-code-hook* 'parse-lisp-or-leap)))

(defun leap-syntax-enable (on)
  (let ((sym (find-symbol "*READ-CODE-HOOK*" :cl-jupyter)))
    (when sym
      (if on
          (set sym 'parse-lisp-or-leap)
          (set sym nil)))))
