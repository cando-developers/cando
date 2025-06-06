
#+tests(in-package #:clasp-tests)

(defparameter *zjs5-from-coords* (multiple-value-call 'chem:cpring-coordinates-to-zjs (chem:cpring-5-sample)))
(defparameter *cp5-from-coords* (multiple-value-list (multiple-value-call 'chem:cpring-5-coordinates-to-cremer-pople (chem:cpring-5-sample))))
(defparameter *cp5-from-coords-zjs* (multiple-value-list (chem:cpring-5-zjs-to-cremer-pople *zjs5-from-coords*)))
(defparameter *zjs5-from-cp5* (apply 'chem:cpring-5-cremer-pople-to-zjs *cp5-from-coords*))
(defparameter *cp5-from-cpzjs* (multiple-value-list (chem:cpring-5-zjs-to-cremer-pople *zjs5-from-cp5*)))

(defvar *5-coords*)
(defvar *5-indexes*)
(multiple-value-setq (*5-coords* *5-indexes*) (chem:cpring-5-sample))
(defparameter *new-5-coords* (chem:nvector-copy *5-coords*))

(defun make-nvector* (len)
 (make-array len :element-type 'double-float :adjustable t))

(defvar *5-dists* (make-nvector* (length *5-indexes*)))
(defvar *5-angles* (make-nvector* (length *5-indexes*)))
(chem:cpring-coordinates-to-distances-and-angles *5-dists* *5-angles* *5-coords* *5-indexes*)
(chem:cpring-generate-coordinates *zjs5-from-cp5* *5-dists* *5-angles* *new-5-coords* *5-indexes*)
(defvar *new-5-dists* (make-nvector* (length *5-indexes*)))
(defvar *new-5-angles* (make-nvector* (length *5-indexes*)))
(chem:cpring-coordinates-to-distances-and-angles *new-5-dists* *new-5-angles* *new-5-coords* *5-indexes*)


(progn
  (format t "cpring-5-sample -> ~a~%" (multiple-value-call 'chem:cpring-mathematica (chem:cpring-5-sample)))
  (format t "cp5-from-coords*     -> ~10,7f~%" *cp5-from-coords*)
  (format t "cp5-from-coords-zjs* -> ~10,7f~%" *cp5-from-coords-zjs*)
  (format t "cp5-from-cpzjs*      -> ~10,7f~%" *cp5-from-cpzjs*)
  (format t "zjs5-from-cp5*    = ( ~{ ~10,7f~^,~} )~%" (coerce *zjs5-from-cp5* 'list))
  (format t "zjs5-from-coords* = ( ~{ ~10,7f~^,~} )~%" (coerce *zjs5-from-coords* 'list))
  (format t "*new-5-coords* = ~a~%" (chem:cpring-mathematica *new-5-coords* *5-indexes*))
  (format t "*5-dists* = ~s~%" *5-dists*)
  (format t "*new-5-dists* = ~s~%" *new-5-dists*)
  (format t "*5-angles* = ~s~%" (mapcar (lambda (aa) (/ aa 0.0174533)) (coerce *5-angles* 'list)))
  (format t "*new-5-angles* = ~s~%" (mapcar (lambda (aa) (/ aa 0.0174533)) (coerce *new-5-angles* 'list)))
  )

#+tests(test-true cp5-first (< (abs (- (first *cp5-from-coords-zjs*) (first *cp5-from-cpzjs*))) 0.0001))
#+tests(test-true cp5-second (< (abs (- (second *cp5-from-coords-zjs*) (second *cp5-from-cpzjs*))) 0.0001))

(defparameter *zjs6-from-coords* (multiple-value-call 'chem:cpring-coordinates-to-zjs (chem:cpring-6-sample)))
(defparameter *cp6-from-coords* (multiple-value-list (multiple-value-call 'chem:cpring-6-coordinates-to-cremer-pople (chem:cpring-6-sample))))
(defparameter *cp6-from-coords-zjs* (multiple-value-list (chem:cpring-6-zjs-to-cremer-pople *zjs6-from-coords*)))
(defparameter *zjs6-from-cp6* (apply 'chem:cpring-6-cremer-pople-to-zjs *cp6-from-coords*))
(defparameter *cp6-from-cpzjs* (multiple-value-list (chem:cpring-6-zjs-to-cremer-pople *zjs6-from-cp6*)))

(defvar *6-coords*)
(defvar *6-indexes*)
(multiple-value-setq (*6-coords* *6-indexes*) (chem:cpring-6-sample))
(defparameter *new-6-coords* (chem:nvector-copy *6-coords*))
(defvar *6-dists* (make-nvector* (length *6-indexes*)))
(defvar *6-angles* (make-nvector* (length *6-indexes*)))
(chem:cpring-coordinates-to-distances-and-angles *6-dists* *6-angles* *6-coords* *6-indexes*)
(chem:cpring-generate-coordinates *zjs6-from-cp6* *6-dists* *6-angles* *new-6-coords* *6-indexes*)
(defvar *new-6-dists* (make-nvector* (length *6-indexes*)))
(defvar *new-6-angles* (make-nvector* (length *6-indexes*)))
(chem:cpring-coordinates-to-distances-and-angles *new-6-dists* *new-6-angles* *new-6-coords* *6-indexes*)

(defparameter *test-cp6-first* (abs (- (first *cp6-from-coords-zjs*) (first *cp6-from-cpzjs*))))
(defparameter *test-cp6-second* (abs (- (second *cp6-from-coords-zjs*) (second *cp6-from-cpzjs*))))
(defparameter *test-cp6-third* (abs (- (third *cp6-from-coords-zjs*) (third *cp6-from-cpzjs*))))
#+tests(test-true cp6-first (< *test-cp6-first* 0.0001))
#+tests(test-true cp6-second (< *test-cp6-second* 0.0001))
#+tests(test-true cp6-third (< *test-cp6-third* 0.0001))
(progn
  (format t "cpring-6-sample -> ~a~%" (multiple-value-call 'chem:cpring-mathematica (chem:cpring-6-sample)))
  (format t "cp6-from-coords*     -> ~10,7f~%" *cp6-from-coords*)
  (format t "cp6-from-coords-zjs* -> ~10,7f~%" *cp6-from-coords-zjs*)
  (format t "cp6-from-cpzjs*      -> ~10,7f~%" *cp6-from-cpzjs*)
  (format t "zjs6-from-cp6*    = ( ~{ ~10,7f~^,~} )~%" (coerce *zjs6-from-cp6* 'list))
  (format t "zjs6-from-coords* = ( ~{ ~10,7f~^,~} )~%" (coerce *zjs6-from-coords* 'list))
  (format t "*new-6-coords* = ~a~%" (chem:cpring-mathematica *new-6-coords* *6-indexes*))
  (format t "*6-dists* = ~s~%" *6-dists*)
  (format t "*new-6-dists* = ~s~%" *new-6-dists*)
  (format t "*6-angles* = ~s~%" (mapcar (lambda (aa) (/ aa 0.0174633)) (coerce *6-angles* 'list)))
  (format t "*new-6-angles* = ~s~%" (mapcar (lambda (aa) (/ aa 0.0174633)) (coerce *new-6-angles* 'list)))
  )





