(in-package #:topology-shape-key-tests)

;;; Policy dispatch is tested without depending on the foldamer package.
(defclass ten-degree-test-foldamer () ())
(defclass twenty-degree-test-foldamer () ())

(defmethod topology:bin-shape-component
    ((foldamer ten-degree-test-foldamer) kind parameters)
  (topology:bin-shape-component-with-widths kind parameters))

(defmethod topology:unbin-shape-component
    ((foldamer ten-degree-test-foldamer) kind bins)
  (topology:unbin-shape-component-with-widths kind bins))

(defmethod topology:bin-shape-component
    ((foldamer twenty-degree-test-foldamer) kind parameters)
  (topology:bin-shape-component-with-widths kind parameters :dihedral-width 20))

(defmethod topology:unbin-shape-component
    ((foldamer twenty-degree-test-foldamer) kind bins)
  (topology:unbin-shape-component-with-widths kind bins :dihedral-width 20))

(defun test-shape-key-binning ()
  (let ((foldamer (make-instance 'ten-degree-test-foldamer)))
    (dolist (pair '((0 0) (360 0) (-360 0) (359 0) (-1 0)
                    (350 35) (-10 35) (180 18) (-180 18)
                    (5 0) (15 2) (355 0)))
      (destructuring-bind (degrees bin) pair
        (let ((angle (* degrees (/ pi 180d0))))
          (assert (equal (list bin)
                         (topology:bin-shape-component foldamer :dihedral (list angle))))
          (assert (equal (list 2 bin)
                         (topology:bin-shape-component foldamer :ring5 (list 0.2d0 angle))))
          (assert (equal (list 2 18 bin)
                         (topology:bin-shape-component foldamer :ring6 (list 0.2d0 pi angle)))))))
    (dotimes (bin 36)
      (dolist (kind '(:dihedral :ring5 :ring6))
        (let ((bins (ecase kind (:dihedral (list bin))
                          (:ring5 (list 3 bin)) (:ring6 (list 3 9 bin)))))
          (assert (equal bins (topology:bin-shape-component
                              foldamer kind
                              (topology:unbin-shape-component foldamer kind bins)))))))
    (dolist (bins '((-1) (36) (12.0) (1 2)))
      (assert (signals-error-p
               (lambda () (topology:unbin-shape-component foldamer :dihedral bins)))))
    (assert (signals-error-p
             (lambda () (topology:unbin-shape-component foldamer :ring6 '(2 19 0)))))
    (assert (signals-error-p
             (lambda () (topology:bin-shape-component foldamer :ring5 '(-0.1d0 0d0))))))
  (let ((foldamer (make-instance 'twenty-degree-test-foldamer)))
    (assert (equal '(6) (topology:bin-shape-component
                        foldamer :dihedral (list (* 120 (/ pi 180d0))))))
    (assert (< (abs (- (first (topology:unbin-shape-component foldamer :dihedral '(6)))
                      (* 120 (/ pi 180d0)))) 1d-12)))
  (format t "~&SHAPE-KEY-BINNING-TESTS-PASSED~%")
  t)
