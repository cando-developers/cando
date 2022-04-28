(ql:quickload :sketch2d)


(defparameter *default-pathname-defaults* (translate-logical-pathname #P"sys:extensions;cando;src;lisp;sketch2d;"))
(defparameter *agg* (cando:load-chem-draw-aggregate "/tmp/methane.cdxml"))
(defparameter *agg* (cando:load-mol2 "examples/017.mol2"))
(defparameter *mol* (chem:content-at *agg* 0))
(cando:save-mol2 *mol* "/tmp/test0.mol2" :use-sybyl-types t)
*mol*
(let ((chem:*verbose* 0))(defparameter *nm* (sketch2d::sketch2d *mol*)))
(cando:save-mol2 *mol* "/tmp/test2.mol2" :use-sybyl-types t)


(progn
  (defparameter *agg* (cando:load-chem-draw-aggregate "/tmp/methane.cdxml"))
  (defparameter *mol* (chem:content-at *agg* 0))
  (defparameter *m* (sketch2d::adjust-structure *mol*))
  (cando:jostle *m* 1.0 t)
  (cando:do-atoms (atom *m*) (chem:set-type atom :sketch))
  (defparameter *sk* (sketch2d::generate-sketch-function *m*))
  (defparameter *v* (make-array (chem:get-nvector-size *sk*) :element-type 'double-float :initial-element 0.0d0))
  (defparameter *f* (make-array (chem:get-nvector-size *sk*) :element-type 'double-float :initial-element 0.0d0))
  *f*
  *v*
  (setf (elt *v* 0) 2.0d0)
  (let ((chem::*verbose* 0))
    (loop for x from 0.5 below 100.0 by 0.5
          for y = (progn
                    (setf (elt *v* 0) x)
                    (chem:evaluate-energy-force *sk* *v* t *f* ))
          do (format t "~a ~a~%" x y)))
  )





(ql:quickload :sketch2d)


(progn
  (defparameter *default-pathname-defaults* (translate-logical-pathname #P"sys:extensions;cando;src;lisp;sketch2d;"))
  (defparameter *agg* (cando:load-chem-draw-aggregate "/tmp/methane.cdxml"))
  (defparameter *agg* (cando:load-mol2 "examples/017.mol2"))
  (defparameter *mol* (chem:content-at *agg* 0))
  (cando:jostle *mol* 0.1 t)
  (cando:save-mol2 *mol* "/tmp/test0.mol2" :use-sybyl-types t)
  (let ((chem:*verbose* 0))(defparameter *nm* (sketch2d::sketch2d *mol*)))
  (cando:save-mol2 *mol* "/tmp/test2.mol2" :use-sybyl-types t)
  )
