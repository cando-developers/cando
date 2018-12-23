(progn
  (load "source-dir:extensions;cando;src;lisp;start-cando.lisp")
  (format t "DONE DONE DONE~%"))

(ext:chdir "~/Development/cando/extensions/cando/src/tests/smarts/" t)

(defparameter *agg1* (cando:load-chem-draw-aggregate "linearAmide.cdxml"))
(defparameter *agg2* (cando:load-chem-draw-aggregate "lactam6.cdxml"))
(defparameter *amide1* (chem:compile-smarts "[N:1][C:2](~[O:3])[C:4]"))
(defparameter *amide2* (chem:compile-smarts "[N:1][C:2]"))
(defparameter *lactam6* (chem:compile-smarts "[N:1]1[C:2](=[O:3])[C:4][C:5][C:6][C:7][N:8]1"))

(defun smarts-find (smarts agg)
  (let (the-match)
    (chem:map-atoms 'nil (lambda (a &aux match)
                           (setf match (chem:matches smarts a))
                           (when match
                             (setf the-match match)))
                    agg)
    the-match))
                             


