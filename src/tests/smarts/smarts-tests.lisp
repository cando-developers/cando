(progn
  (load "sys:extensions;cando;src;lisp;start-cando.lisp")
  (format t "DONE DONE DONE~%"))

(ext:chdir "~/Development/cando/extensions/cando/src/tests/smarts/" t)

(defun <isb> (a) (print a) (eq (chem:get-element a) :C))
(defparameter *agg1* (cando:load-chem-draw-aggregate "linearAmide.cdxml"))
(defparameter *agg2* (cando:load-chem-draw-aggregate "lactam6.cdxml"))
(defparameter *amide1* (chem:compile-smarts "[N:1][C:2](~[O:3])[C:4]"))
(defparameter *amide2* (chem:compile-smarts "[N:1][C:2]"))
(defparameter *lactam6* (chem:compile-smarts "[N:1]1[C:2](=[O:3])[<isb>:4][C:5][C:6][C:7][N:8]1"
                                             :tests (list (cons :ISB #'<ISB>))))
(defparameter *match-lactone-lactam* (chem:compile-smarts "O"))
(defparameter *match-lactone-lactam* (chem:compile-smarts "[$(O=[C:3]1-[O,N]~[*]~[*]~[*]1),$(O=[C:4]1-[O,N]~[*]~[*]~[*]~[*]1),$(O=[C:5]1-[O,N]~[*]~[*]~[*]~[*]~[*]1),$(O=[C:6]1-[O,N]~[*]~[*]~[*]~[*]~[*]~[*]1)]"))


(defun smarts-find (smarts agg)
  (let (the-match)
    (chem:map-atoms 'nil (lambda (a &aux match)
                           (setf match (chem:matches smarts a))
                           (when match
                             (format t "Found a match starting on atom: ~a~%" a)
                             (setf the-match match)))
                    agg)
    the-match))
                             


