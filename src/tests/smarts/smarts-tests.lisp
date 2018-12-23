(progn
  (load "source-dir:extensions;cando;src;lisp;start-cando.lisp")
  (format t "DONE DONE DONE~%"))

(ext:chdir "~/Development/cando/extensions/cando/src/tests/smarts/" t)

(defparameter *agg* (cando:load-chem-draw-aggregate "lactam6.cdxml"))

(defparameter *amide* (chem:compile-smarts "[N:1][C:2]=[O:3]"))

(defparameter *amide* (chem:compile-smarts "[N:1][C:2]"))

(chem:map-atoms 'nil (lambda (a &aux match)
                       (when (eq (chem:get-element a) :N)
                         (smarts:print-smarts *amide*)
                         (format t "Checking ~s ~s -> ~s~%" *amide* a (chem:matches *amide* a))
                         (setf match (chem:matches *amide* a))
                         (when match
                           (print match))))
                *agg*)
                             


