(cl:in-package #:rosetta.params.load)


(defun load-rosetta-params ()
  (format t "rosetta-params-load loading and registering :rosetta force field~%")
  (let* ((ff (rosetta.params:load-lk-solvation-force-field))
         (cff (core:make-cxx-object 'chem:combined-force-field)))
    (chem:add-shadowing-force-field cff ff nil)
    (setf (gethash :rosetta leap.core:*force-fields*) cff)))

(eval-when (:load-toplevel :execute)
  (load-rosetta-params))
