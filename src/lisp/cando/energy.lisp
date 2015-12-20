(in-package :energy)

(defparameter *ff* nil
  "Store the force-field")

(defun setup-amber ()
  (let ((*default-pathname-defaults*
         #P"~/Development/clasp/projects/cando/work/forceField/"))
    (defparameter *parms*
      (let ((parms (chem:make-read-amber-parameters)))
        (with-open-file (fin "ATOMTYPE_GFF.DEF" :direction :input)
          (chem:read-types parms fin))
        (with-open-file (fin "gaff.dat" :direction :input)
          (chem:read-parameters parms fin)
          parms)))
    (defparameter *ff* (chem:get-force-field *parms*))))


(defun minimize (agg &key (restraints-on t) (force-field *ff*))
  "Minimize the conformational energy"
    (format t "Got minimizer~%")
  (let* ((energy-func (chem:make-energy-function agg force-field))
         (minimizer (chem:make-minimizer :energy-function energy-func)))
    (unless restraints-on
      (let ((restraint-term (chem:get-anchor-restraint-component energy-func)))
        (chem:disable restraint-term)))
    (cando:configure-minimizer minimizer
                               :max-steepest-descent-steps 1000
                               :max-conjugate-gradient-steps 50000
                               :max-truncated-newton-steps 0)
    (chem:enable-print-intermediate-results minimizer)
    (chem:set-option energy-func 'chem:nonbond-term nil)
    (cando:minimize-no-fail minimizer)
    (chem:set-option energy-func 'chem:nonbond-term t)
    (cando:minimize-no-fail minimizer)))
