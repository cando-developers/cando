(in-package :smarts)

(defun make-chem-info (&key tests smarts)
  (let ((ci (core:make-cxx-object 'chem:chem-info)))
    (chem:compile-smarts ci smarts)
    ci))
