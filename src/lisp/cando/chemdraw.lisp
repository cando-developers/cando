(in-package :chemdraw)
;;; Load a single aggregate from a cdxml file
;;; Return it
(defun load-only-aggregate (filename)
  (let ((cd-only (with-open-file
                     (fin (probe-file filename) :direction :input)
                   (chem:make-chem-draw fin))))
    (chem:as-aggregate cd-only)))

(apropos "chemdraw")

