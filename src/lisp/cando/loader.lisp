(in-package :cando)
;;; Load a single aggregate from a cdxml file
;;; Return it
(defun load-cdxml-aggregate (filename)
  (let ((cd-only (with-open-file
                     (fin (probe-file filename) :direction :input)
                   (chem:make-chemdraw fin))))
    (chem:as-aggregate cd-only)))

(apropos "chemdraw")

