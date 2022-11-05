(in-package :units)

;;;
;;; Extra functionality for units
;;;

(defmethod print-object ((obj units:quantity) stream)
  (format stream "~a" (core:repr obj)))

