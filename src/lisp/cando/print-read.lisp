(in-package :cando)

;; use
(defmethod cl:print-object ((object core:cxx-object) stream)
  (if (and *print-readably* (core:fieldsp object))
      (core:print-cxx-object object stream)
      (call-next-method)))


(defun read-cxx-object (stream char n)
  (declare (ignore char))
  (let ((description (read stream t nil t)))
    (apply #'core:make-cxx-object (car description) (cdr description))))

(defun as-string (obj)
  (let ((*print-readably* t)
        (*print-pretty* nil)
        (*print-circle* t))
    (with-output-to-string (sout)
      (print obj sout))))


(defun from-string (str)
  (with-input-from-string (sin str)
    (read sin)))


(set-dispatch-macro-character #\# #\I #'read-cxx-object)

