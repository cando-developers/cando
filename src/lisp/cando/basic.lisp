(in-package :cando)

;;; Empty

(defun shallow-copy (object)
  (let* ((class (class-of object))
         (copy (allocate-instance class)))
    (dolist (slotd (clos:class-slots class) copy)
      (when (clos:slot-boundp-using-class class object slotd)
        (setf (clos:slot-value-using-class class copy slotd)
              (clos:slot-value-using-class class object slotd))))))


(defmacro assign (vars expr)
  (if (symbolp vars)
      `(defparameter ,vars ,expr)       ; single variable
      `(progn
         ,@(mapcar (lambda (v) `(defvar ,v)) vars) ; declare each
         (setf (values ,@(loop for v in vars collect `(symbol-value ',v))) ,expr)
         ',vars)))
