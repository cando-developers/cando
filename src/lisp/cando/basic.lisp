(in-package :cando)

;;; Empty

(defun shallow-copy (object)
  (let* ((class (class-of object))
         (copy (allocate-instance class)))
    (dolist (slotd (clos:class-slots class) copy)
      (when (clos:slot-boundp-using-class class object slotd)
        (setf (clos:slot-value-using-class class copy slotd)
              (clos:slot-value-using-class class object slotd))))))
