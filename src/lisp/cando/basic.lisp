(in-package :cando)

(defun lazy-setup-lparallel-kernel ()
  (unless lparallel:*kernel*
    (setf lparallel:*kernel* (lparallel:make-kernel (core:num-logical-processors)))))

(defmacro |.=| (a b)
  `(defparameter ,a ,b))

