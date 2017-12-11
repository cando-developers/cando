(in-package :cando-user)

(defun cl-jupyter-kernel-start ()
  (let ((*readtable* (copy-readtable)))
    (set-macro-character #\~ #'leap.core:leap-lookup-variable-reader-macro)
    (let ((cl-jup (find-symbol "KERNEL-START" "CL-JUPYTER")))
      (if cl-jup
          (funcall cl-jup)
          (error "cl-jupyter is not installed")))))
