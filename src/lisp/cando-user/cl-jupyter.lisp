(in-package :cando-user)

(defun cl-jupyter-kernel-start ()
  (let ((amber-home
          (namestring (or (uiop:ensure-directory-pathname (ext:getenv "AMBERHOME"))
                          (uiop:ensure-directory-pathname "/home/app/amber16-data/")))))
    (setf (logical-pathname-translations "amber")
          (list (list "**;*.*" (concatenate 'string amber-home "/**/*.*"))))
    (format t "Setting amber host pathname translation -> ~a~%" amber-home))
  (let ((*readtable* (copy-readtable)))
    (set-macro-character #\~ #'leap.core:leap-lookup-variable-reader-macro)
    (let ((cl-jup (find-symbol "KERNEL-START" "CL-JUPYTER")))
      (if cl-jup
          (funcall cl-jup)
          (error "cl-jupyter is not installed")))))
