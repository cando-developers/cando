;;;
;;; Build the kernels - use sbcl
;;;


(load "~/quicklisp/setup.lisp")

(asdf:load-asd (probe-file "cando-kernels.asd"))

(setf (uiop:getenv "KERNEL_PATH")
        (namestring (truename (asdf:system-relative-pathname
                               :cando-kernels
                               "../../../include/cando/chem/energyKernels/"))))

(ql:quickload :mathkernel)


(handler-bind ((error (lambda (e)
                          (format *error-output* "~&### ~a~%" e)
                          (sb-debug:print-backtrace :count 60 :stream *error-output*)
                          (finish-output *error-output*))))
  (load "./energy-kernels.lisp"))

