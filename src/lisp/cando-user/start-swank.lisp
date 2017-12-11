(in-package :cando-user)

(progn
  (defparameter *started-swank* nil)
  (defun start-swank ()
    ;; Bad!  This is hard-coded to work with docker
    (if *started-swank*
        (format t "Swank is already running~%")
        (progn
          (let ((swank-loader (probe-file "~/Development/slime/swank-loader.lisp")))
            (format t "swank-loader -> ~a~%" swank-loader)
            (load swank-loader))
          (let ((swank-loader-init (find-symbol "INIT" "SWANK-LOADER")))
            (funcall swank-loader-init :delete nil :reload nil :load-contribs nil))
          (let ((swank-create-server (find-symbol "CREATE-SERVER" "SWANK")))
            (mp:process-run-function 'swank-main
                                     (lambda () (funcall swank-create-server
                                                         :port 4005
                                                         :interface "0.0.0.0")))
            (setf *started-swank* t))))))
