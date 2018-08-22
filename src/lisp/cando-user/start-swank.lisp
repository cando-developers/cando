(in-package :cando-user)

(progn
  (defvar *started-swank* nil)
  (defun start-swank (&optional (port 4005))
    ;; Bad!  This is hard-coded to work with docker
    (if *started-swank*
        (format t "Swank is already running~%")
        (progn
          (funcall (find-symbol "SAVE-JUPYTER-CELL-STATE" :cl-ipywidgets))
          (let* ((slime-home (cond
                               ((and (ext:getenv "SLIME_HOME") (probe-file (pathname (ext:getenv "SLIME_HOME")))))
                               ((probe-file #P"/home/app/slime/"))
                               (t "Could not determine directory for slime - set SLIME_HOME")))
                 (swank-loader (probe-file (merge-pathnames "swank-loader.lisp" slime-home))))
            ;; (format t "swank-loader -> ~a~%" swank-loader)
            (load swank-loader))
          (let ((swank-loader-init (find-symbol "INIT" "SWANK-LOADER")))
            (funcall swank-loader-init :delete nil :reload nil :load-contribs nil))
          (let ((swank-create-server (find-symbol "CREATE-SERVER" "SWANK")))
            (mp:process-run-function 'swank-main
                                     (lambda () (funcall swank-create-server
                                                         :port port
                                                         :interface "0.0.0.0")))
            (setf *started-swank* t))))))
