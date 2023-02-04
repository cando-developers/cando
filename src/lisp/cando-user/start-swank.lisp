(in-package :cando-user)

(defvar *started-swank* nil)

(defvar *swank-home* nil
  "You can set cando-user:*swank-home* in your .clasprc to the pathname where slime is installed to start swank")

(defun start-swank-server (&optional (port 4005))
  (format t "Checking SLIME_HOME and cando-user:*swank-home*~%")
  (let ((slime-home (or (and (ext:getenv "SLIME_HOME") (probe-file (pathname (ext:getenv "SLIME_HOME"))))
                         (and *swank-home* (probe-file (pathname *swank-home*))))))
    (if slime-home
        (let ((swank-loader (probe-file (merge-pathnames "swank-loader.lisp" slime-home))))
          (load swank-loader)
          (let* ((slime-cache (merge-pathnames "fasl/" slime-home))
                 (slime-cache-symbol (find-symbol "*FASL-DIRECTORY*" :SWANK-LOADER)))
            (setf (symbol-value slime-cache-symbol) slime-cache)
            (let ((swank-loader-init (find-symbol "INIT" "SWANK-LOADER")))
              (funcall swank-loader-init :delete nil :reload nil :load-contribs nil))
            (let ((swank-create-server (find-symbol "CREATE-SERVER" "SWANK")))
              (format t "About to start swank server~%")
              (mp:process-run-function 'swank-main
                                       (lambda ()
                                         (format *debug-io* "Starting server thread~%")
                                         (let ((port (funcall swank-create-server
                                                              :port port
                                                              :interface "0.0.0.0")))
                                           (format t "Started swank server on port ~d~%" port))
                                         (format *debug-io* "Leaving server thread~%")))
              (sleep 2)
              (format t "Started swank server~%")
              (setf *started-swank* t))))
        (error "Could not determine directory for slime - set SLIME_HOME"))))

(defun start-swank (&optional (port 4005))
  "Start a swank server to connect slime.
Within emacs use M-x slime-connect to connect."
  ;; Bad!  This is hard-coded to work with docker
  (if *started-swank*
      (format t "Swank is already running~%")
      (if (find-package :cl-ipywidgets)
          (let ((save-jupyter-cell-state-symbol (find-symbol "SAVE-JUPYTER-CELL-STATE" :cl-ipywidgets)))
            (if save-jupyter-cell-state-symbol
                (funcall save-jupyter-cell-state-symbol)
                (error "cl-ipywidgets isn't loaded - can't save cell state"))
            (start-swank-server port))
          (start-swank-server port))))
