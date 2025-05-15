(in-package :cando-user)

(defvar *started-swank* nil)
(defvar *swank-port* nil)

(defvar *slime-home* nil
  "You can set cando-user:*slime-home* in your .clasprc to the pathname where slime is installed to start swank")

(defun start-swank-server (&optional (port 4005))
  (format t "Checking SLIME_HOME and cando-user:*slime-home*~%")
  (let ((slime-home (or (and (ext:getenv "SLIME_HOME") (probe-file (pathname (ext:getenv "SLIME_HOME"))))
                        (and *slime-home* (probe-file (pathname *slime-home*))))))
    (unless slime-home
      (format t "Neither were set~%")
      (let ((slime-dirs (directory "~/.emacs.d/elpa/slime*/")))
        (loop with max-date = 0
              with max-slime-dir = nil
              for slime-dir in slime-dirs
              for date = (file-write-date slime-dir)
              when (> date max-date)
                do (setf max-date date
                         max-slime-dir slime-dir)
              finally (progn
                        (format t "Found slime-home at ~s~%" max-slime-dir)
                        (setf slime-home max-slime-dir)))))
    (if slime-home
        (let ((swank-loader (probe-file (merge-pathnames "swank-loader.lisp" slime-home))))
          (load swank-loader)
          (let* ((slime-cache (merge-pathnames "fasl/" slime-home))
                 (slime-cache-symbol (find-symbol "*FASL-DIRECTORY*" :SWANK-LOADER)))
            (setf (symbol-value slime-cache-symbol) slime-cache)
            (let ((swank-loader-init (find-symbol "INIT" "SWANK-LOADER")))
              (funcall swank-loader-init :delete nil :reload nil :load-contribs nil))
            (let ((swank-create-server (find-symbol "CREATE-SERVER" "SWANK")))
              (mp:process-run-function 'swank-main
                                       (lambda ()
                                         (format *debug-io* "Starting server thread~%")
                                         (let ((port (funcall swank-create-server
                                                              :port port
                                                              :interface "0.0.0.0")))
                                           (format t "Started swank server on port ~d~%" port)
                                           (setf *swank-port* port))
                                         (format *debug-io* "Leaving server thread~%")))
              #+(or)(sleep 1)
              (setf *started-swank* t))))
        (error "Could not determine directory for slime - set SLIME_HOME or cando-user:*slime-home*"))))

(defun start-swank (&optional (port 4005))
  "Start a swank server to connect slime.
Within emacs use M-x slime-connect to connect."
  ;; Bad!  This is hard-coded to work with docker
  (if *started-swank*
      (format t "Swank is already running on port ~a~%" *swank-port*)
      (if (find-package :cl-ipywidgets)
          (let ((save-jupyter-cell-state-symbol (find-symbol "SAVE-JUPYTER-CELL-STATE" :cl-ipywidgets)))
            (if save-jupyter-cell-state-symbol
                (funcall save-jupyter-cell-state-symbol)
                (error "cl-ipywidgets isn't loaded - can't save cell state"))
            (start-swank-server port))
          (start-swank-server port))))


(defvar *started-slynk* nil)
(defvar *slynk-port* nil)

(defvar *slime-home* nil
  "You can set cando-user:*slime-home* in your .clasprc to the pathname where slime is installed to start slynk")

(defun start-slynk-server (&optional (port 4005))
  (format t "Checking SLY_HOME and cando-user:*sly-home*~%")
  (let ((sly-home (or (and (ext:getenv "SLY_HOME") (probe-file (pathname (ext:getenv "SLY_HOME"))))
                        (and *sly-home* (probe-file (pathname *sly-home*))))))
    (unless sly-home
      (format t "Neither were set~%")
      (let ((sly-dirs (directory "~/Development/sly/")))
        (loop with max-date = 0
              with max-sly-dir = nil
              for sly-dir in sly-dirs
              for date = (file-write-date sly-dir)
              when (> date max-date)
                do (setf max-date date
                         max-sly-dir sly-dir)
              finally (progn
                        (format t "Found sly-home at ~s~%" max-sly-dir)
                        (setf sly-home max-sly-dir)))))
    (if sly-home
        (let ((slynk-loader (probe-file (merge-pathnames "slynk/slynk-loader.lisp" sly-home))))
          (load slynk-loader)
          (let* ((sly-cache (merge-pathnames "fasl/" sly-home))
                 (sly-cache-symbol (find-symbol "*FASL-DIRECTORY*" :SLYNK-LOADER)))
            (setf (symbol-value sly-cache-symbol) sly-cache)
            (let ((slynk-loader-init (find-symbol "INIT" "SLYNK-LOADER")))
              (funcall slynk-loader-init :delete nil :reload nil :load-contribs nil))
            (let ((slynk-create-server (find-symbol "CREATE-SERVER" "SLYNK")))
              (mp:process-run-function 'slynk-main
                                       (lambda ()
                                         (format *debug-io* "Starting server thread~%")
                                         (let ((port (funcall slynk-create-server
                                                              :port port
                                                              :interface "0.0.0.0")))
                                           (format t "Started slynk server on port ~d~%" port)
                                           (setf *slynk-port* port))
                                         (format *debug-io* "Leaving server thread~%")))
              #+(or)(sleep 1)
              (setf *started-slynk* t))))
        (error "Could not determine directory for sly - set SLY_HOME or cando-user:*sly-home*"))))

(defun start-slynk (&optional (port 4005))
  "Start a slynk server to connect sly.
Within emacs use M-x sly-connect to connect."
  ;; Bad!  This is hard-coded to work with docker
  (if *started-slynk*
      (format t "Slynk is already running on port ~a~%" *slynk-port*)
      (if (find-package :cl-ipywidgets)
          (let ((save-jupyter-cell-state-symbol (find-symbol "SAVE-JUPYTER-CELL-STATE" :cl-ipywidgets)))
            (if save-jupyter-cell-state-symbol
                (funcall save-jupyter-cell-state-symbol)
                (error "cl-ipywidgets isn't loaded - can't save cell state"))
            (start-slynk-server port))
          (start-slynk-server port))))
