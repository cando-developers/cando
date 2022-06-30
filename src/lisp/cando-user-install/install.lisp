(in-package #:cando-user-install)

(defvar *snapshot-timestamp* nil)

(defun iso8601-now ()
  (multiple-value-bind (second minute hour date month year day daylight-p zone)
      (decode-universal-time (get-universal-time))
    (declare (ignore second day))
    (format nil "~4d-~2,'0d-~2,'0dT~2,'0d:~2,'0d-~2,'0d" year month date hour minute
            (if daylight-p zone (1- zone)))))

(defun install ()
  (let* ((bin-path (merge-pathnames (make-pathname :directory '(:relative ".local" "bin"))
                                    (uiop:getenv-absolute-directory "HOME")))
         (scando (merge-pathnames "scando" bin-path))
         (scleap (merge-pathnames "scleap" bin-path)))
    (cando-jupyter:install :load-system nil :bin-path scando)
    (ensure-directories-exist bin-path)
    (uiop:run-program (format nil "ln -sf ~a ~a" scando scleap)
                      :ignore-error-status t)
    (setf *snapshot-timestamp* (iso8601-now))
    (clos:compile-all-generic-functions)
    (gctools:save-lisp-and-die (namestring scando) :executable t)))

(defun updatep ()
  (or (not (search (software-version)
                   (uiop:run-program "cando --version"
                                     :ignore-error-status t
                                     :output '(:string :stripped t)
                                     :error-output '(:string :stripped t))))
       #+quicklisp (some (lambda (dist)
                           (and (ql:subscribedp dist)
                                (ql:available-update dist)))
                         (ql:all-dists))))

(defun update ()
  (when (and *snapshot-timestamp* (updatep))
    (uiop:run-program "cando-user-install"
                      :ignore-error-status t
                      :output '(:string :stripped t)
                      :error-output '(:string :stripped t))))
