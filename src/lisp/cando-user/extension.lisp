(in-package #:cando-user)

(defun initialize-cando-user ()
  (let ((amber-home (or (uiop:getenv-absolute-directory "AMBERHOME")
                        (probe-file "/opt/amber/")))
        (threads (if (ext:getenv "CANDO_THREADS")
                     (parse-integer (ext:getenv "CANDO_THREADS"))
                     (sys:num-logical-processors))))
    (cond (amber-home
           (unless (or (sys:noinform-p) (sys:noprint-p))
             (format t "Creating logical host AMBER with a root path of ~a~%"
                     amber-home))
           (setf (logical-pathname-translations "amber")
                 (list (list "**;*.*"
                             (concatenate 'string (namestring amber-home) "/**/*.*"))))
           (leap:setup-default-paths))
          ((not (or (sys:noinform-p) (sys:noprint-p)))
           (format t "No AMBERHOME environment variable found!~%")))
    (unless (member :no-auto-lparallel *features*)
      (unless (or (sys:noinform-p) (sys:noprint-p))
        (format t "Creating LPARALLEL kernel with ~a threads.~%" threads))
      (setf lparallel:*kernel* (lparallel:make-kernel threads)))
    (sys:symbol-global-value-set '*package* (find-package :cando-user))
    (setf *package* (find-package :cando-user))))

(push #'initialize-cando-user core:*initialize-hooks*)

(defun terminate-cando-user ()
  (setf (logical-pathname-translations "amber") nil)
  (unless (member :no-auto-lparallel *features*)
    (lparallel:end-kernel :wait t)))

(push #'terminate-cando-user core:*terminate-hooks*)    
