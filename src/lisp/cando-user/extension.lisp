(in-package #:cando-user)

(defun initialize-amber-path ()
  (let ((amber-home (or (uiop:getenv-absolute-directory "AMBERHOME")
                        (uiop:getenv-absolute-directory "AMBER_HOME")
                        (probe-file "/opt/amber/"))))
    (cond (amber-home
           (unless (or (sys:noinform-p) (sys:noprint-p))
             (format t "Creating logical host AMBER with a root path of ~a~%"
                     amber-home))
           (setf (logical-pathname-translations "amber")
                 (list (list "**;*.*"
                             (concatenate 'string (namestring amber-home) "/**/*.*"))))
           (leap:setup-default-paths))
          ((not (or (sys:noinform-p) (sys:noprint-p)))
           (format t "No AMBERHOME or AMBER_HOME environment variable found!~%   Consider installing the latest AMBERTOOLS into /opt/amber/~%")))))

(defun initialize-rosetta-path ()
  (let ((rosetta-home (or (uiop:getenv-absolute-directory "ROSETTA_HOME")
                          (probe-file "/opt/rosetta/"))))
    (cond (rosetta-home
           (unless (or (sys:noinform-p) (sys:noprint-p))
             (format t "Creating logical host ROSETTA with a root path of ~a~%"
                     rosetta-home))
           (setf (logical-pathname-translations "rosetta")
                 (list (list "**;*.*"
                             (concatenate 'string (namestring rosetta-home) "/**/*.*"))))
           (leap:setup-default-paths))
          ((not (or (sys:noinform-p) (sys:noprint-p)))
           (format t "No ROSETTA_HOME environment variable found!~%   Consider cloning the https://github.com/RosettaCommons/rosetta.git into /opt/rosetta/~%")))))

(defun initialize-threading ()
  (let ((threads (if (ext:getenv "CANDO_THREADS")
                     (parse-integer (ext:getenv "CANDO_THREADS"))
                     (sys:num-logical-processors))))
    (unless (member :no-auto-lparallel *features*)
      (unless (or (sys:noinform-p) (sys:noprint-p))
        (format t "Creating LPARALLEL kernel with ~a threads.~%" threads))
      (setf lparallel:*kernel* (lparallel:make-kernel threads)))))

(defun initialize-cando-user-package ()
  (sys:symbol-global-value-set '*package* (find-package :cando-user))
  (setf *package* (find-package :cando-user)))

(defun initialize-cando-user ()
  (initialize-amber-path)
  (initialize-rosetta-path)
  (initialize-threading)
  (initialize-cando-user-package)
  )

(push #'initialize-cando-user core:*initialize-hooks*)

(defun terminate-cando-user ()
  (setf (logical-pathname-translations "amber") nil)
  (unless (member :no-auto-lparallel *features*)
    (lparallel:end-kernel :wait t)))

(push #'terminate-cando-user core:*terminate-hooks*)
