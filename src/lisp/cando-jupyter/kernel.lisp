(in-package :cando-jupyter)

(defvar +display-name+ "Cando")
(defvar +language+ "cando")

(defparameter *leap-syntax* t)

(defun leap-syntax-enable (on)
  (setq *leap-syntax* on))

(defclass kernel (common-lisp-jupyter:kernel)
  ()
  (:default-initargs
    :package (find-package :cando-user)
    :banner "cando-jupyter: a Cando Jupyter kernel
(C) 2020 Christian Schafmeister (LGPL-3.0)"))

; The readtable was copied in the original code. Do we need to do that?
(defmethod jupyter:start :after ((k kernel))
  (let ((amber-home
          (namestring (or (if (ext:getenv "AMBERHOME")
                              (probe-file (ext:getenv "AMBERHOME"))
                              "/usr/local/amber/")
                          (probe-file "/home/app/amber16-data/")
                          "/dev/null"))))
    (setf (logical-pathname-translations "amber")
          (list (list "**;*.*" (concatenate 'string amber-home "/**/*.*"))))
    (jupyter:inform :info k "Setting amber host pathname translation -> ~a" amber-home)))


(defun lisp-code-p (code)
  (do ((index 0 (1+ index)))
      ((>= index (length code)))
    (case (char code index)
      ((#\( #\*)
        (return t))
      ((#\space #\tab #\newline))
      (otherwise
        (return nil)))))

(defun leap-read (code)
  (jupyter:handling-errors
    (architecture.builder-protocol:with-builder ('list)
      (esrap:parse 'leap.parser::leap code))))

(defun leap-eval (ast)
  (jupyter:handling-errors
    ; Not sure if leap needs its own result handler
    (jupyter:make-lisp-result
      (leap.core:evaluate 'list ast leap.core:*leap-env*))))

(defmethod jupyter:evaluate-code ((k kernel) code)
  (if (or (not *leap-syntax*)
          (lisp-code-p code))
    (call-next-method)
    (let ((ast (leap-read code)))
      (if (typep ast 'jupyter:result)
        ast
        (list (leap-eval ast))))))


(defclass cando-installer (jupyter:installer)
  ()
  (:default-initargs
    :class 'kernel
    :language +language+
    :display-name +display-name+
    :kernel-name +language+
    :systems '(:cando-jupyter)))

(defclass system-installer (jupyter:system-installer cando-installer)
  ()
  (:documentation "cando system installer."))

(defclass user-installer (jupyter:user-installer cando-installer)
  ()
  (:documentation "cando user installer."))

(defclass user-image-installer (jupyter:user-image-installer cando-installer)
  ()
  (:documentation "cando user image installer."))

(defmethod jupyter:command-line ((instance user-installer))
  "Get the command line for a user installation."
  (let ((implementation (jupyter:installer-implementation instance)))
    (list
      (or implementation
          (first (uiop:raw-command-line-arguments))
          (format nil "~(~A~)" (uiop:implementation-type)))
      "--eval" "(ql:quickload :cando-jupyter)"
      "--eval" "(jupyter:run-kernel 'cando-jupyter:kernel #\"{connection_file}\")")))

(defmethod jupyter:command-line ((instance system-installer))
  "Get the command line for a system installation."
  (let ((implementation (jupyter:installer-implementation instance)))
    (list
      (or implementation
          (first (uiop:raw-command-line-arguments))
          (format nil "~(~A~)" (uiop:implementation-type)))
      "--load" (namestring (jupyter:installer-path instance :root :program :bundle))
      "--eval" "(asdf:load-system :cando-jupyter)"
      "--eval" "(jupyter:run-kernel 'cando-jupyter:kernel #\"{connection_file}\")")))

(defun install (&key bin-path system local prefix)
  "Install Cando kernel.
- `bin-path` specifies path to LISP binary.
- `system` toggles system versus user installation.
- `local` toggles `/usr/local/share versus` `/usr/share` for system installations.
- `prefix` key specifies directory prefix for packaging.
"
  (jupyter:install
    (make-instance
      (if system
        'system-installer
        'user-installer)
      :implementation bin-path
      :local local
      :prefix prefix)))

(defun install-image (&key prefix)
  "Install Cando kernel based on image.
- `prefix` key specifies directory prefix for packaging."
  (jupyter:install
    (make-instance 'user-image-installer
      :prefix prefix)))


(defun run-kernel-from-slime ()
  (bordeaux-threads:make-thread
   (lambda ()
     (loop for time from 1 below 100
           do (format t "Starting kernel for the ~a time~%" time)           
           do (jupyter:run-kernel 'cando-jupyter:kernel (core:argv (1- (core:argc))))))
   :name :jupyter))
