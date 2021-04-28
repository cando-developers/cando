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
(C) 2020-2021 Christian Schafmeister (LGPL-3.0)"))


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


(defmethod jupyter:start :before ((k kernel))
  (when lparallel:*kernel*
    (error "The lparallel:*kernel* is not NIL and it must be when jupyter starts.   Add -f no-auto-lparallel to cando startup"))
  (setf lparallel:*kernel* (lparallel:make-kernel (core:num-logical-processors))))


(defmethod jupyter:stop :after ((k kernel))
  #+(or)(format t "jupyter:stop :after (mp:all-processes) -> ~s~%" (mp:all-processes))
  (lparallel:end-kernel :wait t))


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
  (architecture.builder-protocol:with-builder ('list)
    (esrap:parse 'leap.parser::leap code)))


(defun leap-eval (ast)
  (jupyter:handling-errors
    (leap.core:evaluate 'list ast leap.core:*leap-env*)))


(defmethod jupyter:evaluate-code ((k kernel) code)
  (if (or (not *leap-syntax*)
          (lisp-code-p code))
    (call-next-method)
    (multiple-value-bind (ast ename evalue traceback)
                         (jupyter:handling-errors (leap-read code))
      (if ename
        (values ename evalue traceback)
        (dolist (expr (cadadr ast) (values))
          (unless (eq :comment (caar expr))
            (multiple-value-bind (result ename evalue traceback)
                                 (leap-eval (list :leap (list :instruction (list expr))))
            (when ename
              (return (values ename evalue traceback)))
            (jupyter:execute-result result))))))))


(defun leap-locate (ast cursor-pos &optional child-pos parents)
  (if (listp (car ast))
    (loop for child in ast
          for index from 0
          for fragment = (leap-locate child cursor-pos index parents)
          when fragment
          return fragment)
    (let ((bounds (getf (cddr ast) :bounds))
          (new-parents (cons ast parents)))
      (unless (and bounds
                   (or (< cursor-pos (car bounds))
                       (>= cursor-pos (cdr bounds))))
        (setf (getf (cddr ast) :position) child-pos)
        (or (leap-locate (second ast) cursor-pos nil new-parents)
            (when bounds
              new-parents))))))


(defun complete-command (match-set partial start end)
  (dolist (pair leap.parser:*function-names/alist*)
    (when (and (<= (length partial) (length (car pair)))
               (string-equal partial (subseq (car pair) 0 (length partial))))
      (jupyter:match-set-add match-set (car pair) start end :type "command"))))


(defmethod jupyter:complete-code ((k kernel) match-set code cursor-pos)
   (if (or (not *leap-syntax*)
           (lisp-code-p code))
    (call-next-method)
    (let* ((ast (ignore-errors (leap-read code)))
           (fragment (unless (typep ast 'jupyter:result)
                       (leap-locate ast (1- cursor-pos)))))
      (when (and fragment
                 (eql :s-expr (caar fragment)))
        (let ((start (car (getf (car fragment) :bounds)))
              (end (cdr (getf (car fragment) :bounds))))
          (call-next-method k (jupyter:make-offset-match-set :parent match-set :offset start)
                            (subseq code start end) (- cursor-pos start))))
      (when (and fragment
                 (eql :s-expr (caar fragment))
                 (symbolp (getf (car fragment) :value))
                 (eql 0 (getf (car fragment) :position))
                 (eql :instruction (caadr fragment)))
        (complete-command match-set (symbol-name (getf (car fragment) :value))
              (car (getf (car fragment) :bounds))
              (cdr (getf (car fragment) :bounds))))
      (values))))


(defclass cando-installer (jupyter:installer)
  ((fork-client
     :reader fork-client
     :initarg :fork-client
     :initform nil)
   (image
     :reader image
     :initarg :image
     :initform nil))
  (:default-initargs
    :class 'kernel
    :language +language+
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


(defmethod jupyter:command-line ((instance cando-installer))
  "Get the command line for a cando installation."
  (let ((implementation (jupyter:installer-implementation instance)))
    (cons (or implementation
              (first (uiop:raw-command-line-arguments))
              (format nil "~(~A~)" (uiop:implementation-type)))
          (if (fork-client instance)
            (list "{connection_file}")
            (append (list "-f" "no-auto-lparallel")
                    (unless (image instance)
                      (list "--eval" "(ql:quickload :cando-jupyter)"))
                    (list "--eval" "(jupyter:run-kernel 'cando-jupyter:kernel)"
                          "--" "{connection_file}"))))))


(defun install (&key bin-path system local prefix root fork image)
  "Install Cando kernel.
- `bin-path` specifies path to LISP binary.
- `system` toggles system versus user installation.
- `local` toggles `/usr/local/share versus` `/usr/share` for system installations.
- `prefix` key specifies directory prefix for packaging.
- `root` key specifies the root under which the Jupyter folder is found. Is automatically determined if not provided.
- `fork` toggles a normal client versus a fork client."
  (jupyter:install
    (make-instance
      (if system
        'system-installer
        'user-installer)
      :implementation bin-path
      :display-name (format nil "~A ~:[~;(Fork)~]" +display-name+ fork)
      :kernel-name (format nil "~A~:[~;_fork~]" +language+ fork)
      :local local
      :prefix prefix
      :image image
      :fork-client fork
      :root root)))


;; cando doesn't create images on the fly yet.
#+(or)(defun install-image (&key prefix)
  "Install Cando kernel based on image.
- `prefix` key specifies directory prefix for packaging."
  (jupyter:install
    (make-instance 'user-image-installer
      :prefix prefix)))


(defun do-run-kernel-from-slime ()
  (jupyter:run-kernel 'cando-jupyter:kernel (clasp-posix:argv (1- (clasp-posix:argc)))))


(defun run-kernel-from-slime ()
  (bordeaux-threads:make-thread
   (lambda ()
     (loop for time from 1 below 100
           do (format t "Starting kernel for the ~a time~%" time)           
           do (jupyter:run-kernel 'cando-jupyter:kernel (clasp-posix:argv (1- (clasp-posix:argc))))))
   :name :jupyter))


(defun jupyterlab ()
  (run-kernel-from-slime))
