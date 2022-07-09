(in-package :cando-jupyter)

(defvar +display-name+ "Cando")
(defvar +language+ "cando")

(defparameter *leap-syntax* t)


(defun leap-syntax-enable (on)
  (setq *leap-syntax* on))


(defclass kernel (jupyter/common-lisp:kernel)
  ()
  (:default-initargs
    :package (find-package :cando-user)
    :banner "cando-jupyter: a Cando Jupyter kernel
(C) 2020-2021 Christian Schafmeister (LGPL-3.0)"))

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
  (leap.core:evaluate 'list ast leap.core:*leap-env*))


(defmethod jupyter:evaluate-code ((k kernel) code &optional source-path breakpoints)
  (if (or (not *leap-syntax*)
          (lisp-code-p code))
    (call-next-method)
    (jupyter:handling-errors
      (dolist (expr (cadadr (leap-read code)) (values))
        (unless (eq :comment (caar expr))
          (jupyter:execute-result (leap-eval (list :leap (list :instruction (list expr))))))))))


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
           (fragment (leap-locate ast (1- cursor-pos))))
      (when (and fragment
                 (eql :s-expr (caar fragment)))
        (let ((start (car (getf (car fragment) :bounds)))
              (end (cdr (getf (car fragment) :bounds))))
          (call-next-method k
                            (j:make-substring-match-set :parent match-set :start start :end end)
                            (concatenate 'string
                                         (substitute-if #\space #'graphic-char-p (subseq code 0 start))
                                         (subseq code start end)
                                         (substitute-if #\space #'graphic-char-p (subseq code end)))
                            cursor-pos)))
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
   (load-system :accessor installer-load-system
                :initarg :load-system
                :type (or null string)))
  (:default-initargs
    :class 'kernel
    :language +language+
    :debugger t
    :systems '(:cando-jupyter)))


(defclass system-installer (jupyter:system-installer cando-installer)
  ()
  (:documentation "cando system installer."))


(defclass user-installer (jupyter:user-installer cando-installer)
  ()
  (:documentation "cando user installer."))


;; cando doesn't create images on the fly yet.
#+(or)(defclass user-image-installer (jupyter:user-image-installer cando-installer)
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
            (append (when (installer-load-system instance)
                      (list "--eval" (installer-load-system instance)))
                    (list "--eval" "(jupyter:run-kernel 'cando-jupyter:kernel)"
                          "--" "{connection_file}"))))))


(defun install (&key bin-path system local prefix jupyter program fork implementation (load-system t))
  "Install Cando kernel.
- `bin-path` specifies path to LISP binary.
- `system` toggles system versus user installation.
- `local` toggles `/usr/local/share versus` `/usr/share` for system installations.
- `prefix` key specifies directory prefix for packaging.
- `jupyter` key specifies the root under which the Jupyter folder is found. Is automatically determined if not provided.
- `program` key specifies the root under which the program folder is found. Is automatically determined if not provided.
- `fork` toggles a normal client versus a fork client."
  (jupyter:install
    (make-instance
      (if system
        'system-installer
        'user-installer)
      :implementation bin-path
      :display-name (format nil "~A~:[~; (Fork)~]~@[ (~A)~]"
                                +display-name+ fork implementation)
      :kernel-name (format nil "~A~:[~;_fork~]~@[_~(~A~)~]"
                               +language+ fork implementation)
      :local local
      :prefix prefix
      :load-system (cond ((or (null load-system)
                              (stringp load-system))
                          load-system)
                         ((find-package :quicklisp)
                          "(ql:quickload :cando-jupyter)")
                         (t
                          "(asdf:load-system :cando-jupyter)"))
      :fork-client fork
      :jupyter-path jupyter
      :program-path program)))


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

(defmethod jupyter:start :after ((kernel kernel))
  (setf cando:*ui-client* (make-instance 'jupyter-ui-client))
  (bordeaux-threads:make-thread (lambda ()
                                  (sleep 600)
                                  (with-output-to-string (*standard-output*)
                                    (let ((*error-output* *standard-output*))
                                      (asdf:load-system :cando-user-install)
                                      (uiop:symbol-call "CANDO-USER-INSTALL" "UPDATE"))))))