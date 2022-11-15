
(in-package :work)

(defclass job ()
  ((scripts :initform nil :initarg :scripts :accessor scripts)
   (inputs :initform nil :initarg :inputs :accessor inputs)
   (outputs :initform nil :initarg :outputs :accessor outputs)
   (parameters :initform nil :initarg :parameters :accessor parameters)
   (makefile-clause :initarg :makefile-clause :accessor makefile-clause)))


(defclass data ()
  ((definers :initarg :definers :initform nil :accessor definers)
   (users :initarg :users :initform nil :accessor users)
   (pathname :initarg :pathname :accessor pathname)))


(defclass job-graph ()
  ((jobs :initarg :jobs :accessor jobs)))


(defclass script (data)
  ())

(defun connect-graph (job)
  (let ((scripts (scripts job))
        (inputs (inputs job))
        (outputs (outputs job)))
    (loop for script in scripts
          do (pushnew job (users script)))
    (loop for input in inputs
          for node-input = (node input)
          if (consp node-input)
            do (mapc (lambda (one) (push job (users one))) node-input)
          else
            do (push job (users node-input)))
    (loop for output in outputs
          do (push job (definers (node output))))
    job))

