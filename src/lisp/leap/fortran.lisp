(in-package :fortran)

(defconstant +fortran-debug-comment-char+ #\=)

(defstruct (fortran-output-file (:conc-name fof-))
  stream format per-line number-on-line wrote-nothing debug)

(defvar *fortran-output-file* nil)

(defmacro with-fortran-output-file ((ff &rest open-args) &body body)
  `(let* ((*fortran-output-file* (make-fortran-output-file
                                  :stream (apply #'cl:open (list ,@open-args))
                                  :format nil
                                  :per-line nil
                                  :number-on-line nil
                                  :wrote-nothing nil))
          (,ff *fortran-output-file*))
     (unwind-protect
          (progn ,@body)
       (close (fof-stream *fortran-output-file*)))))

(defun format (per-line fmt &optional (ff *fortran-output-file*))
  "Define the format for write commands that will follow"
  (setf (fof-format ff) fmt
        (fof-per-line ff) per-line
        (fof-number-on-line ff) 0
        (fof-wrote-nothing ff) t))

(defun write (val &optional (ff *fortran-output-file*))
  "Write _val_ using the current format"
  (core:bformat (fof-stream ff) (fof-format ff) val)
  (incf (fof-number-on-line ff))
  (setf (fof-wrote-nothing ff) nil)
  (when (> (fof-number-on-line ff) (fof-per-line ff))
    (core:bformat (fof-stream ff) "\n")
    (setf (fof-number-on-line ff) 0)))

(defun end-line (&optional (ff *fortran-output-file*))
  (when (or (fof-wrote-nothing ff) (/= (fof-number-on-line ff) 0))
    (core:bformat (fof-stream ff) "\n"))
  (setf (fof-wrote-nothing ff) t
        (fof-number-on-line ff) 0))

(defun debug-on (&optional (ff *fortran-output-file*))
  (setf (fof-debug ff) t))

(defun debug-off (&optional (ff *fortran-output-file*))
  (setf (fof-debug ff) nil))

(defun debug (msg &optional (ff *fortran-output-file*))
  (when (fof-debug ff)
    (core:bformat file "%c%s\n" +fortran-debug-comment-char+ msg)))


;;; ------------------------------------------------------------
;;;
;;; Input routines
;;;

(defvar *fortran-input-file* nil)

(defstruct fortran-input-file
  file need-input)

(defmacro with-fortran-input-file ((ff &rest open-args) &body body)
  `(let* ((*fortran-input-file* (make-fortran-input-file
                           :file (apply #'open open-args)
                           :need-input t))
          (,ff *fortran-input-file*))
     ,@body))

(defun get-input-line (&optional (ff *fortran-input-file*))
  (loop for line = (read-line (fortran-input-file-stream ff) nil :eof)
     until (eq line :eof)
     unless (char= (elt line 0) +fortran-debug-comment-char+)
     do (progn
          (setf (fortran-input-file-input ff) line)
          (return-from get-input-line )))
  :eof)

#|
(defun read-string (&optional (ff *fortran-input-file*))
  (let ((got-one nil))
    (loop until got-one
       do (when (fortran-input-file-need-input ff)
            (when (eq (get-input-line ff) :eof)
              (return-from read-string :eof)))
       do (progn
            (setf (fortran-input-file-needs-input ff) nil)
            (setf (fortran-input-file-input ff)
                  (string-left-trim '(#\space #\tab) (fortran-input-file-input ff)))
            (setf str (remove-first-string
|#
  
