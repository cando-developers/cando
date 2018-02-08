(in-package :fortran)

(defconstant +fortran-debug-comment-char+ #\=)

(defstruct (fortran-output-file (:conc-name fof-))
  stream format per-line number-on-line wrote-nothing
  (debug nil))

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

(defun fformat (per-line fmt &optional (ff *fortran-output-file*))
  "Define the format for write commands that will follow"
  (setf (fof-format ff) fmt
        (fof-per-line ff) per-line
        (fof-number-on-line ff) 0
        (fof-wrote-nothing ff) t))

(defun fwrite (val &optional (ff *fortran-output-file*))
  "Write _val_ using the current format"
  (core:bformat (fof-stream ff) (fof-format ff) val)
  (incf (fof-number-on-line ff))
  (setf (fof-wrote-nothing ff) nil)
  (when (>= (fof-number-on-line ff) (fof-per-line ff))
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
    (core:bformat (fof-stream ff) "%c%s\n" +fortran-debug-comment-char+ msg)))


;;; ------------------------------------------------------------
;;;
;;; Input routines
;;;


(defun format-symbol (format-char)
  (cond
    ((char-equal format-char #\i) :integer)
    ((or (char-equal format-char #\e)
         (char-equal format-char #\f)
         (char-equal format-char #\d)) :double-float)
    ((char-equal format-char #\a) :string)
    (t (error "Handle format-char of ~s~%" format-char))))
  
(defun parse-fortran-format (fortranformat &key (start 0) end)
  (let* ((format-start (position-if #'alpha-char-p fortranformat :start start))
         (format-end (1+ format-start))
         (per-line (parse-integer fortranformat :start start :end format-start))
         (format-char (aref fortranformat format-start))
         (has-period (position-if (lambda (c) (char= #\. c)) fortranformat :start format-end))
         width decimal)
    (if has-period
        (let ((width (parse-integer fortranformat :start format-end :end has-period))
              (decimal (parse-integer fortranformat :start (1+ has-period) :end end)))
          (values per-line format-char width decimal))
        (let ((width (parse-integer fortranformat :start format-end :end end)))
          (values per-line format-char width nil)))))

          
(defconstant %format "%format")
(defun parse-fortran-format-line (line)
  (if (string-equal %format line :end2 (length %format))
      (let ((format-end (position #\) line :start 7)))
        (parse-fortran-format line :start 8 :end format-end))
      (error "Not a %FORMAT line -> ~s" line)))


(defstruct fortran-input-file
  stream look-ahead)

(defun fread-line (fif &optional eof-error-p eof-value)
  (let ((data (read-line (fortran-input-file-stream fif) eof-error-p eof-value)))
    (when (string= "=-" data :start2 0 :end2 2)
      (setf data (read-line (fortran-input-file-stream fif) eof-error-p eof-value)))
    (setf (fortran-input-file-look-ahead fif) data)))

(defun fread-line-or-error (fif)
  (let ((data (fread-line fif nil nil)))
    (if data
        data
        (error "Hit the end of file"))))


(defun parse-integer-line (line result width)
  (loop for start = 0 then end
     for end = (+ start width)
     until (< (length line) end)
     do (let* ((num-pos (position-if (lambda (c) (char/= c #\space)) line :start start :end end))
               (num (parse-integer line :start num-pos :end end)))
          (vector-push-extend num result))))

(defun fread-integer-vector (fif per-line width &optional eof-error-p eof-value)
  (let ((result (make-array 32 :element-type '(signed-byte 32) :fill-pointer 0 :adjustable t)))
    (loop for line = (fortran-input-file-look-ahead fif)
       until (or (eq line eof-value) (eql (aref line 0) #\%))
       do (parse-integer-line line result width)
       do (fread-line fif eof-error-p eof-value))
    result))

(defun parse-string-line (line result width)
  (loop for start = 0 then end
     for end = (+ start width)
     until (< (length line) end)
     do (let* ((label-pos (position-if (lambda (c) (char/= c #\space)) line :start start :end end))
               (label (subseq line label-pos end)))
          (vector-push-extend label result))))

(defun fread-string-vector (fif per-line width &optional eof-error-p eof-value)
  (let ((result (make-array 32 :element-type t :fill-pointer 0 :adjustable t)))
    (loop for line = (fortran-input-file-look-ahead fif)
       until (or (eq line eof-value) (eql (aref line 0) #\%))
       do (parse-string-line line result width)
       do (fread-line fif eof-error-p eof-value))
    result))

;;;parse a double-float in the string starting at start end ending at end.
;;;  Allow whitespace at the start end end of the string
(defun parse-double-float (string &key (start 0) end)
  (let ((float-start (position-if (lambda (c) (and (char/= c #\space) (char/= c #\tab))) string :start start :end end))
        (float-end (or end (length string)))
        (has-period (position #\. string :start start :end end))
        (has-exponent-char (position-if (lambda (c) (or (char-equal #\E c) (char-equal #\D c))) string :start start :end end)))
    (let ((significand-end (or has-exponent-char float-end))
          sig-x
          (sig-y "0000000000")
          significand
          (exponent 0))
      (if has-period
          (setf sig-x (subseq string float-start has-period)
                Sig-y (concatenate 'string (subseq string (1+ has-period) significand-end) sig-y))
          (setf sig-x (subseq string float-start significand-end)))
      (if has-exponent-char
          (setf exponent (parse-integer string :start (1+ has-exponent-char) :end float-end)))
      (setf exponent (- exponent (length sig-y)))
      (let ((significand-str (concatenate 'string sig-x sig-y)))
        (* (float (parse-integer significand-str) 1.0D0) (expt 10.0d0 exponent))))))

(defun parse-double-float-line (line result width)
  (loop for start = 0 then end
     for end = (+ start width)
     until (< (length line) end)
     do (let* ((val-pos (position-if (lambda (c) (char/= c #\space)) line :start start :end end))
               (val (parse-double-float line :start start :end end)))
          (vector-push-extend val result))))

(defun fread-double-float-vector (fif per-line width &optional eof-error-p eof-value)
  (let ((result (make-array 32 :element-type 'double-float :fill-pointer 0 :adjustable t)))
    (loop for line = (fortran-input-file-look-ahead fif)
       until (or (eq line eof-value) (eql (aref line 0) #\%))
       do (parse-double-float-line line result width)
       do (fread-line fif eof-error-p eof-value))
    result))
    
(defun fread-vector (fif per-line format-char width)
  (ecase (format-symbol format-char)
    (:integer (fread-integer-vector fif per-line width))
    (:double-float (fread-double-float-vector fif per-line width))
    (:string (fread-string-vector fif per-line width))))

