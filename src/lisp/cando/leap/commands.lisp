(in-package :leap.commands)

;;; ----------------------------------------------------------------------
;;;
;;; LEaP commands
;;;
    
(defun log-file (filename)
  (let* ((log-stream (open filename :direction :output :if-exists :supersede))
         (broadcast (make-broadcast-stream ext:+process-standard-output+ log-stream)))
    (setf cl:*standard-output* broadcast)))

(defun add-atom-types (info)
  (warn "add-atom-types doesn't do anything"))

(defun load-off (filename)
  (print "In commands:load-off")
  (leap.off:load-off filename))

(defun load-amber-params (filename)
  "Load an AMBER parameter file from filename and return a chem:force-field object."
  (let* ((filename (merge-pathname filename)))
    (let ((parmreader (chem:make-read-amber-parameters)))
      (with-open-file (fin filename :direction :input)
        (chem:read-parameters parms fin)
        (chem:get-force-field parms)))))

(defun source (filename)
  "Load the file of leap commands and execute them one by one.
Nothing is returned."
  (let ((path (leap.core:search-path filename)))
    (with-open-file (stream path :direction :input)
      (let* ((entire-file (make-string (+ (file-length stream) 2)
                                       :initial-element #\newline)))
        (read-sequence entire-file stream)
        (with-input-from-string (sin entire-file)
          (loop for entry = (leap-parse-entry sin nil :eof)
             unless (or (null entry) (eq entry :eof))
             do (interpret-leap-entry entry)
             until (eq entry :eof)))))))

(defun save-amber-parm (aggregate topology-file-name &optional coordinate-file-name (force-field-name :default))
  (let* ((top-pathname (merge-pathnames (pathname topology-file-name)))
         (crd-pathname (if crd-pathname
                           (merge-pathnames crd-pathname)
                           (make-pathname :type "crd" :defaults top-pathname))))
    (leap.topology:save-amber-parm aggregate top-pathname crd-pathname (leap.core:merge-force-fields force-field-name))))

;;; ------------------------------------------------------------
;;;
;;; Commands that work with the leap parser
;;;

(defun leap.log-file (entry)
  (valid-arguments entry 1)
  (let ((filename (pathname (ensure-string (second entry)))))
    (log-file filename)))
  
(defun leap.desc (entry)
  (valid-arguments entry 1)
  (let ((var (second entry)))
    (format t "~a~%" var)))

(defun leap.add-pdb-res-map (entry)
  (valid-arguments entry 1)
  (let ((map (second entry)))
    (leap.pdb:add-pdb-res-map map)))

(defun leap.add-atom-types (entry)
  (valid-arguments entry 1)
  (let ((info (second entry)))
    (add-atom-types info)))

(defun leap.load-off (entry)
  (valid-arguments entry 1)
  (let ((filename (ensure-path (second entry))))
    (load-off filename)))

(defun leap.load-amber-params (entry)
  (valid-arguments entry 1)
  (let* ((filename (ensure-path (second entry))))
    (load-amber-params filename)))

(defun leap.source (entry)
  (valid-arguments entry 1)
  (let* ((filename (ensure-path (second entry))))
    (source filename)))
  
(defun ensure-string (obj)
  (cond
    ((stringp obj) obj)
    ((symbolp obj)
     (let ((val (leap.core:lookup-variable obj nil :not-found)))
       (if (eq val :not-found)
           (string val)
           val)))))

(defun ensure-path (obj)
  "Convert obj into a pathname and search the amber path list for the file and return it"
  (declare (type (or symbol pathname string) obj))
  (let ((pn (cond
              ((symbolp obj) (pathname (string obj)))
              ((pathname obj) obj)
              ((stringp obj) (pathname obj))
              (t (pathname (string obj))))))
    (leap.core:search-path pn)))

(defun valid-arguments (entry num)
  (unless (= (length entry) (1+ num))
    (error "Bad arguments for ~a" (car entry))))

(defparameter *commands*
  '((:logfile . leap.log-file)
    (:desc . leap.desc)
    (:loadoff . leap.load-off)
    (:source . leap.source)
    (:loadamberparams . leap.load-amber-params)
    (:addpdbresmap . leap.add-pdb-res-map)
    (:addatomtypes . leap.add-atom-types)))


;;; ----------------------------------------------------------------------
;;;
;;; Implement a simple leap script parser
;;;
;;; This uses the Common Lisp reader and it would probably be better
;;; to develop a proper parser
;;;
(defun leap-list-reader (stream char)
  "Read a Leap list delimited by #\} and return it."
  (let ((*readtable* (copy-readtable)))
    (set-syntax-from-char #\newline #\space)
    (read-delimited-list #\} stream t)))

(defun leap-parse-entry (stream &optional (eof-error-p t) eof)
  "Parse a leap command from the stream and return it.
A leap command can have the form:
[command] [arg1] [arg2] [arg3] ... #\NEWLINE  - execute command
[variable] = [expression] #\NEWLINE           - set variable
[variable] = #\NEWLINE                        - clear variable
# ... #\NEWLINE                               - comment
An [expression] is either a command or a number or variable."
  (if (eq (peek-char nil stream eof-error-p eof) eof)
      eof
      (let ((*readtable* (copy-readtable))
            (*package* (find-package :keyword)))
        (setf (readtable-case *readtable*) :preserve)
        (set-macro-character #\{ #'leap-list-reader nil)
        (set-syntax-from-char #\# #\;)
        (set-syntax-from-char #\newline #\))
        (read-delimited-list #\newline stream nil))))

(defun interpret-leap-entry (entry)
  "Interpret a leap command or variable assignment. 
Nothing is returned, it's all side effects."
  (cond
    ((eq (second entry) :=)
     (set-variable entry))
    (t (let* ((cmd (intern (string-upcase (string (first entry))) :keyword))
              (cmd-assoc (assoc cmd *commands*))
              (cmd-func (cdr cmd-assoc)))
         (if (null cmd-func)
             (error "Illegal function ~a cmd-assoc: ~a  cmd-func: ~a" cmd cmd-assoc cmd-func)
             (progn
               (funcall cmd-func entry)))))))

(defun leap-interpreter (fin &optional eof-error-p eof)
  "Load the file of leap commands and execute them one by one.
Nothing is returned."
  (loop for entry = (leap-parse-entry fin eof-error-p eof)
     unless (or (null entry) (eq entry eof))
     do (interpret-leap-entry entry)
     until (eq entry eof)))

