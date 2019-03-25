(in-package :leap)


(defun setup-amber-paths ()
  (leap.core::clear-path)
  (leap.core:add-path "amber:dat;leap;prep;")
  (leap.core:add-path "amber:dat;leap;lib;")
  (leap.core:add-path "amber:dat;leap;parm;")
  (leap.core:add-path "amber:dat;leap;cmd;")
  (leap.core:add-path "amber:dat;antechamber;"))

(defun setup-default-paths ()
  (setup-amber-paths)
  (leap.core:add-path *default-pathname-defaults*))

;;; ----------------------------------------------------------------------
;;;
;;; Quick and dirty setup-gaff
;;;

(defun setup-gaff ()
  (let ((*default-pathname-defaults*
          (translate-logical-pathname #P"source-dir:extensions;cando;src;data;force-field;"))
        (parms (chem:make-read-amber-parameters)))
    (with-open-file (fin "ATOMTYPE_GFF.DEF" :direction :input)
      (chem:read-types parms fin))
    (with-open-file (fin "gaff.dat" :direction :input)
      (chem:read-parameters parms fin (symbol-value (find-symbol "*AMBER-SYSTEM*" :leap))))
    (setf energy:*ff* (chem:get-force-field parms))
    energy:*ff*))

;;; ----------------------------------------------------------------------
;;;
;;; LEaP commands
;;;
(defvar *out* t)

(defun set-variable (entry)
  (let ((variable (first entry))
        (value (third entry)))
    (leap.core:register-variable variable value)))


(defun log-file (filename)
  (let* ((log-stream (open filename :direction :output :if-exists :supersede))
         (broadcast (make-broadcast-stream ext:+process-standard-output+ log-stream)))
    (setf cl:*standard-output* broadcast)))
  
(defun leap.log-file (entry)
  (valid-arguments entry 1)
  (let ((filename (pathname (ensure-string (second entry)))))
    (log-file filename)))

(defun desc (name)
  "Print a description of the object with name."
  (let ((val (leap.core:lookup-variable name)))
    (format *out* "~S~%" val)
    val))

(defun object (name)
  "Return the object with name."
  (leap.core:lookup-variable name))

(defun leap.desc (entry)
  (valid-arguments entry 1)
  (let ((var (second entry)))
    (format *out* "~a~%" var)))

(defun leap.add-pdb-res-map (entry)
  (valid-arguments entry 1)
  (let ((map (second entry)))
    (leap.pdb:add-pdb-res-map map)))

(defun add-atom-types (info)
  (warn "add-atom-types doesn't do anything"))

(defun leap.add-atom-types (entry)
  (valid-arguments entry 1)
  (let ((info (second entry)))
    (add-atom-types info)))


(defun load-off (filename)
  "Load an OFF file. OFF files that contain libraries of objects will create those objects in the environment."
  (leap.off:load-off filename))

(defun leap.load-off (entry)
  (valid-arguments entry 1)
  (let ((filename (ensure-path (second entry))))
    (load-off filename)))

(defun load-amber-params (filename &optional (force-field :default))
  (let* ((fn (ensure-path filename))
         (ff (let ((parmreader (chem:make-read-amber-parameters)))
              (with-open-file (fin fn :direction :input)
                (chem:read-parameters parmreader fin leap:*amber-system*)
                (chem:get-force-field parmreader)))))
    (chem:set-title ff filename)
    (format *out* "Adding force field ~a to ~a~%" ff force-field)
    (leap.core:add-force-field-or-modification ff force-field)))

(defun leap.load-amber-params (entry)
  (valid-arguments entry 1)
  (let ((filename (second entry)))
    (load-amber-params filename)))

(defun clear-force-field (&optional (ff :default))
  (leap.core:clear-force-field ff))

(defun list-force-fields ()
  (let ((ffs nil))
    (maphash (lambda (k v) (push (cons k v) ffs)) leap.core:*force-fields*)
    ffs))

(defun load-atom-type-rules (filename &optional (force-field :default))
  (let* ((path (leap.core:search-path filename))
         (fftypedb (with-open-file (fin path)
                       (leap.antechamber-type-definition-parser::read-antechamber-type-rules fin)))
         (ff (core:make-cxx-object 'chem:force-field)))
    (chem:set-title ff filename)
    (chem:set-type-db ff fftypedb)
    (leap.core:add-force-field-or-modification ff force-field)))

(defun assign-atom-types (matter &optional system)
  (chem:map-molecules
   nil
   (lambda (molecule)
     (let* ((force-field (chem:lookup-force-field-for-molecule molecule system))
            (type-rules (chem:get-types force-field)))
       (chem:assign-types type-rules molecule)))
   matter))

(defun source (filename)
  (let* ((path (let ((p (leap.core:search-path filename)))
                 (unless p (error "Could not find file ~a" filename))
                 p))
         (entire-file (alexandria:read-file-into-string path))
         (ast (architecture.builder-protocol:with-builder ('list)
                (esrap:parse 'leap.parser::leap entire-file))))
    (leap.core:evaluate 'list ast leap.core:*leap-env*))
  t)

(defun leap.source (entry)
  (valid-arguments entry 1)
  (let* ((filename (ensure-path (second entry))))
    (source filename)))

(defun save-amber-parm (aggregate topology-file-name &optional crd-pathname (force-field-name :default))
  (let* ((top-pathname (merge-pathnames (pathname topology-file-name)))
         (crd-pathname (if crd-pathname
                           (merge-pathnames crd-pathname)
                           (make-pathname :type "crd" :defaults top-pathname))))
    (leap.topology:save-amber-parm-format aggregate top-pathname crd-pathname :force-field (leap.core:merged-force-field force-field-name) :assign-types t)))

(defun solvate-box (solute solvent buffer &rest iso-closeness)
  (when (numberp buffer)
    (setf buffer (list buffer buffer buffer)))
  ;;; Process the leap syntax for the two optional parameters "iso" and a closeness value
  (let ((iso nil)
        (closeness 1.0))
    (flet ((error-iso-closeness ()
             (error "You must provide either :iso or a number for closeness - you provided ~s" iso-closeness)))
      (cond
        ((= (length iso-closeness) 0))
        ((= (length iso-closeness) 1)
         (cond
           ((eq (first iso-closeness) :iso)
            (setf iso t))
           ((and (stringp (first iso-closeness))
                 (string-equal (first iso-closeness) "iso"))
            (setf iso t))
           ((numberp (first iso-closeness))
            (setf closeness (first iso-closeness)))
           (t (error-iso-closeness))))
        ((and (= (length iso-closeness) 2)
              (or (find :iso iso-closeness)
                  (find "iso" iso-closeness :test #'string-equal))
              (find-if #'numberp iso-closeness))
         (setf iso t
               closeness (find-if #'numberp iso-closeness)))
        (t (error-iso-closeness))))
    (leap.solvate:tool-solvate-and-shell solute solvent buffer :closeness closeness :isotropic iso)))

(defun solvate-oct (solute solvent buffer &rest iso-closeness)
  (when (numberp buffer)
    (setf buffer (list buffer buffer buffer)))
  (let ((iso nil)
        (closeness 1.0))
    (flet ((error-iso-closeness ()
             (error "You must provide either :iso or a number for closeness - you provided ~s" iso-closeness)))
      (cond
        ((= (length iso-closeness) 0))
        ((= (length iso-closeness) 1)
         (cond
           ((eq (first iso-closeness) :iso)
            (setf iso t))
           ((string-equal (first iso-closeness) "iso")
            (setf iso t))
           ((numberp (first iso-closeness))
            (setf closeness (first iso-closeness)))
           (t (error-iso-closeness))))
        ((and (= (length iso-closeness) 2)
              (or (find :iso iso-closeness)
                  (find "iso" iso-closeness :test #'string-equal))
              (find-if #'numberp iso-closeness))
         (setf iso t
               closeness (find-if #'numberp iso-closeness)))
        (t (error-iso-closeness))))
    (leap.solvate:tool-solvate-and-shell solute solvent buffer :closeness closeness :oct t :isotropic iso)))

(defun solvate-shell (solute solvent farness &rest iso-closeness)
  (let ((buffer (list farness farness farness))
        (iso nil)
        (closeness 1.0))
    (flet ((error-iso-closeness ()
             (error "You must provide either :iso or a number for closeness - you provided ~s" iso-closeness)))
      (cond
        ((= (length iso-closeness) 0))
        ((= (length iso-closeness) 1)
         (cond
           ((eq (first iso-closeness) :iso)
            (setf iso t))
           ((string-equal (first iso-closeness) "iso")
            (setf iso t))
           ((numberp (first iso-closeness))
            (setf closeness (first iso-closeness)))
           (t (error-iso-closeness))))
        ((and (= (length iso-closeness) 2)
              (or (find :iso iso-closeness)
                  (find "iso" iso-closeness :test #'string-equal))
              (find-if #'numberp iso-closeness))
         (setf iso t
               closeness (find-if #'numberp iso-closeness)))
        (t (error-iso-closeness))))
    (leap.solvate:tool-solvate-and-shell solute solvent buffer :closeness closeness :farness farness :shell t :isotropic iso)))


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


(defun leap-quit ()
  (throw 'repl-done nil))


(defvar *commands*)
(defun leap-help (&optional arg)
  (if arg
      (format t "There is no specific help for commands yet~%")
      (progn
        (format t "Help will be available on the following subjects: ~%")
        (loop for cmd in *commands*
              for col from 1
              do (format t "~20a" (car cmd))
              when (and (> col 0) (= (rem col 4) 0))
                do (terpri)
              finally (when (/= (rem col 4) 0)
                        (terpri))))))

(eval-when (:load-toplevel :execute)
  (setf *commands*
    '(("logFile" . log-file)
      ("desc" . desc)
      ("loadOff" . load-off)
      ("loadMol2" . load-mol2)
      ("source" . source)
      ("loadAmberParams" . load-amber-params)
      ("addPdbResMap" . leap.pdb:add-pdb-res-map)
      ("addPdbAtomMap" . leap.pdb:add-pdb-atom-map)
      ("addAtomTypes" . add-atom-types)
      ("saveAmberParms" . save-amber-parm)
      ("solvateBox" . solvate-box)
      ("solvateOct" . solvate-oct)
      ("solvateShell" . solvate-shell)
      ("addIons" . leap.add-ions:add-ions)
      ("setBox" . leap.set-box:set-box)
      ("help" . leap-help)
      ("quit" . leap-quit)
      ))
  (dolist (command *commands*)
    (if (fboundp (cdr command))
        (setf (leap.core:function-lookup (car command) leap.core:*leap-env*) (cdr command))
        (error "~a is not a function" (cdr command)))))

  

(defun parse-leap-code (code)
  (let ((ast (architecture.builder-protocol:with-builder ('list)
               (esrap:parse 'leap.parser::leap code))))
    ast))

(defun process-command-line-options ()
  (let ((includes (core:leap-command-line-includes))
        (scripts (core:leap-command-line-scripts)))
    (loop for include-path in includes
          do (format t "Adding include path: ~s~%" include-path)
          do (leap.core:add-path include-path))
    (loop for script in scripts
          do (format t "Sourcing script: ~s~%" script)
          do (source script))
    ))
    
(defun leap-repl ()
  (process-command-line-options)
  (format t "Welcome to Cando-LEaP!~%")
  (clear-input)
  (catch 'repl-done
    (loop for x below 10
          for code = (progn
                       (format t "> ") (finish-output)
                       (read-line))
          do (handler-case
                 (let ((ast (architecture.builder-protocol:with-builder
                                ('list)
                              (handler-bind ((esrap:esrap-parse-error
                                               (lambda (c)
                                                 (format t "Encountered error ~s while parsing ~s~%" c code)
                                                 (break "Encountered error ~s while parsing ~s" c code))))
                                (esrap:parse 'leap.parser::leap code)))))
                   (core:call-with-stack-top-hint
                    (lambda ()
                      (leap.core:evaluate 'list ast leap.core:*leap-env*))))
               (error (var)
                 (format t "Error evaluating ~a~%" code))))))

(defun leap-repl-then-exit ()
  (leap-repl)
  (core:exit 0))



;;(defun solvate-box (solute solvent width-list 
;;
;;)
;;; ----------------------------------------------------------------------
;;;
;;; Implement a simple leap script parser
;;;
;;; This uses the Common Lisp reader and it would probably be better
;;; to develop a proper parser
;;;
#++
(progn
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
           (format *out* "cmd cmd-assoc cmd-func: ~a ~a ~a~%" cmd cmd-assoc cmd-func)
           (if (null cmd-func)
               (error "Illegal function ~a cmd-assoc: ~a  cmd-func: ~a" cmd cmd-assoc cmd-func)
               (progn
                 (funcall cmd-func entry)))))))

  (defun leap-interpreter (fin &optional eof-error-p eof)
    "Load the file of leap commands and execute them one by one.
Nothing is returned."
    (loop for entry = (leap-parse-entry fin eof-error-p eof)
       do (format *out* "entry = ~a~%" entry)
       unless (or (null entry) (eq entry eof))
       do (interpret-leap-entry entry)
       until (eq entry eof)))

  )
