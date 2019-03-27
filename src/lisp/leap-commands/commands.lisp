(in-package :leap.commands)




(defun setup-amber-paths ()
  (leap.core::clear-path)
  (leap.core:add-path "amber:dat;leap;prep;")
  (leap.core:add-path "amber:dat;leap;lib;")
  (leap.core:add-path "amber:dat;leap;parm;")
  (leap.core:add-path "amber:dat;leap;cmd;")
  (leap.core:add-path "amber:dat;antechamber;"))

(defun setup-default-paths ()
  (setup-amber-paths)
  (leap.core:add-path *default-pathname-defaults*)
  (show-paths))

(defun show-paths ()
  (format t "Paths: ~a~%" leap.core::*path*))

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
  "    desc variable
      object                       _variable_

Print a description of the object.
"
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
  "    addPdbAtomMap list
      LIST                         _list_

The atom Name Map is used to try to map atom names read from PDB files
to atoms within residue UNITs when the atom name in the PDB file does
not match an atom in the residue.  This enables PDB files to be read
in without extensive editing of atom names.  The LIST is a LIST of LISTs:
      { {sublist} {sublist} ... }
where each sublist is of the form
      { \"OddAtomName\" \"LibAtomName\" }
Many `odd' atom names can map to one `standard' atom name, but any single
odd atom name maps only to the last standard atom name it was mapped to.
"
  (valid-arguments entry 1)
  (let ((map (second entry)))
    (leap.pdb:add-pdb-res-map map)))

(defun add-atom-types (info)
  (warn "add-atom-types doesn't do anything"))

(defun leap.add-atom-types (entry)
  (valid-arguments entry 1)
  (let ((info (second entry)))
    (add-atom-types info)))

(defun create-atom (name &optional type (charge 0.0))
  "    variable = createAtom name type charge

      ATOM                         _variable_
      STRING                       _name_
      STRING                       _type_
      NUMBER                       _charge_

Return a new ATOM with _name_, _type_, and _charge_.
"
  (let* ((name-sym (intern (string name) :keyword))
         (element (chem:element-from-atom-name-string (string name-sym)))
         (atom (chem:make-atom name-sym element)))
    (chem:set-type atom type)
    (chem:set-charge atom charge)
    atom))

(defun load-off (filename)
"    loadOff filename
      STRING                       _filename_

This command loads the Object File Format library within the file named
_filename_.  All UNITs and PARMSETs within the library will be loaded.
The objects are loaded into LEaP under the variable names the objects
had when they were saved.  Variables already in existence that have the
same names as the objects being loaded will be overwritten.  PARMSETs
loaded using this command are included in LEaP's library of PARMSETs
that is searched whenever parameters are required.
"
  (leap.off:load-off filename))

(defun list-force-fields ()
  (let ((ffs nil))
    (maphash (lambda (k v) (push (cons k v) ffs)) leap.core:*force-fields*)
    ffs))

(defun assign-atom-types (matter &optional system)
  (chem:map-molecules
   nil
   (lambda (molecule)
     (let* ((force-field (chem:lookup-force-field-for-molecule molecule system))
            (type-rules (chem:get-types force-field)))
       (chem:assign-types type-rules molecule)))
   matter))

(defun source (filename)
"    source filename
      STRING                       _filename_

This command executes LEaP commands within a text file.  To display the
commands as they are read, see the verbosity command.  The text within
the source file must be formatted exactly like the text the user types
into LEaP.
"
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
  (let* ((filename (leap.core:ensure-path (second entry))))
    (source filename)))

(defun save-amber-parm (aggregate topology-file-name &optional crd-pathname (force-field-name :default))
  (let* ((top-pathname (merge-pathnames (pathname topology-file-name)))
         (crd-pathname (if crd-pathname
                           (merge-pathnames crd-pathname)
                           (make-pathname :type "crd" :defaults top-pathname))))
    (leap.topology:save-amber-parm-format aggregate top-pathname crd-pathname :force-field (leap.core:merged-force-field force-field-name) :assign-types t)))

(defun solvate-box (solute solvent buffer &rest iso-closeness)
  "    solvateBox solute solvent buffer [ \"iso\" ] [ closeness ]

      UNIT                         _solute_
      UNIT                         _solvent_
      object                       _buffer_
      NUMBER                       _closeness_

The solvateBox command creates a solvent box around the _solute_ UNIT.
The _solute_ UNIT is modified by the addition of _solvent_ RESIDUEs.

The user may want to first align long solutes that are not expected
to tumble using alignAxes, in order to minimize box volume.

The normal choice for a TIP3 _solvent_ UNIT is WATBOX216. Note that
constant pressure equilibration is required to bring the artificial box
to reasonable density, since Van der Waals voids remain due to the
impossibility of natural packing of solvent around the solute and at
the edges of the box.

The solvent box UNIT is copied and repeated in all three spatial directions
to create a box containing the entire solute and a buffer zone defined
by the _buffer_ argument. The _buffer_ argument defines the distance,
in angstroms, between the wall of the box and the closest ATOM in the
solute.

If the buffer argument is a single NUMBER, then the buffer distance is
the same for the x, y, and z directions, unless the \"iso\" option is used
to make the box isometric, with the shortest box clearance = buffer. If
\"iso\" is used, the solute is rotated to orient the principal axes,
otherwise it is just centered on the origin.

If the buffer argument is a LIST of three NUMBERS, then the NUMBERs are
applied to the x, y, and z axes respectively. As the larger box is created
and superimposed on the solute, solvent molecules overlapping the solute
are removed.

The optional _closeness_ parameter can be used to control the extent to
which _solvent_ ATOMs overlap _solute_ ATOMs.  The default value of
the _closeness_ argument is 1.0, which allows no overlap.  Smaller
values allow solvent ATOMs to overlap _solute_ ATOMs by (1 - closeness) *
R*ij, where R*ij is the sum of the Van der Waals radii of solute and
solvent atoms.  Values greater than 1 enforce a minimum gap between
solvent and solute of (closeness - 1) * R*ij.

This command modifies the _solute_ UNIT in several ways.  First, the
coordinates of the ATOMs are modified to move the center of a box
enclosing the Van der Waals radii of the atoms to the origin.  Secondly,
the UNIT is modified by the addition of _solvent_ RESIDUEs copied from
the _solvent_ UNIT. Finally, the box parameter of the new system (still
named for the _solute_) is modified to reflect the fact that a periodic,
rectilinear solvent box has been created around it.
"
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

(defun valid-arguments (entry num)
  (unless (= (length entry) (1+ num))
    (error "Bad arguments for ~a" (car entry))))


(defun leap-quit ()
  (throw 'repl-done nil))


(defun leap-start-swank ()
  (funcall (find-symbol "START-SWANK" :cando-user)))

(defun leap-help (&optional arg)
  "    help [string]

      STRING                       _string_

This command prints a description of the command in _string_.  If
the STRING is not given then a list of legal STRINGs is provided.
"
  (if arg
      (progn
        (let ((found (assoc arg leap.parser:*function-names/alist* :test #'string=)))
          (if found
              (let ((cmd (cdr found)))
                (format t "~%~a~%" (documentation cmd 'function)))
              (format t "No help available for ~a~%" arg))))
      (progn
        (format t "Help is available on the following subjects: ~%")
        (let* ((commands (copy-seq leap.parser:*function-names/alist*))
               (sorted-commands (sort commands #'string< :key #'car)))
          (loop for cmd in sorted-commands
                for col from 1
                do (format t "~20a" (car cmd))
                when (and (> col 0) (= (rem col 4) 0))
                  do (terpri)
                finally (when (/= (rem col 4) 0)
                          (terpri)))))))

(eval-when (:load-toplevel :execute)
  (setf leap.parser:*function-names/alist*
    '(("logFile" . log-file)
      ("desc" . desc)
      ("loadOff" . load-off)
      ("loadMol2" . cando:load-mol2)
      ("source" . source)
      ("loadAmberParams" . leap:load-amber-params)
      ("addPdbResMap" . leap.pdb:add-pdb-res-map)
      ("addPdbAtomMap" . leap.pdb:add-pdb-atom-map)
      ("addAtomTypes" . add-atom-types)
      ("saveAmberParms" . save-amber-parm)
      ("solvateBox" . solvate-box)
      ("solvateOct" . solvate-oct)
      ("solvateShell" . solvate-shell)
      ("setupDefaultPaths" . setup-default-paths )
      ("addIons" . leap.add-ions:add-ions)
      ("setBox" . leap.set-box:set-box)
      ("showPaths" . show-paths)
      ("createAtom" . create-atom )
      ("help" . leap-help)
      ("startSwank" . leap-start-swank)
      ("quit" . leap-quit)
      ))
  (dolist (command leap.parser:*function-names/alist*)
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
          do ;;(handler-case
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
;;               (error (var) (format t "Error ~a - while evaluating: ~a~%" var code)))
    )))

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
                (cmd-assoc (assoc cmd leap.parser:*function-names/alist*))
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
