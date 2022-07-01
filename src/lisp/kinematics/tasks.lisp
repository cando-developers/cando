(in-package :kin.tasks)

(defun sanitize-filename (name)
  (with-output-to-string (sout)
    (loop for ch across name
          if (member ch (list #\( #\) #\[ #\] #\,))
            do (princ #\_ sout)
          else
            do (princ ch sout))))

(defclass node ()
  ((definers :initarg :definers :initform nil :accessor definers)
   (users :initarg :users :initform nil :accessor users)
   (pname :initform nil :initarg :pname :accessor pname)))

(defmethod print-object ((obj node) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (pname obj))))


(defgeneric relative-node-pathname (node))

(defmethod relative-node-pathname ((node node))
  (pname node))

(defclass trainer-node (node)
  ())

(defun make-trainer-node (trainer-name)
  (make-instance 'trainer-node
                 :pname (make-pathname
                            :directory (list :relative (sanitize-filename trainer-name))
                            :name "trainer"
                            :type "trainer")))

(defun project-node-pathname (node project)
  (let ((node-pathname (relative-node-pathname node)))
    (merge-pathnames node-pathname (directory project))))

(defclass grammar-node (node) ())

(defun make-grammar-node ()
  (make-instance 'grammar-node :pname (make-pathname :name "grammar"
                                                        :type "cando")))

(defclass sdf-node (trainer-node) ())

(defun make-sdf-node (trainer-node name)
  (make-instance 'sdf-node
                 :pname (make-pathname :directory (pathname-directory (pname trainer-node))
                                          :name name
                                          :type "sdf")))


(defclass log-node (trainer-node) ())

(defun make-log-node (trainer-node name)
  (make-instance 'log-node
                 :pname (make-pathname :directory (pathname-directory (pname trainer-node))
                                          :name name
                                          :type "log")))

(defclass cando-aggregate-node (trainer-node) ())

(defun make-cando-aggregate-node (trainer-node)
  (make-instance 'cando-aggregate-node
                 :pname (make-pathname :directory (pathname-directory (pname trainer-node))
                                          :name "minimized"
                                          :type "cando-aggregate")))

(defclass script-node (trainer-node)
  ((script :initarg :script :accessor script)))

(defun make-script-node (trainer-node &key script script-name)
  (make-instance 'script-node
                 :script script
                 :pname (make-pathname :directory (pathname-directory (pname trainer-node))
                                          :name script-name
                                          :type "lisp")))

(defclass conformation-node (node)
  ())

(defun make-conformation-node (trainer-node)
  (make-instance 'conformation-node
                 :pname (make-pathname :directory (pathname-directory (pname trainer-node))
                                          :name "conf"
                                          :type "conf")))

(defclass sketch-node (node)
  ())

(defun make-sketch-node (trainer-node)
  (make-instance 'sketch-node
                 :pname (make-pathname :directory (pathname-directory (pname trainer-node))
                                          :name "structure"
                                          :type "svg")))

(defclass temp-conformation-node (node)
  ())

(defun make-temp-node (trainer-node)
  (make-instance 'temp-conformation-node
                 :pname (make-pathname :directory (pathname-directory (pname trainer-node))
                                          :name "conf"
                                          :type "temp-conf")))

(defclass job ()
  ((inputs :initform nil :initarg :inputs :accessor inputs)
   (temps :initform nil :initarg :temps :accessor temps)
   (outputs :initform nil :initarg :outputs :accessor outputs)
   (parameters :initform nil :initarg :parameters :accessor parameters)
   (script :initform nil :initarg :script :accessor script)
   (makefile-clause :initarg :makefile-clause :accessor makefile-clause)))

(defclass build-conformation-job (job)
  ())

(defclass project ()
  ((directory :initarg :directory :initform (pathname "jobs/") :accessor directory)
   (grammar :initarg :grammar :accessor grammar)
   (settings :initform nil :initarg :settings :accessor settings)
   (nodes :initform nil :initarg :nodes :accessor nodes)
   (jobs :initform nil :initarg :jobs :accessor jobs)))

(defun connect-graph (amber-job)
  (let ((script (script amber-job))
        (inputs (inputs amber-job))
        (outputs (outputs amber-job)))
    (pushnew amber-job (users script))
    (loop for input in inputs
          for node-input = (node input)
          if (consp node-input)
            do (mapc (lambda (one) (push amber-job (users one))) node-input)
          else
            do (push amber-job (users node-input)))
    (loop for output in outputs
          do (push amber-job (definers (node output))))
    amber-job))

(defclass argument ()
  ((option :initarg :option :accessor option)
   (node :initarg :node :accessor node)))

(defmethod relative-node-pathname ((argument argument))
  (let ((node (node argument)))
    (relative-node-pathname node)))

(defmethod print-object ((obj argument) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a" (option obj) (node obj))))

(defmethod users ((arg argument))
  (users (node arg)))

(defmethod (setf users) (value (arg argument))
  (setf (users (node arg)) value))

(defmethod definers ((arg argument))
  (definers (node arg)))

(defmethod (setf definers) (value (arg argument))
  (setf (definers (node arg)) value))


(defun single-argument (option node)
  "Handle only a single argument"
  (make-instance 'argument :option option :node node))

(defun arguments (&rest args)
  "Later we will add sanity checks to the args - for now just accumulate a list"
  (loop for cur = args then (cddr cur)
        for option = (first cur)
        for node = (second cur)
        until (null cur)
        while option
        collect (make-instance 'argument :option option :node node)))


(defun replace-all (dict in-script)
  (let (script-result)
    (loop for script = in-script then script-result
          for (match . substitution) in dict
          for match-string = (cond
                               ((stringp match) (format nil ":~a" match))
                               ((symbolp match) (format nil "~s" match))
                               (t (error "Illegal match ~s for regex key" match)))
          do (setf script-result (cl-ppcre:regex-replace-all match-string script substitution))
          finally (return-from replace-all script-result))))


(defun get-makefile-substitutions (project job)
  "This returns an alist of label/string substitutions.
The fancy part is the inputs - inputs that have the form :-xxx are added as option-inputs
If there is one (and only one) input with the argument :. - that is appended to the option-inputs with the
form '--' <list of file node pathnames>.  inputs and outputs with names like :yxxx where 'y' is NOT - or . are
added to inputs and outputs but not option-inputs or option-outputs"
  (let* (dependency-inputs
         inputs
         option-inputs
         outputs
         option-outputs
         temps
         option-temps
         script
         (inputs-job (inputs job))
         (dot-option-arg (find-if (lambda (arg) (eq :. (option arg))) inputs-job)) ; Find argument with :. option
         (argument-inputs-job (remove-if (lambda (arg) (eq :. (option arg))) inputs-job)) ; Remove argument with :. option
         extra-substitutions)
    ;; Add all the scripts as an input
    (setf script (script job))
    (push (namestring (relative-node-pathname script)) dependency-inputs)
    ;; Add all inputs
    (loop for input in argument-inputs-job
          for node-input = (node input)
          do (pushnew (namestring (relative-node-pathname node-input)) inputs :test #'string=)
             (if (char= (char (string (option input)) 0) #\-)
                 (progn
                   (push (string-downcase (option input)) option-inputs)
                   (push (namestring (relative-node-pathname (node input))) option-inputs))
                 (push (cons (intern (string-upcase (option input)) :keyword) (namestring (relative-node-pathname (node input)))) extra-substitutions)))
    (when dot-option-arg
      (mapc (lambda (one) (pushnew (namestring (relative-node-pathname one)) inputs :test #'string=)) (reverse (node dot-option-arg)))
      (push "--" option-inputs)
      (mapc (lambda (one) (pushnew (namestring (relative-node-pathname one)) option-inputs :test #'string=)) (reverse (node dot-option-arg))))
    (loop for output in (outputs job)
          do (pushnew (namestring (relative-node-pathname (node output))) outputs :test #'string=)
             (if (char= (char (string (option output)) 0) #\-)
                 (progn
                   (push (string-downcase (option output)) option-outputs)
                   (push (namestring (relative-node-pathname (node output))) option-outputs))
                 (push (cons (intern (string-upcase (option output)) :keyword) (namestring (relative-node-pathname (node output)))) extra-substitutions)))
    (loop for temp in (temps job)
          do (pushnew (namestring (relative-node-pathname (node temp))) temps :test #'string=)
             (if (char= (char (string (option temp)) 0) #\-)
                 (progn
                   (push (string-downcase (option temp)) option-temps)
                   (push (namestring (relative-node-pathname (node temp))) option-temps))
                 (push (cons (intern (string-upcase (option temp)) :keyword) (namestring (relative-node-pathname (node temp)))) extra-substitutions)))
    (list* (cons :%DEPENDENCY-INPUTS% (format nil "~{~a~^ ~}" (append dependency-inputs (reverse inputs))))
           (cons :%INPUTS% (format nil "~{~a~^ ~}" (reverse inputs)))
           (cons :%DEPENDENCY-OUTPUTS% (format nil "~{~a~^ ~}" (reverse outputs)))
           (cons :%DEPENDENCY-TEMPS% (format nil "~{~a~^ ~}" (reverse temps)))
           (cons :%OUTPUTS% (format nil "~{~a~^ ~}" (reverse outputs)))
           (cons :%OPTION-INPUTS% (format nil "~{~a~^ ~}" (reverse option-inputs)))
           (cons :%OPTION-OUTPUTS% (format nil "~{~a~^ ~}" (reverse option-outputs)))
           (cons :%OPTION-TEMPS% (format nil "~{~a~^ ~}" (reverse option-temps)))
           (cons :%SCRIPT% (format nil "~a" (relative-node-pathname script)))
           (cons :%CANDO% "$(CANDO)")
           extra-substitutions)))

(defun standard-makefile-clause (command)
  (format nil ":%OUTPUTS% &: :%DEPENDENCY-INPUTS%
	$(RUNCMD) -c :%JOBID% -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- :%DEPENDENCY-TEMPS% -- \\
	~a~%" command))

(defparameter *write-files* nil)
(defmacro with-top-directory ((project) &body body)
  (declare (ignore project))
  `(let ((*write-files* t))
     ,@body))

(defun write-file-if-it-has-changed (pathname code)
    (when (probe-file pathname)
      (let ((existing-code (with-open-file (fin pathname :direction :input)
                                    (let ((data (make-string (file-length fin))))
                                      (read-sequence data fin)
                                      data))))
        (when (string= code existing-code)
          (return-from write-file-if-it-has-changed nil))))
    (with-open-file (fout (ensure-jobs-directories-exist pathname) :direction :output :if-exists :supersede)
      (write-string code fout)))

(defun ensure-jobs-directories-exist (pathname)
        (ensure-directories-exist pathname))

(defclass makefile-joblet ()
  ((%job :initarg :%job :accessor %job)
   (jobid :initarg :jobid :accessor jobid)
   (template :initarg :template :accessor template)
   (output-files :initform nil :initarg :output-files :accessor output-files)))

(defclass script-joblet ()
  ((%job :initarg :%job :accessor %job)
   (script :initarg :script :accessor script)
   (target-pathname :initarg :target-pathname :accessor target-pathname)))

(defun job-substitutions (job project)
  (let (substitutions)
    ;; I am removing this because only some jobs
    ;; may have -i <script> as an argument but not all of them
    ;; meister May 2020
    #+(or)(push (cons (format nil "-i")
                      (namestring (project-node-pathname (script job) project)))
                substitutions)
    (loop for input in (inputs job)
          for option = (option input)
          for file = (node input)
          if (eq option :.)
            do (let ((dot-parts (format nil "(list ~{\"~a\" ~})" (reverse (mapcar (lambda (x) (namestring (relative-node-pathname x))) file)))))
                 (push (cons :%.PARTS% dot-parts) substitutions))
          else
            do (push (cons (intern (string-upcase option) :keyword)
                           (namestring (relative-node-pathname file)))
                     substitutions))
    (loop for output in (outputs job)
          for option = (option output)
          for file = (node output)
          do (push (cons (intern (string-upcase option) :keyword)
                         (namestring (relative-node-pathname file)))
                   substitutions))
    substitutions))


(defgeneric substitutions (project job node-file)
  (:documentation "Do the node-file (script?) substitutions for the job"))

(defmethod substitutions (project job (node node))
  (declare (ignore node))
  (append (job-substitutions job project)
          (setting-substitutions project)))

(defvar *testing*)

(defun setting-substitutions (project)
  (let ((settings (settings project)))
    (loop for setting in settings
          for setting-key = (car setting)
          for setting-value = (if (listp (cdr setting))
                                  (if *testing*
                                      (second (cdr setting))
                                      (first (cdr setting)))
                                  (cdr setting))
          for value-string = (cond
                               ((stringp setting-value)
                                setting-value)
                               ((floatp setting-value)
                                (format nil "~f" setting-value))
                               ((integerp setting-value)
                                (format nil "~d" setting-value))
                               (t (error "Add support to stringify setting value ~s:" setting-value)))
          collect (cons setting-key value-string))))

(defun do-joblet (project joblet)
  (let* ((script (script joblet))
         (raw-script (script script))
         (substitutions (substitutions project (%job joblet) script))
         (substituted-script (replace-all substitutions raw-script))
         (pn (target-pathname joblet)))
    (ensure-jobs-directories-exist pn)
    (write-file-if-it-has-changed pn substituted-script)))

(defmethod generate-code (project (job job) script-joblets makefile-joblets visited-nodes)
  ;; Generate script
  (let* ((script (script job))
         (script-joblet (make-instance 'script-joblet :script (node script) :%job job :target-pathname (project-node-pathname script project))))
    (vector-push-extend script-joblet script-joblets))
  (let (makefile-entry
        (raw-makefile-clause (makefile-clause job)))
    (when raw-makefile-clause
      (setf makefile-entry (make-instance 'makefile-joblet
                                          :%job job
                                          :template raw-makefile-clause))
      (vector-push-extend makefile-entry makefile-joblets))
    (loop for output in (outputs job)
          do (when makefile-entry (push output (output-files makefile-entry)))
          do (loop for child in (users output)
                   unless (gethash child visited-nodes)
                     do (setf (gethash child visited-nodes) t)
                        (generate-code project child script-joblets makefile-joblets visited-nodes)))))

(defun generate-all-code (project jobs &key final-outputs progress-callback)
  (with-top-directory (project)
    (let ((visited-nodes (make-hash-table))
          (script-joblets (make-array 256 :adjustable t :fill-pointer 0))
          (makefile-joblets (make-array 256 :adjustable t :fill-pointer 0))
          (makefile-pathname (ensure-jobs-directories-exist (merge-pathnames "makefile" (directory project))))
          (status-pathname (ensure-jobs-directories-exist (merge-pathnames "status" (directory project)))))
      (loop for job in jobs
            with number-of-jobs = (length jobs)
            do (generate-code project job script-joblets makefile-joblets visited-nodes))
      ;;; do (when progress-callback (funcall progress-callback number-of-jobs))
      (let ((number-of-joblets (length script-joblets)))
        (map nil (lambda (joblet)
                   (do-joblet project joblet)
                   (when progress-callback
                     (funcall progress-callback number-of-joblets)))
             script-joblets))
      (format t "Writing makefile to ~a~%" (translate-logical-pathname makefile-pathname))
      (let ((body (with-output-to-string (makefile)
                    (format makefile "RUNCMD ?= ./runcmd_simple~%~%")
                    (loop for entry across makefile-joblets
                          for jobid from 100001
                          for entry-substitutions = (get-makefile-substitutions project (%job entry))
                          for substituted-makefile-clause = (replace-all (list* (cons :%JOBID% (format nil "~a" jobid))
                                                                                entry-substitutions)
                                                                         (template entry))
                          do (setf (jobid entry) jobid)
                          do (write-string substituted-makefile-clause makefile)
                          do (terpri makefile)
                          do (terpri makefile)))))
        (write-file-if-it-has-changed
         makefile-pathname
         (with-output-to-string (makefile)
           (format makefile "# cando version ~a~%" (lisp-implementation-version))
           (format makefile "# task version ~a~%" (tirun:version))
           (format makefile "all : ~{~a~^ ~}~%" (mapcar (lambda (file) (relative-node-pathname file)) final-outputs))
           (format makefile "~aecho DONE~%" #\tab)
           (format makefile "~%")
           (write-string body makefile)
           (terpri makefile))))
      (write-file-if-it-has-changed
       status-pathname
       (with-output-to-string (status)
         (loop for entry across makefile-joblets
               do (loop for output in (output-files entry)
                        for check-pathname = (project-node-pathname (node output) project)
                        for enough-name = (enough-namestring check-pathname *default-pathname-defaults*)
                        do (format status "check-status ~a ~a~%" (jobid entry) enough-name)))))
      (format t "Writing runcmd~%")
      (generate-runcmd project))))

(defun read-file (infile)
  (with-open-file (instream infile :direction :input :if-does-not-exist nil)
    (when instream 
      (let ((string (make-string (file-length instream))))
        (read-sequence string instream)
        string))))

(defun generate-runcmd (project)
  (flet ((top-dir (name)
           (make-pathname :name name
                          :defaults (directory project))))
    (with-open-file (fout (top-dir "runcmd_simple") :direction :output :if-exists :supersede)
      (format fout "#! /bin/sh

try_special ()
{
    local arg
    local count=0
    for arg in \"$@\" ; do
        [ \"$arg\" = -- ] && count=\"$(($count + 1))\"
        shift
        [ \"$count\" = 4 ] && break # more than 4 are not allowed
    done
    if [ \"$count\" = 4 ] ; then
        exec \"$@\"
        # exit # not reachable
    fi
}

try_special \"$@\"
# if try_special didn't like it execute directly
exec \"$@\"
"))
    (clasp-posix:chmod (top-dir "runcmd_simple") #o755)
    (let ((run-in-docker-file #P"source-dir:extensions;cando;src;data;common;run-in-docker"))
      (unless (probe-file run-in-docker-file)
        (error "Could not find file ~s~%" run-in-docker-file))
      (let ((run-in-docker (read-file (probe-file run-in-docker-file))))
        (with-open-file (fout (top-dir "run-in-docker") :direction :output :if-exists :supersede)
          (write-string run-in-docker fout))
        (clasp-posix:chmod (top-dir "run-in-docker") #o755)))
    (let ((runcmd_with_docker-file #P"source-dir:extensions;cando;src;data;common;runcmd_with_docker"))
      (unless (probe-file runcmd_with_docker-file)
        (error "Could not find file ~s~%" runcmd_with_docker-file))
      (let ((runcmd_with_docker (read-file (probe-file runcmd_with_docker-file))))
        (with-open-file (fout (top-dir "runcmd_with_docker") :direction :output :if-exists :supersede)
          (write-string runcmd_with_docker fout))
        (clasp-posix:chmod (top-dir "runcmd_with_docker") #o755)))
    ))

(defun save-project (project grammar trainers)
  (let* ((grammar-node (make-grammar-node))
         (pname (project-node-pathname grammar-node project)))
    (let ((trainer-nodes (loop for trainer in trainers
                               for trainer-node = (make-trainer-node (kin:trainer-filename trainer))
                               for trainer-node-pathname = (kin.tasks:relative-node-pathname trainer-node)
                               do (push trainer-node (nodes project))
                               do (push trainer-node-pathname (kin:trainer-node-pathnames grammar))
                               do (kin:maybe-save-trainer trainer (project-node-pathname trainer-node project))
                               collect trainer-node)))
      ;; Now save the grammar file
      (ensure-directories-exist pname)
      (cando:save-cando (grammar project) pname)
      (multiple-value-bind (roots final-outputs)
          (loop
            for trainer-node in trainer-nodes
            for outputs = (arguments :--output (make-conformation-node trainer-node)
                                     :--sketch (make-sketch-node trainer-node)
                                     :--allsdf (make-sdf-node trainer-node "all")
                                     :--uniquesdf (make-sdf-node trainer-node "unique")
                                     :--log (make-log-node trainer-node "log")
                                     )
            for job = (make-instance 'build-conformation-job
                                     :inputs (arguments :--trainer trainer-node :--grammar grammar-node)
                                     :temps (arguments :--temp (make-temp-node trainer-node))
                                     :outputs outputs
                                     :script (single-argument :--script (make-script-node trainer-node
                                                                                          :script kin.work:*build-conformation*
                                                                                          :script-name "build-conformation"))
                                     :makefile-clause (standard-makefile-clause ":%CANDO% -l :%SCRIPT% -- :%OPTION-INPUTS% :%OPTION-TEMPS% :%OPTION-OUTPUTS%"))
            collect (connect-graph job) into roots
            append outputs into final-outputs
            finally (return (values roots final-outputs)))
        (generate-all-code project roots :final-outputs (mapc #'node final-outputs))
        trainer-nodes))))
      
(defmacro with-project ((project) &body body)
  (let ((dir-gs (gensym)))
    `(let ((,dir-gs (ext:getcwd)))
       (ext:chdir (namestring (merge-pathnames (directory ,project) (pathname ,dir-gs))) t)
       (progn
         ,@body)
       (ext:chdir (namestring (pathname ,dir-gs))))))
