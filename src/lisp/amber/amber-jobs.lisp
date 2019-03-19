(in-package :amber)

(defmacro convert-to-script (&rest code)
  `(cl:format nil "~{~s~%~}" ',code))
             

(defparameter *min-in*
"minimisation
 &cntrl
   imin = 1, ntmin = 2,
   maxcyc = :%MAXCYC%,
   ntpr = 20, ntwe = 20,
   ntb = 1,
   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H=',

   icfe = 1, ifsc = 1, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 0
 /
 &ewald
 / 
")

(defparameter *heat-in*
  "heating
 &cntrl
   imin = 0, nstlim = 10000, irest = 0, ntx = 1, dt = 0.002,
   nmropt = 1,
   ntt = 3, temp0 = 300.0, gamma_ln = 2.0, ig = -1,
   tempi = 5.0, tautp = 1.0,
   vlimit = 20,
   ntb = 1,
   ntc = 2, ntf = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,

   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H=',

 /
 &ewald
 / 

 &wt
   type='TEMP0',
   istep1 = 0, istep2 = 8000,                                      
   value1 = 5.0, value2 = 300.0
 /

 &wt type = 'END'
 /

")

(defparameter *press-in*
  "pressurising
 &cntrl
   imin = 0, nstlim = 10000, irest = 1, ntx = 5, dt = 0.002,
   ntt = 3, temp0 = 300.0, gamma_ln = 2.0, ig = -1,
   tautp = 1.0,
   vlimit = 20,
   ntp = 1, pres0 = 1.0, taup = 2.0,
   ntb = 2,
   ntc = 2, ntf = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,

   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H='
 /
 &ewald
 / 

")

(defparameter *dynamics-in*
  "dynamics
 &cntrl
   imin = 0, nstlim = :%NSTLIM%, irest = 0, ntx = 1, dt = 0.002,
   nmropt = 1,
   ntt = 3, temp0 = :%TEMP0%, gamma_ln = 2.0, ig = -1,
   tempi = 5.0, tautp = 1.0,
   vlimit = 20,
   ntb = 1,
   ntc = 2, ntf = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = :%NTWX%, ntpr = 1000, ntwr = 5000,

   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H=',

 /
 &ewald
 / 

 &wt
   type='TEMP0',
   istep1 = 0, istep2 = 8000,                                      
   value1 = 5.0, value2 = 300.0
 /

 &wt type = 'END'
 /

")


;;; ------------------------------------------------------------
;;;
;;; File node types
;;;
;;; side means either :ligand or :complex
;;; lambda means the lambda value from 0.0 to 1.0
;;; morph refers to the morph between two compounds of a compound graph

(defclass node-file ()
  ((definers :initarg :definers :initform nil :accessor definers)
   (users :initarg :users :initform nil :accessor users)
   (filename :initarg :filename :accessor filename)))

(defun make-node-file (file)
  (make-instance 'node-file :filename file))

(defclass script-mixin ()
  ((script :initarg :script :accessor script)))

(defclass amber-script-file (node-file script-mixin)
  ())



;;; ------------------------------------------------------------
;;;
;;; Return a pathname for the node-file
(defgeneric node-pathname (node))

(defmethod node-pathname :around ((node node-file))
  (ensure-directories-exist (call-next-method)))

(defmethod node-pathname ((node node-file))
  (filename node))

(defclass job ()
  ((script :initform nil :initarg :script :accessor script)
   (inputs :initform nil :initarg :inputs :accessor inputs)
   (outputs :initform nil :initarg :outputs :accessor outputs)
   (parameters :initform nil :initarg :parameters :accessor parameters)
   (makefile-clause :initarg :makefile-clause :accessor makefile-clause)))

(defclass jupyter-job (job)
  ()
  (:default-initargs
   :makefile-clause nil))


(defclass dynamics (job)
  ()
  )

(defun job-substitutions (job)
  (let (substitutions)
    (push (cons (format nil "-i")
                (namestring (node-pathname (script job))))
          substitutions)
    (loop for input in (inputs job)
          for option = (option input)
          for file = (node input)
          if (eq option :.)
            do (let ((dot-parts (format nil "(list ~{\"~a\" ~})" (reverse (mapcar (lambda (x) (namestring (node-pathname x))) file)))))
                 (push (cons :%.PARTS% dot-parts) substitutions))
          else
            do (push (cons (intern (string-upcase option) :keyword)
                           (namestring (node-pathname file)))
                     substitutions))
    (loop for output in (outputs job)
          for option = (option output)
          for file = (node output)
          do (push (cons (intern (string-upcase option) :keyword)
                         (namestring (node-pathname file)))
                   substitutions))
    (loop for parameter in (parameters job)
          do (push (cons (car parameter) (format nil "~a" (cdr parameter)))
                   substitutions))
    substitutions))

(defgeneric substitutions (job node-file))

(defmethod substitutions (job (node-file node-file))
  (job-substitutions job))


(defun makefile-substitutions (job script-code)
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
         (inputs-job (inputs job))
         (dot-option-arg (find-if (lambda (arg) (eq :. (option arg))) inputs-job)) ; Find argument with :. option
         (argument-inputs-job (remove-if (lambda (arg) (eq :. (option arg))) inputs-job)) ; Remove argument with :. option
         extra-substitutions)
    ;; Add the script as an input
    (when (script job)
      (push (namestring (node-pathname (script job))) dependency-inputs))
    ;; Add all inputs
    (loop for input in argument-inputs-job
          for node-input = (node input)
          do (pushnew (namestring (node-pathname node-input)) inputs :test #'string=)
             (if (char= (char (string (option input)) 0) #\-)
                 (progn
                   (push (string-downcase (option input)) option-inputs)
                   (push (namestring (node-pathname (node input))) option-inputs))
                 (push (cons (intern (string-upcase (option input)) :keyword) (namestring (node-pathname (node input)))) extra-substitutions)))
    (when dot-option-arg
      (mapc (lambda (one) (pushnew (namestring (node-pathname one)) inputs :test #'string=)) (reverse (node dot-option-arg)))
      (push "--" option-inputs)
      (mapc (lambda (one) (pushnew (namestring (node-pathname one)) option-inputs :test #'string=)) (reverse (node dot-option-arg))))
    (loop for output in (outputs job)
          do (pushnew (namestring (node-pathname (node output))) outputs :test #'string=)
             (if (char= (char (string (option output)) 0) #\-)
                 (progn
                   (push (string-downcase (option output)) option-outputs)
                   (push (namestring (node-pathname (node output))) option-outputs))
                 (push (cons (intern (string-upcase (option output)) :keyword) (namestring (node-pathname (node output)))) extra-substitutions)))
    (append
     (list* (cons :%DEPENDENCY-INPUTS% (format nil "~{~a ~}" (append dependency-inputs (reverse inputs))))
            (cons :%INPUTS% (format nil "~{~a ~}" (reverse inputs)))
            (cons :%DEPENDENCY-OUTPUTS% (format nil "~{~a ~}" (reverse outputs)))
            (cons :%OUTPUTS% (format nil "~{~a ~}" (reverse outputs)))
            (cons :%OPTION-INPUTS% (format nil "~{~a ~}" (reverse option-inputs)))
            (cons :%OPTION-OUTPUTS% (format nil "~{~a ~}" (reverse option-outputs)))
            (if script-code
                (list (cons :%SCRIPT-CODE% (format nil "~a" script-code)))
                nil))
     extra-substitutions)))



(defun connect-graph (amber-job)
  (let ((script (script amber-job))
        (inputs (inputs amber-job))
        (outputs (outputs amber-job)))
    (when script (push amber-job (users script)))
    (loop for input in inputs
          for node-input = (node input)
          if (consp node-input)
            do (mapc (lambda (one) (push amber-job (users one))) node-input)
          else
            do (push amber-job (users node-input)))
    (loop for output in outputs
          do (push amber-job (definers (node output))))
    amber-job))

(defun job-file (amber-job option)
  (loop for output-arg in (outputs amber-job)
        when (eq option (option output-arg))
          do (return-from job-file (node output-arg)))
  (loop for input-arg in (inputs amber-job)
        when (eq option (option input-arg))
          do (return-from job-file (node input-arg)))
  (error "Could not find option ~a in ~a" option (append (outputs amber-job) (inputs amber-job))))

(defclass argument ()
  ((option :initarg :option :accessor option)
   (node :initarg :node :accessor node)))

(defmethod print-object ((obj argument) stream)
  (print-unreadable-object (obj stream)
    (format stream "~a ~a" (option obj) (node obj))))

(defmethod users ((arg argument))
  (users (node arg)))

(defmethod definers ((arg argument))
  (definers (node arg)))

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

(defmethod generate-code ((job job) makefile visited-nodes &key pathname-defaults)
  ;; Generate script
  (let ((script (script job))
        script-code)
    (when script
      (let* ((raw-script (script script))
             (substituted-script (replace-all (substitutions job script) raw-script)))
        (setf script-code substituted-script)
        (write-file-if-it-has-changed (merge-pathnames (node-pathname script) pathname-defaults) substituted-script)))
    (when (slot-boundp job 'makefile-clause)
      (let ((raw-makefile-clause (makefile-clause job)))
        (when raw-makefile-clause
          (let* ((makefile-substitutions (makefile-substitutions job script-code))
                 (substituted-makefile-clause (replace-all makefile-substitutions raw-makefile-clause)))
            (write-string substituted-makefile-clause makefile))
          (terpri makefile)
          (terpri makefile))))
    (loop for output in (outputs job)
          do (loop for child in (users output)
                   unless (gethash child visited-nodes)
                     do (setf (gethash child visited-nodes) t)
                        (generate-code child makefile visited-nodes :pathname-defaults pathname-defaults)))))

(defun write-file-if-it-has-changed (pathname code)
    (when (probe-file pathname)
      (let ((existing-code (with-open-file (fin pathname :direction :input)
                                    (let ((data (make-string (file-length fin))))
                                      (read-sequence data fin)
                                      data))))
        (when (string= code existing-code)
          (format t "Skipping generation of ~a - it has not changed~%" pathname)
          (return-from write-file-if-it-has-changed nil))))
    (format t "Generating script ~a~%" pathname)
    (with-open-file (fout (ensure-directories-exist pathname) :direction :output :if-exists :supersede)
      (write-string code fout)))

(defun generate-all-code (jupyter-job final-outputs &key pathname-defaults)
  (let ((visited-nodes (make-hash-table))
        (makefile-pathname (ensure-directories-exist (merge-pathnames "makefile" pathname-defaults))))
    (format t "Writing makefile to ~a~%" (translate-logical-pathname makefile-pathname))
    (let ((body (with-output-to-string (makefile)
                  (generate-code jupyter-job makefile visited-nodes :pathname-defaults pathname-defaults))))
      (write-file-if-it-has-changed
       makefile-pathname
       (with-output-to-string (makefile)
         (format makefile "all : ~{~a ~}~%" (mapcar (lambda (file) (node-pathname file)) final-outputs))
         (format makefile "~aecho DONE~%" #\tab)
         (format makefile "~%")
         (write-string body makefile)
         (terpri makefile)))))
  jupyter-job)


(defmacro job-let (argument-pairs)
  

(defun setup-job (&key parameters input-topology-file input-coordinate-file pathname-defaults script makefile-clause)
  (let ((script-file (make-instance 'amber-script-file
                                    :script script
                                    :filename (merge-pathnames (make-pathname :type "in") pathname-defaults))))
    (connect-graph
     (make-instance 'job
                    :script script-file
                    :parameters parameters
                    :inputs (arguments :|-i| script-file
                                       :|-c| input-coordinate-file
                                       :|-ref| input-coordinate-file
                                       :|-p| input-topology-file)
                    :outputs (arguments :|-o| (make-node-file (merge-pathnames (make-pathname :type "out") pathname-defaults))
                                        :|-inf| (make-node-file (merge-pathnames (make-pathname :type "info") pathname-defaults))
                                        :|-e| (make-node-file (merge-pathnames (make-pathname :type "en") pathname-defaults))
                                        :|-r| (make-node-file (merge-pathnames (make-pathname :type "rst7") pathname-defaults))
                                        :|-x| (make-node-file (merge-pathnames (make-pathname :type "nc") pathname-defaults))
                                        :|-l| (make-node-file (merge-pathnames (make-pathname :type "log") pathname-defaults)))
                    :makefile-clause makefile-clause))))


(defun jupyter-job (&key outputs)
  (connect-graph
   (make-instance 'jupyter-job
                  :inputs nil
                  :outputs outputs)))

(defun simple-jupyter-job (&key topology-file coordinate-file)
  (connect-graph
   (make-instance 'jupyter-job
                  :inputs nil
                  :outputs (arguments :|-p| (make-node-file topology-file)
                                      :|-r| (make-node-file coordinate-file)))))
  
(defun minimize (&key previous-job input-topology-file input-coordinate-file (pathname-defaults #P"min" p-d-p) (maxcyc 100))
  (setup-job :input-topology-file (or input-topology-file (job-file previous-job :|-p|))
             :input-coordinate-file (or input-coordinate-file (job-file previous-job :|-r|))
             :script *min-in*
             :parameters (list (cons :%maxcyc% maxcyc))
             :pathname-defaults pathname-defaults
             :makefile-clause ":%OUTPUTS% : :%INPUTS%
	runcmd -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	pmemd -AllowSmallBox :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))

(defun heat (&key previous-job
               input-topology-file
               input-coordinate-file
               (pathname-defaults #P"heat" p-d-p))
  (setup-job :input-topology-file (or input-topology-file (job-file previous-job :|-p|))
             :input-coordinate-file (or input-coordinate-file (job-file previous-job :|-r|))
             :script *heat-in*
             :pathname-defaults pathname-defaults
             :makefile-clause ":%OUTPUTS% : :%INPUTS%
	runcmd -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	pmemd.cuda -AllowSmallBox :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))

(defun pressurize (&key previous-job
                     input-topology-file
                     input-coordinate-file
                     (pathname-defaults #P"pres" p-d-p))
  (setup-job :input-topology-file (or input-topology-file (job-file previous-job :|-p|))
             :input-coordinate-file (or input-coordinate-file (job-file previous-job :|-r|))
             :script *press-in*
             :pathname-defaults pathname-defaults
             :makefile-clause ":%OUTPUTS% : :%INPUTS%
	runcmd -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	pmemd.cuda -AllowSmallBox :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))

(defun dynamics ( &key previous-job
                    input-topology-file input-coordinate-file
                    (pathname-defaults #P"dynamics" p-d-p)
                    (nstlim 10000)
                    (ntwx 100)
                    (temp0 300.0)
                    )
  (format t "pathname-defaults: ~s~%" pathname-defaults)
  (setup-job :input-topology-file (or input-topology-file (job-file previous-job :|-p|))
             :input-coordinate-file (or input-coordinate-file (job-file previous-job :|-r|))
             :script *dynamics-in*
             :parameters (list (cons :%nstlim% nstlim)
                               (cons :%ntwx% ntwx)
                               (cons :%temp0% temp0))
             :pathname-defaults pathname-defaults
             :makefile-clause ":%OUTPUTS% : :%INPUTS%
	runcmd -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	pmemd.cuda -AllowSmallBox :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))

(defun mdcrd (job)
  "Return the mdcrd output for a job"
  (job-file job :|-x|))


(defparameter *replica-exchange-equilibration1*
  "Equilibration
 &cntrl
   irest=0, ntx=1, 
   nstlim=:%NSTLIM%, dt=0.002,
   irest=0, ntt=3, gamma_ln=1.0,
   temp0=:%TEMP0%, ig=:%IG%,
   ntc=2, ntf=2, nscm=1000,
   ntb=0, igb=5,
   cut=999.0, rgbmax=999.0,
   ntpr=500, ntwx=:%NTWX%, ntwr=100000,
   nmropt=1,
 /
 &wt TYPE='END'
 /
")

(defun replica-exchange-equilibrate1 ( &key previous-job
                                         input-topology-file input-coordinate-file
                                         (pathname-defaults #P"equilibrate" p-d-p)
                                         (nstlim 10000)
                                         (ntwx 100)
                                         (temp0 300.0)
                                         (random (random 32768))
                                         )
  (setup-job :input-topology-file (or input-topology-file (job-file previous-job :|-p|))
             :input-coordinate-file (or input-coordinate-file (job-file previous-job :|-r|))
             :script *dynamics-in*
             :parameters (list (cons :%nstlim% nstlim)
                               (cons :%ntwx% ntwx)
                               (cons :%temp0% temp0)
                               (cons :%ig% random))
             :pathname-defaults pathname-defaults
             :makefile-clause ":%OUTPUTS% : :%INPUTS%
	runcmd -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	pmemd.cuda -AllowSmallBox :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))

(defparameter *replica-exchange-run1*
  "Replica exchange
 &cntrl
   irest=0, ntx=1, 
   nstlim=:%NSTLIM%, dt=0.002,
   irest=0, ntt=3, gamma_ln=1.0,
   temp0=:%TEMP0%, ig=:%IG%,
   ntc=2, ntf=2, nscm=1000,
   ntb=0, igb=5,
   cut=999.0, rgbmax=999.0,
   ntpr=100, ntwx=:%NTWX%, ntwr=100000,
   nmropt=1,
   numexchg=1000,
 /
 &wt TYPE='END'
 /
")

(defun replica-exchange1 ( &key previous-job
                            input-topology-file input-coordinate-file
                            (pathname-defaults #P"replica-exchange" p-d-p)
                            (nstlim 10000)
                            (ntwx 100)
                            (temp0 300.0)
                            (random (random 32768))
                            )
  (setup-job :input-topology-file (or input-topology-file (job-file previous-job :|-p|))
             :input-coordinate-file (or input-coordinate-file (job-file previous-job :|-r|))
             :script *replica-exchange-run1*
             :parameters (list (cons :%nstlim% nstlim)
                               (cons :%ntwx% ntwx)
                               (cons :%temp0% temp0)
                               (cons :%ig% random))
             :pathname-defaults pathname-defaults
             :makefile-clause ":%OUTPUTS% : :%INPUTS%
	runcmd -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	pmemd.cuda -AllowSmallBox :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))


