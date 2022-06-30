(in-package :amber)

(defmacro convert-to-script (&rest code)
  `(cl:format nil "~{~s~%~}" ',code))

;;; ------------------------------------------------------------
;;;
;;; job-name is a unique name for every job
;;;
;;; It has two parts, a verb (what is done) and a noun (aggregate name)
;;; names are like MINIMIZE-MOL0 or HEAT-MOL5

(defclass job-name ()
  ((path :initform nil :initarg :path :accessor path)
   (verb :initarg :verb :accessor verb)
   (noun :initarg :noun :accessor noun)))

(defun job-pathname (job-name &key type)
  (make-pathname :name (format nil "~a-~a" (string (verb job-name)) (string (noun job-name)))
                 :type type))



;;; ------------------------------------------------------------
;;;
;;; File node types
;;;

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

(defmethod print-object ((node-file node-file) stream)
  (print-unreadable-object (node-file stream :type t)
    (format stream "~a" (filename node-file))))


;;; ------------------------------------------------------------
;;;
;;; Return a pathname for the node-file
;;;

(defgeneric node-pathname (simulation node &optional relative))

(defmethod node-pathname :around (simulation (node node-file) &optional relative)
  (declare (ignore simulation node))
  (ensure-directories-exist (call-next-method)))

(defmethod node-pathname (simulation (node node-file) &optional relative)
  (let ((pathname (merge-pathnames (merge-pathnames (filename node) (path simulation)))))
    pathname))

(defun node-namestring (simulation node)
  (namestring (node-pathname simulation node)))


(defun in-simulation-node-namestring (simulation node-file)
  (let* ((pathname (node-namestring simulation node-file))
         (enough (enough-namestring (merge-pathnames pathname) (path simulation))))
    enough))

(defclass simulation ()
  ((path :initarg :path :accessor path)
   (latest-jobs :initform nil :accessor latest-jobs)
   (aggregates :initform (make-hash-table) :accessor aggregates)
   (parameters :initform (make-hash-table) :accessor parameters)
   (jobs :initform nil :initarg :jobs :accessor jobs)))

(defmethod print-object ((object simulation) stream)
  (print-unreadable-object (object stream :type t)
    (format stream "#jobs ~a ~s" (length (jobs object)) (sort (alexandria:hash-table-keys (aggregates object)) #'string<))))


(defun simulation ()
  (let ((sim (make-instance 'simulation
                            )))
    sim))

(defclass aggregate-nodes ()
  ((aggregate :initarg :aggregate :accessor aggregate)
   (parm-node :initarg :parm-node :accessor parm-node)
   (coord-node :initarg :coord-node :accessor coord-node)))

(defun addAggregates (simulation aggregate-list)
  (let ((jobs (loop for aggregate in aggregate-list
                    for name = (chem:get-name aggregate)
                    for parm-node = (make-node-file (make-pathname :name (format nil "start-~a" (string name)) :type "parm"))
                    for coord-node = (make-node-file (make-pathname :name (format nil "start-~a" (string name)) :type "crd"))
                    do (setf (gethash name (aggregates simulation)) (make-instance 'aggregate-nodes
                                                                                   :aggregate aggregate
                                                                                   :parm-node parm-node
                                                                                   :coord-node coord-node))
                    collect (connect-graph
                             simulation
                             (make-instance 'jupyter-job
                                            :the-simulation simulation
                                            :aggregate aggregate
                                            :script nil
                                            :job-name (make-instance 'job-name :verb "START" :noun name)
                                            :outputs (arguments :|-p| parm-node
                                                                :|-r| coord-node))))))
    (format t "Setting latest-jobs simulation -> ~a~%" jobs)
    (setf (latest-jobs simulation) jobs)))

(defun addAggregate (simulation aggregate)
  (addAggregates simulation (list aggregate)))

(defun saveSimulation (simulation &key (name "jobs/") (if-exists :overwrite))
  (let ((pathname (uiop:ensure-directory-pathname name)))
    (when (and (probe-file pathname) (not (eq if-exists :overwrite)))
      (format t "A simulation already exists in ~a~%" pathname)
      (return-from saveSimulation nil))
    (setf (path simulation) (merge-pathnames (pathname name)))
    (maphash (lambda (key agg-node)
               (declare (ignore key))
               (let ((aggregate (aggregate agg-node))
                     (top-name (parm-node agg-node))
                     (coord-name (coord-node agg-node)))
                 (ensure-directories-exist (node-pathname simulation top-name))
                 (leap:save-amber-parm aggregate
                                       (node-pathname simulation top-name)
                                       (node-pathname simulation coord-name))))
             (aggregates simulation))
    (generate-runcmd simulation)
    (generate-all-code simulation)))


(defun startSimulation (simulation &key (name "jobs/") (if-exists :overwrite))
  (saveSimulation simulation :name name :if-exists if-exists)
  (let ((cmd (format nil "cd ~a; make &" (namestring (path simulation)))))
    (ext:run-program "/bin/bash" (list "-c" cmd))))

(defclass job ()
  ((the-simulation :initarg :the-simulation :accessor the-simulation)
   (aggregate :initarg :aggregate :accessor aggregate)
   (job-name :initarg :job-name :accessor job-name)
   (script :initform nil :initarg :script :accessor script)
   (inputs :initform nil :initarg :inputs :accessor inputs)
   (outputs :initform nil :initarg :outputs :accessor outputs)
   (parameters :initform nil :initarg :parameters :accessor parameters)
   (makefile-clause :initarg :makefile-clause :accessor makefile-clause)))

(defmethod print-object ((obj job) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (job-pathname (job-name obj)))))

(defclass jupyter-job (job)
  ()
  (:default-initargs
   :makefile-clause nil))


(defclass dynamics (job)
  ()
  )

(defun job-substitutions (job)
  (let (substitutions
        (simulation (the-simulation job)))
    (push (cons (format nil "-i")
                (node-namestring simulation (script job)))
          substitutions)
    (loop for input in (inputs job)
          for option = (option input)
          for file = (node input)
          if (eq option :.)
            do (let ((dot-parts (format nil "(list ~{\"~a\" ~})" (reverse (mapcar (lambda (x) (node-namestring simulation x)) file)))))
                 (push (cons :%.PARTS% dot-parts) substitutions))
          else
            do (push (cons (intern (string-upcase option) :keyword)
                           (node-namestring simulation file))
                     substitutions))
    (loop for output in (outputs job)
          for option = (option output)
          for file = (node output)
          for node-namestring = (node-namestring simulation file)
          do (push (cons (intern (string-upcase option) :keyword)
                         node-namestring)
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
  (let* ((simulation (the-simulation job))
         dependency-inputs
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
      (push (in-simulation-node-namestring simulation (script job)) dependency-inputs))
    ;; Add all inputs
    (loop for input in argument-inputs-job
          for node-input = (node input)
          do (pushnew (in-simulation-node-namestring simulation node-input) inputs :test #'string=)
             (if (char= (char (string (option input)) 0) #\-)
                 (progn
                   (push (string-downcase (option input)) option-inputs)
                   (push (in-simulation-node-namestring simulation (node input)) option-inputs))
                 (push (cons (intern (string-upcase (option input)) :keyword)
                             (in-simulation-node-namestring simulation (node input))) extra-substitutions)))
    (when dot-option-arg
      (mapc (lambda (one) (pushnew (in-simulation-node-namestring simulation one) inputs :test #'string=)) (reverse (node dot-option-arg)))
      (push "--" option-inputs)
      (mapc (lambda (one) (pushnew (in-simulation-node-namestring simulation one) option-inputs :test #'string=)) (reverse (node dot-option-arg))))
    (loop for output in (outputs job)
          do (pushnew (in-simulation-node-namestring simulation (node output)) outputs :test #'string=)
             (if (char= (char (string (option output)) 0) #\-)
                 (progn
                   (push (string-downcase (option output)) option-outputs)
                   (push (in-simulation-node-namestring simulation (node output)) option-outputs))
                 (push (cons (intern (string-upcase (option output)) :keyword) (in-simulation-node-namestring simulation (node output))) extra-substitutions)))
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



(defun connect-graph (simulation job)
  (let ((script (script job))
        (inputs (inputs job))
        (outputs (outputs job)))
    (push job (jobs simulation))
    (when script (push job (users script)))
    (loop for input in inputs
          for node-input = (node input)
          if (consp node-input)
            do (mapc (lambda (one) (push job (users one))) node-input)
          else
            do (push job (users node-input)))
    (loop for output in outputs
          do (push job (definers (node output))))
    job))

(defgeneric job-file (what option))

(defmethod job-file (job option)
  (loop for output-arg in (outputs job)
        when (eq option (option output-arg))
          do (return-from job-file (node output-arg)))
  (loop for input-arg in (inputs job)
        when (eq option (option input-arg))
          do (return-from job-file (node input-arg)))
  (error "Could not find option ~a in ~a" option (append (outputs job) (inputs job))))

(defclass argument ()
  ((option :initarg :option :accessor option)
   (node :initarg :node :accessor node)))

(defmethod print-object ((obj argument) stream)
  (print-unreadable-object (obj stream :type t)
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

(defun generate-code (simulation makefile)
  ;; Generate script
  (loop for job in (jobs simulation)
        do (let ((script (script job))
                 script-code)
             (when script
               (let* ((raw-script (script script))
                      (substituted-script (replace-all (substitutions job script) raw-script)))
                 (setf script-code substituted-script)
                 (write-file-if-it-has-changed (node-pathname simulation script) substituted-script)))
             (when (slot-boundp job 'makefile-clause)
               (let ((raw-makefile-clause (makefile-clause job)))
                 (when raw-makefile-clause
                   (let* ((makefile-substitutions (makefile-substitutions job script-code))
                          (substituted-makefile-clause (replace-all makefile-substitutions raw-makefile-clause)))
                     (write-string substituted-makefile-clause makefile))
                   (terpri makefile)
                   (terpri makefile)))))))

(defun write-file-if-it-has-changed (pathname code)
    (when (probe-file pathname)
      (let ((existing-code (with-open-file (fin pathname :direction :input)
                                    (let ((data (make-string (file-length fin))))
                                      (read-sequence data fin)
                                      data))))
        (when (string= code existing-code)
          (format t "Skipping generation of ~a - it has not changed~%" pathname)
          (return-from write-file-if-it-has-changed nil))))
    (with-open-file (fout (ensure-directories-exist pathname) :direction :output :if-exists :supersede)
      (write-string code fout)))

(defun generate-all-code (simulation)
  (let ((final-outputs (loop for job in (latest-jobs simulation)
                             append (outputs job)))
        (all-outputs (loop for job in (jobs simulation)
                           when (not (typep job 'jupyter-job))
                             append (outputs job)))
        (makefile-pathname (ensure-directories-exist (merge-pathnames "makefile" (path simulation)))))
    (format t "Writing makefile to ~a~%" (translate-logical-pathname makefile-pathname))
    (let ((body (with-output-to-string (makefile)
                  (generate-code simulation makefile))))
      (write-file-if-it-has-changed
       makefile-pathname
       (with-output-to-string (makefile)
         (format makefile "RUNCMD = ./runcmd~%")
         (format makefile "all : ~{~a ~}~%" (mapcar (lambda (argument) (in-simulation-node-namestring simulation (node argument))) final-outputs))
         (format makefile "~aecho all DONE~%" #\tab)
         (format makefile "~%")
         (format makefile "clean :~%~arm -f ~{~a ~}~%" #\tab (mapcar (lambda (arg) (in-simulation-node-namestring simulation (node arg))) all-outputs))
         (format makefile "~aecho clean DONE~%" #\tab)
         (format makefile "~%")
         (write-string body makefile)
         (terpri makefile)))))
  nil
  )


(defparameter *output-types*
  '((:|-o| . "out")
    (:|-inf| . "info")
    (:|-e| . "en")
    (:|-r| . "rst7")
    (:|-x| . "nc")))
(defun setup-job (simulation previous-job &key parameters
                                            job-name
                                            input-topology-file
                                            input-coordinate-file
                                            script
                                            makefile-clause
                                            (output-args '(:|-o| :|-inf| :|-e| :|-r| :|-x|)))
  (let ((script-file (make-instance 'amber-script-file
                                    :script script
                                    :filename (job-pathname job-name :type "in"))))
    (connect-graph
     simulation
     (let ((outputs (loop for out-arg in output-args
                          for type = (cdr (assoc out-arg *output-types*))
                          append (list out-arg (make-node-file (job-pathname job-name :type type))))))
       (make-instance 'job
                      :the-simulation simulation
                      :aggregate (aggregate previous-job)
                      :job-name job-name
                      :script script-file
                      :parameters parameters
                      :inputs (arguments :|-i| script-file
                                         :|-c| input-coordinate-file
                                         :|-ref| input-coordinate-file
                                         :|-p| input-topology-file)
                      :outputs (apply 'arguments outputs)
                      :makefile-clause makefile-clause)))))


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

(defparameter *min-in*
"Minimize
 &cntrl
  imin=1,
  ntx=1,
  ntb=:%BOX%,
  irest=0,
  maxcyc=:%MAXCYC%,
  ncyc=:%SDCYC%,
  ntpr=100,
  ntwx=0,
  cut=:%VDW-CUTOFF%,
 /
 "
)

(defun bounding-box-parameters (aggregate &key constant-pressure)
  (let ((bounding-box 0)
        (cutoff 16.0))
    (when (chem:aggregate-bounding-box-bound-p aggregate)
      (setf bounding-box (if constant-pressure
                             2
                             1)
            cutoff 8.0))
    (list (cons :%BOX% bounding-box)
          (cons :%VDW-CUTOFF% cutoff))))

(defun minimize (simulation &key (name "minimize")
                              (steps 1000)
                              (sd-steps 100 sd-steps-p)
                              )
  (let ((jobs (loop for prev-job in (latest-jobs simulation)
                    for aggregate = (aggregate prev-job)
                    collect (setup-job simulation prev-job
                                       :input-topology-file (job-file prev-job :|-p|)
                                       :input-coordinate-file (job-file prev-job :|-r|)
                                       :script *min-in*
                                       :job-name (make-instance 'job-name :verb name :noun (noun (job-name prev-job)))
                                       :parameters (list* (cons :%MAXCYC% steps)
                                                          (cons :%SDCYC% sd-steps)
                                                          (bounding-box-parameters aggregate))
                                       :output-args '(:|-o| :|-inf| :|-r|)
                                       :makefile-clause ":%OUTPUTS% &: :%INPUTS%
	$(RUNCMD) -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	sander :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))))
    (setf (latest-jobs simulation) jobs))
  simulation)

(defparameter *heat-in*
  "heating
 &cntrl
   imin = 0, nstlim = :%STEPS%, irest = 0, ntx = 1, dt = 0.002,
   nmropt = 1,
   ntt = 3, temp0 = :%TEMPERATURE%, gamma_ln = 2.0, ig = -1,
   tempi = 5.0, tautp = 1.0,
   vlimit = 20,
   ntb = :%BOX%,
   ntc = 2, ntf = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = :%STEPS-PER-WRITE%, ntpr = 1000, ntwr = 5000,
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

(defun heat (simulation &key
                          (name "heat")
                          (temperature 300.0)
                          (steps 10000)
                          (steps-per-write 0)
                          )
  (let ((jobs (loop for prev-job in (latest-jobs simulation)
                    for aggregate = (aggregate prev-job)
                    collect (setup-job simulation prev-job
                                       :input-topology-file (job-file prev-job :|-p|)
                                       :input-coordinate-file (job-file prev-job :|-r|)
                                       :script *heat-in*
                                       :job-name (make-instance 'job-name :verb name :noun (noun (job-name prev-job)))
                                       :parameters (list*
                                                    (cons :%STEPS% steps)
                                                    (cons :%STEPS-PER-WRITE% steps-per-write)
                                                    (cons :%TEMPERATURE% temperature)
                                                    (bounding-box-parameters aggregate))
                                       :makefile-clause ":%OUTPUTS% &: :%INPUTS%
	$(RUNCMD) -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	sander :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))))
    (setf (latest-jobs simulation) jobs))
  simulation)

(defparameter *press-in*
  "pressurising
 &cntrl
   imin = 0, nstlim = :%STEPS%, irest = 1, ntx = 5, dt = 0.002,
   ntt = 3, temp0 = :%TEMPERATURE%, gamma_ln = 2.0, ig = -1,
   tautp = 1.0,
   vlimit = 20,
   ntp = 1, pres0 = 1.0, taup = 2.0,
   ntb = :%BOX%,
   ntc = 2, ntf = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = :%STEPS-PER-WRITE%, ntpr = 1000, ntwr = 5000,

   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H='
 /
 &ewald
 / 

")

(defun pressurize (simulation &key (name "pressurize")
                                (steps 10000)
                                (steps-per-write 0)
                                (temperature 300.0))
  (let ((jobs (loop for prev-job in (latest-jobs simulation)
                    for aggregate = (aggregate prev-job)
                    collect (setup-job simulation prev-job
                                       :input-topology-file (job-file prev-job :|-p|)
                                       :input-coordinate-file (job-file prev-job :|-r|)
                                       :script *press-in*
                                       :job-name (make-instance 'job-name :verb name :noun (noun (job-name prev-job)))
                                       :parameters (list*
                                                    (cons :%STEPS% steps)
                                                    (cons :%STEPS-PER-WRITE% steps-per-write)
                                                    (cons :%TEMPERATURE% temperature)
                                                    (bounding-box-parameters aggregate :constant-pressure t))
                                       :makefile-clause ":%OUTPUTS% &: :%INPUTS%
	$(RUNCMD) -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	sander :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))))
    (setf (latest-jobs simulation) jobs))
  simulation)

(defparameter *dynamics-in*
  "dynamics
 &cntrl
   imin = 0, nstlim = :%STEPS%, irest = 0, ntx = 1, dt = 0.002,
   nmropt = 1,
   ntt = 3, temp0 = :%TEMPERATURE%, gamma_ln = 2.0, ig = -1,
   tempi = 5.0, tautp = 1.0,
   vlimit = 20,
   ntb = :%BOX%,
   ntc = 2, ntf = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = :%STEPS-PER-WRITE%, ntpr = 1000, ntwr = 5000,

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

(defun dynamics (simulation &key (name "dynamics")
                              (steps 10000)
                              (steps-per-write 1000)
                              (temperature 300.0)
                              (constant-pressure nil)
                              )
  (let ((jobs (loop for prev-job in (latest-jobs simulation)
                    for aggregate = (aggregate prev-job)
                    collect (setup-job simulation prev-job
                                       :input-topology-file (job-file prev-job :|-p|)
                                       :input-coordinate-file (job-file prev-job :|-r|)
                                       :script *dynamics-in*
                                       :job-name (make-instance 'job-name :verb name :noun (noun (job-name prev-job)))
                                       :parameters (list*
                                                    (cons :%STEPS% steps)
                                                    (cons :%STEPS-PER-WRITE% steps-per-write)
                                                    (cons :%TEMPERATURE% temperature)
                                                    (bounding-box-parameters aggregate :constant-pressure constant-pressure))
                                       :makefile-clause ":%OUTPUTS% &: :%INPUTS%
	$(RUNCMD) -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	sander :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))))
    (setf (latest-jobs simulation) jobs))
  simulation)

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
             :makefile-clause ":%OUTPUTS% &: :%INPUTS%
	$(RUNCMD) -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
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
             :makefile-clause ":%OUTPUTS% &: :%INPUTS%
	$(RUNCMD) -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	pmemd.cuda -AllowSmallBox :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))

(defun read-file (infile)
  (with-open-file (instream infile :direction :input :if-does-not-exist nil)
    (when instream 
      (let ((string (make-string (file-length instream))))
        (read-sequence string instream)
        string))))


(defun generate-runcmd (simulation)
  (let ((runcmd-namestring (namestring (merge-pathnames #P"runcmd" (path simulation)))))
    (with-open-file (fout runcmd-namestring :direction :output :if-exists :supersede)
      (format fout "#! /bin/sh

try_special ()
{
    local arg
    local count=0
    for arg in \"$@\" ; do
        [ \"$arg\" = -- ] && count=\"$(($count + 1))\"
        shift
        [ \"$count\" = 3 ] && break # more than 3 are allowed
    done
    if [ \"$count\" = 3 ] ; then
        exec \"$@\"
        # exit # not reachable
    fi
}

try_special \"$@\"
# if try_special didn't like it execute directly
exec \"$@\"
"))
    (clasp-posix:chmod runcmd-namestring #o755))
  #+(or)
  (let ((run-in-docker-file #P"sys:extensions;cando;src;data;common;run-in-docker"))
    (unless (probe-file run-in-docker-file)
      (error "Could not find file ~s~%" run-in-docker-file))
    (let ((run-in-docker (read-file (probe-file run-in-docker-file))))
      (with-open-file (fout "run-in-docker" :direction :output :if-exists :supersede)
        (write-string run-in-docker fout))
      (clasp-posix:chmod "run-in-docker" #o755)))
  #+(or)
  (let ((runcmd_with_docker-file #P"sys:extensions;cando;src;data;common;runcmd_with_docker"))
    (unless (probe-file runcmd_with_docker-file)
      (error "Could not find file ~s~%" runcmd_with_docker-file))
    (let ((runcmd_with_docker (read-file (probe-file runcmd_with_docker-file))))
      (with-open-file (fout "runcmd_with_docker" :direction :output :if-exists :supersede)
        (write-string runcmd_with_docker fout))
      (clasp-posix:chmod "runcmd_with_docker" #o755)))
  )
