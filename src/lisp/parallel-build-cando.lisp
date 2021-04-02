(load "quicklisp:setup.lisp")

(defparameter *number-of-jobs* 8)

;;;#+(or)
(progn
  (defparameter *log* (open "/tmp/clasp-builder-log.txt" :direction :output :if-exists :supersede))
  (si:fset 'mmsg #'(lambda (whole env)
                           (let ((fmt (cadr whole))
                                 (args (cddr whole)))
                             `(progn
                                (core:bformat *log* ,fmt ,@args)
                                (finish-output *log*))))
           t))

#+(or)
(si:fset 'mmsg #'(lambda (whole env)
                         nil)
         t)


(defun optimize-tree (levels)
  (labels ((earlier-p (system level)
             (loop for index from level below (length levels)
                   thereis (member system (elt levels index) :test #'string=))))
    (let ((new-levels (make-array (length levels))))
      (loop for index from 0 below (length levels)
            do (loop for system in (copy-list (elt levels index))
                     do (pop (elt levels index))
 ;;                  do (format t "Looking at ~a on level ~a~%" system index)
                     unless (earlier-p system index)
                       do (push system (elt new-levels index))))
             ;; Now sort the names
      (loop for index from 0 below (length new-levels)
            for copy = (copy-list (elt new-levels index))
            do (setf (elt new-levels index) (sort copy #'string<)))
      new-levels)))

(defun walk-tree (root &optional (depth 0) (levels (make-array 1 :adjustable t)))
  (labels ((seen (name depth levels)
             ;;(format t "seen lookin for ~s from ~s in ~s~%" name depth levels)
             (let ((is-seen (loop for idx from depth below (length levels)
                                    thereis (let ((level (elt levels idx)))
                                              (format t "looking for ~s level ~s list ~s~%" name idx level)
                                              (member name level :test #'string=)))))
               (format t "    -----> ~s~%" is-seen)
               is-seen))
           (walk-impl (root depth levels)
             (when (>= depth (length levels))
               (adjust-array levels (1+ depth)))
             (let* ((root-system (asdf:find-system root))
                    (children (mapcar #'coerce-system-name (asdf:system-depends-on root-system)))
                    (depth (1+ depth)))
               (loop for child-name in children
                     when child-name
                       do (walk-impl child-name depth levels)
                          (pushnew child-name (elt levels depth) :test #'string=))))
           (coerce-system-name (system-designator)
             (cond
               ((or (symbolp system-designator) (stringp system-designator))
                (string system-designator))
               ((consp system-designator)
                (case (first system-designator)
                  (:version (second system-designator))
                  (:feature nil)
                  (otherwise nil)))
               (t (error "Handle system-designator ~s" system-designator)))))
    ;;  (format t "walk-tree-impl ~s~%" root)
    (when root
      (walk-impl root depth levels)
      (push root (elt levels 0)))
    levels))



(defstruct pjob done pid signals entry counter child-stdout child-stderr)


(defun wait-for-child-to-exit (jobs)
  (mmsg "About to waitpid sigchld-count: %s%N" (clasp-posix:sigchld-count))
  (multiple-value-bind (wpid status)
      (clasp-posix:wait)
    (mmsg "Returned from waitpid with wpid: %s status:%s%N" wpid status)
    (if (not (= wpid 0))
        (progn
          (unless (= status 0)
            (core:bformat t "wpid -> %s  status -> %s%N" wpid status))
          (when (clasp-posix:wifexited status)
            (mmsg "A child exited wpid: %s  status: %s%N" wpid status)
            (return-from wait-for-child-to-exit (values wpid status nil)))
          (when (clasp-posix:wifsignaled status)
            (let ((signal (clasp-posix:wtermsig status)))
              (mmsg "Child process with pid %s got signal %s%N" wpid signal)
              (warn "Child process with pid ~a got signal ~a" wpid signal)))))
    ;; If we drop through to here the child died for some reason - return and inform the parent
    (values wpid status t)))


(defun quickload-parallel (systems level &key (parallel-jobs *number-of-jobs*))
  #+dbg-print(bformat t "DBG-PRINT compile-system files: %s%N" files)
  (let ((total (length systems))
        (counter 0)
        (job-counter 0)
        (child-count 0)
        (jobs (make-hash-table :test #'eql)))
    (labels ((started-quickload (entry counter child-pid)
               (let ((filename entry))
                 (format t "quickloading [level ~d/~d of ~d (child-pid: ~d)] ~a~%~%" level counter total child-pid filename)))
             (report-result-stream (result-stream child-pid)
               (when result-stream ;;;can be empty
                 (with-open-stream (sin result-stream)
                   (do ((line (read-line sin nil nil) (read-line sin nil nil)))
                       ((null line))
                     (write-string "[")
                     (princ child-pid)
                     (write-string "]--> ")
                     (princ line)
                     (terpri)))))
             (finished-report-one (job entry child-pid)
               (let ((filename entry))
                 (format t "Finished [~a of ~d (child pid: ~d)] ~a output follows...~%" (pjob-counter job) total child-pid filename)))
             (read-fd-into-buffer (fd)
               (let ((buffer (make-array 1024 :element-type 'base-char :adjustable nil))
                     (sout (make-string-output-stream)))
                 (clasp-posix:lseek fd 0 :seek-set)
                 (block readloop
                   (tagbody
                    top
                      (mmsg "About to read-fd%N")
                      (multiple-value-bind (num-read errno)
                          (clasp-posix:read-fd fd buffer)
                        (when (= num-read 0)
                          (return-from readloop))
                        (when (< num-read 0)
                          (mmsg "Three was an error reading the stream errno %d%N" errno)
                          (core:bformat t "There was an error reading the stream errno %d%N" errno))
                        (when (> num-read 0)
                          (write-sequence buffer sout :start 0 :end num-read)
                          (mmsg "Wrote <%s>%N" (subseq buffer 0 num-read))))
                      (go top)))
                 (mmsg "Returning with buffer%N")
                 (clasp-posix:close-fd fd)
                 (get-output-stream-string sout)))
             (report-stream (pid fd name)
               (mmsg "In report-stream%N")
               (let* ((buffer (read-fd-into-buffer fd))
                      (sin (make-string-input-stream buffer)))
                 (do ((line (read-line sin nil nil) (read-line sin nil nil)))
                     ((null line))
                   (core:bformat t "--%d(%s)--> %s%N" pid name line)
                   (finish-output))))
             (report-child-exited (child-pid child-stdout child-stderr)
               (mmsg "In report-child-exited: %s %s%N" child-stdout child-stderr)
               (report-stream child-pid child-stdout "out")
               (report-stream child-pid child-stderr "err"))
             (finished-some (child-pid pjob)
               (let* ((entry (pjob-entry pjob))
                      (child-stdout (pjob-child-stdout pjob))
                      (child-stderr (pjob-child-stderr pjob)))
                 (setf (pjob-done pjob) t)
                 (finished-report-one pjob entry child-pid)
                 (report-child-exited child-pid child-stdout child-stderr)))
             (one-quickload (entry counter)
               (ql:quickload entry :verbose t))) ;; Yeah silent and verbose
      (let (entry wpid status)
        (tagbody
         top
           (when (or (and (null systems) (> child-count 0))
                     (> job-counter parallel-jobs))
             (multiple-value-setq (wpid status child-died)
               (wait-for-child-to-exit jobs))
             (if (null child-died)
                 (if (and (numberp wpid) (>= wpid 0))
                     (let* ((pjob (gethash wpid jobs))
                            (finished-entry (pjob-entry pjob)))
                       (finished-some wpid pjob)
                       (when (and (numberp status)
                                  (not (zerop status)))
                         (format *error-output* "wait returned for process ~d status ~d: exiting compile-system~%" wpid status)
                         (ext:quit 1))
                       (decf child-count))
                     (error "wait returned ~d  status ~d~%" wpid status))
                 (if (and (numberp wpid) (>= wpid 0))
                     (progn
                       (finished-some wpid (gethash wpid jobs))
                       (error "The child with wpid ~a died with status ~a - terminating build" wpid status)
                       (ext:quit 1))
                     (progn
                       (error "A child died with wpid ~a status ~a" wpid status)
                       (ext:quit 1)))))
           (when systems
             (incf counter)
             (setq entry (pop systems))
             (incf job-counter)
             (let ((child-stdout (clasp-posix:mkstemp-fd "cando-build-stdout"))
                   (child-stderr (clasp-posix:mkstemp-fd "cando-build-stderr")))
               (multiple-value-bind (maybe-error pid-or-error result-stream)
                   (clasp-posix:fork-redirect child-stdout child-stderr)
                 (if maybe-error
                     (error "Could not fork when trying to build ~a" entry)
                     (let ((pid pid-or-error))
                       (if (= pid 0)
                           (progn
                             (ext:disable-debugger)
                             (let ((new-sigset (core:make-cxx-object 'clasp-posix:sigset))
                                   (old-sigset (core:make-cxx-object 'clasp-posix:sigset)))
                               (clasp-posix:sigset-sigaddset new-sigset 'clasp-posix:signal-sigint)
                               (clasp-posix:sigset-sigaddset new-sigset 'clasp-posix:signal-sigchld)
                               (multiple-value-bind (fail errno)
                                   (core:sigthreadmask :sig-setmask new-sigset old-sigset)
                                 (one-quickload entry counter)
                                 (core:sigthreadmask :sig-setmask old-sigset nil)
                                 (when fail
                                   (error "sigthreadmask has an error errno = ~a" errno))
                                 (finish-output)
                                 (sleep 1)
                                 (ext:quit))))
                           (let ((pjob (make-pjob
                                        :done nil
                                        :pid pid
                                        :signals nil
                                        :entry entry
                                        :counter counter
                                        :child-stdout child-stdout
                                        :child-stderr child-stderr)))
                             (started-quickload entry counter pid)
                             (setf (gethash pid jobs) pjob)
                             (incf child-count))))))))
           (when (> child-count 0) (go top))))))
  (format t "Leaving quickload-parallel~%"))



;;; ------------------------------------------------------------
;;;
;;; Setup asd files for cando
;;;
;;;
(defun all-subdirs (dir)
  (let (dirs)
    (labels ((trav (d)
               (dolist (d (uiop:subdirectories d))
                 (push d dirs)
                 (trav d))))
      (trav dir))
    dirs))

;;; Add the cando hostname
(format t "Setting CANDO hostname~%")
(progn
  (setf (logical-pathname-translations "cando")
        '(("**;*.*" "source-dir:extensions;cando;src;**;*.*"))))

;;; Add directories for ASDF to search for systems
(let* ((topdir (translate-logical-pathname #P"cando:lisp;"))
       (dirs (all-subdirs topdir)))
  (push topdir asdf:*central-registry*)
  (dolist (dir dirs)
    (format t "Pushing dir: ~a~%" dir)
    (push dir asdf:*central-registry*)))

(progn
  (format t "Pushing dir: ~a~%" *default-pathname-defaults*)
  (push *default-pathname-defaults* asdf:*central-registry*))

;;;(make-package :cando)

(let ((amber-home
        (namestring (uiop:ensure-directory-pathname (or (ext:getenv "AMBERHOME") "/amber/")))))
  (setf (logical-pathname-translations "amber")
        (list (list "**;*.*" (concatenate 'string amber-home "/**/*.*"))))
  (format t "Setting *amber-home* -> ~a~%" amber-home))


;;; ------------------------------------------------------------
;;;
;;; Start the parallel build
;;;


(defparameter *load-systems*
  #(("cando-user") ("leap-commands" "leap" "design" "charges" "cando")
    ("architecture.builder-protocol" "parser.common-rules.operators"
     "parser.common-rules" "esrap" "alexandria" "smirnoff" "amber" "charges"
     "leap" "cando" "aromaticity" "smarts" "inet" "utility")
    ("let-plus" "split-sequence" "plump" "cl-ppcre" "aromaticity" "alexandria"
     "smirnoff" "cando" "inet" "utility" "smarts" "architecture.builder-protocol"
     "parser.common-rules.operators" "parser.common-rules" "language.smarts"
     "esrap" "trivial-http")
    ("anaphora" "documentation-utils" "array-utils" "plump" "cando" "aromaticity"
     "smarts" "inet" "utility" "trivial-http" "language.smarts" "let-plus"
     "split-sequence" "language.smiles" "architecture.builder-protocol"
     "parser.common-rules.operators" "parser.common-rules" "esrap" "alexandria"
     "usocket")
    ("trivial-indent" "documentation-utils" "array-utils" "inet" "utility"
     "smarts" "language.smarts" "trivial-http" "usocket" "language.smiles"
     "parser.common-rules.operators" "anaphora" "parser.common-rules"
     "architecture.builder-protocol" "esrap" "let-plus" "alexandria"
     "split-sequence")
    ("trivial-indent" "trivial-http" "language.smarts" "language.smiles"
     "parser.common-rules.operators" "usocket" "parser.common-rules"
     "architecture.builder-protocol" "esrap" "let-plus" "split-sequence"
     "anaphora" "alexandria")
    ("usocket" "language.smiles" "parser.common-rules.operators"
     "parser.common-rules" "architecture.builder-protocol" "esrap" "let-plus"
     "split-sequence" "anaphora" "alexandria")
    ("parser.common-rules" "architecture.builder-protocol" "esrap" "let-plus"
     "split-sequence" "anaphora" "alexandria")
    ("esrap" "let-plus" "split-sequence" "anaphora" "alexandria")
    ("anaphora" "alexandria")))

(defun quickload-calculate-order ()
  "Calculate the order of systems to build"
  (walk-tree "cando-user"))


(defun quickload-all () 
  (loop for level downfrom (1- (length *load-systems*)) to 1
        for systems = (elt *load-systems* level)
        do (when systems
             (format t "Loading systems: ~a~%" systems)
             (quickload-parallel systems level))))
 

(quickload-all)

