(load "~/quicklisp/setup.lisp")

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


(defun quickload-parallel (systems level &key (parallel-jobs *number-of-jobs*))
  #+dbg-print(bformat t "DBG-PRINT compile-system files: %s\n" files)
  (let ((total (length systems))
        (counter 0)
        (job-counter 0)
        (child-count 0)
        (jobs (make-hash-table :test #'eql)))
    (labels ((started-quickload (entry counter child-pid)
               (let* ((filename entry))
                 (format t "Compiling [level ~d/~d of ~d (child-pid: ~d)] ~a~%~%" level counter total child-pid filename)))
             (finished-report-one (entry counter child-pid)
               (let* ((filename entry))
                 (format t "Finished [level ~d/~d of ~d (child pid: ~d)] system '~a' - output follows...~%" level counter total child-pid filename)))
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
             (finished-one (entry counter child-pid result-stream)
               (finished-report-one entry counter child-pid)
               (report-result-stream result-stream child-pid))
             (one-quickload (entry counter)
               (ql:quickload entry :verbose t)))  ;; Yeah silent and verbose
      (let (entry wpid status)
        (tagbody
         top
           (incf counter)
           (setq entry (pop systems))
           (incf job-counter)
           (when (> job-counter parallel-jobs)
             (block wait-loop
               (tagbody
                top
                  (multiple-value-setq (wpid status) (core:wait))
                  (core:bformat t "wpid -> %s  status -> %s\n" wpid status)
                  (when (core:wifexited status) (go done))
                  (when (core:wifsignaled status)
                    (let ((signal (core:wtermsig status)))
                      (warn "Child process with pid ~a got signal ~a" signal)))
                  (go top)
                done
                  ))
             (if (>= wpid 0)
                 (let* ((finished-entries-triplet (gethash wpid jobs))
                        (finished-entries (first finished-entries-triplet))
                        (finished-counter (second finished-entries-triplet))
                        (result-stream (third finished-entries-triplet)))
                   (finished-one finished-entries finished-counter wpid result-stream)
                   (decf child-count))
                 (error "wait returned ~d  status ~d~%" wpid status)))
           (when entry
             (multiple-value-bind (maybe-error pid-or-error result-stream)
                 (core:fork t)
               (if maybe-error
                   (error "Could not fork when trying to build ~a" entry)
                   (let ((pid pid-or-error))
                     (if (= pid 0)
                         (progn
                           ;; Turn off interactive mode so that errors cause clasp to die with backtrace
                           (core:set-interactive-lisp nil)
                           (let ((new-sigset (core:make-cxx-object 'core:sigset))
                                 (old-sigset (core:make-cxx-object 'core:sigset)))
                             (core:sigset-sigaddset new-sigset 'core:signal-sigint)
                             (core:sigset-sigaddset new-sigset 'core:signal-sigchld)
                             (multiple-value-bind (fail errno)
                                 (core:sigthreadmask :sig-setmask new-sigset old-sigset)
                               (one-quickload entry counter)
                               (core:sigthreadmask :sig-setmask old-sigset nil)
                               (when fail
                                 (error "sigthreadmask has an error errno = ~a" errno))
                               (core:exit))))
                         (progn
                           (started-quickload entry counter pid)
                           (setf (gethash pid jobs) (list entry counter result-stream))
                           (incf child-count)))))))
           (when (> child-count 0) (go top))))))
  (format t "Leaving compile-system-parallel~%"))


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


(defparameter *number-of-jobs* 8)

(defparameter *load-systems*
  #(("cando-user")
    ("cando-jupyter" "charges" "design")
    ("aromaticity" "cl-nglview" "leap")
    ("cando" "cl-ipywidgets" "smarts")
    ("cl-ipykernel" "inet" "language.smarts" "utility")
    ("cl-jupyter" "closer-mop" "language.smiles" "parser.common-rules.operators"
     "trivial-http")
    ("architecture.builder-protocol" "bordeaux-threads" "cl-base64"
     "parser.common-rules" "pzmq" "trivial-arguments" "trivial-backtrace"
     "trivial-gray-streams" "usocket" "uuid")
    ("cffi" "esrap" "let-plus" "split-sequence" "trivial-utf-8")
    ("anaphora" "babel" "uiop") ("alexandria" "trivial-features")))


(defun quickload-all () 
  (loop for level downfrom (1- (length *load-systems*)) to 1
        for systems = (elt *load-systems* level)
        do (when systems
             (format t "Loading systems: ~a~%" systems)
             (quickload-parallel systems level))))
 

(quickload-all)

