(in-package :amber)

(defclass info-file ()
  ((job-name :initarg :job-name :accessor job-name)
   (micro-name :initarg :micro-name :accessor micro-name)
   (filename :initarg :filename :accessor filename))
  )

(defmethod print-object ((obj info-file) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a:~a" (action (job-name obj)) (target (job-name obj)))))


(defclass dynamics-info-file (info-file)
  ())

(defclass minimize-info-file (info-file)
  ((job-steps :initarg :job-steps :accessor job-steps)))

(defclass panel ()
  ((info-files :initform nil :initarg :info-files :accessor info-files)
   (summary :initform nil :accessor summary)
   (state :initform 0 :accessor state)))

(defun current-info-file (panel)
  (if (< (state panel) (length (info-files panel)))
      (elt (info-files panel) (state panel))
      (elt (info-files panel) (1- (length (info-files panel))))))

(defclass dashboard ()
  ((panels :initarg :panels :accessor panels)
   ))


(defun make-dashboard (simulation)
  (let ((db (make-instance
             'dashboard
             :panels (loop for job in (latest-jobs simulation)
                           with info-files = nil
                           do (loop named inner
                                    for pjob = job then (previous-job pjob)
                                    when (null pjob)
                                      do (return-from inner nil)
                                    do (let* ((info-file-node (job-file pjob :|-inf| nil))
                                              (info-file (when info-file-node
                                                           (if (eq (job-kind pjob) :minimize)
                                                               (make-instance 'minimize-info-file
                                                                              :job-name (job-name pjob)
                                                                              :micro-name (micro-name pjob)
                                                                              :filename (node-pathname simulation info-file-node )
                                                                              :job-steps (minimize-job-steps pjob))
                                                               (make-instance 'dynamics-info-file
                                                                              :job-name (job-name pjob)
                                                                              :micro-name (micro-name pjob)
                                                                              :filename (node-pathname simulation info-file-node))))))
                                         (when info-file
                                           (push info-file info-files))))
                           collect (make-instance 'panel
                                                  :info-files info-files)))))
    db))



(defstruct panel-summary
  total-steps
  completed-steps
  etot
  temp
  press
  eamber
  ewald-error-estimate
  estimated-time-remaining
  )

(defun panel-summary-percent-done (summary)
  (ceiling (* 100.0 (/ (panel-summary-completed-steps summary)
                       (panel-summary-total-steps summary)))))

(defun parse-running-info (lines)
  (let* ((progress-line (elt lines 13))
         (total-steps (parse-integer progress-line :start 14 :junk-allowed t))
         (completed-steps (parse-integer progress-line :start 37 :junk-allowed t))
         (energy-line (elt lines 2))
         (etot (fortran:parse-double-float energy-line :start 10 :end 24))
         (estimated-time-remaining (when (> (length lines) 23)
                                         (let ((time-line (elt lines 23)))
                                           (subseq time-line 30)))))
    (make-panel-summary
     :total-steps total-steps
     :completed-steps completed-steps
     :etot etot
     :estimated-time-remaining estimated-time-remaining)))

(defun read-info-file (info-file)
  (let ((filename (filename info-file)))
    (when (probe-file filename)
      (let ((lines (with-open-file (fin filename :direction :input)
                     (loop named readlines
                           for line = (read-line fin nil :eof)
                           when (eq line :eof)
                             do (return-from readlines llines)
                           collect line into llines))))
        lines))))

(defmethod parse-info-file ((info-file minimize-info-file))
  (let ((lines (read-info-file info-file)))
    (if (< (length lines) 5)
        (values nil nil)
        (let* ((nstep-line (elt lines 3))
               (nstep (parse-integer nstep-line :junk-allowed t))
               (etot (fortran:parse-double-float nstep-line :start 10 :end 24))
               )
          (values (make-panel-summary
                   :total-steps (job-steps info-file)
                   :completed-steps nstep
                   :etot etot)
                  (= nstep (job-steps info-file)))))))


(defmethod parse-info-file ((info-file dynamics-info-file))
  (let ((lines (read-info-file info-file)))
    (if (< (length lines) 13)
        (values nil nil)
        (let ((overall-line (elt lines 11)))
          (cond
            ((string= overall-line "| Final Performance Info:")
             (values nil t))
            ((string= overall-line "| Current Timing Info")
             (values (parse-running-info lines) nil))
            (t (values nil nil)))))))

(defun update-panel (panel)
  (with-slots (name info-files summary state) panel
    (unless (>= state (length info-files))
      (loop named update
            for cur-state from state below (length info-files)
            for info-file = (elt info-files cur-state)
            do (multiple-value-bind (parsed-summary finalp)
                   (parse-info-file info-file)
                 (when (not finalp)
                   (setf summary parsed-summary
                         state cur-state)
                   (return-from update nil)))
            finally (setf state cur-state)
            ))
    (< state (length info-files))))

(defun update (dashboard)
  (some #'identity
        (loop for panel in (panels dashboard)
              collect (update-panel panel))))

(defun panel-line (panel)
  (cond
    ((>= (state panel) (length (info-files panel)))
     (let ((job-name (job-name (current-info-file panel))))
       (format nil "~10a Finished." (format nil "~a:~a" (target job-name) (action job-name)))))
    (t
     (let ((summary (summary panel)))
       (when summary
         (let ((job-name (job-name (current-info-file panel))))
           (format nil "~10a ~a%  ~a" (format nil "~a:~a" (target job-name) (action job-name))
                   (ceiling (* 100.0 (/ (panel-summary-completed-steps summary)
                                        (panel-summary-total-steps summary))))
                   (if (panel-summary-estimated-time-remaining summary)
                       (format nil "~a" (panel-summary-estimated-time-remaining summary))))))))))

(defun panel-lines (dashboard)
  (loop for panel in (panels dashboard)
        do (format t "~a~%" (panel-line panel))))

(defmethod watch (dashboard)
  (loop named outer
        unless (update dashboard)
          do (return-from outer)
        do (finish-output t)
        do (panel-lines dashboard)
        do (sleep 10)))
