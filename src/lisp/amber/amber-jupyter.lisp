(in-package :amber.jupyter)


(defun monitor (simulation)
  (let* ((dashboard (amber:make-dashboard simulation)))
    (amber:update dashboard)
    (flet ((progress-info (panel)
             (let* ((info-file (amber:current-info-file panel))
                    (job-name (amber:job-name info-file))
                    (summary (amber:summary panel)))
               (values (amber:target job-name)
                       (amber:action job-name)
                       (when summary (or (amber:panel-summary-percent-done summary) 0))
                       (when summary (or (amber:panel-summary-estimated-time-remaining summary) "Starting")))))
           (progress-label (action percent time)
             (format nil "~a ~a% - remaining ~a" action percent time)))
      (let* ((jobs (loop for panel in (amber:panels dashboard)
                         for info-file = (amber:current-info-file panel)
                         for job-name = (amber:job-name info-file)
                         for summary = (amber:summary panel)
                         append (multiple-value-bind (target action percent time)
                                    (progress-info panel)
                                  (list (jw:make-int-progress :description target
                                                              :value percent)
                                        (jw:make-label :value (progress-label action
                                                                              percent
                                                                              time))))))
             (bars (jw:make-grid-box :children jobs
                                     :layout (make-instance 'jw:layout
                                                            :width "100%"
                                                            :grid-gap "0.1em 1em"
                                                            :grid-template-columns "auto 1fr"))))
        (mp:process-run-function
         "monitor"
         (lambda ()
           (loop named outer
                 do (unless (amber:update dashboard)
                      (loop named monitor
                            for panel in (amber:panels dashboard)
                            for index from 0
                            for progress-index = (* 2 index)
                            for label-index = (1+ progress-index)
                            for progress-widget = (elt jobs progress-index)
                            for label-widget = (elt jobs label-index)
                            do (setf (jw:widget-value progress-widget) 100
                                     (jw:widget-value label-widget) "Finished"))
                      (return-from outer nil))
                 do (loop named monitor
                          for panel in (amber:panels dashboard)
                          for index from 0
                          for progress-index = (* 2 index)
                          for label-index = (1+ progress-index)
                          for progress-widget = (elt jobs progress-index)
                          for label-widget = (elt jobs label-index)
                          do (multiple-value-bind (target action percent time)
                                 (progress-info panel)
                               (declare (ignore target))
                               (setf (jw:widget-value progress-widget) percent
                                     (jw:widget-value label-widget) (progress-label action percent time))))
                 do (sleep 2))))
        bars))))
