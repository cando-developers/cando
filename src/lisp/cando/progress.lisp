(in-package :cando)

(defstruct (progress-bar
            (:constructor make-progress-bar
                (&key (total 100.0) (divisions 100.0)
                   (message "Completed")
                   (message-width 10)
                   (bar-character #\*)
                   (width 50)
                   (style :time)
                   (advance-every-n-calls 1000) ;; Only advance every 1000 calls to progress-advance
                   (on t)
                 &aux (spacing (if (>= divisions 0)
                                   (/ total divisions)
                                   (error "divisions must be >0"))))
                ))
  on total (advance-counter 1) advance-every-n-calls
  (style :time)
  (bar-character #\*)
  (width 50)
  (start-time (get-universal-time))
  divisions message message-width spacing (next 0))

(defun progress-convenient-time (time)
  (block nil
    (when (= time 0) (return "--"))
    (when (< time 60) (return (format nil "~5,1f seconds" time)))
    (setf time (/ time 60.0))
    (when (< time 60) (return (format nil "~5,1f minutes" time)))
    (setf time (/ time 60.0))
    (when (< time 24) (return (format nil "~5,2f hours" time)))
    (setf time (/ time 24.0))
    (format nil "~5,2f days" time)))

(defun progress-display (bar counter)
  (let* ((fraction (/ (float counter) (progress-bar-total bar)))
         (elapsed-time (- (get-universal-time) (progress-bar-start-time bar)))
         (expected-time (* (/ 1.0 fraction) elapsed-time))
         (remaining-time (max 0.0 (* (- 1.0 fraction) expected-time))))
    (case (progress-bar-style bar)
      (:bar
       (let* ((star-width (floor (* fraction (progress-bar-width bar))))
              (stars (make-string (max 0 star-width) :initial-element (progress-bar-bar-character bar)))
              (lines (make-string (max 0 (- (progress-bar-width bar) star-width)) :initial-element #\.)))
         (format t "~C~va [~a~a] ETC: ~a"
                 (if (core:is-interactive-lisp) #\return #\space)
                 (progress-bar-message-width bar)
                 (progress-bar-message bar)
                 stars
                 lines
                 (progress-convenient-time remaining-time))))
      (otherwise
       (format t "~C~va ~5,1f%~@[ elapsed: ~a~]~@[ ETC: ~a~]   "
               (if (core:is-interactive-lisp) #\return #\space)
               (progress-bar-message-width bar)
               (progress-bar-message bar)
               (floor (* fraction 100.0))
               (progress-convenient-time elapsed-time)
               (progress-convenient-time remaining-time))))
    (finish-output)
    (let* ((new-next (* (1+ (floor (/ counter (progress-bar-spacing bar)))) (progress-bar-spacing bar))))
      (setf (progress-bar-next bar) new-next))))

(defun progress-advance (bar counter &optional message)
  (when message
    (setf (progress-bar-message bar) message))
  (if (<= (decf (progress-bar-advance-counter bar)) 0)
      (setf (progress-bar-advance-counter bar) (progress-bar-advance-every-n-calls bar))
      (when (and (progress-bar-on bar) (>= counter (progress-bar-next bar)))
        (progress-display bar counter))))

(defun progress-done (bar)
  (progress-advance bar (progress-bar-total bar))
  (setf (progress-bar-advance-counter bar) 0)
  (terpri)
  (finish-output))
