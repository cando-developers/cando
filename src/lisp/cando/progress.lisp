(in-package :cando)

(defstruct (progress-bar
             (:constructor make-progress-bar
                           (&key total (divisions 10.0)
                                 (message "Completed")
                                 (advance-every-n-calls 1000) ;; Only advance every 1000 calls to progress-advance
                                 (on t)
                                 &aux (spacing (if (>= divisions 0)
                                                   (/ total divisions)
                                                   (error "divisions must be >0"))))))
  on total (advance-counter 1) advance-every-n-calls
  (start-time (get-universal-time))
  divisions message spacing (next 0))

(defun progress-convenient-time (time)
  (block nil
    (when (= time 0) (return nil))
    (when (< time 60) (return (format nil "~5,2f seconds" time)))
    (setf time (/ time 60.0))
    (when (< time 60) (return (format nil "~5,2f minutes" time)))
    (setf time (/ time 60.0))
    (when (< time 24) (return (format nil "~5,2f hours" time)))
    (setf time (/ time 24.0))
    (format nil "~5,2f days" time)))

(defun progress-advance (bar counter)
  (if (<= (decf (progress-bar-advance-counter bar)) 0)
    (setf (progress-bar-advance-counter bar) (progress-bar-advance-every-n-calls bar))
    (when (and (progress-bar-on bar) (> counter (progress-bar-next bar)))
      (let* ((fraction (/ counter (progress-bar-total bar)))
             (elapsed-time (- (get-universal-time) (progress-bar-start-time bar)))
             (expected-time (* (/ 1.0 fraction) elapsed-time))
             (remaining-time (* (- 1.0 fraction) expected-time)))
        (format t "~&~a ~3,1f%~@[ elapsed: ~a~]~@[ est.remaining: ~a~]~%"
                (progress-bar-message bar)
                (floor (* fraction 100.0))
                (progress-convenient-time elapsed-time)
                (progress-convenient-time remaining-time))
        (finish-output t)
        (let* ((new-next (* (1+ (floor (/ counter (progress-bar-spacing bar)))) (progress-bar-spacing bar))))
          (setf (progress-bar-next bar) new-next))))))

(defun progress-done (bar)
  (setf (progress-bar-advance-counter bar) 0)
  (progress-advance bar (progress-bar-total bar)))
