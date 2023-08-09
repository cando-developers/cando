(ql:quickload :yampi)

;;; Start up several terminals. Then in one terminal do
;;; cando --load example-2.lisp --eval "(defvar ch (make-server \"~/connect.sexp\"))"
;;; In each of the other terminals do
;;; cando --load example-2.lisp --eval "(defvar ch (make-client \"~/connect.sexp\"))"

(defclass server (yampi:server)
  ((remaining-jobs :accessor remaining-jobs
                   :initform (loop for i below 30
                                   collect i))
   (pending-jobs :accessor pending-jobs
                 :initform (make-hash-table))))

(defmethod yampi:receive ((channel server) identity (code (eql :request)) &rest parts)
  (declare (ignore parts))
  (with-accessors ((remaining-jobs remaining-jobs)
                   (pending-jobs pending-jobs))
      channel
    (cond (remaining-jobs
           (let ((job (pop remaining-jobs)))
             (setf (gethash job pending-jobs) nil)
             (yampi:send channel identity :assign
                         (format nil "~a" job))))
          ((zerop (hash-table-count pending-jobs))
                                        ; no pending jobs so tell everyone to quit
           (yampi:send channel nil :stop)
           (sleep 10)
           (yampi:stop channel)))))

(defmethod yampi:receive ((channel server) identity (code (eql :done)) &rest parts)
  (with-accessors ((remaining-jobs remaining-jobs)
                   (pending-jobs pending-jobs))
      channel
    (remhash (car parts) pending-jobs)
    (format *debug-io* "Job done ~2d, ~2d remaining, ~2d pending.~%"
            (car parts) (length remaining-jobs) (hash-table-count pending-jobs))
    (when (and (null remaining-jobs)
               (zerop (hash-table-count pending-jobs)))
      (yampi:send channel nil :stop)
      (sleep 10)
      (yampi:stop channel))))

(defmethod yampi:receive ((channel server) identity (code (eql :fail)) &rest parts)
  (with-accessors ((remaining-jobs remaining-jobs)
                   (pending-jobs pending-jobs))
      channel
    (remhash (car parts) pending-jobs)
    (format *debug-io* "Job fail ~2d, ~2d remaining, ~2d pending.~%"
            (car parts) (length remaining-jobs) (hash-table-count pending-jobs))
    (push (car parts) remaining-jobs)))

(defun make-server (connection-path)
  (let ((channel (make-instance 'server))
        (endpoint "tcp://0.0.0.0:*" (machine-instance)))
    (yampi:start channel connection-path
                 :threaded nil
                 :control-endpoint endpoint
                 :broadcast-endpoint endpoint)
    channel))

(defclass client (yampi:client)
  ((kernel-channel :accessor kernel-channel
                   :initform (lparallel:make-channel))))

(defmethod yampi:receive ((channel client) identity (code (eql :assign)) &rest parts)
  (lparallel:submit-task (kernel-channel channel)
                         (lambda (id &aux (status :fail))
                           (sleep (random 20))
                           (list (if (zerop (random 10)) :fail :done) id))
                         (car parts)))

(defmethod yampi:receive ((channel client) identity (code (eql :stop)) &rest parts)
  (declare (ignore parts))
  (yampi:stop channel))

(defmethod yampi:idle ((channel client))
  (multiple-value-bind (result validp)
      (lparallel:try-receive-result (kernel-channel channel) :timeout .1)
    (when validp
      (yampi:send channel nil (first result) (second result))
      (when (eql (car result) :done)
        (yampi:send channel nil :request)))))

(defun make-client (connection-path)
  (let ((channel (make-instance 'client)))
    (yampi:start channel connection-path :threaded nil)
    channel))

(defmethod yampi:initialize ((channel client))
  (yampi:subscribe channel :stop)
  (loop repeat (lparallel:kernel-worker-count)
        do (yampi:send channel nil :request)))

(defmethod yampi:receive (channel identity code &rest parts)
  (format *debug-io* "~&No handler for code ~a~%" code))

(defmethod yampi:deserialize (channel (position (eql 0)) (part string))
  (or (find-symbol part :keyword)
      part))

(defmethod yampi:deserialize (channel (position (eql 1)) (part string))
  (read-from-string part nil nil))

(defmethod yampi:serialize (channel (position (eql 1)) (part number))
  (write-to-string part))
