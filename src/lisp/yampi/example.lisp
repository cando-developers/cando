(ql:quickload :yampi)

;;; Start up several terminals. Then in one terminal do
;;; cando --load example.lisp --eval "(defvar ch (make-server \"~/connect.sexp\"))"
;;; In each of the other terminals do
;;; cando --load example.lisp --eval "(defvar ch (make-client \"~/connect.sexp\"))"
;;; In the clients you can do things like (yampi:send ch nil :job-request)
;;; The server will respond with a job id.
;;; In the server you can also broadcast a stop message (yampi:send ch nil :stop)

(defclass server (yampi:server)
  ((next-job-id :accessor next-job-id
                :initform 0)))

(defmethod yampi:receive ((channel server) identity (code (eql :job-request)) &rest parts)
  (declare (ignore parts))
  (format *debug-io* "~&Received job-request from ~s~%" identity)
  (yampi:send channel identity code
              (format nil "~a" (incf (next-job-id channel)))))

(defmethod yampi:receive ((channel server) identity (code (eql :job-done)) &rest parts)
  (declare (ignore parts))
  (format *debug-io* "~&Received job-done from ~s~%" identity))

(defmethod yampi:receive ((channel server) identity (code (eql :job-abort)) &rest parts)
  (declare (ignore parts))
  (format *debug-io* "~&Received job-abort from ~s~%" identity))

(defun make-server (connection-path)
  (let ((channel (make-instance 'server))
        (endpoint "tcp://0.0.0.0:*" (machine-instance)))
    (yampi:start channel connection-path
                 :threaded t
                 :control-endpoint endpoint
                 :broadcast-endpoint endpoint)
    channel))

(defclass client (yampi:client) ())

(defmethod yampi:receive ((channel client) identity (code (eql :job-request)) &rest parts)
  (format *debug-io* "~&Received job-request reply ~s~%" parts))

(defmethod yampi:receive ((channel client) identity (code (eql :stop)) &rest parts)
  (declare (ignore parts))
  (format *debug-io* "~&Received stop~%")
  (yampi:stop channel))

(defun make-client (connection-path)
  (let ((channel (make-instance 'client)))
    (yampi:start channel connection-path :threaded t)
    (yampi:subscribe channel :stop)
    channel))

(defmethod yampi:receive (channel identity code &rest parts)
  (format *debug-io* "~&No handler for code ~a~%" code))

(defmethod yampi:deserialize (channel (position (eql 0)) part)
  (find part '(:job-request :job-done :job-abort :stop)
        :test #'equalp
        :key #'symbol-name))
