(in-package #:yampi)

(defconstant +zmq-poll-timeout+ 500)

(defun more-parts-p (msg)
  (not (zerop (pzmq::%msg-more msg))))

(defgeneric start (channel connection-path &rest initargs &key))

(defgeneric stop (channel))

(defgeneric send (channel identity code &rest parts))

(defgeneric receive (channel identity code &rest parts))

(defgeneric idle (channel)
  (:method (channel)
    (declare (ignore channel))))

(defgeneric initialize (channel)
  (:method (channel)
    (declare (ignore channel))))

(defgeneric serialize (channel position part)
  (:method (channel position part)
    (declare (ignore channel position))
    part)
  (:method (channel position (part symbol))
    (declare (ignore channel position))
    (symbol-name part)))

(defgeneric translate (channel position part length)
  (:method (channel position part length)
    (declare (ignore channel position))
    (handler-case
        (cffi:foreign-string-to-lisp part
                                     :count length
                                     :encoding :utf-8)
      (babel-encodings:character-decoding-error ()
        "")))
  (:method (channel (position (eql nil)) part length)
    (cffi:foreign-array-to-lisp part
                                (list :array :uint8 length)
                                :element-type '(unsigned-byte 8))))

(defgeneric deserialize (channel position part)
  (:method (channel position part)
    (declare (ignore channel position))
    part))

(defun read-part (socket msg channel position)
  (pzmq:msg-recv msg socket)
  (deserialize channel position
               (translate channel position
                          (pzmq:msg-data msg)
                          (pzmq:msg-size msg))))

(defclass channel ()
  ((context :reader context
            :initform (pzmq:ctx-new))
   (control :accessor control)
   (control-send-lock :reader control-send-lock
                      :initform (bordeaux-threads:make-lock))
   (broadcast :accessor broadcast)
   (broadcast-send-lock :reader broadcast-send-lock
                        :initform (bordeaux-threads:make-lock))
   (thread :accessor thread
           :initform nil)))

(defmethod stop ((channel channel))
  (if (thread channel)
      (bordeaux-threads:interrupt-thread
       (thread channel)
       (lambda ()
         (throw 'shutdown nil)))
      (throw 'shutdown nil)))

(defclass server (channel) ())

(defmethod initialize-instance :after ((instance server) &rest initargs &key)
  (declare (ignore initargs))
  (setf (control instance) (pzmq:socket (context instance) :router)
        (broadcast instance) (pzmq:socket (context instance) :pub))
  instance)

(defun server-message-loop (channel)
  (catch 'shutdown
    (with-accessors ((control control)
                     (broadcast broadcast))
        channel
      (unwind-protect
           (pzmq:with-poll-items items ((control :pollin))
             (loop for poll = (pzmq:poll items +zmq-poll-timeout+)
                   unless (zerop poll)
                     do (pzmq:with-message msg
                          (apply #'receive channel (read-part control msg channel nil)
                                 (loop for position from 0
                                       while (more-parts-p msg)
                                       collect (read-part control msg channel position))))
                   do (idle channel)))
        (pzmq:close control)
        (pzmq:close broadcast)
        (pzmq:ctx-destroy (context channel))))))

(defmethod start ((channel server) connection-path &rest initargs
                  &key threaded control-endpoint broadcast-endpoint)
  (declare (ignore initargs))
  (with-accessors ((control control)
                   (broadcast broadcast)
                   (thread thread))
      channel
    (pzmq:bind control control-endpoint)
    (pzmq:bind broadcast broadcast-endpoint)
    (with-open-file (stream connection-path :direction :output
                                            :if-does-not-exist :create
                                            :if-exists :supersede)
      (with-standard-io-syntax
        (write `(:control ,(pzmq:getsockopt control :last-endpoint)
                 :broadcast ,(pzmq:getsockopt broadcast :last-endpoint))
               :stream stream)))
    (initialize channel)
    (if threaded
        (setf thread
              (bordeaux-threads:make-thread
               (lambda ()
                 (handler-case
                     (server-message-loop channel)
                   (error (condition)
                     (format t "~a" condition))))))
        (server-message-loop channel))))

(defmethod send ((channel server) (identity null) code &rest parts)
  (bordeaux-threads:with-lock-held ((broadcast-send-lock channel))
    (let ((broadcast (broadcast channel)))
      (pzmq:send broadcast (serialize channel 0 code) :sndmore (and parts t))
      (loop for (part . remaining) on parts
            for position from 1
            do (pzmq:send broadcast (serialize channel position part)
                          :sndmore (and remaining t))))))

(defmethod send ((channel server) identity code &rest parts)
  (bordeaux-threads:with-lock-held ((control-send-lock channel))
    (let ((control (control channel)))
      (pzmq:send control (serialize channel nil identity) :sndmore t)
      (pzmq:send control (serialize channel 0 code) :sndmore (and parts t))
      (loop for (part . remaining) on parts
            for position from 1
            do (pzmq:send control (serialize channel position part)
                          :sndmore (and remaining t))))))

(defclass client (channel) ())

(defmethod initialize-instance :after ((instance client) &rest initargs &key)
  (declare (ignore initargs))
  (setf (control instance) (pzmq:socket (context instance) :dealer)
        (broadcast instance) (pzmq:socket (context instance) :sub))
  instance)

(defun client-message-loop (channel)
  (catch 'shutdown
    (with-accessors ((control control)
                     (broadcast broadcast))
        channel
      (unwind-protect
           (flet ((recv (socket)
                    (pzmq:with-message msg
                      (apply #'receive channel nil
                             (loop for position from 0
                                   while (or (zerop position) (more-parts-p msg))
                                   collect (read-part socket msg channel position))))))
             (pzmq:with-poll-items items ((control :pollin) (broadcast :pollin))
               (loop for poll = (pzmq:poll items +zmq-poll-timeout+)
                     when (pzmq:revents items 0)
                       do (recv control)
                     when (pzmq:revents items 1)
                       do (recv broadcast)
                     do (idle channel))))
        (pzmq:close control)
        (pzmq:close broadcast)
        (pzmq:ctx-destroy (context channel))))))

(defmethod start ((channel client) connection-path &rest initargs &key threaded)
  (declare (ignore initargs))
  (with-accessors ((control control)
                   (broadcast broadcast)
                   (thread thread))
      channel
    (with-open-file (stream connection-path)
      (let ((data (with-standard-io-syntax (read stream nil nil))))
        (pzmq:connect control (getf data :control))
        (pzmq:connect broadcast (getf data :broadcast))))
    (initialize channel)
    (if threaded
        (setf thread
              (bordeaux-threads:make-thread
               (lambda ()
                 (handler-case
                     (client-message-loop channel)
                   (error (condition)
                     (format t "~a" condition))))))
        (client-message-loop channel))))

(defmethod send ((channel client) (identity null) code &rest parts)
  (bordeaux-threads:with-lock-held ((control-send-lock channel))
    (let ((control (control channel)))
      (pzmq:send control (serialize channel 0 code) :sndmore (and parts t))
      (loop for (part . remaining) on parts
            for position from 1
            do (pzmq:send control (serialize channel position part)
                          :sndmore (and remaining t))))))

(defun subscribe (channel code)
  (pzmq:setsockopt (broadcast channel)
                   :subscribe
                   (typecase code
                     (symbol (symbol-name code))
                     (string code)
                     (otherwise (prin1-to-string code)))))

(defun unsubscribe (channel code)
  (pzmq:setsockopt (broadcast channel)
                   :unsubscribe
                   (typecase code
                     (symbol (symbol-name code))
                     (string code)
                     (otherwise (prin1-to-string code)))))
