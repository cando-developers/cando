(in-package #:yampi)

(defconstant +zmq-poll-timeout+ 500)

(defun read-binary-part (socket msg)
  (pzmq:msg-recv msg socket)
  (cffi:foreign-array-to-lisp (pzmq:msg-data msg)
                              (list :array :uint8 (pzmq:msg-size msg))
                              :element-type '(unsigned-byte 8)))

(defun read-string-part (socket msg)
  (pzmq:msg-recv msg socket)
  (handler-case
      (cffi:foreign-string-to-lisp (pzmq:msg-data msg)
                                   :count (pzmq:msg-size msg)
                                   :encoding :utf-8)
    (babel-encodings:character-decoding-error ()
      "")))

(defun more-parts-p (msg)
  (not (zerop (pzmq::%msg-more msg))))

(defgeneric start (channel connection-path &rest initargs &key))

(defgeneric stop (channel))

(defgeneric send (channel identity code &rest parts))

(defgeneric receive (channel identity code &rest parts))

(defclass channel ()
  ((context :reader context
            :initform (pzmq:ctx-new))
   (control :accessor control)
   (control-send-lock :reader control-send-lock
                      :initform (bordeaux-threads:make-lock))
   (broadcast :accessor broadcast)
   (broadcast-send-lock :reader broadcast-send-lock
                        :initform (bordeaux-threads:make-lock))
   (thread :accessor thread)))

(defmethod stop ((channel channel))
  (bordeaux-threads:interrupt-thread (thread channel)
                                     (lambda ()
                                       (throw 'shutdown nil)))
  (pzmq:close (control channel))
  (pzmq:close (broadcast channel))
  (pzmq:ctx-destroy (context channel)))

(defclass server (channel) ())

(defmethod initialize-instance :after ((instance server) &rest initargs &key)
  (declare (ignore initargs))
  (setf (control instance) (pzmq:socket (context instance) :router)
        (broadcast instance) (pzmq:socket (context instance) :pub))
  instance)

(defmethod start ((channel server) connection-path &rest initargs
                  &key control-endpoint broadcast-endpoint)
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
    (setf thread
          (bordeaux-threads:make-thread
           (lambda ()
             (pzmq:with-poll-items items ((control :pollin))
               (catch 'shutdown
                 (loop with identity
                       with code
                       with parts
                       for poll = (pzmq:poll items +zmq-poll-timeout+)
                       unless (zerop poll)
                         do (pzmq:with-message msg
                              (setf identity (read-binary-part control msg))
                              (unless (more-parts-p msg)
                                (error "Need more parts"))
                              (setf code (read-binary-part control msg))
                              (unless (= (length code) 1)
                                (error "Code should only have length 1"))
                              (setf parts (loop while (more-parts-p msg)
                                                collect (read-string-part control msg)))
                              (apply #'receive channel identity (aref code 0) parts))))))))))

(defmethod send ((channel server) (identity null) code &rest parts)
  (bordeaux-threads:with-lock-held ((broadcast-send-lock channel))
    (let ((broadcast (broadcast channel)))
      (pzmq:send broadcast
                 (make-array 1 :initial-element code :element-type '(unsigned-byte 8))
                 :sndmore (and parts t))
      (loop for (part . remaining) on parts
            do (pzmq:send broadcast part
                          :sndmore (and remaining t))))))

(defmethod send ((channel server) identity code &rest parts)
  (bordeaux-threads:with-lock-held ((control-send-lock channel))
    (let ((control (control channel)))
      (pzmq:send control identity :sndmore t)
      (pzmq:send control
                 (make-array 1 :initial-element code :element-type '(unsigned-byte 8))
                 :sndmore (and parts t))
      (loop for (part . remaining) on parts
            do (pzmq:send control part
                          :sndmore (and remaining t))))))

(defclass client (channel) ())

(defmethod initialize-instance :after ((instance client) &rest initargs &key)
  (declare (ignore initargs))
  (setf (control instance) (pzmq:socket (context instance) :dealer)
        (broadcast instance) (pzmq:socket (context instance) :sub))
  instance)

(defmethod start ((channel client) connection-path &rest initargs &key)
  (declare (ignore initargs))
  (with-accessors ((control control)
                   (broadcast broadcast)
                   (thread thread))
      channel
    (with-open-file (stream connection-path)
      (let ((data (with-standard-io-syntax (read stream nil nil))))
        (pzmq:connect control (getf data :control))
        (pzmq:connect broadcast (getf data :broadcast))))
    (setf thread
          (bordeaux-threads:make-thread
           (lambda ()
             (flet ((recv (socket)
                      (pzmq:with-message msg
                        (setf code (read-binary-part socket msg))
                        (unless (= (length code) 1)
                          (error "Code should only have length 1"))
                        (setf parts (loop while (more-parts-p msg)
                                          collect (read-string-part socket msg)))
                        (apply #'receive channel nil (aref code 0) parts))))
               (pzmq:with-poll-items items ((control :pollin) (broadcast :pollin))
                 (catch 'shutdown
                   (loop with code
                         with parts
                         for poll = (pzmq:poll items +zmq-poll-timeout+)
                         when (pzmq:revents items 0)
                           do (recv control)
                         when (pzmq:revents items 1)
                           do (recv broadcast))))))))))

(defmethod send ((channel client) (identity null) code &rest parts)
  (bordeaux-threads:with-lock-held ((control-send-lock channel))
    (let ((control (control channel)))
      (pzmq:send control
                 (make-array 1 :initial-element code :element-type '(unsigned-byte 8))
                 :sndmore (and parts t))
      (loop for (part . remaining) on parts
            do (pzmq:send control part
                          :sndmore (and remaining t))))))

(defun subscribe (channel code)
  (pzmq:setsockopt (broadcast channel)
                   :subscribe
                   (make-array 1 :initial-element code :element-type '(unsigned-byte 8))))

(defun unsubscribe (channel code)
  (pzmq:setsockopt (broadcast channel)
                   :unsubscribe
                   (make-array 1 :initial-element code :element-type '(unsigned-byte 8))))
