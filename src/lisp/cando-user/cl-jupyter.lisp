(in-package :cando-user)


(defun jupyterlab-fork-server (&optional (server-info-directory #P"/tmp/clasp-fork-server/"))
  (let ((server-info-directory (pathname server-info-directory))
        (listen (make-instance 'sb-bsd-sockets:inet-socket
                               :type :stream
                               :protocol :tcp)))
    (setf (sb-bsd-sockets:sockopt-reuse-address listen) t)
    (sb-bsd-sockets:socket-bind listen
                                (sb-bsd-sockets:make-inet-address "127.0.0.1")
                                0)
    (multiple-value-bind (address port)
        (sb-bsd-sockets:socket-name listen)
      (ensure-directories-exist server-info-directory)
      (let ((pidfile (make-pathname :name "pidfile":defaults server-info-directory))
            (portfile (make-pathname :name "portfile" :defaults server-info-directory)))
        (with-open-file (fout pidfile :direction :output :if-exists :supersede)
          (format fout "~d~%" (core:getpid)))
        (with-open-file (fout portfile :direction :output :if-exists :supersede)
          (format fout "~d~%" port))))
    ;; Create a connection queue, ignore child exit details and wait for clients.
    (sb-bsd-sockets:socket-listen listen 5) ; listen(server_sockfd, 5);
    (unwind-protect
         (loop
           (format t "Server waiting - server-info-directory: ~s~%" server-info-directory)
           ;; Accept connection.
           (let* ((socket (sb-bsd-sockets:socket-accept listen))
                  (stream (sb-bsd-sockets:socket-make-stream socket
                                                             :input t
                                                             :output t
                                                             :buffering :none
                                                             :external-format :default
                                                             :element-type :default)))
             (let ((connection-filename (read-line stream)))
               (format t "Received connection-filename: ~s~%" connection-filename)
               (finish-output)
               (multiple-value-bind (maybe-error pid child-stream)
                   (core:fork)
                 (if (= pid 0)
                     (progn
                       (format t "child starting  connection-filename: ~s~%" connection-filename)
                       (format t "child server-info-directory: ~s~%" server-info-directory)
                       (let* ((child-directory (merge-pathnames (make-pathname :directory (list :relative (format nil "~d" (core:getpid)))) server-info-directory))
                              (child-connection-filename (make-pathname :name "child-connection-filename" #|connection-filename|# :defaults child-directory)))
                         (format t "child-directory: ~a~%" child-directory)
                         (format t "child-connection-filename: ~a~%" child-connection-filename)
                         (format stream "~d~%" (core:getpid))
                         (finish-output stream)
                         (ensure-directories-exist child-connection-filename)
                         (with-open-file (fout child-connection-filename :direction :output :if-exists :supersede)
                           (format fout "~a~%" connection-filename))
                         (with-open-file (fout (make-pathname :name "pidfile" :defaults child-directory) :direction :output :if-exists :supersede)
                           (format fout "~d~%" (core:getpid)))
                         (let* ((client-current-working-dir (read-line stream))
                                (_ (format t "Read client-current-working-dir: ~a~%" client-current-working-dir))
                                (client-pid-str (read-line stream))
                                (_ (format t "Read client-pid-str: ~a~%" client-pid-str))
                                (client-pid (parse-integer client-pid-str)))
                           (format t "Read client-pid: ~a~%" client-pid)
                           (with-open-file (fout (make-pathname :name "client-pidfile" :defaults child-directory) :direction :output :if-exists :supersede)
                             (format fout "~d~%" client-pid))
                           (close stream)
                           (format t "Starting jupyter:run-kernel~%")
                           (let ((jupyter-symbol (find-symbol "RUN-KERNEL" (find-package "JUPYTER")))
                                 (kernel (find-symbol "KERNEL" (find-package "CANDO-JUPYTER"))))
                             (funcall jupyter-symbol kernel connection-filename))
                           (sb-bsd-sockets:socket-close listen)
                           (format t "Child ~a is exiting - client ~a should be killed as well...~%" (core:getpid) client-pid)
                           (core:c_exit))))
                     (progn             ; parent
                       (close stream)))))))
      (sb-bsd-sockets:socket-close listen))))

