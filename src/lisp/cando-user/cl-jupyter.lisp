(in-package :cando-user)

(defun show (matter)
  (cond
    ((typep matter 'chem:aggregate))
    ((typep matter 'chem:molecule)
     (let ((agg (chem:make-aggregate nil)))
       (chem:add-matter agg matter)
       (setf matter agg)))
    ((typep matter 'chem:residue)
     (let ((agg (chem:make-aggregate nil))
           (mol (chem:make-molecule nil)))
       (chem:add-matter agg mol)
       (chem:add-matter mol matter)
       (setf matter agg)))
    (t (error "You cannot show a ~a" matter)))
  (cando:show-aggregate matter))

(defun repr (widget representation &optional (selection "all"))
  (funcall (find-symbol "ADD-REPRESENTATION" :nglv) widget representation :selection selection))

(defun cl-jupyter-kernel-start (&optional connection-file-name)
  (let ((amber-home
          (namestring (or (if (ext:getenv "AMBERHOME")
                              (probe-file (ext:getenv "AMBERHOME"))
                              "/usr/local/amber/")
                          (probe-file "/home/app/amber16-data/")
                          "/dev/null"))))
    (setf (logical-pathname-translations "amber")
          (list (list "**;*.*" (concatenate 'string amber-home "/**/*.*"))))
    (format t "Setting amber host pathname translation -> ~a~%" amber-home))
  (let ((*readtable* (copy-readtable)))
    (set-macro-character #\~ #'leap.core:leap-lookup-variable-reader-macro)
    (let ((cl-jup (find-symbol "KERNEL-START" "CL-JUPYTER")))
      (if cl-jup
          (if connection-file-name
              (funcall cl-jup connection-file-name)
              (funcall cl-jup))
          (error "cl-jupyter is not installed")))))

(defun jupyterlab-fork-server (&optional (server-info-directory #P"/tmp/clasp-fork-server/"))
  (let ((listen (make-instance 'sb-bsd-sockets:inet-socket
                               :type :stream
                               :protocol :tcp)))
    (setf (sb-bsd-sockets:sockopt-reuse-address listen) t)
    (sb-bsd-sockets:socket-bind listen
                                (sb-bsd-sockets:make-inet-address "127.0.0.1")
                                0)
    (multiple-value-bind (address port)
        (sb-bsd-sockets:socket-name listen)
      (ensure-directories-exist server-info-directory)
      (let ((pidfile (make-pathname :name "pidfile" :defaults server-info-directory))
            (portfile (make-pathname :name "portfile" :defaults server-info-directory)))
        (with-open-file (fout pidfile :direction :output :if-exists :supersede)
          (format fout "~d~%" (core:getpid)))
        (with-open-file (fout portfile :direction :output :if-exists :supersede)
          (format fout "~d~%" port))))
    ;; Create a connection queue, ignore child exit details and wait for clients.
    (sb-bsd-sockets:socket-listen listen 5) ; listen(server_sockfd, 5);
    (unwind-protect
         (loop
           (format t "Server waiting~%")
           ;; Accept connection.
           (let* ((socket (sb-bsd-sockets:socket-accept listen))
                  (stream (sb-bsd-sockets:socket-make-stream socket
                                                             :input t
                                                             :output t
                                                             :buffering :none
                                                             :external-format :default
                                                             :element-type :default)))
             (let ((connection-file-name (read-line stream)))
               (format t "Received: ~s~%" connection-file-name)
               (finish-output)
               (multiple-value-bind (maybe-error pid child-stream)
                   (core:fork)
                 (if (= pid 0)
                     (let* ((child-directory (make-pathname :directory (list :relative (format nil "~d" (core:getpid))) :defaults server-info-directory))
                            (child-connection-filename (make-pathname :name connection-file-name :defaults child-directory)))
                       (format stream "~d~%" (core:getpid))
                       (finish-output stream)
                       (with-open-file (fout child-connection-filename :direction :output :if-exists :supersede)
                         (format fout "~a~%" connection-file-name))
                       (with-open-file (fout (make-pathname :name "pidfile" :defaults child-directory) :direction :output :if-exists :supersede)
                         (format fout "~d~%" (core:getpid)))
                       (let ((client-pid (parse-integer (read-line stream))))
                         (with-open-file (fout (make-pathname :name "client-pidfile" :defaults child-directory) :direction :output :if-exists :supersede)
                           (format fout "~d~%" client-pid))
                         (close stream)
                         (format t "Starting cl-jupyter-kernel-start~%")
                         (sleep 10)
                         ;; (cl-jupyter-kernel-start connection-file-name)
                         (sb-bsd-sockets:socket-close listen)
                         (format t "Child ~a is exiting - client ~a should be killed as well...~%" (core:getpid) client-pid)
                         (core:exit)))
                     (progn             ; parent
                       (close stream)))))))
      (sb-bsd-sockets:socket-close listen))))
