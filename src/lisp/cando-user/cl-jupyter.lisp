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

(defun jupyterlab-fork-server (&optional (port 9734))
  (let ((listen (make-instance 'sb-bsd-sockets:inet-socket
                               :type :stream
                               :protocol :tcp)))
    (setf (sb-bsd-sockets:sockopt-reuse-address listen) t)
    (sb-bsd-sockets:socket-bind listen
                                (sb-bsd-sockets:make-inet-address "127.0.0.1")
                                port)
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
                     (progn ; child
                       (prin1 (core:getpid) stream)
                       (terpri stream)
                       (finish-output stream)
                       (let ((client-pid (parse-integer (read-line stream))))
                         (format t "client-pid: ~s~%" client-pid)
                         (close stream)
                         (format t "Starting cl-jupyter-kernel-start~%")
                         (cl-jupyter-kernel-start connection-file-name))
                       (sb-bsd-sockets:socket-close listen)
                       (format t "Child ~a is exiting - client ~a should be killed as well...~%" (core:getpid) client-pid)
                       (core:exit))
                     (progn ; parent
                       (close stream)))))))
      (sb-bsd-sockets:socket-close listen))))
