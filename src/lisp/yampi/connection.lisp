(in-package #:yampi)

(defgeneric write-connection-file (connection-path control broadcast backend)
  (:documentation "Publish bound CONTROL and BROADCAST endpoint strings.
BACKEND selects storage; T uses a local pathname. Return the configuration plist.
One server owns each connection path; stop it before reusing that path."))

(defgeneric wait-for-connection-file (connection-path backend &key poll-interval timeout)
  (:documentation "Wait for connection information and return CONTROL and BROADCAST
as two values. T selects filesystem storage. POLL-INTERVAL defaults to 10 seconds;
TIMEOUT defaults to NIL (no limit), or is a nonnegative number of seconds.
Only absence is retried; malformed configurations and storage errors propagate."))

(defgeneric delete-connection-file (connection-path backend)
  (:documentation "Remove server discovery information from BACKEND.
T uses a local pathname. A missing configuration is harmless. Called on shutdown;
custom storage implementations must supply this method along with read/write."))

(defun connection-endpoints (data)
  "Validate the endpoint plist without traversing potentially circular input."
  (unless (and (consp data) (consp (cdr data)) (consp (cddr data))
               (consp (cdddr data)) (null (cddddr data)))
    (error "Invalid connection configuration"))
  (destructuring-bind (key1 value1 key2 value2) data
    (unless (and (or (and (eq key1 :control) (eq key2 :broadcast))
                     (and (eq key1 :broadcast) (eq key2 :control)))
                 (stringp value1) (plusp (length value1))
                 (stringp value2) (plusp (length value2)))
      (error "Connection configuration requires :CONTROL and :BROADCAST endpoint strings")))
  (values (getf data :control) (getf data :broadcast)))

(defun call-with-connection-wait (reader connection-path poll-interval timeout)
  "READER returns a configuration or NIL only when it is absent."
  (check-type poll-interval (real (0)))
  (check-type timeout (or null (real 0)))
  (let ((start (get-internal-real-time)))
    (loop
      (let ((data (funcall reader)))
        (when data (return (connection-endpoints data))))
      (let ((remaining (when timeout
                         (- timeout (/ (- (get-internal-real-time) start)
                                       internal-time-units-per-second)))))
        (when (and remaining (<= remaining 0))
          (error "Timed out waiting for connection configuration ~s" connection-path))
        (format t "Waiting for ~s to appear~%" connection-path)
        (sleep (if remaining (min remaining poll-interval) poll-interval))))))

(defmethod write-connection-file (connection-path control broadcast (backend (eql t)))
  (let ((data (list :control control :broadcast broadcast))
        (destination (merge-pathnames connection-path)))
    (connection-endpoints data)
    (when (probe-file destination)
      (error "Connection file ~s already exists; stop the old server and remove it" destination))
    (ensure-directories-exist destination)
    (uiop:with-temporary-file
        (:pathname path :directory (uiop:pathname-directory-pathname destination)
         :prefix "yampi-connection-" :type "config")
      (with-open-file (stream path :direction :output :if-exists :supersede
                                  :external-format :utf-8)
        (with-standard-io-syntax (write data :stream stream)))
      (rename-file path destination))
    data))

(defmethod wait-for-connection-file (connection-path (backend (eql t))
                                    &key (poll-interval 10) timeout)
  (call-with-connection-wait
   (lambda ()
     (when (probe-file connection-path)
       (with-open-file (stream connection-path :external-format :utf-8)
         (with-standard-io-syntax
           (let ((*read-eval* nil) (eof (gensym "EOF")))
             (let ((data (read stream nil eof)))
               (when (eq data eof) (error "Empty connection configuration"))
               (unless (eq (read stream nil eof) eof)
                 (error "Trailing forms in connection configuration"))
               ;; NIL in an existing file is invalid, not absence.
               (connection-endpoints data)
               data))))))
   connection-path poll-interval timeout))

(defmethod delete-connection-file (connection-path (backend (eql t)))
  (when (probe-file connection-path) (delete-file connection-path)))
