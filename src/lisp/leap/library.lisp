
(in-package :leap.off)


(defclass off-header ()
  ((name :initarg :name :accessor name)
   (type :initarg :type :accessor type)
   (fields :initarg :fields :accessor fields)))

(define-condition bad-header ()
  ((line :initarg :line :accessor line)))

(defun parse-off-header (line)
  (unless (char= (elt line 0) #\!)
    (error 'bad-header :line line))
  (let* ((parts (split-sequence:split-sequence-if (lambda (ch)
                                                    (or (char= #\Space ch)
                                                        (char= #\Tab ch)))
                                                  (subseq line 1)
                                                  :remove-empty-subseqs t))
         (name (elt parts 0))
         (type (intern (string-upcase (nth 1 parts)) :keyword))
         (fields (cond
                   ((eq type :table)
                    (loop for index from 2 below (length parts) by 2
                       collect (cons (intern (string-upcase (elt parts index)) :keyword)
                                     (intern (string-upcase (elt parts (1+ index))) :keyword))))
                   ((eq type :array)
                    (list (cons (intern (string-upcase (elt parts 2)) :keyword)
                                :only)))
                   ((eq type :single)
                    (list (cons (intern (string-upcase (elt parts 2)) :keyword)
                                :only)))
                   (t (error "Unknown OFF file entry type: ~a" type)))))
    (make-instance 'off-header
                   :name name
                   :type type
                   :fields fields)))

(defun parse-off-data-line (line header)
  (let* ((data (make-array (length (fields header)) :element-type t)))
    (with-input-from-string (sin line)
      (loop for x from 0 below (length (fields header))
         do (setf (elt data x) (read sin))))
    data))


(defun read-off-data-line (fin header error-eof-p eof)
  (let ((first-char (read-char fin error-eof-p eof)))
    (cond
      ((eq first-char eof) eof)
      ((eq first-char #\!)
       (unread-char first-char fin)
       nil)
      (t (parse-off-data-line (read-line fin error-eof-p eof) header)))))

(defun read-off-data-block (fin &optional (eof-error-p t) eof)
  (let ((header-line (read-line fin eof-error-p eof)))
    (if (eq header-line eof)
        eof
        (let ((header (parse-off-header header-line)))
          (when (eq header eof)
            (return-from read-off-data-block eof))
          (let ((data (make-array 0 :element-type t :fill-pointer t)))
            (loop for x = (read-off-data-line fin header nil :eof)
               until (or (null x) (eq x :eof))
               do (vector-push-extend x data))
            (values header data))))))

      
