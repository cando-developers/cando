;;;; xtal.lisp

(in-package #:xtal)

(defclass spacegroup ()
  ((num :initarg :num :reader num)
   (name :initarg :name :reader name)
   (equivalent-positions :initarg :equivalent-positions :reader equivalent-positions)
   (equivalent-position-strings :initarg :equivalent-position-strings :reader equivalent-position-strings)
   ))

(defmethod print-object ((spacegroup spacegroup) stream)
  (print-unreadable-object (spacegroup stream :type t)
    (format stream ":num ~d :name ~s" (num spacegroup) (name spacegroup))))

(defclass spacegroups ()
  ((spacegroup-by-num :initform (make-hash-table :test 'eql) :reader spacegroup-by-num)
   (spacegroup-by-name :initform (make-hash-table :test 'equal) :reader spacegroup-by-name)))

(defun parse-spacegroup-header (line)
  "Return two values: (values number name-string). Signals an error if malformed."
  (labels ((skip-ws (s i)
             (loop while (and (< i (length s))
                              (find (char s i) " \t"))
                   do (incf i)
                   finally (return i))))
    (let* ((s (string-trim " \t" line))
           (prefix "#")
           (plen (length prefix)))
      (unless (and (>= (length s) plen)
                   (string-equal prefix s :end2 plen))
        (error "Header must start with '#': ~A" line))
      (let ((i (skip-ws s plen)))
        (multiple-value-bind (num end)
            (parse-integer s :start i :junk-allowed t)
          (unless num
            (error "Missing spacegroup number in: ~A" line))
          (setf end (skip-ws s end))
          (let ((name (string-trim " \t" (subseq s end))))
            (when (zerop (length name))
              (error "Missing spacegroup name in: ~A" line))
            (values num name)))))))

(defun add-spacegroup (spacegroups header equivalent-position-strings equivalent-positions)
  (multiple-value-bind (num name)
      (parse-spacegroup-header header)
    (let ((sg (make-instance 'spacegroup
                             :num num
                             :name name
                             :equivalent-position-strings equivalent-position-strings
                             :equivalent-positions equivalent-positions)))
      (when (gethash name (spacegroup-by-name spacegroups))
        (error "The spacegroup name ~s has already been defined" name))
      (when (gethash num (spacegroup-by-num spacegroups))
        (error "The spacegroup num ~s has already been defined" num))
      (setf (gethash num (spacegroup-by-num spacegroups)) sg
            (gethash name (spacegroup-by-name spacegroups)) sg))))


(defun read-spacegroup-file (pathname)
  "Read PATHNAME and, for each block, call (foo header positions).

HEADER is the literal header line (e.g., \"spacegroup 5 C2\").
POSITIONS is a list of lines between the header and 'end'.
If POSITION-TRANSFORM is provided, it is mapped over each position
string before calling FOO (useful to parse to symbols/matrices)."
  (let ((spacegroups (make-instance 'spacegroups)))
    (labels ((trim (s)
               (string-trim '(#\Space #\Tab #\Return #\Linefeed) s))
             (starts-with (s prefix)
               (and (<= (length prefix) (length s))
                    (string-equal prefix s :end2 (length prefix)))))
      (with-open-file (in pathname :direction :input)
        (loop with header = nil
              with acc = nil
              for raw = (read-line in nil :eof)
              while (not (eq raw :eof))
              for line = (trim raw)
              do (let ((line (trim raw)))
                   (cond
                     ;; start of a block
                     ((starts-with line "#")
                      (setf header line
                            acc '()))
                     ;; end of a block
                     ((string-equal line "")
                      (when header
                        (let* ((equivalent-position-strings (nreverse acc))
                               (equivalent-positions (mapcar (lambda (str) (esrap:parse 'xtal.symop-parser:symop str)) equivalent-position-strings)))
                          (add-spacegroup spacegroups header equivalent-position-strings equivalent-positions)))
                      (setf header nil acc nil))
                     ;; accumulate positions within a block
                     (t
                      (when header
                        (push line acc)))))
              finally
                 (when header
                   (error "File ended before 'end' for header: ~A" header)))))
    spacegroups))
