(in-package #:rosetta.params)

(defclass rosetta-atom-type-properties ()
  ((name :initarg :name :accessor rosetta-atom-type-name)
   (element :initarg :element :accessor rosetta-atom-type-element)
   (lj-radius :initarg :lj-radius :accessor rosetta-atom-type-lj-radius)
   (lj-wdepth :initarg :lj-wdepth :accessor rosetta-atom-type-lj-wdepth)
   (lk-dgfree :initarg :lk-dgfree :accessor rosetta-atom-type-lk-dgfree)
   (lk-lambda :initarg :lk-lambda :accessor rosetta-atom-type-lk-lambda)
   (lk-volume :initarg :lk-volume :accessor rosetta-atom-type-lk-volume)
   (flags :initarg :flags :accessor rosetta-atom-type-flags)))

(defmethod print-object ((obj rosetta-atom-type-properties) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a lj(~,3f/~,3f) lk(~,3f/~,3f/~,3f)"
            (rosetta-atom-type-name obj)
            (rosetta-atom-type-element obj)
            (rosetta-atom-type-lj-radius obj)
            (rosetta-atom-type-lj-wdepth obj)
            (rosetta-atom-type-lk-dgfree obj)
            (rosetta-atom-type-lk-lambda obj)
            (rosetta-atom-type-lk-volume obj))
    (when (rosetta-atom-type-flags obj)
      (format stream " ~{~a~^ ~}" (rosetta-atom-type-flags obj)))))

(defclass rosetta-atom-type-table ()
  ((entries :initarg :entries :accessor rosetta-atom-type-table-entries)
   (by-name :initarg :by-name :accessor rosetta-atom-type-table-by-name)
   (source-path :initarg :source-path :accessor rosetta-atom-type-table-source-path)))

(defmethod print-object ((obj rosetta-atom-type-table) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~d entries" (length (rosetta-atom-type-table-entries obj)))
    (when (rosetta-atom-type-table-source-path obj)
      (format stream " ~a" (rosetta-atom-type-table-source-path obj)))))

(defun default-fa-standard-atom-properties-path ()
  (let ((home (probe-file #P"rosetta:")))
    (when home
      (merge-pathnames
       "database/chemical/atom_type_sets/fa_standard/atom_properties.txt"
       home))))

(defun scan-floats-in-token (token)
  (let ((len (length token))
        (pos 0)
        (results nil))
    (labels ((digitp (ch)
               (and ch (<= (char-code #\0) (char-code ch) (char-code #\9)))))
      (loop while (< pos len) do
        (let ((ch (char token pos)))
          (if (or (digitp ch) (char= ch #\.) (char= ch #\+) (char= ch #\-))
              (let ((start pos)
                    (seen-digit nil)
                    (seen-exp nil))
                (when (or (char= ch #\+) (char= ch #\-))
                  (incf pos))
                (loop while (< pos len) do
                  (let ((c (char token pos)))
                    (cond
                      ((digitp c)
                       (setf seen-digit t)
                       (incf pos))
                      ((char= c #\.)
                       (incf pos))
                      ((and (not seen-exp) (or (char= c #\e) (char= c #\E)))
                       (setf seen-exp t)
                       (incf pos)
                       (when (< pos len)
                         (let ((next (char token pos)))
                           (when (or (char= next #\+) (char= next #\-))
                             (incf pos)))))
                      (t (return)))))
                (let ((frag (subseq token start pos)))
                  (when seen-digit
                    (push (parse-float frag) results))))
              (incf pos)))))
    (nreverse results)))

(defun parse-atom-properties-line (line)
  (let* ((clean (string-trim '(#\Space #\Tab #\Newline #\Return)
                             (strip-comment line))))
    (when (> (length clean) 0)
      (let ((fields (split-fields clean)))
        (when (and fields (not (string= (string-upcase (first fields)) "NAME")))
          (let* ((name (make-keyword (first fields)))
                 (element (make-keyword (second fields)))
                 (numbers nil)
                 (num-count 0)
                 (flags nil))
            (dolist (tok (cddr fields))
              (if (< num-count 5)
                  (let ((vals (scan-floats-in-token tok)))
                    (if vals
                        (dolist (v vals)
                          (when (< num-count 5)
                            (push v numbers)
                            (incf num-count)))
                        (push tok flags)))
                  (push tok flags)))
            (setf numbers (nreverse numbers))
            (setf flags (nreverse flags))
            (when (< (length numbers) 5)
              (error "Could not parse 5 numeric fields from atom_properties line: ~a" line))
            (make-instance 'rosetta-atom-type-properties
                           :name name
                           :element element
                           :lj-radius (first numbers)
                           :lj-wdepth (second numbers)
                           :lk-dgfree (third numbers)
                           :lk-lambda (fourth numbers)
                           :lk-volume (fifth numbers)
                           :flags (mapcar #'make-keyword flags))))))))

(defun parse-atom-properties-file (pathname)
  (let ((entries nil)
        (by-name (make-hash-table :test #'eq)))
    (with-open-file (fin pathname)
      (loop for line = (read-line fin nil nil)
            while line do
              (let ((entry (parse-atom-properties-line line)))
                (when entry
                  (push entry entries)
                  (setf (gethash (rosetta-atom-type-name entry) by-name) entry)))))
    (make-instance 'rosetta-atom-type-table
                   :entries (nreverse entries)
                   :by-name by-name
                   :source-path pathname)))

(defun load-lk-solvation-force-field (&key (pathname (default-fa-standard-atom-properties-path))
                                        (type-prefix nil))
  "Return a ForceField_sp with LK solvation params loaded from atom_properties.txt.
  :type-prefix when non-nil is prepended to the atom type name (string), e.g. \"lk-solvation-\"."
  (unless pathname
    (error "No atom_properties.txt pathname available"))
  (let* ((table (parse-atom-properties-file pathname))
         (ff (chem:force-field/make))
         (lkdb (chem:get-lksolvation-db ff)))
    (dolist (entry (rosetta-atom-type-table-entries table))
      (let* ((raw-name (rosetta-atom-type-name entry))
             (type (if type-prefix
                       (make-keyword (format nil "~a~a" type-prefix (symbol-name raw-name)))
                       raw-name))
             (term (chem:make-fflksolvation
                    type
                    :lj-radius (rosetta-atom-type-lj-radius entry)
                    :lj-wdepth (rosetta-atom-type-lj-wdepth entry)
                    :lk-dgfree (rosetta-atom-type-lk-dgfree entry)
                    :lk-lambda (rosetta-atom-type-lk-lambda entry)
                    :lk-volume (rosetta-atom-type-lk-volume entry))))
        (chem:fflksolvation-db-add lkdb term)))
    ff))
