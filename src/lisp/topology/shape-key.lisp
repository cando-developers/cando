(in-package #:topology)

(defun binned-shape-key-proper-list-p (value)
  (and (listp value)
       (handler-case (integerp (list-length value))
         (type-error () nil))))

(defun shape-key-name-less-p (left right)
  "Canonical order: symbol name, then home package name."
  (or (string< (symbol-name left) (symbol-name right))
      (and (string= (symbol-name left) (symbol-name right))
           (string< (package-name (symbol-package left))
                    (package-name (symbol-package right))))))

(defun canonical-shape-key-entries (entries)
  "Validate, copy, and sort entries without modifying the caller's data."
  (unless (binned-shape-key-proper-list-p entries)
    (error "Shape-key entries must be a proper list"))
  (let ((seen (make-hash-table :test #'eq)))
    (dolist (entry entries)
      (unless (and (consp entry) (car entry) (symbolp (car entry))
                   (symbol-package (car entry))
                   (binned-shape-key-proper-list-p (cdr entry))
                   (<= 1 (length (cdr entry)) 3)
                   (every #'integerp (cdr entry)))
        (let ((*print-circle* t))
          (error "Invalid shape-key entry ~s" entry)))
      (multiple-value-bind (previous presentp) (gethash (car entry) seen)
        (when presentp
          (error "Duplicate shape-key name ~s: existing ~s, incoming ~s"
                 (car entry) previous (cdr entry))))
      (setf (gethash (car entry) seen) (cdr entry))))
  (sort (copy-tree entries) #'shape-key-name-less-p :key #'car))

(defclass binned-shape-key-cache (cando.serialize:serializable)
  ((entries :initarg :entries :initform nil))
  (:documentation "Canonical alist of names and integer-bin lists.
Use ADD-BINNED-SHAPE-KEY-ENTRY to modify; readers return copies."))

(defclass binned-shape-key (binned-shape-key-cache) ()
  (:documentation "Selected sidechain shape, distinct from a backbone cache."))

(defmethod initialize-instance :after ((key binned-shape-key-cache) &key)
  (setf (slot-value key 'entries)
        (canonical-shape-key-entries (slot-value key 'entries))))

(defun make-binned-shape-key-cache (entries)
  (make-instance 'binned-shape-key-cache :entries entries))

(defun make-binned-shape-key (entries)
  (make-instance 'binned-shape-key :entries entries))

(defmethod binned-shape-key-entries ((key binned-shape-key-cache))
  (copy-tree (slot-value key 'entries)))

(defmethod binned-shape-key-names ((key binned-shape-key-cache))
  (mapcar #'car (slot-value key 'entries)))

(defmethod binned-shape-key-components ((key binned-shape-key-cache))
  (mapcar (lambda (entry) (copy-list (cdr entry))) (slot-value key 'entries)))

(defun add-binned-shape-key-entry (key name bins)
  "Reject duplicates, including identical values, before changing KEY."
  (check-type key binned-shape-key-cache)
  (setf (slot-value key 'entries)
        (canonical-shape-key-entries
         (cons (cons name bins) (slot-value key 'entries))))
  key)

(defun binned-shape-key-component (key name)
  (let ((entry (assoc name (slot-value key 'entries) :test #'eq)))
    (unless entry
      (error "No shape-key component ~s; available names: ~s"
             name (binned-shape-key-names key)))
    (copy-list (cdr entry))))

(defun binned-shape-key= (left right)
  "Selected keys compare by their named bins, not object identity."
  (and (typep left 'binned-shape-key) (typep right 'binned-shape-key)
       (equal (slot-value left 'entries) (slot-value right 'entries))))

(defun binned-shape-key-pathname-part (key)
  "Encode canonical contents without depending on printer settings or package.
Names include their home package. Escape punctuation so names cannot collide
with the field, package, or bin separators. The result is one pathname component."
  (check-type key binned-shape-key)
  (with-output-to-string (stream)
    (labels ((write-name (name)
               (loop for char across name
                     do (if (or (char<= #\a char #\z) (char<= #\A char #\Z)
                                (char<= #\0 char #\9))
                            (write-char char stream)
                            (format stream "_~x_" (char-code char))))))
      (write-string "shape-key" stream)
      (dolist (entry (slot-value key 'entries))
        (write-string "--" stream)
        (write-name (package-name (symbol-package (car entry))))
        (write-char #\~ stream)
        (write-name (symbol-name (car entry)))
        (write-char #\= stream)
        (format stream "~{~d~^,~}" (cdr entry))))))

(defmethod print-object ((key binned-shape-key-cache) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (key stream :type t)
        (write (slot-value key 'entries) :stream stream))))
