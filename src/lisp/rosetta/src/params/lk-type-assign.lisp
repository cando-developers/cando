(in-package #:rosetta.params)

(defclass rosetta-lk-type-assignment ()
  ((lk-type :initarg :lk-type :accessor rosetta-lk-type-assignment-type)
   (smarts :initarg :smarts :accessor rosetta-lk-type-assignment-smarts)
   (compiled-smarts :initarg :compiled-smarts :accessor rosetta-lk-type-assignment-compiled-smarts)))

(defmethod print-object ((obj rosetta-lk-type-assignment) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a"
            (rosetta-lk-type-assignment-type obj)
            (rosetta-lk-type-assignment-smarts obj))))

(defun default-lk-type-assignments-path ()
  (let* ((base (probe-file (or (asdf:system-source-directory :rosetta-params)
               *load-pathname* *compile-file-pathname*)))
         (dir (and base (uiop:pathname-directory-pathname base))))
    (when dir
      (merge-pathnames "data/types.lisp" dir))))

(defun parse-lk-type-assignments-form (form)
  (unless (listp form)
    (error "Expected a list for LK type assignment form, got ~s" form))
  (let ((table (make-hash-table :test #'eq)))
    (dolist (group form)
      (unless (and (consp group) (symbolp (first group)))
        (error "Invalid LK type assignment group: ~s" group))
      (let ((element (first group))
            (assignments nil))
        (dolist (entry (rest group))
          (unless (and (consp entry) (= (length entry) 2))
            (error "Invalid LK type assignment entry: ~s" entry))
          (let ((lk-type (first entry))
                (smarts (second entry)))
            (unless (symbolp lk-type)
              (error "LK type must be a symbol, got ~s" lk-type))
            (unless (stringp smarts)
              (error "SMARTS must be a string, got ~s" smarts))
            (push (make-instance 'rosetta-lk-type-assignment
                                 :lk-type lk-type
                                 :smarts smarts
                                 :compiled-smarts (chem:compile-smarts smarts))
                  assignments)))
        (let ((existing (gethash element table)))
          (setf (gethash element table)
                (if existing
                    (nconc existing (nreverse assignments))
                    (nreverse assignments))))))
    table))

(defun load-lk-type-assignments (&key (pathname (default-lk-type-assignments-path)))
  (unless pathname
    (error "No LK type assignments path available"))
  (with-open-file (fin pathname)
    (let ((form (read fin nil :eof)))
      (when (eq form :eof)
        (error "LK type assignments file ~a is empty" pathname))
      (parse-lk-type-assignments-form form))))

(defun lookup-lk-type-assignments (element lk-type-assignments)
  (or (gethash element lk-type-assignments)
      (let ((upper (intern (string-upcase (symbol-name element)) :keyword)))
        (gethash upper lk-type-assignments))
      (let ((target (symbol-name element))
            (found nil))
        (maphash (lambda (key val)
                   (when (string-equal (symbol-name key) target)
                     (setf found val)))
                 lk-type-assignments)
        found)))

(define-condition missing-lk-type (chem:parameterization-error)
  ((molecule :initarg :molecule :accessor molecule)
   (residue :initarg :residue :accessor residue)
   (a1-name :initarg :a1-name :accessor a1-name)
   (a1-element :initarg :a1-element :accessor a1-element))
  (:report (lambda (obj stream)
             (let* ((name (chem:get-name (molecule obj)))
                    (properties (chem:matter/properties (molecule obj)))
                    (maybe-description (getf properties :description)))
               (format stream "Molecule(res) name: ~s(~s) cannot identify lk-solvation type for atom(element)s: (~a(~a)); ~a"
                       (chem:get-name (molecule obj))
                       (chem:get-name (residue obj))
                       (a1-name obj)
                       (a1-element obj)
                       (if maybe-description
                           (format nil "description: ~a" maybe-description)
                           "")
                       )))))

(defun missing-lk-type-error (a1 residue molecule)
  (let* ((a1-name (chem:get-name a1))
         (a1-element (chem:get-element a1)))
    (error 'missing-lk-type
           :molecule molecule
           :residue residue
           :a1-name a1-name
           :a1-element a1-element
           )))

(defun assign-lk-types (molecule lk-type-assignments &key (error-on-missing t) verbose)
  "Return a hash-table mapping atoms to LK type symbols."
  (let ((results (make-hash-table :test #'eq)))
    (chem:do-residues (res molecule)
      (chem:do-atoms (atom res)
        (let* ((element (chem:get-element atom))
               (assignments (lookup-lk-type-assignments element lk-type-assignments)))
          (if (not assignments)
              (if error-on-missing
                  (error "No LK type assignments for element ~s" element)
                  (setf (gethash atom results) nil))
              (let ((assigned nil))
                (dolist (assignment assignments)
                  ;;(format t "Looking at assignment: ~s~%" assignment)
                  (let* ((compiled (rosetta-lk-type-assignment-compiled-smarts assignment)))
                    (when (chem:matches compiled atom)
                      (setf assigned (rosetta-lk-type-assignment-type assignment)))))
                (if assigned
                    (progn
                      (when verbose (format t "Assigned ~s(~s) -> ~s~%" (chem:get-name res) (chem:get-name atom) assigned))
                      (setf (gethash atom results) assigned))
                    (if error-on-missing
                        (progn
                          (restart-case
                              (missing-lk-type-error atom res molecule)
                            (chem:skip-missing-parameters ()
                              :report "skip this missing type and continue"
                              nil)))
                        (setf (gethash atom results) nil))))))))
    results))
