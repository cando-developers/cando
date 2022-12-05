(in-package :monomer-context)


(defclass monomer-match-node ()
  ((names :initarg :names :accessor names)))

(cando:make-class-save-load
 monomer-match-node
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (names obj)))))

(defclass plug-to-monomer-node ()
  ((plug-name :initarg :plug-name :accessor plug-name)
   (monomer-match-node :initarg :monomer-match-node :accessor monomer-match-node)))

(cando:make-class-save-load
 plug-to-monomer-node
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a" (plug-name obj) (monomer-match-node obj)))))

(defclass chain-node ()
  ((head :initarg :head :accessor head)
   (tail :initarg :tail :accessor tail)))

(cando:make-class-save-load
 chain-node
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a" (head obj) (tail obj)))))


(defclass branch-node ()
  ((left :initarg :left :accessor left)
   (right :initarg :right :accessor right)))

(cando:make-class-save-load
 branch-node
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a" (left obj) (right obj)))))

(defclass match ()
  ((parts :initform (make-array 16 :adjustable t :fill-pointer 0) :accessor parts)))

(defun match-as-string (match)
  (with-output-to-string (sout)
    (format sout "~{~a~^_~}" (coerce (parts match) 'list))))

(defun cursor (match)
  (fill-pointer (parts match)))

(defun add-match (match thing)
  (vector-push-extend thing (parts match)))

(defun unwind-cursor (match cursor)
  (setf (fill-pointer (parts match)) cursor))

(defun parse-monomer (name)
  (cond
    ((and (consp name) (eq (car name) :cap))
     nil)
    ((symbolp name)
     (let ((names (gethash name topology::*topology-groups*)))
       (if names
           (make-instance 'monomer-match-node
                          :names names)
           (error "Unrecognizable topology group name ~a" name))))))

(defun parse-plug-to-something (sexp)
  (when sexp
    (let ((maybe-plug-name (car sexp)))
      (cond
        ((consp maybe-plug-name)
         (let ((left (parse-plug-to-something (car sexp)))
               (right (parse-plug-to-something (cdr sexp))))
           (when (or left right)
             (if (and left right)
                 (make-instance 'branch-node
                                :left left
                                :right right)
                 (or left right)))))
        ((symbolp maybe-plug-name)
         (let ((head (let ((monomer (parse-monomer (cadr sexp))))
                                (when monomer
                                  (make-instance 'plug-to-monomer-node
                                                 :plug-name maybe-plug-name
                                                 :monomer-match-node monomer))))
               (tail (parse-plug-to-something (cddr sexp))))
           (when (or head tail)
             (make-instance 'chain-node
                            :head head
                            :tail tail))))))))

(defun parse (sexp)
  (let ((head (car sexp)))
    (make-instance 'chain-node
                   :head (parse-monomer head)
                   :tail (parse-plug-to-something (cdr sexp)))))


(defgeneric matches-impl (pattern monomer oligomer match))

(defun ensure-monomer-or-nil (monomer)
  (when monomer
    (if (typep monomer 'topology:monomer)
        monomer
        (error "monomer ~a must be a topology:monomer" monomer))))

(defmethod matches-impl ((pattern null) (monomer topology:monomer) oligomer match)
  (ensure-monomer-or-nil monomer))

(defmethod matches-impl ((pattern monomer-match-node) (monomer topology:monomer) (oligomer-space topology:oligomer-space) match)
  (ensure-monomer-or-nil
   (progn
     (loop for monomer-name in (topology:monomers monomer)
           unless (member monomer-name (names pattern))
             do (return-from matches-impl nil))
     (add-match match (topology:monomers monomer))
     monomer)))

(defmethod matches-impl ((pattern monomer-match-node) (monomer topology:monomer) (oligomer topology:oligomer) match)
  (ensure-monomer-or-nil
   (let ((monomer-name (topology:current-stereoisomer-name monomer oligomer)))
     (when (member monomer-name (names pattern))
       (add-match match monomer-name)
       monomer))))

(defmethod matches-impl ((pattern plug-to-monomer-node) (monomer topology:monomer) oligomer match)
  (ensure-monomer-or-nil
   (let* ((plug-name (plug-name pattern))
          (coupling (gethash plug-name (topology:couplings monomer))))
     (when coupling
       (add-match match plug-name)
       (let ((other-monomer (topology:other-monomer coupling monomer)))
         (when (matches-impl (monomer-match-node pattern) other-monomer oligomer match)
           other-monomer))))))

(defmethod matches-impl ((pattern chain-node) (monomer topology:monomer) oligomer match)
  (ensure-monomer-or-nil
   (let ((head-monomer (matches-impl (head pattern) monomer oligomer match)))
    (when head-monomer
      (matches-impl (tail pattern) head-monomer oligomer match)))))

(defmethod matches-impl ((pattern branch-node) (monomer topology:monomer) oligomer match)
  (ensure-monomer-or-nil
   (let ((cursor (cursor match)))
     (add-match match #\{)
     (if (matches-impl (left pattern) monomer oligomer match)
         (progn
           (add-match match #\})
           (matches-impl (right pattern) monomer oligomer match))
         (progn
           (unwind-cursor match cursor)
           nil)))))

(defun match (pattern monomer oligomer)
  (let ((match (make-instance 'match)))
    (when (matches-impl pattern monomer oligomer match)
      match)))
