(in-package :monomer-context)


(defgeneric %copy-specialize (obj)
  (:documentation "Deep copy the node except for some slots that I know are immutable"))

(defmethod %copy-specialize ((obj null))
  nil)

(defclass monomer-match-node ()
  ((names :initarg :names :reader names)))

(defclass wild-card-monomer-match-node (cando.serialize:serializable) ())

(defclass untrained-monomer-match-node ()
  ())

(defclass trained-monomer-match-node (untrained-monomer-match-node)
  ((name :initarg :name :accessor name)))

(defmethod print-object ((obj trained-monomer-match-node) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~s" (name obj))))

(defmethod %copy-specialize ((obj monomer-match-node))
  (make-instance 'untrained-monomer-match-node))

(cando.serialize:make-class-save-load monomer-match-node
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (names obj)))))

(defclass plug-to-monomer-node ()
  ((plug-name :initarg :plug-name :reader plug-name)
   (monomer-match-node :initarg :monomer-match-node :accessor monomer-match-node)))

(defmethod %copy-specialize ((obj (eql :***)))
  (make-instance 'wild-card-monomer-match-node))

(defmethod %copy-specialize ((obj plug-to-monomer-node))
  (make-instance 'plug-to-monomer-node
                 :plug-name (plug-name obj)
                 :monomer-match-node (%copy-specialize (monomer-match-node obj))))

(cando.serialize:make-class-save-load plug-to-monomer-node
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a" (plug-name obj) (monomer-match-node obj)))))

(defclass context-node ()
  ((node-count :initarg :node-count :reader node-count)))

(defclass chain-node (context-node)
  ((head :initarg :head :reader head)
   (tail :initarg :tail :reader tail)))

(defmethod %copy-specialize ((obj chain-node))
  (make-instance 'chain-node
                 :head (%copy-specialize (head obj))
                 :tail (%copy-specialize (tail obj))
                 :node-count (1+ (if (tail obj)
                                     (node-count (tail obj))
                                     0))))

(cando.serialize:make-class-save-load chain-node
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a ~a" (head obj) (tail obj) (node-count obj)))))

(defclass branch-node (context-node)
  ((left :initarg :left :reader left)
   (right :initarg :right :reader right)))

(defmethod %copy-specialize ((obj branch-node))
  (make-instance 'branch-node
                 :left (%copy-specialize (left obj))
                 :right (%copy-specialize (right obj))
                 :node-count (+ (node-count (left obj))
                                (node-count (right obj)))))

(cando.serialize:make-class-save-load branch-node
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a ~a" (left obj) (right obj) (node-count obj)))))

(defclass match ()
  ((parts :initform (make-array 16 :adjustable t :fill-pointer 0) :accessor parts)
   (maybe-monomers :initform (make-array 16 :adjustable t :fill-pointer 0) :accessor maybe-monomers)))

(defun match-iterator (match)
  "Return an iterator that will return successive monomer-contexts and then nil"
  (let* ((name-indexes (coerce (loop for index below (length (parts match))
                                     for maybe-names = (elt (parts match) index)
                                     when (consp maybe-names)
                                       collect (cons index (length maybe-names)))
                               'vector))
         (counters (make-array (length name-indexes) :element-type 'ext:byte32 :initial-element 0))
         (done nil))
    (lambda ()
      (let (result
            focus-monomer-name)
        (when (not done)
          (setf focus-monomer-name (elt (elt (parts match) 0) (elt counters 0))
                result (intern (format nil "~{~a~^_~}" (loop with counter-index = 0
                                                             for index below (length (parts match))
                                                             for maybe-name = (elt (parts match) index)
                                                             if (consp maybe-name)
                                                               collect (let* ((name (elt maybe-name (elt counters counter-index))))
                                                                         (incf counter-index)
                                                                         name)
                                                             else
                                                               collect maybe-name)) :keyword))
          (loop named advance
                for counter-index below (length counters)
                do (incf (elt counters counter-index))
                do (when (< (elt counters counter-index) (cdr (elt name-indexes counter-index)))
                     (return-from advance nil))
                do (setf (elt counters counter-index) 0)
                finally (setf done t)))
        (values result focus-monomer-name counters)))))

(defun match-as-symbol (match)
  (intern (with-output-to-string (sout)
            (format sout "~{~a~^_~}" (coerce (parts match) 'list))) :keyword))

(defun cursor (match)
  (fill-pointer (parts match)))

(defun add-match (match thing &optional monomer)
  (vector-push-extend thing (parts match))
  (vector-push-extend monomer (maybe-monomers match)))

(defun unwind-cursor (match cursor)
  (setf (fill-pointer (parts match)) cursor))

(defun parse-monomer (name)
  (cond
    ((and (consp name) (eq (car name) :cap))
     nil)
    ((and (consp name) (eq (car name) :ring-cap))
     :*** ; This is a wildcard
     )
    ((symbolp name)
     (let ((names (gethash name topology::*topology-groups*)))
       (if names
           (make-instance 'monomer-match-node
                          :names names)
           (error "Unrecognizable topology group name ~a" name))))))


(defun parse-plug-to-something (sexp)
  (if sexp
      (let ((maybe-plug-name (car sexp)))
        (cond
          ((consp maybe-plug-name)
           (multiple-value-bind (left left-node-count)
               (parse-plug-to-something (car sexp))
             (multiple-value-bind (right right-node-count)
                 (parse-plug-to-something (cdr sexp))
               (cond
                 ((and left right)
                  (values
                   (make-instance 'branch-node
                                  :left left
                                  :right right
                                  :node-count (+ left-node-count right-node-count))
                   (+ left-node-count right-node-count)))
                 (left
                  (values left left-node-count))
                 (right
                  (values right right-node-count))
                 (t
                  (values nil 0))))))
          ((symbolp maybe-plug-name)
           (let ((head (let ((monomer (parse-monomer (cadr sexp))))
                         (when monomer
                           (make-instance 'plug-to-monomer-node
                                          :plug-name maybe-plug-name
                                          :monomer-match-node monomer)))))
             (if head
                 (multiple-value-bind (tail node-count)
                     (parse-plug-to-something (cddr sexp))
                   (values (make-instance 'chain-node
                                          :head head
                                          :tail tail
                                          :node-count (1+ node-count))
                           (1+ node-count)))
                 (values nil 0))))))
      (values nil 0) ;; no plug to anything
      ))

(defun parse-sexp (sexp)
  (let ((head (car sexp)))
    (multiple-value-bind (tail node-count)
        (parse-plug-to-something (cdr sexp))
      (values (make-instance 'chain-node
                             :head (parse-monomer head)
                             :tail tail
                             :node-count (1+ node-count))
              (1+ node-count)))))

(defgeneric matches-impl (pattern monomer oligomer match cache))

(defclass wild-card-monomer () ())

(defparameter *wild-card-monomer* (make-instance 'wild-card-monomer))

(defun ensure-monomer-or-nil (monomer)
  (when monomer
    (cond
      ((typep monomer 'wild-card-monomer) monomer)
      ((typep monomer 'topology:monomer)
       monomer)
      (t (error "monomer ~a must be a topology:monomer" monomer)))))

(defmethod matches-impl ((pattern null) (monomer topology:monomer) oligomer match cache)
  (declare (ignore oligomer match))
  (ensure-monomer-or-nil monomer))

(defmethod matches-impl ((pattern monomer-match-node) (monomer topology:monomer) (oligomer-space topology:oligomer-space) match cache)
  (ensure-monomer-or-nil
   (progn
     (loop for monomer-name in (topology:monomers monomer)
           unless (member monomer-name (names pattern))
             do (return-from matches-impl nil))
     (add-match match (topology:monomers monomer) monomer)
     monomer)))

(defstruct match-cache monomer-names)

(defun cached-stereoisomer-name (cache monomer oligomer)
  (let ((cached-name (gethash monomer (match-cache-monomer-names cache))))
    (if cached-name
        cached-name
        (let ((name (topology:current-stereoisomer-name monomer oligomer)))
          (setf monomer (match-cache-monomer-names cache))
          name))))

(defmethod matches-impl ((pattern monomer-match-node) (monomer topology:monomer) (oligomer topology:oligomer) match cache)
  (ensure-monomer-or-nil
   (let ((monomer-name (cached-stereoisomer-name cache monomer oligomer)))
     (when (member monomer-name (names pattern))
       (add-match match monomer-name monomer)
       monomer))))

(defmethod matches-impl ((pattern trained-monomer-match-node) (monomer topology:monomer) (oligomer topology:oligomer) match cache)
   (let ((monomer-name (cached-stereoisomer-name cache monomer oligomer)))
     (when (eq monomer-name (name pattern))
       (add-match match monomer-name monomer)
       monomer)))

(defmethod matches-impl ((pattern untrained-monomer-match-node) (monomer topology:monomer) (oligomer topology:oligomer) match cache)
  (let ((monomer-name (cached-stereoisomer-name cache monomer oligomer)))
    (change-class pattern 'trained-monomer-match-node
                  :name monomer-name)
    monomer))

(defmethod matches-impl ((pattern wild-card-monomer-match-node) (monomer topology:monomer) (oligomer topology:oligomer) match cache)
  (warn "Returning monomer ~s for wild-card-monomer-match-node" monomer)
  monomer)

(defmethod matches-impl ((pattern plug-to-monomer-node) (monomer topology:monomer) oligomer match cache)
  (ensure-monomer-or-nil
   (let* ((plug-name (plug-name pattern))
          (coupling (gethash plug-name (topology:couplings monomer))))
     (when coupling
       (add-match match plug-name)
       (let ((other-monomer (topology:other-monomer coupling monomer)))
         (when (matches-impl (monomer-match-node pattern) other-monomer oligomer match cache)
           other-monomer))))))

(defmethod matches-impl ((pattern (eql :***)) (monomer topology:monomer) oligomer match cache)
  (ensure-monomer-or-nil
   *wild-card-monomer*))

(defmethod matches-impl ((pattern chain-node) (monomer topology:monomer) oligomer match cache)
  (ensure-monomer-or-nil
   (let ((head-monomer (matches-impl (head pattern) monomer oligomer match cache)))
    (when head-monomer
      (matches-impl (tail pattern) head-monomer oligomer match cache)))))

(defmethod matches-impl ((pattern branch-node) (monomer topology:monomer) oligomer match cache)
  (ensure-monomer-or-nil
   (let ((cursor (cursor match)))
     (add-match match #\{)
     (if (matches-impl (left pattern) monomer oligomer match cache)
         (progn
           (add-match match #\})
           (matches-impl (right pattern) monomer oligomer match cache))
         (progn
           (unwind-cursor match cursor)
           nil)))))

(defparameter *match-cache* nil)

(defmacro with-match-cache (&body body)
  `(let ((*match-cache* (make-match-cache :monomer-names (make-hash-table))))
     ,@body))

(defun match (pattern monomer oligomer-or-space)
  "Try to match the pattern starting on the monomer in the oligomer-space.
If oligomer-or-space is an oligomer-space then the match can be iterated over to generate
all monomer-contexts that match.  Use (match-iterator match) to generate that iterator."
  (if *match-cache*
      (let ((match (make-instance 'match)))
        (when (matches-impl pattern monomer oligomer-or-space match *match-cache*)
          match))
      (with-match-cache
          (let ((match (make-instance 'match)))
            (when (matches-impl pattern monomer oligomer-or-space match *match-cache*)
              match))
        )))


(defun copy-specialized (matcher focus-monomer oligomer)
  (let ((untrained (%copy-specialize matcher)))
    (with-match-cache
        (monomer-context:match untrained focus-monomer oligomer))
    untrained))
