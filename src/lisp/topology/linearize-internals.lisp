(in-package :topology)

(defclass linearized-matched-fragment-conformations-holder ()
  ((internals-values :initarg :internals-values :accessor internals-values)
   (internals-types :initarg :internals-types :accessor internals-types)
   (internals-index-vector :initarg :internals-index-vector :accessor internals-index-vector)
   (internals-count-vector :initarg :internals-count-vector :accessor internals-count-vector)
   (probabilities-vector :initarg :probabilities-vector :accessor probabilities-vector)
   (energies-vector :initarg :energies-vector :accessor energies-vector)
   (trainer-index-vector :initarg :trainer-index-vector :accessor trainer-index-vector)
   (fragment-conformations-start-vector :initarg :fragment-conformations-start-vector :accessor fragment-conformations-start-vector)
   (fragment-conformations-count-vector :initarg :fragment-conformations-count-vector :accessor fragment-conformations-count-vector)
   (monomer-context-names-vector :initarg :monomer-context-names-vector :accessor monomer-context-names-vector)
   (monomer-context-names-start :initarg :monomer-context-names-start :accessor monomer-context-names-start)
   (monomer-context-names-count :initarg :monomer-context-names-count :accessor monomer-context-names-count)
   (focus-monomer-names-vector :initarg :focus-monomer-names-vector :accessor focus-monomer-names-vector)
   (focus-monomer-names-start :initarg :focus-monomer-names-start :accessor focus-monomer-names-start)
   (focus-monomer-names-count :initarg :focus-monomer-names-count :accessor focus-monomer-names-count)
   (fragment-match-fragment-conformations-index-vector  :initarg :fragment-match-fragment-conformations-index-vector  :accessor fragment-match-fragment-conformations-index-vector )
   (fragment-match-fragment-conformations-start-vector  :initarg :fragment-match-fragment-conformations-start-vector  :accessor fragment-match-fragment-conformations-start-vector )
   (fragment-match-fragment-conformations-count-vector  :initarg :fragment-match-fragment-conformations-count-vector  :accessor fragment-match-fragment-conformations-count-vector )
   (fragment-match-start-vector  :initarg :fragment-match-start-vector  :accessor fragment-match-start-vector )
   (fragment-match-count-vector  :initarg :fragment-match-count-vector  :accessor fragment-match-count-vector )
   (fragment-match-key-from-vector  :initarg :fragment-match-key-from-vector  :accessor fragment-match-key-from-vector )
   (fragment-match-key-to-vector :initarg :fragment-match-key-to-vector :accessor fragment-match-key-to-vector)))

#|
monomer-context-to-index maps monomer-context to a monomer-context-index
monomer-context-index indexes into ...
... fragment-conformations-start-vector
... fragment-conformations-count-vector
...... these index into ...
...... probabilities-vector
...... energies-vector
...... trainer-index-vector
...... internals-index-vector
...... internals-count-vector
......... these index into ...
......... internals-values
......... internals-types


|#
(defclass linearized-matched-fragment-conformations-map ()
  ((linearized-matched-fragment-conformations-holder
    :accessor linearized-matched-fragment-conformations-holder
    :initarg :linearized-matched-fragment-conformations-holder )
   (monomer-context-to-index :accessor monomer-context-to-index
                             :initarg :monomer-context-to-index)
   (monomer-context-to-fragment-conformations
    :accessor monomer-context-to-fragment-conformations
    :initarg :monomer-context-to-fragment-conformations)
   (focus-monomer-names :accessor focus-monomer-names
                        :initarg :focus-monomer-names)
   (fragment-matches :accessor fragment-matches
                                :initarg :fragment-matches)
   ))

(defclass linearized-fragment-conformations ()
  ((fragments :initarg :fragments :accessor fragments)))

(defclass linearized-fragment-internals ()
  ((energy :initarg :energy :accessor energy)
   (probability :initarg :probability :accessor probability)
   (trainer-index :initarg :trainer-index :accessor trainer-index)
   (internals-types :initarg :internals-types :accessor internals-types)
   (internals-values :initarg :internals-values :accessor internals-values)))

(defgeneric compact-fragment-internal-type (internal))
(defmethod compact-fragment-internal-type ((internal jump-internal)) #\j)
(defmethod compact-fragment-internal-type ((internal bonded-internal)) #\b)
(defmethod compact-fragment-internal-type ((internal complex-bonded-internal)) #\c)

(defgeneric fill-vector-values (internal internals-values))
(defmethod fill-vector-values ((internal jump-internal) internals-values)
  (vector-push-extend 0.0s0 internals-values)
  (vector-push-extend 0.0s0 internals-values)
  (vector-push-extend 0.0s0 internals-values))

(defmethod fill-vector-values ((internal bonded-internal) internals-values)
  (vector-push-extend (bond internal) internals-values)
  (vector-push-extend (angle internal) internals-values)
  (vector-push-extend (dihedral internal) internals-values))

(defun fill-vectors (internal internals-types internals-values)
  (vector-push-extend (compact-fragment-internal-type internal) internals-types)
  (fill-vector-values internal internals-values))

(defun vector-push-cstring (string byte-vector)
  (loop for cc across string
        do (vector-push-extend cc byte-vector))
  (vector-push-extend (code-char 0) byte-vector)
  (length string))

(defmacro debug-linearize (limit format &rest args)
  `(when (and ,limit (> ,limit 0)) (let ((*print-pretty* nil)) (format t ,format ,@args))))

(defmacro debug-linearize-decf (limit)
  `(decf ,limit))

(defun create-linearized-matched-fragment-conformations-holder (matched-fragment-conformations-map &key (debug-limit 0))
  (let ((internals-values (make-array (* 1024 1024 16) :element-type 'single-float :adjustable t :fill-pointer 0))
        (internals-types (make-array (* 1024 1024 16) :element-type 'simple-char :adjustable t :fill-pointer 0))
        (monomer-context-to-index (make-hash-table))
        (monomer-context-index -1)
        (probabilities-vector (make-array (* 1024 1024 2) :element-type 'single-float :adjustable t :fill-pointer 0))
        (energies-vector (make-array (* 1024 1024 2) :element-type 'single-float :adjustable t :fill-pointer 0))
        (trainer-index-vector (make-array (* 1024 1024 2) :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (internals-index-vector (make-array (* 1024 1024 2) :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (internals-count-vector (make-array (* 1024 1024 2) :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (fragment-conformations-start-vector (make-array 1024 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (fragment-conformations-count-vector (make-array 1024 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (monomer-context-names-vector (make-array 65536 :element-type 'base-char :adjustable t :fill-pointer 0))
        (monomer-context-names-start (make-array 65536 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (monomer-context-names-count (make-array 65536 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (focus-monomer-names-vector (make-array 65536 :element-type 'base-char :adjustable t :fill-pointer 0))
        (focus-monomer-names-start (make-array 65536 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (focus-monomer-names-count (make-array 65536 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (dl debug-limit))
    (maphash (lambda (key fragment-conformations)
               (debug-linearize-decf dl)
               (debug-linearize dl "Looking at monomer-context ~a~%" key)
               (let ((compact-fragments-start-index (length internals-index-vector)))
                 (loop for fragment-internals across (fragments fragment-conformations)
                       for internals-index = (fill-pointer internals-values)
                       do (loop for internal across (internals fragment-internals)
                                do (fill-vectors internal internals-types internals-values))
                       do (vector-push-extend (if (probability fragment-internals)
                                                  (float (probability fragment-internals) 1.0s0)
                                                  -1.0s0)
                                              probabilities-vector)
                       do (vector-push-extend (if (energy fragment-internals)
                                                  (float (energy fragment-internals) 1.0s0)
                                                  0.0s0)
                                              energies-vector)
                       do (vector-push-extend (trainer-index fragment-internals) trainer-index-vector)
                       do (vector-push-extend internals-index internals-index-vector)
                       do (vector-push-extend (length (internals fragment-internals)) internals-count-vector))
                 (let ((compact-fragments-end-index (length internals-index-vector)))
                   (vector-push-extend compact-fragments-start-index fragment-conformations-start-vector)
                   (vector-push-extend (- compact-fragments-end-index compact-fragments-start-index) fragment-conformations-count-vector))
                 (setf (gethash key monomer-context-to-index) (incf monomer-context-index))
                 (vector-push-extend (length monomer-context-names-vector) monomer-context-names-start)
                 (let ((count (vector-push-cstring (string key) monomer-context-names-vector)))
                   (vector-push-extend count monomer-context-names-count))
                 (vector-push-extend (length focus-monomer-names-vector) focus-monomer-names-start)
                 (let ((count (vector-push-cstring (string (focus-monomer-name fragment-conformations)) focus-monomer-names-vector)))
                   (vector-push-extend count focus-monomer-names-count))))
             (monomer-context-to-fragment-conformations matched-fragment-conformations-map))
    (let ((fragment-match-fragment-conformations-index-vector (make-array (* 1024 1024 2) :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
          (fragment-match-fragment-conformations-start-vector (make-array 65536 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
          (fragment-match-fragment-conformations-count-vector (make-array 65536 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
          (fragment-match-start-vector (make-array 1024 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
          (fragment-match-count-vector (make-array 1024 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
          (fragment-match-key-from-vector (make-array 1024 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
          (fragment-match-key-to-vector (make-array 1024 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
          (dl debug-limit))
      (maphash (lambda (key fragment-conformation-indices-vectors)
                 (let* ((key-from (car key))
                        (key-from-index (gethash key-from monomer-context-to-index))
                        (key-to (cdr key))
                        (key-to-index (gethash key-to monomer-context-to-index)))
                   (debug-linearize-decf dl)
                   (debug-linearize dl "fragment-match key ~s~%" key)
                   (vector-push-extend key-from-index fragment-match-key-from-vector)
                   (vector-push-extend key-to-index fragment-match-key-to-vector)
                   (let ((fm-start (length fragment-match-fragment-conformations-start-vector))) ; was fragment-match-start-vector)))
                     (debug-linearize dl "fm-start ~a~%" fm-start)
                     (vector-push-extend fm-start fragment-match-start-vector)
                     (debug-linearize dl "fragment-match-start-vector ~a~%" fm-start)
                     (loop for fragment-conformation-indices across fragment-conformation-indices-vectors
                           do (let ((fc-start (length fragment-match-fragment-conformations-index-vector)))
                                (debug-linearize dl "fc-start = ~a~%" fc-start) 
                                (vector-push-extend fc-start fragment-match-fragment-conformations-start-vector)
                                (debug-linearize dl "fragment-match-fragment-conformations-start-vector ~a~%" fc-start)
                                (when fragment-conformation-indices
                                  (loop for fragment-conformation-index across fragment-conformation-indices
                                        do (progn
                                             (vector-push-extend fragment-conformation-index fragment-match-fragment-conformations-index-vector)
                                             (debug-linearize dl "fragment-match-fragment-conformations-index-vector ~a~%" fragment-conformation-index)
                                        )))
                                (let ((count (- (length fragment-match-fragment-conformations-index-vector) fc-start)))
                                  (vector-push-extend count fragment-match-fragment-conformations-count-vector)
                                  (debug-linearize dl "fragment-match-fragment-conformations-count-vector ~a~%" count)
                                  )))
                     (let ((count (- (length fragment-match-fragment-conformations-start-vector) fm-start)))
                       (vector-push-extend count fragment-match-count-vector)
                       (debug-linearize dl "fragment-match-count-vector ~a~%" count))
                     )))
               (fragment-matches matched-fragment-conformations-map))
      (make-instance 'linearized-matched-fragment-conformations-holder
                     :internals-values (copy-seq internals-values)
                     :internals-types (copy-seq internals-types)
                     :internals-index-vector (copy-seq internals-index-vector)
                     :internals-count-vector (copy-seq internals-count-vector)
                     :fragment-conformations-start-vector (copy-seq fragment-conformations-start-vector)
                     :fragment-conformations-count-vector (copy-seq fragment-conformations-count-vector)
                     :monomer-context-names-vector (copy-seq monomer-context-names-vector)
                     :monomer-context-names-start (copy-seq monomer-context-names-start)
                     :monomer-context-names-count (copy-seq monomer-context-names-count)
                     :focus-monomer-names-vector (copy-seq focus-monomer-names-vector)
                     :focus-monomer-names-start (copy-seq focus-monomer-names-start)
                     :focus-monomer-names-count (copy-seq focus-monomer-names-count)
                     :probabilities-vector probabilities-vector
                     :energies-vector energies-vector
                     :trainer-index-vector trainer-index-vector
                     :fragment-match-fragment-conformations-index-vector  (copy-seq fragment-match-fragment-conformations-index-vector)
                     :fragment-match-fragment-conformations-start-vector  (copy-seq fragment-match-fragment-conformations-start-vector)
                     :fragment-match-fragment-conformations-count-vector  (copy-seq fragment-match-fragment-conformations-count-vector)
                     :fragment-match-start-vector  (copy-seq fragment-match-start-vector)
                     :fragment-match-count-vector  (copy-seq fragment-match-count-vector)
                     :fragment-match-key-from-vector  (copy-seq fragment-match-key-from-vector)
                     :fragment-match-key-to-vector (copy-seq fragment-match-key-to-vector))
      )))

(defmacro maybe-describe (sym)
  `(progn
     #+(or)(format t "Initializing ~a~%" ',sym)
     ,sym))

(defun save-linearized-matched-fragment-conformations-holder (linearized-fragment-conformations-holder filename)
  (format t "save-compact-matched-fragment-conformations-holder~%")
  (with-accessors ((internals-values internals-values)
                   (internals-types internals-types)
                   (monomer-context-names-vector monomer-context-names-vector)
                   (monomer-context-names-start monomer-context-names-start)
                   (monomer-context-names-count monomer-context-names-count)
                   (focus-monomer-names-vector focus-monomer-names-vector)
                   (focus-monomer-names-start focus-monomer-names-start)
                   (focus-monomer-names-count focus-monomer-names-count)
                   (internals-index-vector internals-index-vector)
                   (internals-count-vector internals-count-vector)
                   (probabilities-vector probabilities-vector)
                   (energies-vector energies-vector)
                   (trainer-index-vector trainer-index-vector)
                   (fragment-conformations-start-vector fragment-conformations-start-vector)
                   (fragment-conformations-count-vector fragment-conformations-count-vector)
                   (fragment-match-fragment-conformations-index-vector fragment-match-fragment-conformations-index-vector )
                   (fragment-match-fragment-conformations-start-vector fragment-match-fragment-conformations-start-vector )
                   (fragment-match-fragment-conformations-count-vector fragment-match-fragment-conformations-count-vector )
                   (fragment-match-start-vector  fragment-match-start-vector )
                   (fragment-match-count-vector  fragment-match-count-vector )
                   (fragment-match-key-from-vector  fragment-match-key-from-vector )
                   (fragment-match-key-to-vector fragment-match-key-to-vector)
                   )
      linearized-fragment-conformations-holder
    (static-vectors:with-static-vectors ((static-internals-values (length internals-values)
                                                                  :element-type 'single-float
                                                                  :initial-contents (maybe-describe internals-values))
                                         (static-internals-types (length internals-types) :element-type 'base-char :initial-contents (maybe-describe internals-types))
                                         (static-internals-index-vector (length internals-index-vector) :element-type 'ext:byte32 :initial-contents (maybe-describe internals-index-vector))
                                         (static-internals-count-vector (length internals-count-vector) :element-type 'ext:byte32 :initial-contents (maybe-describe internals-count-vector))
                                         (static-probabilities-vector (length probabilities-vector) :element-type 'single-float :initial-contents (maybe-describe probabilities-vector))
                                         (static-energies-vector (length energies-vector) :element-type 'single-float :initial-contents (maybe-describe energies-vector))
                                         (static-trainer-index-vector (length trainer-index-vector) :element-type 'ext:byte32 :initial-contents (maybe-describe trainer-index-vector))
                                         (static-fragment-conformations-start-vector (length fragment-conformations-start-vector) :element-type 'ext:byte32 :initial-contents (maybe-describe fragment-conformations-start-vector))
                                         (static-fragment-conformations-count-vector (length fragment-conformations-count-vector) :element-type 'ext:byte32 :initial-contents (maybe-describe fragment-conformations-count-vector))
                                         (static-monomer-context-names-vector (length monomer-context-names-vector) :element-type 'base-char :initial-contents (maybe-describe monomer-context-names-vector))
                                         (static-monomer-context-names-start (length monomer-context-names-start) :element-type 'ext:byte32 :initial-contents (maybe-describe monomer-context-names-start))
                                         (static-monomer-context-names-count (length monomer-context-names-count) :element-type 'ext:byte32 :initial-contents (maybe-describe monomer-context-names-count))
                                         (static-focus-monomer-names-vector (length focus-monomer-names-vector) :element-type 'base-char :initial-contents (maybe-describe focus-monomer-names-vector))
                                         (static-focus-monomer-names-start (length focus-monomer-names-start) :element-type 'ext:byte32 :initial-contents (maybe-describe focus-monomer-names-start))
                                         (static-focus-monomer-names-count (length focus-monomer-names-count) :element-type 'ext:byte32 :initial-contents (maybe-describe focus-monomer-names-count))
                                         (static-fragment-match-fragment-conformations-index-vector  (length fragment-match-fragment-conformations-index-vector ) :element-type 'ext:byte32 :initial-contents (maybe-describe fragment-match-fragment-conformations-index-vector ))
                                         (static-fragment-match-fragment-conformations-start-vector  (length fragment-match-fragment-conformations-start-vector ) :element-type 'ext:byte32 :initial-contents (maybe-describe fragment-match-fragment-conformations-start-vector ))
                                         (static-fragment-match-fragment-conformations-count-vector  (length fragment-match-fragment-conformations-count-vector ) :element-type 'ext:byte32 :initial-contents (maybe-describe fragment-match-fragment-conformations-count-vector ))
                                         (static-fragment-match-start-vector  (length fragment-match-start-vector ) :element-type 'ext:byte32 :initial-contents (maybe-describe fragment-match-start-vector ))
                                         (static-fragment-match-count-vector  (length fragment-match-count-vector ) :element-type 'ext:byte32 :initial-contents (maybe-describe fragment-match-count-vector ))
                                         (static-fragment-match-key-from-vector  (length fragment-match-key-from-vector ) :element-type 'ext:byte32 :initial-contents (maybe-describe fragment-match-key-from-vector ))
                                         (static-fragment-match-key-to-vector (length fragment-match-key-to-vector) :element-type 'ext:byte32 :initial-contents (maybe-describe fragment-match-key-to-vector)))
      (let ((cdf (nc:create (namestring (merge-pathnames filename)))))
        (nc:def-dim cdf "probabilities-vector" (length static-probabilities-vector))
        (nc:def-var cdf "probabilities-vector" nc-c:+float+ '("probabilities-vector"))
        (nc:def-dim cdf "energies-vector" (length static-energies-vector))
        (nc:def-var cdf "energies-vector" nc-c:+float+ '("energies-vector"))
        (nc:def-dim cdf "trainer-index-vector" (length static-trainer-index-vector))
        (nc:def-var cdf "trainer-index-vector" nc-c:+int+ '("trainer-index-vector"))
        (nc:def-dim cdf "internals-values" (length static-internals-values))
        (nc:def-var cdf "internals-values" nc-c:+float+ '("internals-values"))
        (nc:def-dim cdf "internals-types" (length internals-types))
        (nc:def-var cdf "internals-types" nc-c:+byte+ '("internals-types"))
        (nc:def-dim cdf "monomer-context-names-vector" (length static-monomer-context-names-vector))
        (nc:def-var cdf "monomer-context-names-vector" nc-c:+byte+ '("monomer-context-names-vector"))
        (nc:def-dim cdf "monomer-context-names-start" (length static-monomer-context-names-start))
        (nc:def-var cdf "monomer-context-names-start" nc-c:+int+ '("monomer-context-names-start"))
        (nc:def-dim cdf "monomer-context-names-count" (length static-monomer-context-names-count))
        (nc:def-var cdf "monomer-context-names-count" nc-c:+int+ '("monomer-context-names-count"))
        (nc:def-dim cdf "focus-monomer-names-vector" (length static-focus-monomer-names-vector))
        (nc:def-var cdf "focus-monomer-names-vector" nc-c:+byte+ '("focus-monomer-names-vector"))
        (nc:def-dim cdf "focus-monomer-names-start" (length static-focus-monomer-names-start))
        (nc:def-var cdf "focus-monomer-names-start" nc-c:+int+ '("focus-monomer-names-start"))
        (nc:def-dim cdf "focus-monomer-names-count" (length static-focus-monomer-names-count))
        (nc:def-var cdf "focus-monomer-names-count" nc-c:+int+ '("focus-monomer-names-count"))
        (nc:def-dim cdf "internals-index-vector" (length static-internals-index-vector))
        (nc:def-var cdf "internals-index-vector" nc-c:+int+ '("internals-index-vector"))
        (nc:def-dim cdf "internals-count-vector" (length static-internals-count-vector))
        (nc:def-var cdf "internals-count-vector" nc-c:+int+ '("internals-count-vector"))
        (nc:def-dim cdf "fragment-conformations-start-vector" (length static-fragment-conformations-start-vector))
        (nc:def-var cdf "fragment-conformations-start-vector" nc-c:+int+ '("fragment-conformations-start-vector"))
        (nc:def-dim cdf "fragment-conformations-count-vector" (length static-fragment-conformations-count-vector))
        (nc:def-var cdf "fragment-conformations-count-vector" nc-c:+int+ '("fragment-conformations-count-vector"))
        (nc:def-dim cdf "fragment-match-fragment-conformations-index-vector" (length static-fragment-match-fragment-conformations-index-vector ))
        (nc:def-var cdf "fragment-match-fragment-conformations-index-vector" nc-c:+int+ '("fragment-match-fragment-conformations-index-vector"))
        (nc:def-dim cdf "fragment-match-fragment-conformations-start-vector" (length static-fragment-match-fragment-conformations-start-vector ))
        (nc:def-var cdf "fragment-match-fragment-conformations-start-vector" nc-c:+int+ '("fragment-match-fragment-conformations-start-vector"))
        (nc:def-dim cdf "fragment-match-fragment-conformations-count-vector" (length static-fragment-match-fragment-conformations-count-vector ))
        (nc:def-var cdf "fragment-match-fragment-conformations-count-vector" nc-c:+int+ '("fragment-match-fragment-conformations-count-vector"))
        (nc:def-dim cdf "fragment-match-start-vector" (length static-fragment-match-start-vector ))
        (nc:def-var cdf "fragment-match-start-vector" nc-c:+int+ '("fragment-match-start-vector"))
        (nc:def-dim cdf "fragment-match-count-vector" (length static-fragment-match-count-vector ))
        (nc:def-var cdf "fragment-match-count-vector" nc-c:+int+ '("fragment-match-count-vector"))
        (nc:def-dim cdf "fragment-match-key-from-vector" (length static-fragment-match-key-from-vector ))
        (nc:def-var cdf "fragment-match-key-from-vector" nc-c:+int+ '("fragment-match-key-from-vector"))
        (nc:def-dim cdf "fragment-match-key-to-vector" (length static-fragment-match-key-to-vector))
        (nc:def-var cdf "fragment-match-key-to-vector" nc-c:+int+ '("fragment-match-key-to-vector"))
        (nc:enddef cdf)
        (nc:put-var-static-vector-float cdf "probabilities-vector" static-probabilities-vector)
        (nc:put-var-static-vector-float cdf "energies-vector" static-energies-vector)
        (nc:put-var-static-vector-uint cdf "trainer-index-vector" static-trainer-index-vector)
        (nc:put-var-static-vector-float cdf "internals-values" static-internals-values)
        (nc:put-var-static-vector-ubyte cdf "internals-types" static-internals-types)
        (nc:put-var-static-vector-ubyte cdf "monomer-context-names-vector" static-monomer-context-names-vector)
        (nc:put-var-static-vector-uint cdf "monomer-context-names-start" static-monomer-context-names-start)
        (nc:put-var-static-vector-uint cdf "monomer-context-names-count" static-monomer-context-names-count)
        (nc:put-var-static-vector-ubyte cdf "focus-monomer-names-vector" static-focus-monomer-names-vector)
        (nc:put-var-static-vector-uint cdf "focus-monomer-names-start" static-focus-monomer-names-start)
        (nc:put-var-static-vector-uint cdf "focus-monomer-names-count" static-focus-monomer-names-count)
        (nc:put-var-static-vector-uint cdf "internals-index-vector" static-internals-index-vector)
        (nc:put-var-static-vector-uint cdf "internals-count-vector" static-internals-count-vector)
        (nc:put-var-static-vector-uint cdf "fragment-conformations-start-vector" static-fragment-conformations-start-vector)
        (nc:put-var-static-vector-uint cdf "fragment-conformations-count-vector" static-fragment-conformations-count-vector)
        (nc:put-var-static-vector-uint cdf "fragment-match-fragment-conformations-index-vector" static-fragment-match-fragment-conformations-index-vector )
        (nc:put-var-static-vector-uint cdf "fragment-match-fragment-conformations-start-vector" static-fragment-match-fragment-conformations-start-vector )
        (nc:put-var-static-vector-uint cdf "fragment-match-fragment-conformations-count-vector" static-fragment-match-fragment-conformations-count-vector )
        (nc:put-var-static-vector-uint cdf "fragment-match-start-vector" static-fragment-match-start-vector )
        (nc:put-var-static-vector-uint cdf "fragment-match-count-vector" static-fragment-match-count-vector )
        (nc:put-var-static-vector-uint cdf "fragment-match-key-from-vector" static-fragment-match-key-from-vector )
        (nc:put-var-static-vector-uint cdf "fragment-match-key-to-vector" static-fragment-match-key-to-vector)
        (nc:nc-close cdf)))))



(defun read-linearized-matched-fragment-conformations-holder (filename)
  (let* ((cdf (nc:nc-open (namestring (merge-pathnames filename))))
         (dim-internals-values (nc:len (nc:get-dimension cdf "internals-values")))
         (dim-probabilities-vector (nc:len (nc:get-dimension cdf "probabilities-vector")))
         (dim-energies-vector (nc:len (nc:get-dimension cdf "energies-vector")))
         (dim-trainer-index-vector (nc:len (nc:get-dimension cdf "trainer-index-vector")))
         (dim-internals-types (nc:len (nc:get-dimension cdf "internals-types")))
         (dim-internals-index-vector (nc:len (nc:get-dimension cdf "internals-index-vector")))
         (dim-internals-count-vector (nc:len (nc:get-dimension cdf "internals-count-vector")))
         (dim-fragment-conformations-start-vector (nc:len (nc:get-dimension cdf "fragment-conformations-start-vector")))
         (dim-fragment-conformations-count-vector (nc:len (nc:get-dimension cdf "fragment-conformations-count-vector")))
         (dim-monomer-context-names-vector (nc:len (nc:get-dimension cdf "monomer-context-names-vector")))
         (dim-monomer-context-names-start (nc:len (nc:get-dimension cdf "monomer-context-names-start")))
         (dim-monomer-context-names-count (nc:len (nc:get-dimension cdf "monomer-context-names-count")))
         (dim-focus-monomer-names-vector (nc:len (nc:get-dimension cdf "focus-monomer-names-vector")))
         (dim-focus-monomer-names-start (nc:len (nc:get-dimension cdf "focus-monomer-names-start")))
         (dim-focus-monomer-names-count (nc:len (nc:get-dimension cdf "focus-monomer-names-count")))
         (dim-fragment-match-fragment-conformations-index-vector (nc:len (nc:get-dimension cdf "fragment-match-fragment-conformations-index-vector")))
         (dim-fragment-match-fragment-conformations-start-vector (nc:len (nc:get-dimension cdf "fragment-match-fragment-conformations-start-vector")))
         (dim-fragment-match-fragment-conformations-count-vector (nc:len (nc:get-dimension cdf "fragment-match-fragment-conformations-count-vector")))
         (dim-fragment-match-start-vector (nc:len (nc:get-dimension cdf "fragment-match-start-vector")))
         (dim-fragment-match-count-vector (nc:len (nc:get-dimension cdf "fragment-match-count-vector")))
         (dim-fragment-match-key-from-vector (nc:len (nc:get-dimension cdf "fragment-match-key-from-vector")))
         (dim-fragment-match-key-to-vector (nc:len (nc:get-dimension cdf "fragment-match-key-to-vector"))))
    (static-vectors:with-static-vectors (
                                         (static-internals-values dim-internals-values :element-type 'single-float)
                                         (static-probabilities-vector dim-probabilities-vector :element-type 'single-float)
                                         (static-energies-vector dim-energies-vector :element-type 'single-float)
                                         (static-trainer-index-vector dim-trainer-index-vector :element-type 'ext:byte32)
                                         (static-internals-types dim-internals-types :element-type 'base-char)
                                         (static-monomer-context-names-vector dim-monomer-context-names-vector :element-type 'base-char)
                                         (static-monomer-context-names-start dim-monomer-context-names-start :element-type 'ext:byte32)
                                         (static-monomer-context-names-count dim-monomer-context-names-count :element-type 'ext:byte32)
                                         (static-focus-monomer-names-vector dim-focus-monomer-names-vector :element-type 'base-char)
                                         (static-focus-monomer-names-start dim-focus-monomer-names-start :element-type 'ext:byte32)
                                         (static-focus-monomer-names-count dim-focus-monomer-names-count :element-type 'ext:byte32)
                                         (static-internals-index-vector dim-internals-index-vector :element-type 'ext:byte32)
                                         (static-internals-count-vector dim-internals-count-vector :element-type 'ext:byte32)
                                         (static-fragment-conformations-start-vector dim-fragment-conformations-start-vector :element-type 'ext:byte32)
                                         (static-fragment-conformations-count-vector dim-fragment-conformations-count-vector :element-type 'ext:byte32)
                                         (static-fragment-match-fragment-conformations-index-vector dim-fragment-match-fragment-conformations-index-vector :element-type 'ext:byte32)
                                         (static-fragment-match-fragment-conformations-start-vector dim-fragment-match-fragment-conformations-start-vector :element-type 'ext:byte32)
                                         (static-fragment-match-fragment-conformations-count-vector dim-fragment-match-fragment-conformations-count-vector :element-type 'ext:byte32)
                                         (static-fragment-match-start-vector dim-fragment-match-start-vector :element-type 'ext:byte32)
                                         (static-fragment-match-count-vector dim-fragment-match-count-vector :element-type 'ext:byte32)
                                         (static-fragment-match-key-from-vector dim-fragment-match-key-from-vector :element-type 'ext:byte32)
                                         (static-fragment-match-key-to-vector dim-fragment-match-key-to-vector :element-type 'ext:byte32))
      (nc:get-var-static-vector-float cdf "internals-values" static-internals-values)
      (nc:get-var-static-vector-float cdf "probabilities-vector" static-probabilities-vector)
      (nc:get-var-static-vector-float cdf "energies-vector" static-energies-vector)
      (nc:get-var-static-vector-uint cdf "trainer-index-vector" static-trainer-index-vector)
      (nc:get-var-static-vector-uchar cdf "internals-types" static-internals-types)
      (nc:get-var-static-vector-uchar cdf "monomer-context-names-vector" static-monomer-context-names-vector)
      (nc:get-var-static-vector-uint cdf "monomer-context-names-start" static-monomer-context-names-start)
      (nc:get-var-static-vector-uint cdf "monomer-context-names-count" static-monomer-context-names-count)
      (nc:get-var-static-vector-uchar cdf "focus-monomer-names-vector" static-focus-monomer-names-vector)
      (nc:get-var-static-vector-uint cdf "focus-monomer-names-start" static-focus-monomer-names-start)
      (nc:get-var-static-vector-uint cdf "focus-monomer-names-count" static-focus-monomer-names-count)
      (nc:get-var-static-vector-uint cdf "internals-index-vector" static-internals-index-vector)
      (nc:get-var-static-vector-uint cdf "internals-count-vector" static-internals-count-vector)
      (nc:get-var-static-vector-uint cdf "fragment-conformations-start-vector" static-fragment-conformations-start-vector)
      (nc:get-var-static-vector-uint cdf "fragment-conformations-count-vector" static-fragment-conformations-count-vector)
      (nc:get-var-static-vector-uint cdf "fragment-match-fragment-conformations-index-vector" static-fragment-match-fragment-conformations-index-vector )
      (nc:get-var-static-vector-uint cdf "fragment-match-fragment-conformations-start-vector" static-fragment-match-fragment-conformations-start-vector )
      (nc:get-var-static-vector-uint cdf "fragment-match-fragment-conformations-count-vector" static-fragment-match-fragment-conformations-count-vector )
      (nc:get-var-static-vector-uint cdf "fragment-match-start-vector" static-fragment-match-start-vector )
      (nc:get-var-static-vector-uint cdf "fragment-match-count-vector" static-fragment-match-count-vector )
      (nc:get-var-static-vector-uint cdf "fragment-match-key-from-vector" static-fragment-match-key-from-vector )
      (nc:get-var-static-vector-uint cdf "fragment-match-key-to-vector" static-fragment-match-key-to-vector)
      (nc:nc-close cdf)
      (let (
            (probabilities-vector (make-array (length static-probabilities-vector) :element-type 'single-float :initial-contents static-probabilities-vector))
            (energies-vector (make-array (length static-energies-vector) :element-type 'single-float :initial-contents static-energies-vector))
            (trainer-index-vector (make-array (length static-trainer-index-vector) :element-type 'ext:byte32 :initial-contents static-trainer-index-vector))
            (internals-values (make-array (length static-internals-values) :element-type 'single-float :initial-contents static-internals-values))
            (internals-types (make-array (length static-internals-types) :element-type 'base-char :initial-contents static-internals-types))
            (monomer-context-names-vector (make-array (length static-monomer-context-names-vector) :element-type 'base-char :initial-contents static-monomer-context-names-vector))
            (monomer-context-names-start (make-array (length static-monomer-context-names-start) :element-type 'ext:byte32 :initial-contents static-monomer-context-names-start))
            (monomer-context-names-count (make-array (length static-monomer-context-names-count) :element-type 'ext:byte32 :initial-contents static-monomer-context-names-count))
            (focus-monomer-names-vector (make-array (length static-focus-monomer-names-vector) :element-type 'base-char :initial-contents static-focus-monomer-names-vector))
            (focus-monomer-names-start (make-array (length static-focus-monomer-names-start) :element-type 'ext:byte32 :initial-contents static-focus-monomer-names-start))
            (focus-monomer-names-count (make-array (length static-focus-monomer-names-count) :element-type 'ext:byte32 :initial-contents static-focus-monomer-names-count))
            (internals-index-vector (make-array (length static-internals-index-vector) :element-type 'ext:byte32 :initial-contents static-internals-index-vector))
            (internals-count-vector (make-array (length static-internals-count-vector) :element-type 'ext:byte32 :initial-contents static-internals-count-vector))
            (fragment-conformations-start-vector (make-array (length static-fragment-conformations-start-vector) :element-type 'ext:byte32 :initial-contents static-fragment-conformations-start-vector :element-type 'ex))
            (fragment-conformations-count-vector (make-array (length static-fragment-conformations-count-vector) :element-type 'ext:byte32 :initial-contents static-fragment-conformations-count-vector))
            (fragment-match-fragment-conformations-index-vector (make-array (length static-fragment-match-fragment-conformations-index-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-fragment-conformations-index-vector))
            (fragment-match-fragment-conformations-start-vector (make-array (length static-fragment-match-fragment-conformations-start-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-fragment-conformations-start-vector))
            (fragment-match-fragment-conformations-count-vector (make-array (length static-fragment-match-fragment-conformations-count-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-fragment-conformations-count-vector))
            (fragment-match-start-vector (make-array (length static-fragment-match-start-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-start-vector))
            (fragment-match-count-vector (make-array (length static-fragment-match-count-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-count-vector))
            (fragment-match-key-from-vector (make-array (length static-fragment-match-key-from-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-key-from-vector))
            (fragment-match-key-to-vector (make-array (length static-fragment-match-key-to-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-key-to-vector)))
        (make-instance 'linearized-matched-fragment-conformations-holder
                       :internals-values internals-values
                       :internals-types internals-types
                       :probabilities-vector probabilities-vector
                       :energies-vector energies-vector
                       :trainer-index-vector trainer-index-vector
                       :monomer-context-names-vector monomer-context-names-vector
                       :monomer-context-names-start monomer-context-names-start
                       :monomer-context-names-count monomer-context-names-count
                       :focus-monomer-names-vector focus-monomer-names-vector
                       :focus-monomer-names-start focus-monomer-names-start
                       :focus-monomer-names-count focus-monomer-names-count
                       :internals-index-vector internals-index-vector
                       :internals-count-vector internals-count-vector
                       :fragment-conformations-start-vector fragment-conformations-start-vector
                       :fragment-conformations-count-vector fragment-conformations-count-vector
                       :fragment-match-fragment-conformations-index-vector fragment-match-fragment-conformations-index-vector
                       :fragment-match-fragment-conformations-start-vector fragment-match-fragment-conformations-start-vector
                       :fragment-match-fragment-conformations-count-vector fragment-match-fragment-conformations-count-vector
                       :fragment-match-start-vector fragment-match-start-vector
                       :fragment-match-count-vector fragment-match-count-vector
                       :fragment-match-key-from-vector fragment-match-key-from-vector
                       :fragment-match-key-to-vector fragment-match-key-to-vector
                       )))))

(defun mimic-matched-fragment-conformations-map (linearized-matched-fragment-conformations-holder)
  (let* ((monomer-context-to-index (make-hash-table))
         (monomer-context-to-fragment-conformations (make-hash-table :test 'equal))
         (monomer-contexts (make-array (length (focus-monomer-names-start linearized-matched-fragment-conformations-holder))))
         (focus-monomer-names (make-array (length (focus-monomer-names-start linearized-matched-fragment-conformations-holder)))))
    (with-slots (
                 fragment-conformations-start-vector
                 fragment-conformations-count-vector
                 internals-index-vector
                 internals-count-vector
                 internals-values
                 internals-types
                 probabilities-vector
                 energies-vector
                 trainer-index-vector
                 )
        linearized-matched-fragment-conformations-holder
      (loop with monomer-context-names-vector = (monomer-context-names-vector linearized-matched-fragment-conformations-holder)
            with monomer-context-names-start = (monomer-context-names-start linearized-matched-fragment-conformations-holder)
            with monomer-context-names-count = (monomer-context-names-count linearized-matched-fragment-conformations-holder)
            for index below (length monomer-context-names-start)
            for start across monomer-context-names-start
            for count across monomer-context-names-count
            for monomer-context-string = (make-array count
                                                     :element-type 'base-char
                                                     :displaced-to monomer-context-names-vector
                                                     :displaced-index-offset start)
            for monomer-context = (intern monomer-context-string :keyword)
            do (setf (elt monomer-contexts index) monomer-context)
            do (setf (gethash monomer-context monomer-context-to-index) index)
            do (let* ((fragment-conformations-start (elt fragment-conformations-start-vector index))
                      (fragment-conformations-count (elt fragment-conformations-count-vector index))
                      (linearized-fragment-conformations (make-array fragment-conformations-count)))
                 (loop for fci from 0 below fragment-conformations-count
                       for fcio = (+ fci fragment-conformations-start)
                       for internals-index = (elt internals-index-vector fcio)
                       for internals-count = (elt internals-count-vector fcio)
                       for the-internals-values = (make-array (* 3 internals-count)
                                                              :element-type 'single-float
                                                              :displaced-to internals-values
                                                              :displaced-index-offset internals-index)
                       for the-internals-types = (make-array internals-count
                                                            :element-type 'base-char
                                                            :displaced-to internals-types
                                                            :displaced-index-offset (/ internals-index 3))
                       for linearized-fragment-conformation = (make-instance 'linearized-fragment-internals
                                                                             :probability (elt probabilities-vector fcio)
                                                                             :energy (elt energies-vector fcio)
                                                                             :trainer-index (elt trainer-index-vector fcio)
                                                                             :internals-types the-internals-types
                                                                             :internals-values the-internals-values)
                       do (setf (elt linearized-fragment-conformations fci) linearized-fragment-conformation))
                 (setf (gethash monomer-context monomer-context-to-fragment-conformations)
                       (make-instance 'linearized-fragment-conformations
                                      :fragments linearized-fragment-conformations)))
            ))
    (loop with focus-monomer-names-vector = (focus-monomer-names-vector linearized-matched-fragment-conformations-holder)
          with focus-monomer-names-start = (focus-monomer-names-start linearized-matched-fragment-conformations-holder)
          with focus-monomer-names-count = (focus-monomer-names-count linearized-matched-fragment-conformations-holder)
          for index below (length focus-monomer-names-start)
          for start across focus-monomer-names-start
          for count across focus-monomer-names-count
          for focus-monomer-string = (make-array count
                                                 :element-type 'base-char
                                                 :displaced-to focus-monomer-names-vector
                                                 :displaced-index-offset start)
          for focus-monomer-name = (intern focus-monomer-string :keyword)
          do (setf (aref focus-monomer-names index) focus-monomer-name))
    (with-slots (fragment-match-key-from-vector
                 fragment-match-key-to-vector
                 fragment-match-start-vector
                 fragment-match-count-vector
                 fragment-match-fragment-conformations-start-vector
                 fragment-match-fragment-conformations-count-vector
                 fragment-match-fragment-conformations-index-vector
                 )
        linearized-matched-fragment-conformations-holder
      (let ((fragment-matches (make-hash-table :test 'equal)))
        (loop for ii below (length fragment-match-key-from-vector)
          for from-index = (elt fragment-match-key-from-vector ii)
          for to-index = (elt fragment-match-key-to-vector ii)
          for from-key = (elt monomer-contexts from-index)
              for to-key = (elt monomer-contexts to-index)
              for fragment-match-start = (elt fragment-match-start-vector ii)
              for fragment-match-count = (elt fragment-match-count-vector ii)
              for key = (cons from-key to-key)
              do (let ((vec (make-array fragment-match-count)))
                   (loop for oii from 0 below fragment-match-count
                         for fragment-match-fragment-conformations-start = (elt fragment-match-fragment-conformations-start-vector (+ oii fragment-match-start))
                         for fragment-match-fragment-conformations-count = (elt fragment-match-fragment-conformations-count-vector (+ oii fragment-match-start))
                         for disp = (make-array fragment-match-fragment-conformations-count
                                                :element-type 'ext:byte32
                                                :displaced-to fragment-match-fragment-conformations-index-vector
                                                :displaced-index-offset fragment-match-fragment-conformations-start)
                         do (setf (elt vec oii) disp))
                   (setf (gethash key fragment-matches) vec)))
    (make-instance 'linearized-matched-fragment-conformations-map
                   :linearized-matched-fragment-conformations-holder linearized-matched-fragment-conformations-holder
                   :monomer-context-to-fragment-conformations monomer-context-to-fragment-conformations
                   :monomer-context-to-index monomer-context-to-index
                   :focus-monomer-names focus-monomer-names
                   :fragment-matches fragment-matches)))))

(defun load-linearized-matched-fragment-conformations-map (filename)
  (let ((linearized-matched-fragment-conformations-holder (read-linearized-matched-fragment-conformations-holder filename)))
    (mimic-matched-fragment-conformations-map linearized-matched-fragment-conformations-holder)))

(defun linearize-and-save-matched-fragment-conformations-map (map filename)
  (let ((linearized (create-linearized-matched-fragment-conformations-holder map)))
    (save-linearized-matched-fragment-conformations-holder linearized filename)))

(defmethod apply-fragment-internals-to-atresidue ((fragment-internals linearized-fragment-internals) atresidue)
  (loop for joint across (joints atresidue)
        for index from 0 by 3
        for bond = (elt (internals-values fragment-internals) index)
        for angle-rad = (elt (internals-values fragment-internals) (+ 1 index))
        for dihedral-rad = (elt (internals-values fragment-internals) (+ 2 index))
        do (fill-joint-internals joint bond angle-rad dihedral-rad)))
