(in-package :topology)

#|

Linearize a data structure that looks like the following JSON file ...

{
    "foldamer-name" : "spiroligomers",
    "rotamers" : {
        [
            "csar_pro4ss" : {
                "focus-monomer-name" : "csar",
                "rotamers" : [ ; <- (A)
                               {
                                   "index" : 123,
                                   "energy" : 1.0,
                                   "probability" : 2.0,
                                   "internals" : [ 1.09, 120.0, 60.0, ... ]
                               },
                               ...
                             ]
            }
            ...
        ],
        ...
    },
    "backbone-dependent-rotamer-indices" : {
        "csar_pro4ss" : [
            [ 180, 150, 1, 2, 9, 20, 34, ... ], ; <- first two numbers are dihedral angle key, remaining are indices into (A)
            [ 180, 160, 1, 3, 5, ... ],
            [ 170, 160, 0, 3, 5, ... ],
            ...
        ],
        ...
    }
}


|#



(defclass linearized-matched-context-rotamers-holder ()
  ((internals-values :initarg :internals-values :accessor internals-values)
   (internals-types :initarg :internals-types :accessor internals-types)
   (internals-index-vector :initarg :internals-index-vector :accessor internals-index-vector)
   (internals-count-vector :initarg :internals-count-vector :accessor internals-count-vector)
   (probabilities-vector :initarg :probabilities-vector :accessor probabilities-vector)
   (energies-vector :initarg :energies-vector :accessor energies-vector)
   (trainer-index-vector :initarg :trainer-index-vector :accessor trainer-index-vector)
   (context-rotamers-start-vector :initarg :context-rotamers-start-vector :accessor context-rotamers-start-vector)
   (context-rotamers-count-vector :initarg :context-rotamers-count-vector :accessor context-rotamers-count-vector)
   (monomer-context-names-vector :initarg :monomer-context-names-vector :accessor monomer-context-names-vector)
   (monomer-context-names-start :initarg :monomer-context-names-start :accessor monomer-context-names-start)
   (monomer-context-names-count :initarg :monomer-context-names-count :accessor monomer-context-names-count)
   (focus-monomer-names-vector :initarg :focus-monomer-names-vector :accessor focus-monomer-names-vector)
   (focus-monomer-names-start :initarg :focus-monomer-names-start :accessor focus-monomer-names-start)
   (focus-monomer-names-count :initarg :focus-monomer-names-count :accessor focus-monomer-names-count)
   (fragment-match-context-rotamers-index-vector  :initarg :fragment-match-context-rotamers-index-vector  :accessor fragment-match-context-rotamers-index-vector )
   (fragment-match-context-rotamers-start-vector  :initarg :fragment-match-context-rotamers-start-vector  :accessor fragment-match-context-rotamers-start-vector )
   (fragment-match-context-rotamers-count-vector  :initarg :fragment-match-context-rotamers-count-vector  :accessor fragment-match-context-rotamers-count-vector )
   (fragment-match-start-vector  :initarg :fragment-match-start-vector  :accessor fragment-match-start-vector )
   (fragment-match-count-vector  :initarg :fragment-match-count-vector  :accessor fragment-match-count-vector )
   (fragment-match-key-from-vector  :initarg :fragment-match-key-from-vector  :accessor fragment-match-key-from-vector )
   (fragment-match-key-to-vector :initarg :fragment-match-key-to-vector :accessor fragment-match-key-to-vector)))

#|
monomer-context-to-index maps monomer-context to a monomer-context-index
monomer-context-index indexes into ...
... context-rotamers-start-vector
... context-rotamers-count-vector
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
(defclass linearized-connected-rotamers-map ()
  ((linearized-matched-context-rotamers-holder
    :accessor linearized-matched-context-rotamers-holder
    :initarg :linearized-matched-context-rotamers-holder )
   (monomer-context-to-index :accessor monomer-context-to-index
                             :initarg :monomer-context-to-index)
   (monomer-context-to-context-rotamers
    :accessor monomer-context-to-context-rotamers
    :initarg :monomer-context-to-context-rotamers)
   (focus-monomer-names :accessor focus-monomer-names
                        :initarg :focus-monomer-names)
   (rotamer-context-connections :accessor rotamer-context-onnections
                         :initarg :rotamer-context-connections)
   ))

(defclass linearized-context-rotamers ()
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
  (vector-push-extend (angle-rad internal) internals-values)
  (vector-push-extend (dihedral-rad internal) internals-values))

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

(defun create-linearized-matched-context-rotamers-holder (connected-rotamers-map &key (debug-limit 0))
  (let ((internals-values (make-array (* 1024 1024 16) :element-type 'single-float :adjustable t :fill-pointer 0))
        (internals-types (make-array (* 1024 1024 16) :element-type 'simple-char :adjustable t :fill-pointer 0))
        (monomer-context-to-index (make-hash-table))
        (monomer-context-index -1)
        (probabilities-vector (make-array (* 1024 1024 2) :element-type 'single-float :adjustable t :fill-pointer 0))
        (energies-vector (make-array (* 1024 1024 2) :element-type 'single-float :adjustable t :fill-pointer 0))
        (trainer-index-vector (make-array (* 1024 1024 2) :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (internals-index-vector (make-array (* 1024 1024 2) :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (internals-count-vector (make-array (* 1024 1024 2) :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (context-rotamers-start-vector (make-array 1024 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (context-rotamers-count-vector (make-array 1024 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (monomer-context-names-vector (make-array 65536 :element-type 'base-char :adjustable t :fill-pointer 0))
        (monomer-context-names-start (make-array 65536 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (monomer-context-names-count (make-array 65536 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (focus-monomer-names-vector (make-array 65536 :element-type 'base-char :adjustable t :fill-pointer 0))
        (focus-monomer-names-start (make-array 65536 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (focus-monomer-names-count (make-array 65536 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (dl debug-limit))
    (maphash (lambda (key context-rotamers)
               (debug-linearize-decf dl)
               (debug-linearize dl "Looking at monomer-context ~a~%" key)
               (let ((compact-fragments-start-index (length internals-index-vector)))
                 (loop for fragment-internals across (fragments context-rotamers)
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
                   (vector-push-extend compact-fragments-start-index context-rotamers-start-vector)
                   (vector-push-extend (- compact-fragments-end-index compact-fragments-start-index) context-rotamers-count-vector))
                 (setf (gethash key monomer-context-to-index) (incf monomer-context-index))
                 (vector-push-extend (length monomer-context-names-vector) monomer-context-names-start)
                 (let ((count (vector-push-cstring (string key) monomer-context-names-vector)))
                   (vector-push-extend count monomer-context-names-count))
                 (vector-push-extend (length focus-monomer-names-vector) focus-monomer-names-start)
                 (let ((count (vector-push-cstring (string (focus-monomer-name context-rotamers)) focus-monomer-names-vector)))
                   (vector-push-extend count focus-monomer-names-count))))
             (monomer-context-to-context-rotamers connected-rotamers-map))
    (let ((fragment-match-context-rotamers-index-vector (make-array (* 1024 1024 2) :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
          (fragment-match-context-rotamers-start-vector (make-array 65536 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
          (fragment-match-context-rotamers-count-vector (make-array 65536 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
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
                   (let ((fm-start (length fragment-match-context-rotamers-start-vector))) ; was fragment-match-start-vector)))
                     (debug-linearize dl "fm-start ~a~%" fm-start)
                     (vector-push-extend fm-start fragment-match-start-vector)
                     (debug-linearize dl "fragment-match-start-vector ~a~%" fm-start)
                     (loop for fragment-conformation-indices across fragment-conformation-indices-vectors
                           do (let ((fc-start (length fragment-match-context-rotamers-index-vector)))
                                (debug-linearize dl "fc-start = ~a~%" fc-start) 
                                (vector-push-extend fc-start fragment-match-context-rotamers-start-vector)
                                (debug-linearize dl "fragment-match-context-rotamers-start-vector ~a~%" fc-start)
                                (when fragment-conformation-indices
                                  (loop for rotamer-index across fragment-conformation-indices
                                        do (progn
                                             (vector-push-extend rotamer-index fragment-match-context-rotamers-index-vector)
                                             (debug-linearize dl "fragment-match-context-rotamers-index-vector ~a~%" rotamer-index)
                                        )))
                                (let ((count (- (length fragment-match-context-rotamers-index-vector) fc-start)))
                                  (vector-push-extend count fragment-match-context-rotamers-count-vector)
                                  (debug-linearize dl "fragment-match-context-rotamers-count-vector ~a~%" count)
                                  )))
                     (let ((count (- (length fragment-match-context-rotamers-start-vector) fm-start)))
                       (vector-push-extend count fragment-match-count-vector)
                       (debug-linearize dl "fragment-match-count-vector ~a~%" count))
                     )))
               (rotamer-context-connections connected-rotamers-map))
      (make-instance 'linearized-matched-context-rotamers-holder
                     :internals-values (copy-seq internals-values)
                     :internals-types (copy-seq internals-types)
                     :internals-index-vector (copy-seq internals-index-vector)
                     :internals-count-vector (copy-seq internals-count-vector)
                     :context-rotamers-start-vector (copy-seq context-rotamers-start-vector)
                     :context-rotamers-count-vector (copy-seq context-rotamers-count-vector)
                     :monomer-context-names-vector (copy-seq monomer-context-names-vector)
                     :monomer-context-names-start (copy-seq monomer-context-names-start)
                     :monomer-context-names-count (copy-seq monomer-context-names-count)
                     :focus-monomer-names-vector (copy-seq focus-monomer-names-vector)
                     :focus-monomer-names-start (copy-seq focus-monomer-names-start)
                     :focus-monomer-names-count (copy-seq focus-monomer-names-count)
                     :probabilities-vector probabilities-vector
                     :energies-vector energies-vector
                     :trainer-index-vector trainer-index-vector
                     :fragment-match-context-rotamers-index-vector  (copy-seq fragment-match-context-rotamers-index-vector)
                     :fragment-match-context-rotamers-start-vector  (copy-seq fragment-match-context-rotamers-start-vector)
                     :fragment-match-context-rotamers-count-vector  (copy-seq fragment-match-context-rotamers-count-vector)
                     :fragment-match-start-vector  (copy-seq fragment-match-start-vector)
                     :fragment-match-count-vector  (copy-seq fragment-match-count-vector)
                     :fragment-match-key-from-vector  (copy-seq fragment-match-key-from-vector)
                     :fragment-match-key-to-vector (copy-seq fragment-match-key-to-vector))
      )))

(defmacro maybe-describe (sym)
  `(progn
     #+(or)(format t "Initializing ~a~%" ',sym)
     ,sym))

(defun save-linearized-matched-context-rotamers-holder (linearized-context-rotamers-holder filename)
  (format t "save-compact-matched-context-rotamers-holder~%")
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
                   (context-rotamers-start-vector context-rotamers-start-vector)
                   (context-rotamers-count-vector context-rotamers-count-vector)
                   (fragment-match-context-rotamers-index-vector fragment-match-context-rotamers-index-vector )
                   (fragment-match-context-rotamers-start-vector fragment-match-context-rotamers-start-vector )
                   (fragment-match-context-rotamers-count-vector fragment-match-context-rotamers-count-vector )
                   (fragment-match-start-vector  fragment-match-start-vector )
                   (fragment-match-count-vector  fragment-match-count-vector )
                   (fragment-match-key-from-vector  fragment-match-key-from-vector )  ; changed
                   (fragment-match-key-to-vector fragment-match-key-to-vector)        ; changed
                   )
      linearized-context-rotamers-holder
    (static-vectors:with-static-vectors ((static-internals-values (length internals-values)
                                                                  :element-type 'single-float
                                                                  :initial-contents (maybe-describe internals-values))
                                         (static-internals-types (length internals-types) :element-type 'base-char :initial-contents (maybe-describe internals-types))
                                         (static-internals-index-vector (length internals-index-vector) :element-type 'ext:byte32 :initial-contents (maybe-describe internals-index-vector))
                                         (static-internals-count-vector (length internals-count-vector) :element-type 'ext:byte32 :initial-contents (maybe-describe internals-count-vector))
                                         (static-probabilities-vector (length probabilities-vector) :element-type 'single-float :initial-contents (maybe-describe probabilities-vector))
                                         (static-energies-vector (length energies-vector) :element-type 'single-float :initial-contents (maybe-describe energies-vector))
                                         (static-trainer-index-vector (length trainer-index-vector) :element-type 'ext:byte32 :initial-contents (maybe-describe trainer-index-vector))
                                         (static-context-rotamers-start-vector (length context-rotamers-start-vector) :element-type 'ext:byte32 :initial-contents (maybe-describe context-rotamers-start-vector))
                                         (static-context-rotamers-count-vector (length context-rotamers-count-vector) :element-type 'ext:byte32 :initial-contents (maybe-describe context-rotamers-count-vector))
                                         (static-monomer-context-names-vector (length monomer-context-names-vector) :element-type 'base-char :initial-contents (maybe-describe monomer-context-names-vector))
                                         (static-monomer-context-names-start (length monomer-context-names-start) :element-type 'ext:byte32 :initial-contents (maybe-describe monomer-context-names-start))
                                         (static-monomer-context-names-count (length monomer-context-names-count) :element-type 'ext:byte32 :initial-contents (maybe-describe monomer-context-names-count))
                                         (static-focus-monomer-names-vector (length focus-monomer-names-vector) :element-type 'base-char :initial-contents (maybe-describe focus-monomer-names-vector))
                                         (static-focus-monomer-names-start (length focus-monomer-names-start) :element-type 'ext:byte32 :initial-contents (maybe-describe focus-monomer-names-start))
                                         (static-focus-monomer-names-count (length focus-monomer-names-count) :element-type 'ext:byte32 :initial-contents (maybe-describe focus-monomer-names-count))
                                         (static-fragment-match-context-rotamers-index-vector  (length fragment-match-context-rotamers-index-vector ) :element-type 'ext:byte32 :initial-contents (maybe-describe fragment-match-context-rotamers-index-vector ))
                                         (static-fragment-match-context-rotamers-start-vector  (length fragment-match-context-rotamers-start-vector ) :element-type 'ext:byte32 :initial-contents (maybe-describe fragment-match-context-rotamers-start-vector ))
                                         (static-fragment-match-context-rotamers-count-vector  (length fragment-match-context-rotamers-count-vector ) :element-type 'ext:byte32 :initial-contents (maybe-describe fragment-match-context-rotamers-count-vector ))
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
        (nc:def-dim cdf "context-rotamers-start-vector" (length static-context-rotamers-start-vector))
        (nc:def-var cdf "context-rotamers-start-vector" nc-c:+int+ '("context-rotamers-start-vector"))
        (nc:def-dim cdf "context-rotamers-count-vector" (length static-context-rotamers-count-vector))
        (nc:def-var cdf "context-rotamers-count-vector" nc-c:+int+ '("context-rotamers-count-vector"))
        (nc:def-dim cdf "fragment-match-context-rotamers-index-vector" (length static-fragment-match-context-rotamers-index-vector ))
        (nc:def-var cdf "fragment-match-context-rotamers-index-vector" nc-c:+int+ '("fragment-match-context-rotamers-index-vector"))
        (nc:def-dim cdf "fragment-match-context-rotamers-start-vector" (length static-fragment-match-context-rotamers-start-vector ))
        (nc:def-var cdf "fragment-match-context-rotamers-start-vector" nc-c:+int+ '("fragment-match-context-rotamers-start-vector"))
        (nc:def-dim cdf "fragment-match-context-rotamers-count-vector" (length static-fragment-match-context-rotamers-count-vector ))
        (nc:def-var cdf "fragment-match-context-rotamers-count-vector" nc-c:+int+ '("fragment-match-context-rotamers-count-vector"))
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
        (nc:put-var-static-vector-uint cdf "context-rotamers-start-vector" static-context-rotamers-start-vector)
        (nc:put-var-static-vector-uint cdf "context-rotamers-count-vector" static-context-rotamers-count-vector)
        (nc:put-var-static-vector-uint cdf "fragment-match-context-rotamers-index-vector" static-fragment-match-context-rotamers-index-vector )
        (nc:put-var-static-vector-uint cdf "fragment-match-context-rotamers-start-vector" static-fragment-match-context-rotamers-start-vector )
        (nc:put-var-static-vector-uint cdf "fragment-match-context-rotamers-count-vector" static-fragment-match-context-rotamers-count-vector )
        (nc:put-var-static-vector-uint cdf "fragment-match-start-vector" static-fragment-match-start-vector )
        (nc:put-var-static-vector-uint cdf "fragment-match-count-vector" static-fragment-match-count-vector )
        (nc:put-var-static-vector-uint cdf "fragment-match-key-from-vector" static-fragment-match-key-from-vector )
        (nc:put-var-static-vector-uint cdf "fragment-match-key-to-vector" static-fragment-match-key-to-vector)
        (nc:nc-close cdf)))))



(defun read-linearized-matched-context-rotamers-holder (filename)
  (let* ((cdf (nc:nc-open (namestring (merge-pathnames filename))))
         (dim-internals-values (nc:len (nc:get-dimension cdf "internals-values")))
         (dim-probabilities-vector (nc:len (nc:get-dimension cdf "probabilities-vector")))
         (dim-energies-vector (nc:len (nc:get-dimension cdf "energies-vector")))
         (dim-trainer-index-vector (nc:len (nc:get-dimension cdf "trainer-index-vector")))
         (dim-internals-types (nc:len (nc:get-dimension cdf "internals-types")))
         (dim-internals-index-vector (nc:len (nc:get-dimension cdf "internals-index-vector")))
         (dim-internals-count-vector (nc:len (nc:get-dimension cdf "internals-count-vector")))
         (dim-context-rotamers-start-vector (nc:len (nc:get-dimension cdf "context-rotamers-start-vector")))
         (dim-context-rotamers-count-vector (nc:len (nc:get-dimension cdf "context-rotamers-count-vector")))
         (dim-monomer-context-names-vector (nc:len (nc:get-dimension cdf "monomer-context-names-vector")))
         (dim-monomer-context-names-start (nc:len (nc:get-dimension cdf "monomer-context-names-start")))
         (dim-monomer-context-names-count (nc:len (nc:get-dimension cdf "monomer-context-names-count")))
         (dim-focus-monomer-names-vector (nc:len (nc:get-dimension cdf "focus-monomer-names-vector")))
         (dim-focus-monomer-names-start (nc:len (nc:get-dimension cdf "focus-monomer-names-start")))
         (dim-focus-monomer-names-count (nc:len (nc:get-dimension cdf "focus-monomer-names-count")))
         (dim-fragment-match-context-rotamers-index-vector (nc:len (nc:get-dimension cdf "fragment-match-context-rotamers-index-vector")))
         (dim-fragment-match-context-rotamers-start-vector (nc:len (nc:get-dimension cdf "fragment-match-context-rotamers-start-vector")))
         (dim-fragment-match-context-rotamers-count-vector (nc:len (nc:get-dimension cdf "fragment-match-context-rotamers-count-vector")))
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
                                         (static-context-rotamers-start-vector dim-context-rotamers-start-vector :element-type 'ext:byte32)
                                         (static-context-rotamers-count-vector dim-context-rotamers-count-vector :element-type 'ext:byte32)
                                         (static-fragment-match-context-rotamers-index-vector dim-fragment-match-context-rotamers-index-vector :element-type 'ext:byte32)
                                         (static-fragment-match-context-rotamers-start-vector dim-fragment-match-context-rotamers-start-vector :element-type 'ext:byte32)
                                         (static-fragment-match-context-rotamers-count-vector dim-fragment-match-context-rotamers-count-vector :element-type 'ext:byte32)
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
      (nc:get-var-static-vector-uint cdf "context-rotamers-start-vector" static-context-rotamers-start-vector)
      (nc:get-var-static-vector-uint cdf "context-rotamers-count-vector" static-context-rotamers-count-vector)
      (nc:get-var-static-vector-uint cdf "fragment-match-context-rotamers-index-vector" static-fragment-match-context-rotamers-index-vector )
      (nc:get-var-static-vector-uint cdf "fragment-match-context-rotamers-start-vector" static-fragment-match-context-rotamers-start-vector )
      (nc:get-var-static-vector-uint cdf "fragment-match-context-rotamers-count-vector" static-fragment-match-context-rotamers-count-vector )
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
            (context-rotamers-start-vector (make-array (length static-context-rotamers-start-vector) :element-type 'ext:byte32 :initial-contents static-context-rotamers-start-vector :element-type 'ex))
            (context-rotamers-count-vector (make-array (length static-context-rotamers-count-vector) :element-type 'ext:byte32 :initial-contents static-context-rotamers-count-vector))
            (fragment-match-context-rotamers-index-vector (make-array (length static-fragment-match-context-rotamers-index-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-context-rotamers-index-vector))
            (fragment-match-context-rotamers-start-vector (make-array (length static-fragment-match-context-rotamers-start-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-context-rotamers-start-vector))
            (fragment-match-context-rotamers-count-vector (make-array (length static-fragment-match-context-rotamers-count-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-context-rotamers-count-vector))
            (fragment-match-start-vector (make-array (length static-fragment-match-start-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-start-vector))
            (fragment-match-count-vector (make-array (length static-fragment-match-count-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-count-vector))
            (fragment-match-key-from-vector (make-array (length static-fragment-match-key-from-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-key-from-vector))
            (fragment-match-key-to-vector (make-array (length static-fragment-match-key-to-vector) :element-type 'ext:byte32 :initial-contents static-fragment-match-key-to-vector)))
        (make-instance 'linearized-matched-context-rotamers-holder
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
                       :context-rotamers-start-vector context-rotamers-start-vector
                       :context-rotamers-count-vector context-rotamers-count-vector
                       :fragment-match-context-rotamers-index-vector fragment-match-context-rotamers-index-vector
                       :fragment-match-context-rotamers-start-vector fragment-match-context-rotamers-start-vector
                       :fragment-match-context-rotamers-count-vector fragment-match-context-rotamers-count-vector
                       :fragment-match-start-vector fragment-match-start-vector
                       :fragment-match-count-vector fragment-match-count-vector
                       :fragment-match-key-from-vector fragment-match-key-from-vector
                       :fragment-match-key-to-vector fragment-match-key-to-vector
                       )))))

(defun mimic-connected-rotamers-map (linearized-matched-context-rotamers-holder)
  (let* ((monomer-context-to-index (make-hash-table))
         (monomer-context-to-context-rotamers (make-hash-table :test 'equal))
         (monomer-contexts (make-array (length (focus-monomer-names-start linearized-matched-context-rotamers-holder))))
         (focus-monomer-names (make-array (length (focus-monomer-names-start linearized-matched-context-rotamers-holder)))))
    (with-slots (
                 context-rotamers-start-vector
                 context-rotamers-count-vector
                 internals-index-vector
                 internals-count-vector
                 internals-values
                 internals-types
                 probabilities-vector
                 energies-vector
                 trainer-index-vector
                 )
        linearized-matched-context-rotamers-holder
      (loop with monomer-context-names-vector = (monomer-context-names-vector linearized-matched-context-rotamers-holder)
            with monomer-context-names-start = (monomer-context-names-start linearized-matched-context-rotamers-holder)
            with monomer-context-names-count = (monomer-context-names-count linearized-matched-context-rotamers-holder)
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
            do (let* ((context-rotamers-start (elt context-rotamers-start-vector index))
                      (context-rotamers-count (elt context-rotamers-count-vector index))
                      (linearized-context-rotamers (make-array context-rotamers-count)))
                 (loop for fci from 0 below context-rotamers-count
                       for fcio = (+ fci context-rotamers-start)
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
                       do (setf (elt linearized-context-rotamers fci) linearized-fragment-conformation))
                 (setf (gethash monomer-context monomer-context-to-context-rotamers)
                       (make-instance 'linearized-context-rotamers
                                      :fragments linearized-context-rotamers)))
            ))
    (loop with focus-monomer-names-vector = (focus-monomer-names-vector linearized-matched-context-rotamers-holder)
          with focus-monomer-names-start = (focus-monomer-names-start linearized-matched-context-rotamers-holder)
          with focus-monomer-names-count = (focus-monomer-names-count linearized-matched-context-rotamers-holder)
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
                 fragment-match-context-rotamers-start-vector
                 fragment-match-context-rotamers-count-vector
                 fragment-match-context-rotamers-index-vector
                 )
        linearized-matched-context-rotamers-holder
      (let ((rotamer-context-connections (make-rotamer-context-connections)))
        (loop for ii below (length fragment-match-key-from-vector)
              for from-index = (elt fragment-match-key-from-vector ii)
              for to-index = (elt fragment-match-key-to-vector ii)
              for from-key = (elt monomer-contexts from-index)
              for to-key = (elt monomer-contexts to-index)
              for fragment-match-start = (elt fragment-match-start-vector ii)
              for fragment-match-count = (elt fragment-match-count-vector ii)
              do (let ((vec (progn
                              (error "the following make-array needs to be a make-rotamer-shape-connections ")
                              (make-array fragment-match-count))))
                   (loop for oii from 0 below fragment-match-count
                         for fragment-match-context-rotamers-start = (elt fragment-match-context-rotamers-start-vector (+ oii fragment-match-start))
                         for fragment-match-context-rotamers-count = (elt fragment-match-context-rotamers-count-vector (+ oii fragment-match-start))
                         for disp = (make-array fragment-match-context-rotamers-count
                                                :element-type 'ext:byte32
                                                :displaced-to fragment-match-context-rotamers-index-vector
                                                :displaced-index-offset fragment-match-context-rotamers-start)
                         do (setf (elt vec oii) disp))
                   (set-rotamer-context-connections rotamer-context-connections from-key to-key vec)
                   ))
    (make-instance 'linearized-connected-rotamers-map
                   :linearized-matched-context-rotamers-holder linearized-matched-context-rotamers-holder
                   :monomer-context-to-context-rotamers monomer-context-to-context-rotamers
                   :monomer-context-to-index monomer-context-to-index
                   :focus-monomer-names focus-monomer-names
                   :rotamer-context-connections rotamer-context-connections)))))

(defun load-linearized-connected-rotamers-map (filename)
  (let ((linearized-matched-context-rotamers-holder (read-linearized-matched-context-rotamers-holder filename)))
    (mimic-connected-rotamers-map linearized-matched-context-rotamers-holder)))

(defun linearize-and-save-connected-rotamers-map (map filename)
  (let ((linearized (create-linearized-matched-context-rotamers-holder map)))
    (save-linearized-matched-context-rotamers-holder linearized filename)))

(defmethod apply-fragment-internals-to-atresidue ((fragment-internals linearized-fragment-internals) atresidue)
  (loop for joint across (joints atresidue)
        for index from 0 by 3
        for bond = (elt (internals-values fragment-internals) index)
        for angle-rad = (elt (internals-values fragment-internals) (+ 1 index))
        for dihedral-rad = (elt (internals-values fragment-internals) (+ 2 index))
        do (fill-joint-internals joint bond angle-rad dihedral-rad)))
