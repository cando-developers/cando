(in-package :topology)

(defconstant +phi+   :phi)
(defconstant +phi-1+ :phi-1)
(defconstant +psi+   :psi)
(defconstant +psi-1+ :psi-1)
(defun shape-key-dihedral-name-p (name)
  (member name (list +phi+ +phi-1+ +psi+ +psi-1+)))

(defclass monomer-shape-info ()
  ((monomer :initarg :monomer :accessor monomer)
   (monomer-context :initarg :monomer-context :accessor monomer-context)
   (monomer-shape-kind :initarg :monomer-shape-kind :accessor monomer-shape-kind)
   ))

(defmethod print-object ((object monomer-shape-info) stream)
  (print-unreadable-object (object stream :type t)
    (format stream "~s" (monomer-context object))))


(defclass monomer-shape ()
  (()))

(defgeneric copy-monomer-shape (monomer-shape))

(defclass rotamer-shape (monomer-shape)
  ((rotamer-index :initarg :rotamer-index :accessor rotamer-index)
   (rotamer-shape-context-rotamers :initarg :rotamer-shape-context-rotamers :accessor rotamer-shape-context-rotamers))
  (:documentation "Monomer-shape that chooses between different rotamers"))

(defmethod print-object ((obj rotamer-shape) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a"
            (if (slot-boundp obj 'rotamer-index)
                (rotamer-index obj)
                "#<unbound>")
            )))

(defmethod copy-monomer-shape ((monomer-shape rotamer-shape))
  (if (slot-boundp monomer-shape 'rotamer-index)
      (make-instance 'rotamer-shape :rotamer-index (rotamer-index monomer-shape)
                     :rotamer-shape-context-rotamers (rotamer-shape-context-rotamers monomer-shape))
      (make-instance 'rotamer-shape
                     :rotamer-shape-context-rotamers (rotamer-shape-context-rotamers monomer-shape))))

(defun make-rotamer-shape (monomer oligomer &key ignore-sidechains)
  "Convenience function to make a rotamer-shape for a MONOMER in an OLIGOMER.
If ORIGINAL-ROTAMER-SHAPE is defined then it must be a ROTAMER-SHAPE and we copy its ROTAMER-INDEX."
  (let* ((foldamer (foldamer (oligomer-space oligomer)))
         (monomer-context (foldamer-monomer-context monomer oligomer foldamer))
         (rotamers-database (foldamer-rotamers-database foldamer))
         (context-rotamers (topology:lookup-rotamers-for-context rotamers-database monomer-context :ignore-sidechains ignore-sidechains :errorp nil :error-value :no-context-rotamers)))
    (if (eq context-rotamers :no-context-rotamers)
        :no-rotamer-shape
        (make-instance 'rotamer-shape :rotamer-shape-context-rotamers context-rotamers))))

(defmethod apply-monomer-shape-to-atresidue-internals (assembler oligomer-shape (rotamer-shape rotamer-shape) monomer-context atresidue &key verbose)
  (when verbose
    (let ((*print-pretty* nil))
      (format t "apply-monomer-shape-to-atresidue-internals ~%   oligomer-shape ~s~%   rotamer-shape ~s~%  monomer-context ~s~%"
              oligomer-shape rotamer-shape monomer-context)))
  (let* ((rotamers-database (rotamers-database oligomer-shape))
         (context-rotamers (let ((rots (gethash monomer-context (monomer-context-to-context-rotamers rotamers-database))))
                             (unless rots (error "Could not find monomer-context ~a in rotamers-database" monomer-context))
                             rots))
         (rotamer-index (rotamer-index rotamer-shape))
         (rotamers (let ((rots (rotamer-vector context-rotamers)))
                     (unless (> (length rots) 0)
                       (error "rots is empty for context ~a" monomer-context))
                     (unless (< rotamer-index (length rots))
                       (error "Check rots ~a because rotamer-index ~a is out of bounds" rots rotamer-index))
                     (elt rots rotamer-index))))
    (setf (rotamer-index atresidue) rotamer-index)
    (apply-fragment-internals-to-atresidue assembler rotamers rotamer-index atresidue)))

(defmethod apply-monomer-shape-to-atresidue-internals (assembler oligomer-shape (no-rotamer-shape symbol) monomer-context atresidue &key verbose)
  (format t "Entered apply-monomer-shape-to-atresidue-internals with the no-rotamer-shape symbol")
  )

(defclass residue-shape (monomer-shape)
  ((residue :initarg :residue :reader residue)
   (molecule :initarg :molecule :reader molecule)))

(defclass backbone-residue-shape (residue-shape)
  ((backbone-dihedral-cache-deg :initarg :backbone-dihedral-cache-deg
                                :accessor backbone-dihedral-cache-deg)))

(defmethod copy-monomer-shape ((residue-shape residue-shape))
  ;;Don't copy anything - it is immutable
  residue-shape)

(defmethod print-object ((obj residue-shape) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~s" (chem:get-name (residue obj)))))


(defmethod apply-monomer-shape-to-atresidue-internals (assembler oligomer-shape (residue-shape residue-shape) monomer-context atresidue &key verbose)
  (when verbose
    (let ((*print-pretty* nil))
      (format t "apply-monomer-shape-to-atresidue-internals ~%   oligomer-shape ~s~%   residue-shape ~s~%  monomer-context ~s~%" oligomer-shape residue-shape monomer-context)))
  (let ((residue (residue residue-shape)))
    (loop for joint across (joints atresidue)
          for name = (kin:joint/name joint)
          for atm = (chem:first-atom-with-name residue name)
          do (kin:xyz-joint/set-atom joint atm))))

(defclass backbone-monomer-shape (residue-shape)
  ()
  (:documentation "monomer-shape for protein backbone monomers"))

(defmethod copy-monomer-shape ((backbone-monomer-shape backbone-monomer-shape))
  (if (slot-boundp backbone-monomer-shape 'topology:rotamer-index)
      (make-instance 'backbone-monomer-shape
                     :rotamer-index (rotamer-index backbone-monomer-shape)
                     :residue (residue backbone-monomer-shape))
      (make-instance 'backbone-monomer-shape
                     :residue (residue backbone-monomer-shape))))



(defclass oligomer-shape ()
  ((name :initarg :name :reader name)
   (rotamers-state :initarg :rotamers-state :accessor rotamers-state
                   :type (member :undefined :incomplete-no-rotamers :incomplete-backbone-rotamers :complete-sidechain-and-backbone-rotamers)
                   :documentation "Mutable slot")
   (oligomer :initarg :oligomer :reader oligomer)
   (rotamers-database :initarg :rotamers-database :reader rotamers-database)
   (monomer-shape-info-vector :initarg :monomer-shape-info-vector :reader monomer-shape-info-vector)
   (the-root-monomer :initarg :the-root-monomer :reader the-root-monomer)
   (in-monomers :initarg :in-monomers :reader in-monomers)
   (out-monomers :initarg :out-monomers :reader out-monomers)
   (monomer-shape-vector :initarg :monomer-shape-vector :accessor monomer-shape-vector
                         :documentation "This is mutable")
   (monomer-shape-to-index :initarg :monomer-shape-to-index :reader monomer-shape-to-index)
   (monomer-shape-map :initarg :monomer-shape-map :reader monomer-shape-map)
   ))

(defmacro do-oligomer-shape ((monomer-shape monomer monomer-context monomer-shape-kind) oligomer-shape &body body)
  "Walk the components of the OLIGOMER-SHAPE.  Successively bind MONOMER-SHAPE, MONOMER, MONOMER-CONTEXT, MONOMER-SHAPE-KIND
to the components and evaluate the BODY in that scope."
  (let ((monomer-shape-info (gensym)))
  `(loop for ,monomer-shape across (monomer-shape-vector ,oligomer-shape)
         for ,monomer-shape-info across (monomer-shape-info-vector ,oligomer-shape)
         for ,monomer = (monomer ,monomer-shape-info)
         for ,monomer-context = (monomer-context ,monomer-shape-info)
         for ,monomer-shape-kind = (monomer-shape-kind ,monomer-shape-info)
         do (progn
              ,@body))))

(defun ensure-valid-rotamers-state (state)
  (unless (member state
                  '(:incomplete-no-rotamers :incomplete-backbone-rotamers :complete-sidechain-and-backbone-rotamers))
    (error "Illegal rotamers state ~s - must be one of ~s"
           state
           '(:incomplete-no-rotamers :incomplete-backbone-rotamers :complete-sidechain-and-backbone-rotamers)
           ))
  state)

(defun set-rotamers-state (oligomer-shape state)
  (ensure-valid-rotamers-state state)
  (setf (rotamers-state oligomer-shape) state)
  state)

(defgeneric copy-oligomer-shape (oligomer-shape))

(defmethod copy-oligomer-shape ((oligomer-shape oligomer-shape))
  (let ((monomer-shape-vector (make-array (length (monomer-shape-vector oligomer-shape))))
        (monomer-shape-map (make-hash-table))
        (monomer-shape-to-index (make-hash-table)))
    (loop for monomer-shape across (monomer-shape-vector oligomer-shape)
          for index from 0
          for monomer-shape-info across (monomer-shape-info-vector oligomer-shape)
          for monomer = (monomer monomer-shape-info)
          for new-monomer-shape = (copy-monomer-shape monomer-shape)
          do (setf (aref monomer-shape-vector index) new-monomer-shape
                   (gethash monomer monomer-shape-map) new-monomer-shape
                   (gethash monomer-shape monomer-shape-to-index) index))
    (make-instance 'oligomer-shape
                   :name (name oligomer-shape)
                   :rotamers-state (rotamers-state oligomer-shape)
                   :oligomer (oligomer oligomer-shape)
                   :rotamers-database (rotamers-database oligomer-shape)
                   :monomer-shape-info-vector (monomer-shape-info-vector oligomer-shape)
                   :the-root-monomer (the-root-monomer oligomer-shape)
                   :monomer-shape-to-index monomer-shape-to-index
                   :in-monomers  (in-monomers oligomer-shape)
                   :out-monomers (out-monomers oligomer-shape)
                   :monomer-shape-vector monomer-shape-vector
                   :monomer-shape-map monomer-shape-map)))

(defmethod oligomer-space ((oligomer-shape oligomer-shape))
  "Return the oligomer-space for the oligomer of the OLIGOMER-SHAPE"
  (oligomer-space (oligomer oligomer-shape)))

(defmethod print-object ((obj oligomer-shape) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~s" (name obj))))

(defmethod orientation ((oligomer-shape oligomer-shape) (orientations orientations))
  (gethash oligomer-shape (orientations orientations)))

(defmethod orientation ((oligomer-shape oligomer-shape) (assembler assembler))
  (gethash oligomer-shape (orientations (orientations assembler))))

(defvar *orientation*)

(defmethod kin:orientation-transform ((jump-joint kin:jump-joint))
  (unless (boundp '*orientation*)
    (error "~s must be bound for kin:orientation-transform to work" '*orientation*))
  (kin:orientation-transform *orientation*))

(defmethod kin:orientation-transform ((joint kin:xyz-joint))
  (unless (boundp '*orientation*)
    (error "For joint ~s ~s must be bound for kin:orientation-transform to work" joint '*orientation*))
  (kin:orientation-transform *orientation*))

(defmethod kin:orientation-transform ((joint kin:stub-joint))
  (unless (boundp '*orientation*)
    (error "For joint ~s ~s must be bound for kin:orientation-transform to work" joint '*orientation*))
  (kin:orientation-transform *orientation*))

(defmacro with-orientation (orientation &body body)
  "Set the *orientation* dynamic variable so that external coordinate building can function"
  `(progn
     (if (boundp '*orientation*)
         (error "*orientation* is already bound - with-orientation cannot be nested")
         (let ((*orientation* ,orientation))
           (progn
             ,@body)))))

(defclass receptor-shape (oligomer-shape)
  ((aggregate :initarg :aggregate :accessor aggregate)
    ))

(defgeneric monomer-kind (foldamer monomer oligomer-space))

(defun sorted-build-order (foldamer oligomer kind-order)
  (let* ((number-of-monomers (length (topology:monomers oligomer)))
         (canonical-sequence (topology:ordered-monomers oligomer))
         (oligomer-space (oligomer-space oligomer))
         (foldamer (foldamer oligomer-space))
         (indexed-sequence (loop for monomer in canonical-sequence
                                 for index from 0
                                 for monomer-context-kind = (shape-kind foldamer monomer oligomer-space)
                                 for monomer-kind-pos = (position monomer-context-kind kind-order)
                                 collect (cons (+ (* number-of-monomers monomer-kind-pos) index) monomer))))
    (let ((sorted-monomers (sort indexed-sequence #'< :key #'car)))
      (mapcar #'cdr sorted-monomers))))

(defun needs-build-order (oligomer)
  (let ((max-out-plugs (loop for monomer across (monomers oligomer)
                             maximize (let ((out-couplings 0))
                                        (maphash (lambda (name coupling)
                                                   (unless (is-in-plug-name name)
                                                     (incf out-couplings)))
                                                 (couplings monomer))
                                        out-couplings))))
    (> max-out-plugs 1)))

(defclass kind-key ()
  ((kind :initarg :kind :accessor kind)
   (keys :initarg :keys :accessor keys)
   (properties :initarg :properties :accessor properties)))

(cando.serialize:make-class-save-load kind-key)

(defclass shape-info ()
  ((kind-keys :initarg :kind-keys :accessor kind-keys)
   (in-plug-names :initarg :in-plug-names :accessor in-plug-names)))

(cando.serialize:make-class-save-load shape-info)

(defun make-shape-info (info)
  (let ((kind-keys (mapcar (lambda (kind-key-list)
                                   (make-instance 'kind-key
                                                  :kind (car kind-key-list)
                                                  :keys (cadr kind-key-list)
                                                  :properties (cddr kind-key-list)))
                                 info)))
    (make-instance 'shape-info
                   :kind-keys kind-keys)))

(defgeneric shape-info (foldamer))

(defgeneric shape-kind (foldamer monomer oligomer))

(defun build-shape-key (shape-info monomer-shape-map shape-kind monomer oligomer in-monomers out-monomers)
  (let ((key-names (loop for kind-keys in (kind-keys shape-info)
                         when (eq shape-kind (kind kind-keys))
                           do (return (car (keys kind-keys))))))
    (cond
      ((null key-names) nil)
      ((and (= (length key-names) 1))
       (let* ((in-plug-name (car key-names))
              (in-plug-monomer (monomer-on-other-side monomer in-plug-name))
              (in-plug-monomer-shape (gethash in-plug-monomer monomer-shape-map)))
         in-plug-monomer-shape))
      ((and (= (length key-names) 2))
       (let* ((in-plug-name (car key-names))
              (out-plug-name (cadr key-names))
              (in-plug-monomer (monomer-on-other-side monomer in-plug-name))
              (other-monomer (monomer-on-other-side in-plug-monomer out-plug-name))
              (in-plug-monomer-shape (gethash in-plug-monomer monomer-shape-map))
              (other-monomer-shape (gethash other-monomer monomer-shape-map)))
         (cons in-plug-monomer-shape other-monomer-shape)))
      (t (error "Unrecognized shape-kind ~a" shape-kind)))))

(defgeneric make-oligomer-shape (oligomer-or-space &key &allow-other-keys))

(defmethod make-oligomer-shape ((oligomer oligomer)
                                &key
                                  (oligomer-shape-class-name 'oligomer-shape)
                                  monomer-to-monomer-shape-map
                                  monomer-shape-factory
                                  callback-backbone-rotamer-indexes
                                  callback-sidechain-rotamer-indexes
                                  extra-arguments
                                  name)
  "Make an oligomer-shape for the OLIGOMER.  Use the MONOMER-TO-MONOMER-SHAPE-MAP to lookup monomer-shapes for each monomer.
Use the CALLBACK-BACKBONE-ROTAMER-INDEXES and CALLBACK-SIDECHAIN-ROTAMER-INDEXES to set rotamer-indexes of any rotamer-shapes."
  (check-type callback-backbone-rotamer-indexes (or null function))
  (check-type callback-sidechain-rotamer-indexes (or null function))
  (let* ((foldamer (topology:foldamer (topology:oligomer-space oligomer)))
         (rotamers-database (topology:foldamer-rotamers-database foldamer))
         (shape-info (shape-info foldamer))
         #+(or)(kind-order (loop for kind-keys in (kind-keys shape-info)
                                 collect (kind kind-keys))))
    (multiple-value-bind (monomer-shape-vector monomer-shape-info-vector monomer-shape-to-index
                                               the-root-monomer in-monomers out-monomers monomer-shape-map)
        (loop with monomer-shape-vector = (make-array (length (monomers oligomer)))
              with monomer-shape-info-vector = (make-array (length (monomers oligomer)))
              with monomer-shape-to-index = (make-hash-table)
              with in-monomers = (make-hash-table)
              with out-monomers = (make-hash-table)
              with the-root-monomer = nil
              with monomer-shape-map = (make-hash-table)
              for index from 0
              for monomer across (monomers (oligomer-space oligomer)) ; (ordered-monomers oligomer)
              for monomer-context = (topology:foldamer-monomer-context monomer oligomer foldamer)
              ;; The following is not good - to not define a slot when rotamers-database is NIL
              ;; for context-rotamers = (topology:lookup-rotamers-for-context rotamers-database monomer-context)
              for shape-kind = (topology:shape-kind foldamer monomer oligomer)
              for couplings = (couplings monomer)
              for in-monomer = (let (in-monomer)
                                 (maphash (lambda (key coupling)
                                            (when (in-plug-name-p key)
                                              (setf in-monomer (topology:source-monomer coupling))
                                              (setf (gethash monomer in-monomers) (topology:source-monomer coupling))))
                                          couplings)
                                 in-monomer)
              for out-mons = (let (out-monomers)
                               (maphash (lambda (key coupling)
                                          (cond
                                            ((in-plug-name-p key)
                                             )
                                            ((typep coupling 'topology::ring-coupling)
                                             (push (topology:other-monomer coupling monomer) out-monomers))
                                            (t 
                                             (push (topology:target-monomer coupling) out-monomers)
                                             #+(or)(format t "Out plug coupling ~a ~a~%" key coupling))))
                                        couplings)
                               out-monomers)
              for in-monomer-context = (if in-monomer
                                           (topology:foldamer-monomer-context in-monomer oligomer foldamer)
                                           nil)
              for monomer-shape-info = (make-instance 'monomer-shape-info
                                                      :monomer-context monomer-context
                                                      :monomer-shape-kind shape-kind
                                                      :monomer monomer
                                                      )
              for monomer-shape = (if monomer-to-monomer-shape-map
                                      (gethash monomer monomer-to-monomer-shape-map)
                                      (make-rotamer-shape monomer oligomer))
              do (setf (gethash monomer monomer-shape-map) monomer-shape)
              do (unless in-monomer (setf the-root-monomer monomer))
              do (setf (gethash monomer out-monomers) out-mons)
              do (setf (aref monomer-shape-vector index) monomer-shape)
              do (setf (gethash monomer-shape monomer-shape-to-index) index)
              do (setf (aref monomer-shape-info-vector index) monomer-shape-info)
                 ;;            do (format t "monomer-context ~a~%" monomer-context)
              finally (return (values monomer-shape-vector monomer-shape-info-vector monomer-shape-to-index the-root-monomer in-monomers out-monomers monomer-shape-map)))
      (let* ((os (apply 'make-instance
                        oligomer-shape-class-name
                        :oligomer oligomer
                        :rotamers-database rotamers-database
                        :monomer-shape-info-vector monomer-shape-info-vector
                        :monomer-shape-vector monomer-shape-vector
                        :monomer-shape-to-index monomer-shape-to-index
                        :monomer-shape-map monomer-shape-map
                        :the-root-monomer the-root-monomer
                        :in-monomers in-monomers
                        :out-monomers out-monomers
                        :name name
                        :rotamers-state (ensure-valid-rotamers-state :incomplete-no-rotamers)
                        extra-arguments)))
        (set-rotamers os :callback-backbone-rotamer-indexes callback-backbone-rotamer-indexes
                         :callback-sidechain-rotamer-indexes callback-sidechain-rotamer-indexes)
        os))))


(defmethod make-oligomer-shape-of-rotamer-shapes ((oligomer oligomer)
                                                  &rest args
                                                  &key (name :default)
                                                    uninitialized
                                                    ignore-sidechains
                                                    )
  "Build an oligomer-shape using only rotamer-shapes for OLIGOMER.
By default initialize the rotamer-indexes to random allowed values.
If UNINITIALIZED then leave them unbound."
  (let ((monomer-to-monomer-shape-map (loop with map = (make-hash-table)
                                            for monomer across (monomers (oligomer-space oligomer))
                                            for monomer-shape = (make-rotamer-shape monomer oligomer :ignore-sidechains ignore-sidechains)
                                            do (setf (gethash monomer map) monomer-shape)
                                            finally (return map))))
    (if uninitialized
        (let ((args (remf args :uninitialized)))
          (apply 'make-oligomer-shape oligomer :monomer-to-monomer-shape-map monomer-to-monomer-shape-map
                 args))
        (apply 'make-oligomer-shape oligomer :monomer-to-monomer-shape-map monomer-to-monomer-shape-map
                                             :callback-backbone-rotamer-indexes #'random-rotamers
                                             :callback-sidechain-rotamer-indexes #'random-rotamers
                                             args)
        )))

(defun set-rotamers (oligomer-shape
                     &key
                       callback-backbone-rotamer-indexes
                       callback-sidechain-rotamer-indexes
                       )
  (let* ((permissible-backbone-rotamers (make-permissible-backbone-rotamers oligomer-shape))
         (backbone-rotamer-indexes (if callback-backbone-rotamer-indexes
                                       (funcall callback-backbone-rotamer-indexes
                                                permissible-backbone-rotamers))))
    (write-rotamers oligomer-shape permissible-backbone-rotamers backbone-rotamer-indexes)
    (let* ((permissible-sidechain-rotamers (make-permissible-sidechain-rotamers oligomer-shape))
           (sidechain-rotamer-indexes (if callback-sidechain-rotamer-indexes
                                          (funcall callback-sidechain-rotamer-indexes
                                                   permissible-sidechain-rotamers))))
      (write-rotamers oligomer-shape permissible-sidechain-rotamers sidechain-rotamer-indexes)
      (values backbone-rotamer-indexes sidechain-rotamer-indexes))))



#+(or)
(defmethod make-oligomer-shape ((oligomer-space oligomer-space)
                                &rest args
                                &key (oligomer-index 0)
                                  monomer-to-monomer-shape-map
                                  (oligomer-shape-class-name 'oligomer-shape)
                                  monomer-shape-factory
                                  callback-backbone-rotamer-indexes
                                  callback-sidechain-rotamer-indexes
                                  extra-arguments
                                  name)
  (let ((oligomer (make-oligomer oligomer-space oligomer-index)))
    (remf args :oligomer-index)
    (apply 'make-oligomer-shape oligomer rotamers-database args)))


(defgeneric read-oligomer-shape-rotamers (oligomer-shape)
  (:documentation "Read the rotamers from the object"))

(defgeneric write-oligomer-shape-rotamers (oligomer-shape vec)
  (:documentation "Write the rotamers into the object"))

(defmethod read-oligomer-shape-rotamers ((obj oligomer-shape))
  (let ((vec (make-array (length (monomer-shape-vector obj)) :element-type 'ext:byte32)))
    (map-into vec #'rotamer-index (monomer-shape-vector obj))
    vec))

(defmethod write-oligomer-shape-rotamers ((obj oligomer-shape) vec)
  (loop for rotamer-index across vec
        for monomer-shape across (monomer-shape-vector obj)
        do (setf (rotamer-index monomer-shape) rotamer-index)))

(defun all-monomers-impl (root shape)
  #+(or)(format t "monomer ~a in: ~a~%" root (gethash root (in-monomers shape)))
  (let ((out-monomers (gethash root (out-monomers shape))))
    (loop for out-monomer in out-monomers
          do (all-monomers-impl out-monomer shape))))

(defun all-monomers (shape)
  (let ((root (the-root-monomer shape)))
    (all-monomers-impl root shape)))

(defun random-rotamer-index-impl (root-monomer-shape oligomer-shape)
  (let ((out-monomers (gethash (monomer root-monomer-shape) (out-monomers oligomer-shape))))
    (loop for out-monomer in out-monomers
          for out-monomer-shape = (gethash out-monomer (monomer-shape-map oligomer-shape))
          for fragment-match-key = (cons (monomer-context root-monomer-shape) (monomer-context out-monomer-shape))
          for allowed-fragment-vec = (let* ((thing (rotamers-database oligomer-shape))
                                            (ht (topology:context-to-rotamers thing))
                                            (val (gethash fragment-match-key ht)))
                                       (unless val (break "Could not find value for key ~a in ht: ~a" fragment-match-key ht))
                                       val)
          for allowed-fragment-indexes = (progn
                                           #+(or)(format t "alowed-fragment-vec ~s (rotamer-index root-monomer-shape) -> ~s~%" allowed-fragment-vec (rotamer-index root-monomer-shape))
                                           (elt allowed-fragment-vec (rotamer-index root-monomer-shape)))
          for rotamer-index = (if allowed-fragment-indexes
                                                (elt allowed-fragment-indexes (random (length allowed-fragment-indexes)))
                                                :BADBADBAD)
          do (setf (rotamer-index out-monomer-shape) rotamer-index)
          do (format t "rotamer-index ~a for monomer-shape ~a~%"
                     rotamer-index out-monomer-shape)
          do (unless (< (rotamer-index out-monomer-shape)
                        (length (topology:fragments (context-rotamers out-monomer-shape))))
               (error "rotamer-index ~a is out of bounds ~a"
                      (rotamer-index root-monomer-shape)
                      (length (topology:fragments (context-rotamers out-monomer-shape)))))
          do (random-rotamer-index-impl out-monomer-shape oligomer-shape))))

(defun random-rotamer-index (oligomer-shape)
  (let* ((root (the-root-monomer oligomer-shape))
         (root-monomer-shape (gethash root (monomer-shape-map oligomer-shape))))
    (format t "root-monomer-shape ~a~%" root-monomer-shape)
    (let* ((context-rotamers (context-rotamers root-monomer-shape)))
      (format t "context-rotamers: ~a~%" context-rotamers)
      (setf (rotamer-index root-monomer-shape)
            (random (length (topology:rotamers context-rotamers))))
      (unless (< (rotamer-index root-monomer-shape)
                 (length (topology:rotamers context-rotamers)))
        (error "rotamer-index ~a is out of bounds ~a"
               (rotamer-index root-monomer-shape)
               (length (topology:rotamers context-rotamers))))
      (random-rotamer-index-impl root-monomer-shape oligomer-shape))))


#+(or)
(defun build-shapes (oligomer-shapes assembler &key monomer-order)
  (let ((coordinates (chem:make-coordinates (topology:energy-function assembler))))
    (chem:energy-function/load-coordinates-into-vector (topology:energy-function assembler) coordinates)
    (loop for oligomer-shape in oligomer-shapes
          for oligomer = (oligomer oligomer-shape)
          do (topology:fill-internals-from-oligomer-shape assembler oligomer-shape)
             (if monomer-order
                 (loop for monomer in monomer-order
                       for monomer-position = (gethash monomer (monomer-positions assembler))
                       for molecule-index = (molecule-index monomer-position)
                       for residue-index = (residue-index monomer-position)
                       for atmolecule = (aref (ataggregate assembler) molecule-index)
                       for atresidue = (aref (atresidues atmolecule) residue-index)
                       do (build-atresidue-atom-tree-external-coordinates atresidue coordinates))
                 (topology:build-all-atom-tree-external-coordinates assembler oligomer coordinates)))
    (chem:energy-function/save-coordinates-from-vector (topology:energy-function assembler) coordinates)
    (topology:aggregate assembler)))


(defun extract-shape-key (monomer-shape)
  (let ((shape-key (shape-key monomer-shape)))
    (cond
      ((null shape-key)
       nil)
      ((eq shape-key :phi/psi)
       (cons (rotamer-index (car shape-key))
             (rotamer-index (cdr shape-key))))
      (t (rotamer-index shape-key)))))


(defun lookup-dihedral-joint-impl (assembler oligomer monomer dihedral-name)
  (let* ((monomer-index (position monomer (monomers oligomer)))
         (monomer-name (oligomer-monomer-name oligomer monomer))
         (topology (chem:find-topology monomer-name))
         (constitution (constitution topology))
         (residue-properties (residue-properties constitution))
         (dihedrals (getf residue-properties :dihedrals)))
    (cond
      (dihedrals
       (let* ((dihedral-info (find dihedral-name dihedrals :key #'name)))
         (etypecase dihedral-info
           (dihedral-info-external
            (loop for plug-name in (plug-path dihedral-info)
                  do (setf monomer (monomer-on-other-side monomer plug-name)))
            (lookup-dihedral-joint-impl assembler oligomer monomer (external-dihedral-name dihedral-info)))
           (dihedral-info-atom
            (let* ((monomer-pos (gethash monomer (monomer-positions assembler)))
                   (atagg (ataggregate assembler))
                   (atmol (aref (atmolecules atagg) (molecule-index monomer-pos)))
                   (atres (aref (atresidues atmol) (residue-index monomer-pos)))
                   (atom-name (atom-name dihedral-info))
                   (joint (joint-with-name atres atom-name)))
              joint)))))
      (t (error "there was no :dihedrals property in ~s" topology)))))

(defun lookup-dihedral-joint (assembler oligomer monomer dihedral-name)
  "Lookup the joint in the assembler corresponding to the DIHEDRAL-NAME"
  (lookup-dihedral-joint-impl assembler oligomer monomer dihedral-name))


(defun lookup-dihedral-cache-impl (oligomer-shape monomer-shape dihedral-name &key ignore-degrees)
  (let* ((monomer-shape-index (position monomer-shape (monomer-shape-vector oligomer-shape)))
         (monomer-shape-info (aref (monomer-shape-info-vector oligomer-shape) monomer-shape-index))
         (monomer (monomer monomer-shape-info))
         (oligomer (oligomer oligomer-shape))
         (monomer-name (oligomer-monomer-name oligomer monomer))
         (topology (chem:find-topology monomer-name))
         (constitution (constitution topology))
         (residue-properties (residue-properties constitution))
         (dihedrals (getf residue-properties :dihedrals)))
    (cond
      (dihedrals
       (let* ((dihedral-info (find dihedral-name dihedrals :key #'name)))
         (etypecase dihedral-info
           (dihedral-info-external
            (loop for plug-name in (plug-path dihedral-info)
                  do (setf monomer (monomer-on-other-side monomer plug-name)))
            (let ((other-monomer-shape (gethash monomer (monomer-shape-map oligomer-shape))))
              (lookup-dihedral-cache-impl oligomer-shape other-monomer-shape (external-dihedral-name dihedral-info)
                                          :ignore-degrees ignore-degrees)))
           (dihedral-info-atom
            (let* ((monomer-context (monomer-context monomer-shape-info))
                   (val (gethash monomer-context (monomer-context-to-context-rotamers (rotamers-database oligomer-shape))))
                   (rotamers (rotamer-vector val)))
              (if ignore-degrees
                  (values :ignore-degrees monomer-shape)
                  (etypecase monomer-shape
                    (rotamer-shape
                     (let* ((rotamer-index (rotamer-index monomer-shape))
                            (rotamer (aref rotamers rotamer-index))
                            (dihedral-cache-deg (backbone-dihedral-cache-deg rotamer))
                            (deg (find-in-cache dihedral-cache-deg dihedral-name)))
                       (unless deg
                         (error "Could not find ~s in ~s" dihedral-name dihedral-cache-deg))
                       (values deg monomer-shape)))
                    (backbone-residue-shape
                     (let ((deg (find-in-cache (backbone-dihedral-cache-deg monomer-shape) dihedral-name)))
                       (values deg monomer-shape))))))))))
      (t (error "there was no :dihedrals property in ~s" topology)))))

(defun lookup-dihedral-cache (oligomer-shape monomer-shape dihedral-name &key ignore-degrees)
  (lookup-dihedral-cache-impl oligomer-shape monomer-shape dihedral-name :ignore-degrees ignore-degrees))

#+(or)
(defun build-externals-from-internals (assembler &key oligomer-shape into-coords)
  (let ((coords (make-coordinates-for-assembler assembler)))
    (build-atom-tree-external-coordinates-and-adjust assembler coords oligomer-shape)
    (copy-joint-positions-into-atoms assembler coords oligomer-shape)))

#+(or)
(defun build-oligomer-shape-in-aggregate (assembler oligomer-shape)
  (let ((coords (make-coordinates-for-assembler assembler)))
    (update-internals assembler oligomer-shape)
    (build-atom-tree-external-coordinates-and-adjust assembler coords oligomer-shape)
    (copy-joint-positions-into-atoms assembler coords oligomer-shape)))

#+(or)
(defgeneric build-initial-oligomer-shape-externals (assembler oligomer-shape &key orientation coords)
  (:documentation "Build the oligomer-shape in the coordinates by generating the initial externals"))

#+(or)
(defmethod build-initial-oligomer-shape-externals (assembler oligomer-shape &key (coords (make-coordinates-for-assembler assembler)) orientation)
  (build-atom-tree-external-coordinates-and-adjust assembler coords oligomer-shape orientation)
  coords)

#+(or)
(defun build-all-externals (assembler &key (coords (make-coordinates-for-assembler assembler)))
  (loop for oligomer-shape in (oligomer-shapes assembler)
        do (build-atom-tree-external-coordinates-and-adjust assembler coords oligomer-shape oligomer-shape))
  coords)

(defun build-oligomer-shape-keep-internals-in-coordinates (assembler oligomer-shape coords)
  (build-atom-tree-external-coordinates-and-adjust assembler coords oligomer-shape))

(defmethod aggregate (oligomer-shape)
  "Generate a random oligomer-shape and return the aggregate"
  (let* ((bs (make-permissible-backbone-rotamers oligomer-shape)))
    (write-rotamers oligomer-shape bs (random-rotamers bs))
    (let ((ss (make-permissible-sidechain-rotamers oligomer-shape)))
      (write-rotamers oligomer-shape ss (random-rotamers ss))
      (let* ((ass (make-assembler (list oligomer-shape)))
             (coords (make-coordinates-for-assembler ass))
             )
        (update-internals ass oligomer-shape)
        (update-externals ass :oligomer-shape oligomer-shape
                              :orientation oligomer-shape
                              :coords coords)
        (copy-all-joint-positions-into-atoms ass coords)
        (aggregate ass)))))


#+(or)
(defun assembler-aggregate (assembler oligomer-shapes)
  (let ((coords (make-coordinates-for-assembler assembler)))
    (loop for oligomer-shape in oligomer-shapes
          do (progn
               (update-internals assembler oligomer-shape)
               (build-all-atom-tree-external-coordinates-and-adjust assembler coords)
               (copy-all-joint-positions-into-atoms assembler coords)))
    (aggregate assembler)))

(defmethod aggregate ((oligomer-shape oligomer-shape))
  "Generate an aggregate for the OLIGOMER-SHAPE"
  (let* ((ass (make-assembler (list oligomer-shape)))
         (coords (make-coordinates-for-assembler ass)))
    (update-internals ass oligomer-shape)
    (build-all-atom-tree-external-coordinates-and-adjust ass coords)
    (copy-all-joint-positions-into-atoms ass coords)
    (aggregate ass)))

(defmethod aggregate ((oligomer oligomer))
  "Generate an aggregate for the OLIGOMER"
  (let* ((oligomer-shape (make-oligomer-shape-of-rotamer-shapes oligomer)))
    (aggregate oligomer-shape)))

#+(or)
(defun analyze-oligomer-shape (oligomer-shape)
  "Print an analysis of the oligomer-shape. Print the number of backbone and sidechain monomer-shapes that have defined rotamer-index's "
  (let ((backbone-count 0)
        (backbone-defined 0)
        (backbone-residue-shapes 0)
        (backbone-rotamer-shapes 0)
        (backbone-zero 0)
        (sidechain-count 0)
        (sidechain-defined 0)
        (sidechain-residue-shapes 0)
        (sidechain-rotamer-shapes 0)
        (sidechain-zero 0))
    (loop for monomer-shape across (monomer-shape-vector oligomer-shape)
          for monomer-shape-info across (monomer-shape-info-vector oligomer-shape)
          if (eq (monomer-shape-kind monomer-shape-info) :backbone)
            do (progn
                 (incf backbone-count)
                 (etypecase monomer-shape
                   (rotamer-shape
                    (incf backbone-rotamer-shapes)
                    (when (slot-boundp monomer-shape 'rotamer-index)
                      (incf backbone-defined)
                      (when (/= 0 (rotamer-index monomer-shape))
                        (incf backbone-zero))))
                   (residue-shape
                    (incf backbone-residue-shapes)
                    (incf backbone-defined))
                   ))
          else
            do (progn
                 (incf sidechain-count)
                 (etypecase monomer-shape
                   (rotamer-shape
                    (incf sidechain-rotamer-shapes)
                    (when (slot-boundp monomer-shape 'rotamer-index)
                      (incf sidechain-defined)
                      (when (/= 0 (rotamer-index monomer-shape))
                        (incf sidechain-zero))))
                   (residue-shape
                    (incf sidechain-residue-shapes)
                    (incf sidechain-defined)))))
    (format t "ANALYZE oligomer-shape ~s monomer-shapes ~d~%" (name oligomer-shape) (length (monomer-shape-vector oligomer-shape)))
    (format t "  backbones ~d  rotamer-index defined ~d  = zero ~d~%" backbone-count backbone-defined backbone-zero)
    (format t "       residue-shapes ~d      rotamer-shapes ~d~%" backbone-residue-shapes backbone-rotamer-shapes)
    (format t "  sidechains ~d  rotamer-index defined ~d  = zero ~d~%" sidechain-count sidechain-defined sidechain-zero)
    (format t "       residue-shapes ~d      rotamer-shapes ~d~%" sidechain-residue-shapes sidechain-rotamer-shapes)
    ))

(defun random-oligomer-shape-aggregate (oligomer-shape &key prepare-assembler)
  "Generate a random oligomer-shape and return the aggregate"
  (let* ((bs (make-permissible-backbone-rotamers oligomer-shape)))
    (write-rotamers oligomer-shape bs (random-rotamers bs))
    (let ((ss (make-permissible-sidechain-rotamers oligomer-shape)))
      (write-rotamers oligomer-shape ss (random-rotamers ss))
      (let* ((ass (make-assembler (list oligomer-shape)))
             (coords (make-coordinates-for-assembler ass)))
        (when prepare-assembler
          (funcall prepare-assembler ass))
        (update-internals ass oligomer-shape)
        (update-externals ass :oligomer-shape oligomer-shape
                              :orientation oligomer-shape
                              :coords coords)
        (copy-all-joint-positions-into-atoms ass coords)
        (aggregate ass)))))



(defun build-rotamer-shape-map (oligomer)
  "Build a monomer-to-monomer-shape-map of only rotamer-shapes.
The ROTAMER-INDEX slots of the ROTAMER-SHAPEs will be unbound on return."
  (let ((new-monomer-shape-map (make-hash-table)))
    (loop for monomer across (monomers oligomer)
          for rotamer-shape = (make-rotamer-shape monomer oligomer)
          do (setf (gethash monomer new-monomer-shape-map) rotamer-shape))
    new-monomer-shape-map))

(defun build-mutant-monomer-shape-map (original-oligomer-shape oligomer)
  "Build a monomer-to-monomer-shape-map using an ORIGINAL-OLIGOMER-SHAPE as a template and an OLIGOMER for the rotamer-shapes.
TOPOLOGY:RESIDUE-SHAPEs are used directly from the ORIGINAL-OLIGOMER-SHAPE because they are immutable and fresh TOPOLOGY:ROTAMER-SHAPEs
are generated for each TOPOLOGY:ROTAMER-SHAPE in the ORIGINAL-OLIGOMER-SHAPE.

The ROTAMER-INDEX slots of the ROTAMER-SHAPEs will be unbound on return."
  (unless (eq (topology:oligomer-space (topology:oligomer original-oligomer-shape))
              (topology:oligomer-space oligomer))
    (error "The oligomer spaces of the two arguments MUST match"))
  (let ((original-monomer-to-monomer-shape-map (topology:monomer-shape-map original-oligomer-shape))
        (new-monomer-shape-map (make-hash-table)))
    (maphash (lambda (monomer original-monomer-shape)
               (let ((monomer-shape (etypecase original-monomer-shape
                                      (topology:residue-shape original-monomer-shape)
                                      (topology:rotamer-shape (topology:make-rotamer-shape monomer oligomer)))))
                 (setf (gethash monomer new-monomer-shape-map) monomer-shape)))
             original-monomer-to-monomer-shape-map)
    new-monomer-shape-map))

(defun mutate-oligomer-shape (original-oligomer-shape ligand-oligomer &key name sidechain-rotamer-indexes)
  "Generate a new oligomer-shape based on ORIGINAL-OLIGOMER-SHAPE with the sequence LIGAND-OLIGOMER.
If the LIGAND-OLIGOMER mutates any residues in the backbone then random rotamers are chosen for the new oligomer-shape.
Assign the SIDECHAIN-ROTAMER-INDEXES to the shape of the sidechains."
  (check-type original-oligomer-shape topology:oligomer-shape)
  (let ((new-monomer-shape-map (build-mutant-monomer-shape-map original-oligomer-shape ligand-oligomer)))
    ;; Now we have a new-monomer-shape-map with residue-shape and rotamer-shapes
    ;; but the rotamer-shapes have unbound rotamer-index
    ;; Fix up the backbones by copying rotamer-index if the rotamer-shape-context-rotamers is unchanged
    ;;   or assigning a random one if it has changed
    (let ((callback-backbone-rotamer-indexes
            (lambda (permissible-backbone-rotamers)
              (loop with result-vec = (make-array 16 :adjustable t :fill-pointer 0)
                    for permissible-rotamer across (topology:permissible-rotamer-vector permissible-backbone-rotamers)
                    for monomer-shape-locus = (topology:monomer-shape-locus permissible-rotamer)
                    for allowed-rotamer-indexes = (topology:allowed-rotamer-indexes permissible-rotamer)
                    for monomer-shape = (aref (topology:monomer-shape-vector original-oligomer-shape) monomer-shape-locus)
                    for monomer-shape-info = (aref (topology:monomer-shape-info-vector original-oligomer-shape) monomer-shape-locus)
                    for monomer = (topology:monomer monomer-shape-info)
                    for new-monomer-shape = (gethash monomer new-monomer-shape-map)
                    for new-index = (if (eq (topology:rotamer-shape-context-rotamers new-monomer-shape)
                                            (topology:rotamer-shape-context-rotamers monomer-shape))
                                        (topology:rotamer-index monomer-shape)
                                        (progn
                                          (let ((len (length allowed-rotamer-indexes)))
                                            (aref allowed-rotamer-indexes (random len)))))
                    do (vector-push-extend new-index result-vec)
                    finally (return (make-instance 'topology:rotamer-indexes :rotamer-indexes result-vec)))))
          (callback-sidechain-rotamer-indexes
            (lambda (permissible-sidechain-rotamers)
              sidechain-rotamer-indexes)))
      (topology:make-oligomer-shape ligand-oligomer
                                    :monomer-to-monomer-shape-map new-monomer-shape-map
                                    :name name
                                    :callback-backbone-rotamer-indexes callback-backbone-rotamer-indexes
                                    :callback-sidechain-rotamer-indexes callback-sidechain-rotamer-indexes))))
