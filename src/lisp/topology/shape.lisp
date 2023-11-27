(in-package :topology)

(defgeneric set-to-closest-matching-fragment-conformation-index (monomer-shape oligomer info)
  (:documentation "Set the fragment-conformation-index to the closest matching context-rotamers index"))

(defconstant +phi+   :phi)
(defconstant +phi-1+ :phi-1)
(defconstant +psi+   :psi)
(defconstant +psi-1+ :psi-1)
(defun shape-key-dihedral-name-p (name)
  (member name (list +phi+ +phi-1+ +psi+ +psi-1+)))

(defun bin-dihedral-deg (dih-deg)
  "Take a dihedral in degrees and round it to the nearest 10 degrees"
  (floor (degrees-limit (floor (* 10 (round dih-deg 10.0))))))

(defclass monomer-shape ()
  ((fragment-conformation-index :initform 0 :initarg :fragment-conformation-index :accessor fragment-conformation-index)
   (monomer :initarg :monomer :accessor monomer)
   (monomer-context-key :initarg :monomer-context-key :accessor monomer-context-key)
   (monomer-context :initarg :monomer-context :accessor monomer-context)
   (monomer-shape-kind :initarg :monomer-shape-kind :accessor monomer-shape-kind)
   (shape-key :initarg :shape-key :accessor shape-key)
   (context-rotamers :initarg :context-rotamers :accessor context-rotamers)
   (keys :initarg :keys :accessor keys)
   ))

(defmethod print-object ((obj monomer-shape) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a"
            (monomer obj)
            (if (slot-boundp obj 'fragment-conformation-index)
                (fragment-conformation-index obj)
                "#<unbound>"))))

(defclass oligomer-shape ()
  ((orientation :initarg :orientation :accessor orientation)
   (oligomer :initarg :oligomer :accessor oligomer)
   (rotamers-map :initarg :rotamers-map :accessor rotamers-map)
   (monomer-shape-vector :initarg :monomer-shape-vector :accessor monomer-shape-vector)
   (monomer-shape-map :initarg :monomer-shape-map :accessor monomer-shape-map)
   (the-root-monomer :initarg :the-root-monomer :accessor the-root-monomer)
   (in-monomers :initarg :in-monomers :accessor in-monomers)
   (out-monomers :initarg :out-monomers :accessor out-monomers)
   (rotor3 :initarg :rotor3 :initform #(1.0 0.0 0.0 0.0) :accessor rotor3)
   (translation :initarg :translation :initform (geom:vec 0 0 0) :accessor translation)
   (monomer-shape-build-order :initarg :monomer-shape-build-order
                              :accessor monomer-shape-build-order)
   ))

(defmethod print-object ((obj oligomer-shape) stream)
  (let ((indices (mapcar (lambda (ms) (fragment-conformation-index ms)) (coerce (monomer-shape-vector obj) 'list))))
    (print-unreadable-object (obj stream :type t)
      (let ((*print-length* 5)
            (*print-pretty* nil))
        (format stream "~s" indices)))))

(defmethod kin:orientation-transform ((oligomer-shape oligomer-shape))
  (kin:orientation-transform (orientation oligomer-shape)))

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

(defgeneric make-oligomer-shape (oligomer-or-space rotamers-db &key))

(defmethod make-oligomer-shape ((oligomer oligomer) rotamers-db &key monomer-shape-factory (orientation (make-orientation)))
  (let* ((foldamer (topology:foldamer (topology:oligomer-space oligomer)))
         (shape-info (shape-info foldamer))
         (kind-order (loop for kind-keys in (kind-keys shape-info)
                           collect (kind kind-keys))))
    (multiple-value-bind (monomer-shape-vector the-root-monomer in-monomers out-monomers monomer-shape-map)
        (loop with monomer-shape-vector = (make-array (length (monomers oligomer)))
              with in-monomers = (make-hash-table)
              with out-monomers = (make-hash-table)
              with the-root-monomer = nil
              with foldamer = (foldamer (oligomer-space oligomer))
              with monomer-shape-map = (make-hash-table)
              for index from 0
              for monomer in (ordered-monomers oligomer)
              for monomer-context = (topology:foldamer-monomer-context monomer oligomer foldamer)
              ;; The following is not good - to not define a slot when rotamers-db is NIL
              for context-rotamers = (when rotamers-db 
                                       (topology:lookup-rotamers-for-context rotamers-db monomer-context))
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
              for monomer-context-key = (cons in-monomer-context monomer-context)
              for monomer-shape = (if monomer-shape-factory
                                      (funcall monomer-shape-factory
                                               monomer
                                               :monomer-context monomer-context
                                               :monomer-context-key monomer-context-key
                                               :monomer-shape-kind shape-kind
                                               :context-rotamers context-rotamers)
                                      (make-instance 'monomer-shape
                                                     :monomer monomer
                                                     :monomer-context monomer-context
                                                     :monomer-context-key monomer-context-key
                                                     :monomer-shape-kind shape-kind
                                                     :context-rotamers context-rotamers))
              ;;            do (format t "monomer = ~a~%" monomer)
              do (setf (gethash monomer monomer-shape-map) monomer-shape)
              do (unless in-monomer (setf the-root-monomer monomer))
              do (setf (gethash monomer out-monomers) out-mons)
              do (setf (aref monomer-shape-vector index) monomer-shape)
                 ;;            do (format t "monomer-context ~a~%" monomer-context)
              finally (return (values monomer-shape-vector the-root-monomer in-monomers out-monomers monomer-shape-map)))
      (let* ((os (make-instance 'oligomer-shape
                                :orientation orientation
                                :oligomer oligomer
                                :rotamers-map rotamers-db
                                :monomer-shape-vector monomer-shape-vector
                                :monomer-shape-map monomer-shape-map
                                :the-root-monomer the-root-monomer
                                :in-monomers in-monomers
                                :out-monomers out-monomers))
             (ss (make-permissible-sidechain-rotamers os)))
        (write-rotamers ss (first-rotamers ss))
        os))))

(defmethod make-oligomer-shape ((oligomer-space oligomer-space) rotamers-db &key (oligomer-index 0) (orientation (make-orientation)))
  (let ((oligomer (make-oligomer oligomer-space oligomer-index)))
    (make-oligomer-shape oligomer rotamers-db :orientation orientation)))

(defgeneric read-rotamers (obj)
  (:documentation "Read the rotamers from the object"))

(defgeneric write-rotamers (obj vec)
  (:documentation "Write the rotamers into the object"))

(defmethod read-rotamers ((obj oligomer-shape))
  (let ((vec (make-array (length (monomer-shape-vector obj)) :element-type 'ext:byte32)))
    (map-into vec #'fragment-conformation-index (monomer-shape-vector obj))
    vec))

(defmethod write-rotamers ((obj oligomer-shape) vec)
  (loop for rotamer-index across vec
        for monomer-shape across (monomer-shape-vector obj)
        do (setf (fragment-conformation-index monomer-shape) rotamer-index)))

(defun all-monomers-impl (root shape)
  #+(or)(format t "monomer ~a in: ~a~%" root (gethash root (in-monomers shape)))
  (let ((out-monomers (gethash root (out-monomers shape))))
    (loop for out-monomer in out-monomers
          do (all-monomers-impl out-monomer shape))))

(defun all-monomers (shape)
  (let ((root (the-root-monomer shape)))
    (all-monomers-impl root shape)))

#|
(defun random-fragment-conformation-index-impl (root-monomer-shape oligomer-shape)
  (let ((out-monomers (gethash (monomer root-monomer-shape) (out-monomers oligomer-shape))))
    (loop for out-monomer in out-monomers
          for out-monomer-shape = (gethash out-monomer (monomer-shape-map oligomer-shape))
          for fragment-match-key = (cons (monomer-context root-monomer-shape) (monomer-context out-monomer-shape))
          for allowed-fragment-vec = (let* ((thing (rotamers-map oligomer-shape))
                                            (ht (topology:context-to-rotamers thing))
                                            (val (gethash fragment-match-key ht)))
                                       (unless val (break "Could not find value for key ~a in ht: ~a" fragment-match-key ht))
                                       val)
          for allowed-fragment-indices = (progn
                                           #+(or)(format t "alowed-fragment-vec ~s (fragment-conformation-index root-monomer-shape) -> ~s~%" allowed-fragment-vec (fragment-conformation-index root-monomer-shape))
                                           (elt allowed-fragment-vec (fragment-conformation-index root-monomer-shape)))
          for fragment-conformation-index = (if allowed-fragment-indices
                                                (elt allowed-fragment-indices (random (length allowed-fragment-indices)))
                                                :BADBADBAD)
          do (setf (fragment-conformation-index out-monomer-shape) fragment-conformation-index)
          do (format t "fragment-conformation-index ~a for monomer-shape ~a~%"
                     fragment-conformation-index out-monomer-shape)
          do (unless (< (fragment-conformation-index out-monomer-shape)
                        (length (topology:fragments (context-rotamers out-monomer-shape))))
               (error "fragment-conformation-index ~a is out of bounds ~a"
                      (fragment-conformation-index root-monomer-shape)
                      (length (topology:fragments (context-rotamers out-monomer-shape)))))
          do (random-fragment-conformation-index-impl out-monomer-shape oligomer-shape))))

(defun random-fragment-conformation-index (oligomer-shape)
  (let* ((root (the-root-monomer oligomer-shape))
         (root-monomer-shape (gethash root (monomer-shape-map oligomer-shape))))
    (format t "root-monomer-shape ~a~%" root-monomer-shape)
    (let* ((context-rotamers (context-rotamers root-monomer-shape)))
      (format t "context-rotamers: ~a~%" context-rotamers)
      (setf (fragment-conformation-index root-monomer-shape)
            (random (length (topology:rotamers context-rotamers))))
      (unless (< (fragment-conformation-index root-monomer-shape)
                 (length (topology:rotamers context-rotamers)))
        (error "fragment-conformation-index ~a is out of bounds ~a"
               (fragment-conformation-index root-monomer-shape)
               (length (topology:rotamers context-rotamers))))
      (random-fragment-conformation-index-impl root-monomer-shape oligomer-shape))))
|#


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
       (cons (fragment-conformation-index (car shape-key))
             (fragment-conformation-index (cdr shape-key))))
      (t (fragment-conformation-index shape-key)))))


(defun lookup-dihedral-cache-impl (oligomer-shape monomer-shape dihedral-name)
  (let* ((monomer (monomer monomer-shape))
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
              (lookup-dihedral-cache-impl oligomer-shape other-monomer-shape (external-dihedral-name dihedral-info))))
           (dihedral-info-atom
            (let* ((monomer-context (monomer-context monomer-shape))
                   (rotamers (rotamers (gethash monomer-context (monomer-context-to-context-rotamers (rotamers-map oligomer-shape)))))
                   (rotamer-index (fragment-conformation-index monomer-shape))
                   (rotamer (aref rotamers rotamer-index))
                   (dihedral-cache-deg (backbone-dihedral-cache-deg rotamer))
                   (deg (getf dihedral-cache-deg dihedral-name)))
              (values deg monomer-shape))))))
      (t (error "there was no :dihedrals property in ~s" topology)))))

(defun lookup-dihedral-cache (oligomer-shape monomer-shape dihedral-name)
  (lookup-dihedral-cache-impl oligomer-shape monomer-shape dihedral-name))


(defun build-externals-from-internals (assembler &key oligomer-shape into-coords)
  (let ((coords (make-coordinates-for-assembler assembler)))
    (build-atom-tree-external-coordinates-and-adjust assembler coords oligomer-shape)
    (copy-joint-positions-into-atoms assembler coords oligomer-shape)))

(defun build-oligomer-shape-in-aggregate (assembler oligomer-shape)
  (let ((coords (make-coordinates-for-assembler assembler)))
    (fill-internals-from-oligomer-shape-and-adjust assembler oligomer-shape)
    (build-atom-tree-external-coordinates-and-adjust assembler coords oligomer-shape)
    (copy-joint-positions-into-atoms assembler coords oligomer-shape)))

(defun build-oligomer-shape-in-coordinates (assembler oligomer-shape coords)
    (fill-internals-from-oligomer-shape-and-adjust assembler oligomer-shape)
    (build-atom-tree-external-coordinates-and-adjust assembler coords oligomer-shape))

(defun build-oligomer-shape-keep-internals-in-coordinates (assembler oligomer-shape coords)
  (build-atom-tree-external-coordinates-and-adjust assembler coords oligomer-shape))

(defun build-all-oligomer-shapes-in-coordinates (assembler coords)
  (loop for oligomer-shape in (topology:oligomer-shapes assembler)
        do (fill-internals-from-oligomer-shape-and-adjust assembler oligomer-shape)
        do (build-atom-tree-external-coordinates-and-adjust assembler coords oligomer-shape)))

(defun build-all-oligomer-shapes-from-internals-in-coordinates (assembler coords)
  (loop for oligomer-shape in (topology:oligomer-shapes assembler)
        do (build-atom-tree-external-coordinates-and-adjust assembler coords oligomer-shape)))

(defun random-oligomer-shape-aggregate (oligomer-shape)
  "Generate a random oligomer-shape and return the aggregate"
  (let* ((bs (make-permissible-backbone-rotamers oligomer-shape)))
    (write-rotamers bs (random-rotamers bs))
    (let ((ss (make-permissible-sidechain-rotamers oligomer-shape)))
      (write-rotamers ss (random-rotamers ss))
      (let* ((ass (make-assembler (list oligomer-shape)))
             (coords (make-coordinates-for-assembler ass))
             )
        (fill-internals-from-oligomer-shape-and-adjust ass oligomer-shape)
        (build-all-atom-tree-external-coordinates-and-adjust ass coords)
        (copy-all-joint-positions-into-atoms ass coords)
        (aggregate ass)))))


(defun assembler-aggregate (assembler oligomer-shapes)
  (let ((coords (make-coordinates-for-assembler assembler)))
    (loop for oligomer-shape in oligomer-shapes
          do (progn
               (fill-internals-from-oligomer-shape-and-adjust assembler oligomer-shape)
               (build-all-atom-tree-external-coordinates-and-adjust assembler coords)
               (copy-all-joint-positions-into-atoms assembler coords)))
    (aggregate assembler)))

(defmethod aggregate ((oligomer-shape oligomer-shape))
  "Generate an aggregate for the oligomer-shape"
  (let* ((ass (make-assembler (list oligomer-shape)))
         (coords (make-coordinates-for-assembler ass))
         )
    (fill-internals-from-oligomer-shape-and-adjust ass oligomer-shape)
    (build-all-atom-tree-external-coordinates-and-adjust ass coords)
    (copy-all-joint-positions-into-atoms ass coords)
    (aggregate ass)))
