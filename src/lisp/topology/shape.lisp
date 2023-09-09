(in-package :topology)

(defgeneric set-to-closest-matching-fragment-conformation-index (monomer-shape oligomer)
  (:documentation "Set the fragment-conformation-index to the closest matching context-rotamers index"))

(defconstant +phi+   :phi)
(defconstant +phi-1+ :phi-1)
(defconstant +psi+   :psi)
(defconstant +psi-1+ :psi-1)
(defun shape-key-dihedral-name-p (name)
  (member name (list +phi+ +phi-1+ +psi+ +psi-1+)))

(defun bin-dihedral (dih)
  "Take a dihedral in degrees and round it to the nearest 10 degrees"
  (floor (* 10 (round dih 10.0))))

(defclass monomer-shape ()
  ((fragment-conformation-index :initarg :fragment-conformation-index :accessor fragment-conformation-index)
   (monomer :initarg :monomer :accessor monomer)
   (monomer-context-key :initarg :monomer-context-key :accessor monomer-context-key)
   (monomer-context :initarg :monomer-context :accessor monomer-context)
   (monomer-shape-kind :initarg :monomer-shape-kind :accessor monomer-shape-kind)
   (shape-key :initarg :shape-key :accessor shape-key)
   (rotamers :initarg :rotamers :accessor rotamers)
   (allowed-rotamer-indices :initarg :allowed-rotamer-indices :accessor allowed-rotamer-indices)
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
  ((oligomer :initarg :oligomer :accessor oligomer)
   (connected-rotamers-map :initarg :connected-rotamers-map :accessor connected-rotamers-map)
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
   (keys :initarg :keys :accessor keys)))

(cando.serialize:make-class-save-load kind-key)

(defclass shape-info ()
  ((kind-keys :initarg :kind-keys :accessor kind-keys)))

(cando.serialize:make-class-save-load shape-info)

(defun make-shape-info (info)
  (let ((kind-keys (mapcar (lambda (kind-key-list)
                                   (make-instance 'kind-key
                                                  :kind (car kind-key-list)
                                                  :keys (cadr kind-key-list)))
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

(defun make-oligomer-shape (oligomer connected-rotamers-map &key monomer-shape-factory)
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
              for context-rotamers = (let ((fc (gethash monomer-context (topology:monomer-context-to-context-rotamers connected-rotamers-map))))
                                             (unless fc
                                               (error "Could not find monomer-context ~s" monomer-context))
                                             fc)
              for shape-kind = (topology:shape-kind foldamer monomer oligomer)
              for shape-key = (build-shape-key shape-info monomer-shape-map shape-kind monomer oligomer in-monomers out-monomers)
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
                                            (unless (in-plug-name-p key)
                                              (push (topology:target-monomer coupling) out-monomers)
                                              #+(or)(format t "Out plug coupling ~a ~a~%" key coupling)))
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
                                               :shape-key shape-key
                                               :context-rotamers context-rotamers)
                                      (make-instance 'monomer-shape
                                                     :monomer monomer
                                                     :monomer-context monomer-context
                                                     :monomer-context-key monomer-context-key
                                                     :shape-key shape-key 
                                                     :context-rotamers context-rotamers))
              ;;            do (format t "monomer = ~a~%" monomer)
              do (setf (gethash monomer monomer-shape-map) monomer-shape)
              do (unless in-monomer (setf the-root-monomer monomer))
              do (setf (gethash monomer out-monomers) out-mons)
              do (setf (aref monomer-shape-vector index) monomer-shape)
                 ;;            do (format t "monomer-context ~a~%" monomer-context)
              finally (return (values monomer-shape-vector the-root-monomer in-monomers out-monomers monomer-shape-map)))
      (let ((monomer-shape-build-order (sorted-build-order foldamer oligomer kind-order)))
        (make-instance 'oligomer-shape
                       :oligomer oligomer
                       :connected-rotamers-map connected-rotamers-map
                       :monomer-shape-vector monomer-shape-vector
                       :monomer-shape-map monomer-shape-map
                       :monomer-shape-build-order monomer-shape-build-order
                       :the-root-monomer the-root-monomer
                       :in-monomers in-monomers
                       :out-monomers out-monomers)))))


(defun all-monomers-impl (root shape)
  #+(or)(format t "monomer ~a in: ~a~%" root (gethash root (in-monomers shape)))
  (let ((out-monomers (gethash root (out-monomers shape))))
    (loop for out-monomer in out-monomers
          do (all-monomers-impl out-monomer shape))))

(defun all-monomers (shape)
  (let ((root (the-root-monomer shape)))
    (all-monomers-impl root shape)))




(defun random-fragment-conformation-index-impl (root-monomer-shape oligomer-shape)
  (let ((out-monomers (gethash (monomer root-monomer-shape) (out-monomers oligomer-shape))))
    (loop for out-monomer in out-monomers
          for out-monomer-shape = (gethash out-monomer (monomer-shape-map oligomer-shape))
          for fragment-match-key = (cons (monomer-context root-monomer-shape) (monomer-context out-monomer-shape))
          for allowed-fragment-vec = (let* ((ht (topology:rotamer-context-connections (topology:connected-rotamers-map oligomer-shape)))
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
            (random (length (topology:fragments context-rotamers))))
      (unless (< (fragment-conformation-index root-monomer-shape)
                 (length (topology:fragments context-rotamers)))
        (error "fragment-conformation-index ~a is out of bounds ~a"
               (fragment-conformation-index root-monomer-shape)
               (length (topology:fragments context-rotamers))))
      (random-fragment-conformation-index-impl root-monomer-shape oligomer-shape))))

(defun build-shapes (oligomer-shapes conf &key monomer-order)
  (let ((coordinates (chem:make-coordinates (topology:energy-function conf))))
    (chem:energy-function/load-coordinates-into-vector (topology:energy-function conf) coordinates)
    (loop for oligomer-shape in oligomer-shapes
          for oligomer = (oligomer oligomer-shape)
          do (topology:fill-internals-from-oligomer-shape conf oligomer-shape)
             (if monomer-order
                 (loop for monomer in monomer-order
                       for monomer-position = (gethash monomer (monomer-positions conf))
                       for molecule-index = (molecule-index monomer-position)
                       for residue-index = (residue-index monomer-position)
                       for atmolecule = (aref (ataggregate conf) molecule-index)
                       for atresidue = (aref (atresidues atmolecule) residue-index)
                       do (build-atresidue-atom-tree-external-coordinates atresidue coordinates))
                 (topology:build-all-atom-tree-external-coordinates conf oligomer coordinates)))
    (chem:energy-function/save-coordinates-from-vector (topology:energy-function conf) coordinates)
    (topology:aggregate conf)))


(defun extract-shape-key (monomer-shape)
  (let ((shape-key (shape-key monomer-shape)))
    (cond
      ((null shape-key)
       nil)
      ((eq shape-key :phi/psi)
       (cons (fragment-conformation-index (car shape-key))
             (fragment-conformation-index (cdr shape-key))))
      (t (fragment-conformation-index shape-key)))))

(defun lookup-allowed-rotamers (monomer-shape oligomer-shape)
  (let* ((connected-rotamers-map (connected-rotamers-map oligomer-shape))
         (monomer-context-key (monomer-context-key monomer-shape))
         (rotamer-shape-connections (lookup-rotamer-context-connections (rotamer-context-connections connected-rotamers-map) monomer-context-key))
         (shape-key (extract-shape-key monomer-shape))
         (allowed-rotamers (lookup-rotamer-shape-connections rotamer-shape-connections shape-key)))
    allowed-rotamers))

(defun maybe-set-shape-key (monomer-shape oligomer-shape connected-rotamers-map rotamer-index)
  "If the monomer-shape is a shape-kind :backbone then look for the :phi,:psi,:phi-1,:psi-1 dihedrals and store them
in the shape-key slot of the monomer-shape corresponding to the sidechain"
  (cond
    ((eq :backbone (shape-kind monomer-shape))
     (error "Do what you need for a backbone"))
    ((eq :sidechain (shape-kind monomer-shape))
     (error "Do what you need for a side-chain"))
    (t (error "Unknown shape-kind ~s" (shape-kind monomer-shape)))))

(defun assign-backbone-rotamer-indices (oligomer-shape
                                        &key
                                          (monomer-shape-build-order (monomer-context-to-context-rotamers rotamers-database)))
  "For each backbone monomer-shape set the rotamer-indices."
    (loop with rotamers-database = (connected-rotamers-map oligomer-shape)
          for monomer-shape in monomer-shape-build-order
          for monomer-context = (monomer-context monomer-shape)
          for shape-kind = (shape-kind monomer-shape)
          if (eq (shape-kind monomer-shape) :backbone)
            do (let ((rotamer-indices (gethash monomer-context (monomer-context-to-context-rotamers rotamers-database))))
                 (setf (allowed-rotamer-indices monomer-shape) rotamer-indices))
          else ; once we hit a sidechain we stop
          do (return nil)))

(defun assign-backbone-random-rotamer-indices (oligomer-shape
                                               &key (monomer-shape-build-order (monomer-shape-build-order oligomer-shape)))
  (loop with rotamers-database = (connected-rotamers-map oligomer-shape)
        for monomer-shape in monomer-shape-build-order
        if (eq (shape-kind monomer-shape) :backbone)
          do (let ((rnd (random (length (allowed-rotamers monomer-shape)))))
               (setf (fragment-conformation-index monomer-shape) (aref (allowed-rotamers monomer-shape) rnd)))
        else ; stop once we hit a sidechain
          do (return nil)))

(defun assign-backbone-shape-keys (oligomer-shape &key (monomer-shape-build-order (monomer-shape-build-order oligomer-shape)))
  (loop with rotamers-database = (connected-rotamers-map oligomer-shape)
        for monomer-shape in monomer-shape-build-order
        if (eq (shape-kind monomer-shape) :backbone)
          do (let* ((topology (monomer-topology (monomer monomer-shape) (oligomer oligomer-shape)))
                    (constitution (topology:constitution topology))
                    (residue-properties (topology:residue-properties constitution))
                    (_ (error "At this point get the shape-key values out of the rotamer fragments and put it into the monomer-shape"))
                    )
               (format t "monomer: ~a   topology: ~a~%" monomer topology))
        else
          do (return nil)))
  
