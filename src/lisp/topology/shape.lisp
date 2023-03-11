(in-package :topology)


(defclass monomer-shape ()
  ((fragment-conformation-index :initform nil :initarg :fragment-conformation-index :accessor fragment-conformation-index)
   (monomer :initarg :monomer :accessor monomer)
   (monomer-context :initarg :monomer-context :accessor monomer-context)
   (fragment-conformations :initarg :fragment-conformations :accessor fragment-conformations)
   ))

(defmethod print-object ((obj monomer-shape) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a" (monomer obj) (fragment-conformation-index obj))))

(defclass oligomer-shape ()
  ((oligomer :initarg :oligomer :accessor oligomer)
   (matched-fragment-conformations-map :initarg :matched-fragment-conformations-map :accessor matched-fragment-conformations-map)
   (foldamer :initarg :foldamer :accessor foldamer)
   (monomer-shape-vector :initarg :monomer-shape-vector :accessor monomer-shape-vector)
   (monomer-shape-map :initarg :monomer-shape-map :accessor monomer-shape-map)
   (the-root-monomer :initarg :the-root-monomer :accessor the-root-monomer)
   (in-monomers :initarg :in-monomers :accessor in-monomers)
   (out-monomers :initarg :out-monomers :accessor out-monomers)
   ))


(defun make-oligomer-shape (oligomer matched-fragment-conformations-map foldamer)
  (multiple-value-bind (monomer-shape-vector the-root-monomer in-monomers out-monomers monomer-shape-map)
      (loop with monomer-shape-vector = (make-array (length (monomers oligomer)))
            with in-monomers = (make-hash-table)
            with out-monomers = (make-hash-table)
            with the-root-monomer = nil
            with monomer-shape-map = (make-hash-table)
            for index from 0
            for monomer across (monomers oligomer)
            for monomer-context = (topology:foldamer-monomer-context monomer oligomer foldamer)
            for fragment-conformations = (let ((fc (gethash monomer-context (topology:monomer-context-to-fragment-conformations matched-fragment-conformations-map))))
                                           (unless fc
                                             (error "Could not find monomer-context ~s" monomer-context))
                                           fc)
            for monomer-shape = (make-instance 'monomer-shape
                                               :monomer monomer
                                               :monomer-context monomer-context
                                               :fragment-conformations fragment-conformations)
            for couplings = (couplings monomer)
            for in-monomer = nil
            for out-mons = nil
;;            do (format t "monomer = ~a~%" monomer)
            do (setf (gethash monomer monomer-shape-map) monomer-shape)
            do (maphash (lambda (key coupling)
                          (if (in-plug-name-p key)
                              (progn
                                (setf in-monomer (topology:source-monomer coupling))
                                (setf (gethash monomer in-monomers) (topology:source-monomer coupling))
                                #+(or)(format t "In plug coupling ~a ~a~%" key coupling))
                              (progn
                                (push (topology:target-monomer coupling) out-mons)
                                #+(or)(format t "Out plug coupling ~a ~a~%" key coupling))))
                        couplings)
            do (unless in-monomer
                 (setf the-root-monomer monomer))
            do (setf (gethash monomer out-monomers) out-mons)
            do (setf (aref monomer-shape-vector index) monomer-shape)
;;            do (format t "monomer-context ~a~%" monomer-context)
            finally (return (values monomer-shape-vector the-root-monomer in-monomers out-monomers monomer-shape-map)))
    (make-instance 'oligomer-shape
                   :oligomer oligomer
                   :matched-fragment-conformations-map matched-fragment-conformations-map
                   :foldamer foldamer
                   :monomer-shape-vector monomer-shape-vector
                   :monomer-shape-map monomer-shape-map
                   :the-root-monomer the-root-monomer
                   :in-monomers in-monomers
                   :out-monomers out-monomers)))


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
          for allowed-fragment-vec = (let* ((ht(topology:fragment-matches (topology:matched-fragment-conformations-map oligomer-shape)))
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
          do (random-fragment-conformation-index-impl out-monomer-shape oligomer-shape))))

(defun random-fragment-conformation-index (oligomer-shape)
  (let* ((root (the-root-monomer oligomer-shape))
         (root-monomer-shape (gethash root (monomer-shape-map oligomer-shape))))
    (format t "root-monomer-shape ~a~%" root-monomer-shape)
    (let* ((fragment-conformations (fragment-conformations root-monomer-shape)))
      (format t "fragment-conformations: ~a~%" fragment-conformations)
      (setf (fragment-conformation-index root-monomer-shape)
            (random (length (topology:fragments fragment-conformations))))
      (random-fragment-conformation-index-impl root-monomer-shape oligomer-shape))))


(defun build-shape (oligomer-shape)
  (let* ((oligomer (oligomer oligomer-shape))
         (conf (topology:make-conformation oligomer))
         (fragment-conformations (matched-fragment-conformations-map oligomer-shape)))
    (topology::fill-internals-from-oligomer-shape conf fragment-conformations oligomer-shape)
    #+(or)(topology:zero-all-atom-tree-external-coordinates conf)
    (let ((coordinates (chem:make-coordinates (topology:energy-function conf))))
      (topology:build-all-atom-tree-external-coordinates conf coordinates)
      (warn "Handle copying vector into atoms")
      (chem:energy-function/save-coordinates-from-vector (topology:energy-function conf) coordinates)
      (topology:aggregate conf))))

