(in-package :topology)

#+(or)
(eval-when (:compile-toplevel :load-toplevel :execute)
  (pushnew :debug-mover *features*))

(defgeneric lookup-backbone-shape-key (oligomer-shape monomer-shape errorp errorval))

(defmethod lookup-backbone-shape-key (oligomer-shape monomer-shape errorp errorval)
  (declare (ignore oligomer-shape monomer-shape errorp errorval))
  (error "Specialize this on the arguments"))

(defun make-phi-psi (phi-deg psi-deg)
  (cons (topology:bin-dihedral-deg phi-deg)
        (topology:bin-dihedral-deg psi-deg)))

(defun lookup-phi-psi-monomer-shape (oligomer-shape monomer-shape)
  "Lookup the phi,psi values in the dihedral cache and return them as a (phi . psi) pair"
  (let ((phi (lookup-dihedral-cache oligomer-shape monomer-shape :phi))
        (psi (lookup-dihedral-cache oligomer-shape monomer-shape :psi)))
    (cons phi psi)))

(define-condition no-rotamers (error)
  ((message :initarg :message :reader message))
  (:report (lambda (condition stream)
             (format stream (message condition)))))

(defun calculate-allowed-rotamers (monomer-shape-index sidechain-monomer-shape sidechain-monomer-shape-info oligomer-shape shape-key-to-allowed-rotamers)
  (declare (ignorable monomer-shape-index))
    #+debug-mover
    (let ((*print-pretty* nil))
      (format t "calculate-allowed-rotamers monomer-shape-index ~a~%" monomer-shape-index))
    (let* ((phi-deg (topology:lookup-dihedral-cache oligomer-shape sidechain-monomer-shape :phi))
           (psi-deg (topology:lookup-dihedral-cache oligomer-shape sidechain-monomer-shape :psi))
           (shape-key nil)
           (result (if (and phi-deg psi-deg)
                       (progn
                         (setf shape-key (make-phi-psi phi-deg psi-deg))
                         (gethash shape-key shape-key-to-allowed-rotamers))
                       nil)))
      (unless result
        (let* ((*print-pretty* nil)
               (monomer (monomer sidechain-monomer-shape-info)))
          (error 'no-rotamers :message (format nil "Found no allowed-rotamers for ~s - phi-deg: ~s  psi-deg: ~s - shape-key ~s -  ~s in ~s - this shouldn't be possible if a backbone-stepper was created and applied to the oligomer-shape" monomer phi-deg psi-deg shape-key sidechain-monomer-shape oligomer-shape))))
      #+debug-mover
      (let ((*print-pretty* nil))
        (format t "   calculate-allowed-rotamers  phi-deg ~d  psi-deg ~d    allowed-rotamers -> ~s~%" phi-deg psi-deg result))
      result
      ))

(defclass permissible-rotamer ()
  ((monomer-shape-locus :initarg :monomer-shape-locus :reader monomer-shape-locus)
   (allowed-rotamer-indexes :initarg :allowed-rotamer-indexes :reader allowed-rotamer-indexes))
  )

(defmethod print-object ((obj permissible-rotamer) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "@~a/~s" (monomer-shape-locus obj) (allowed-rotamer-indexes obj))))

(defclass permissible-rotamers ()
  ((permissible-rotamer-vector :initarg :permissible-rotamer-vector :reader permissible-rotamer-vector)))

(defmethod print-object ((obj permissible-rotamers) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (mapcar (lambda (permissible-rotamer) (format nil "#~s/~s"
                                                                      (monomer-shape-locus permissible-rotamer)
                                                                      (allowed-rotamer-indexes permissible-rotamer)))
                                (coerce (permissible-rotamer-vector obj) 'list)))))

(defclass permissible-backbone-rotamers (permissible-rotamers)
  ())

(defclass permissible-sidechain-rotamers (permissible-rotamers)
  ((required-backbone-rotamer-indexes :initarg :required-backbone-rotamer-indexes
                                      :accessor required-backbone-rotamer-indexes
                                   :documentation "Store a vector of CONS cells (locus . backbone-rotamer-index). This will be used to ensure that the backbone rotamer-indexes are consistent with this permissible-sidechain-rotamers object")))

(defgeneric ensure-oligomer-shape-is-consistent-with-permissible-rotamers (oligomer-shape permissible-rotamers))

(defmethod ensure-oligomer-shape-is-consistent-with-permissible-rotamers (oligomer-shape (permissible-rotamers permissible-rotamers))
  (loop for permissible-rotamer across (permissible-rotamer-vector permissible-rotamers)
        for locus = (monomer-shape-locus permissible-rotamer)
        for monomer-shape = (aref (monomer-shape-vector oligomer-shape) locus)
        do (when (slot-boundp monomer-shape 'rotamer-index)
             (let ((read-rotamer-index (rotamer-index monomer-shape)))
               (unless (position read-rotamer-index (allowed-rotamer-indexes permissible-rotamer))
                 (error "The oligomer-shape is out of sync with the permissible-rotamers~% The read-rotamer-index at index ~a was rotamer-index ~a and must be one of ~s"
                        locus read-rotamer-index (allowed-rotamer-indexes permissible-rotamer)))))))

(defmethod ensure-oligomer-shape-is-consistent-with-permissible-rotamers (oligomer-shape (permissible-rotamers permissible-sidechain-rotamers))
  (loop for pair across (required-backbone-rotamer-indexes permissible-rotamers)
        for locus = (car pair)
        for required-rotamer-index = (cdr pair)
        for monomer-shape = (aref (monomer-shape-vector oligomer-shape) locus)
        for read-rotamer-index = (rotamer-index monomer-shape)
        unless (= required-rotamer-index read-rotamer-index)
          do (error "The required-rotamer-index ~a does not match the read-rotamer-index ~a when we are trying to ensure that the backbone rotamer-indexes are consistent with the permissible-sidechain-rotamers, this will happen if the oligomer-shape has had its backbone rotamer-indexes altered and a stale permissible-sidechain-rotamers is used to write sidechain rotamer-indexes.  Make sure the permissible-sidechain-rotamers is updated to reflect the backbone rotamer-indexes" required-rotamer-index read-rotamer-index))
  (call-next-method))



(defun one-backbone-permissible-rotamer (monomer-context rotamers-db rotamers index sidechain sidechain-monomer-context)
  (declare (ignore sidechain))
  (let ((allowed-rotamer-indexes (make-array 16 :adjustable t :fill-pointer 0)))
    (loop :with backbone-rotamers = rotamers
          :with search-key-history = nil
          :with sidechain-rotamers = (gethash sidechain-monomer-context (topology:context-to-rotamers rotamers-db))
          :with shape-key-to-index = (topology:shape-key-to-index sidechain-rotamers)
          :for ii :below (length backbone-rotamers)
          :for backbone-rotamer = (aref backbone-rotamers ii)
          :for backbone-dihedral-cache-deg = (topology:backbone-dihedral-cache-deg backbone-rotamer)
          :for phi-1 = (getf backbone-dihedral-cache-deg :phi-1)
          :for psi-1 = (getf backbone-dihedral-cache-deg :psi-1)
          :for binned-phi-1 = (topology:bin-dihedral-deg phi-1)
          :for binned-psi-1 = (topology:bin-dihedral-deg psi-1)
          ;; Get the shape-key for the backbone
          :for shape-key = (make-phi-psi binned-phi-1 binned-psi-1)
          :for has-allowed-sidechain-rotamer-indices = (gethash shape-key shape-key-to-index)
          :do (push shape-key search-key-history)
          ;; If the sidechain has a shape-key that matches this backbone then
          ;;    add the backbone rotamer index to monomer-allowed-rotamer-indexes
          :do (progn
;;;                #+debug-mover
                (let ((*print-pretty* nil))
                  (format t "For backbone-index ~a rotamer-index-index ~a   shape-key ~a matches-sidechain -> ~a~%" index ii shape-key has-allowed-sidechain-rotamer-indices)))
          :when has-allowed-sidechain-rotamer-indices
            :do (let ((seen (position ii allowed-rotamer-indexes)))
                  (unless seen (vector-push-extend ii allowed-rotamer-indexes))
                  )
          :finally (when (= 0 (length allowed-rotamer-indexes))
                     (error "There are no allowed-rotamers for monomer-shape at ~a~%The backbone monomer-context is ~s~%The sidechain monomer-context is ~s~%search-key-history ~s~%keys of shape-key-to-index: ~s" index monomer-context sidechain-monomer-context
                            search-key-history
                            (alexandria:hash-table-keys shape-key-to-index)))
          )
    (let ((permissible-rotamer (make-instance 'permissible-rotamer
                                              :allowed-rotamer-indexes allowed-rotamer-indexes
                                              :monomer-shape-locus index)))
      permissible-rotamer)))

(defun make-permissible-backbone-rotamers (oligomer-shape)
  "All backbone rotamers are allowed since we enforced that every backbone rotamer has side-chain rotamers"
  (let ((rotamers-db (topology:rotamers-database oligomer-shape))
        (permissible-rotamer-vector (make-array 16 :adjustable t :fill-pointer 0)))
    (with-slots (topology:monomer-shape-vector
                 topology:monomer-shape-info-vector
                 topology:oligomer
                 topology:rotamers-database
                 topology:monomer-shape-map) oligomer-shape
      #+debug-mover(format t "-------- Next stage~%")
      (loop :for index below (length topology:monomer-shape-vector)
            :for monomer-shape = (aref topology:monomer-shape-vector index)
            :for monomer-shape-info = (aref topology:monomer-shape-info-vector index)
            :for monomer-shape-kind = (topology:monomer-shape-kind monomer-shape-info)
            :for monomer-context = (topology:monomer-context monomer-shape-info)
            :for orotamers = (gethash monomer-context (topology:context-to-rotamers rotamers-db))
            :for rotamers = (rotamer-vector orotamers)
            :do (cond
                  ((eq :backbone monomer-shape-kind)
                   ;; We have a backbone that doesn't influence any sidechain
                   ;; - so all backbone rotamers are accessible
                   (let ((allowed-rotamer-indexes (make-array 16 :adjustable t :fill-pointer 0)))
                     (loop :with backbone-rotamers = rotamers
                           :for ii :below (length backbone-rotamers)
                           :do (vector-push-extend ii allowed-rotamer-indexes))
                     (let ((permissible-rotamer (make-instance 'permissible-rotamer
                                                               :allowed-rotamer-indexes allowed-rotamer-indexes
                                                               :monomer-shape-locus index)))
                       (vector-push-extend permissible-rotamer permissible-rotamer-vector))))
                  ((eq :sidechain monomer-shape-kind)
                   ;; Do nothing
                   )
                  (t (error "How did we get here?")))))
    (make-instance 'permissible-backbone-rotamers
                   :permissible-rotamer-vector (copy-seq permissible-rotamer-vector))))


#+(or) (defun make-permissible-backbone-rotamers (oligomer-shape)
  "Limit the backbone rotamer indices to the ones that have sidechain rotamers with matching shape-keys"
  (let ((rotamers-db (topology:rotamers-database oligomer-shape))
        (permissible-rotamer-vector (make-array 16 :adjustable t :fill-pointer 0)))
    (with-slots (topology:monomer-shape-vector
                 topology:monomer-shape-info-vector
                 topology:oligomer
                 topology:rotamers-database
                 topology:monomer-shape-map) oligomer-shape
      (let ((special-sidechain-monomer-shapes (make-hash-table)))
        (loop :for index below (length topology:monomer-shape-vector)
                    :for monomer-shape = (aref topology:monomer-shape-vector index)
                    :for monomer-shape-info = (aref topology:monomer-shape-info-vector index)
                    :for monomer = (topology:monomer monomer-shape-info)
                    :for monomer-shape-kind = (topology:monomer-shape-kind monomer-shape-info)
                    :for monomer-context = (topology:monomer-context monomer-shape-info)
                    :do (progn
                          #+debug-mover
                          (let* ((out-monomers (gethash monomer (topology:out-monomers oligomer-shape))))
                            (format t "monomer-shape-index ~a monomer-shape-kind: ~s ~s~%" index monomer-shape-kind monomer-context)
                            (loop :for out-monomer in out-monomers
                                  :for out-monomer-shape = (gethash out-monomer topology:monomer-shape-map )
                                  :for out-monomer-shape-index = (position out-monomer-shape topology:monomer-shape-vector)
                                  :do (format t "   out: ~s ~s~%" out-monomer-shape-index out-monomer))
                            )
                          )
                    :when (eq monomer-shape-kind :sidechain)
                      :do (multiple-value-bind (deg backbone1-monomer-shape)
                              (topology:lookup-dihedral-cache oligomer-shape monomer-shape :psi :ignore-degrees t)
                            #+debug-mover
                            (let ((backbone1-monomer-shape-index (position backbone1-monomer-shape topology:monomer-shape-vector)))
                              (format t "     sidechain-monomer-shape-index ~a   backbone1-monomer-shape-index ~a~%" index backbone1-monomer-shape-index))
                            (setf (gethash backbone1-monomer-shape special-sidechain-monomer-shapes) (cons monomer-shape monomer-context))))
        #+debug-mover(format t "-------- Next stage~%")
        (loop :for index below (length topology:monomer-shape-vector)
              :for monomer-shape = (aref topology:monomer-shape-vector index)
              :for monomer-shape-info = (aref topology:monomer-shape-info-vector index)
              :for monomer-shape-kind = (topology:monomer-shape-kind monomer-shape-info)
              :for monomer-context = (topology:monomer-context monomer-shape-info)
              :for orotamers = (gethash monomer-context (topology:context-to-rotamers rotamers-db))
              :for rotamers = (topology:rotamers orotamers)
              :for maybe-sidechain-info = (gethash monomer-shape special-sidechain-monomer-shapes)
              :for maybe-sidechain = (car maybe-sidechain-info)
              :for maybe-sidechain-monomer-context = (cdr maybe-sidechain-info)
              #+debug-mover :do #+debug-mover(format t " monomer-context ~s ~s~%   ~s ~s~%"
                                                     monomer-context
                                                     monomer-shape-kind
                                                     maybe-sidechain
                                                     maybe-sidechain-monomer-context
                                                     )
              :do (cond
                    ((and maybe-sidechain
                          (eq :backbone monomer-shape-kind))
                     ;; We have a backbone that influences a sidechain
                     ;; identify what backbone rotamers are allowed given the sidechain that it effects
                     (let ((permissible-rotamer (one-backbone-permissible-rotamer
                                                 monomer-context
                                                 rotamers-db
                                                 rotamers
                                                 index
                                                 maybe-sidechain
                                                 maybe-sidechain-monomer-context)))
                       (vector-push-extend permissible-rotamer permissible-rotamer-vector)))
                    ((eq :backbone monomer-shape-kind)
                     ;; We have a backbone that doesn't influence any sidechain
                     ;; - so all backbone rotamers are accessible
                     (let ((allowed-rotamer-indexes (make-array 16 :adjustable t :fill-pointer 0)))
                       (loop :with backbone-rotamers = rotamers
                             :for ii :below (length backbone-rotamers)
                             :do (vector-push-extend ii allowed-rotamer-indexes))
                       (let ((permissible-rotamer (make-instance 'permissible-rotamer
                                                                 :allowed-rotamer-indexes allowed-rotamer-indexes
                                                                 :monomer-shape-locus index)))
                         (vector-push-extend permissible-rotamer permissible-rotamer-vector))))
                    ((eq :sidechain monomer-shape-kind)
                     ;; Do nothing
                     )
                    (t (error "How did we get here?")))))
      (make-instance 'permissible-backbone-rotamers
                     :permissible-rotamer-vector (copy-seq permissible-rotamer-vector)))))

(defun make-permissible-rotamer-for-monomer* (oligomer-shape index monomer monomer-context monomer-shape monomer-shape-info rotamer-index rotamers-db
                                              )
  (declare (ignore monomer))
  (let* ((rotamers (gethash monomer-context (topology:context-to-rotamers rotamers-db)))
         (allowed-rotamers (progn
                             (calculate-allowed-rotamers index monomer-shape monomer-shape-info oligomer-shape (topology:shape-key-to-index rotamers)))))
    (when (= 0 (length allowed-rotamers))
      (error "There are no allowed-rotamers for monomer-shape at ~a" index))
    (when (or (and rotamer-index (null (position rotamer-index allowed-rotamers)))
              (null rotamer-index))
      ;; The current sidechain rotamer-index is not in the new permissible rotamers
      ;; make it one of them
      (let* ((rnd (random (length allowed-rotamers)))
             (new-rotamer-index (aref allowed-rotamers rnd)))
        (setf (topology:rotamer-index monomer-shape) new-rotamer-index)))
    (unless allowed-rotamers
      (error "Could not find any allowed-rotamers for ~s" monomer-shape))
    (let ((permissible-rotamer (make-instance 'permissible-rotamer
                                              :monomer-shape-locus index
                                              :allowed-rotamer-indexes allowed-rotamers)))
      permissible-rotamer)))

(defun make-permissible-sidechain-rotamer-for-monomer (oligomer-shape desired-monomer)
  (let ((rotamers-db (topology:rotamers-database oligomer-shape)))
    ;; Loop over all monomer-shape in the oligomer-shape and downselect the ones that
    ;; are sidechains and if monomers is not nil then only the ones in monomers
    (loop :named main
          :for index :from 0
          :for monomer-shape :across (topology:monomer-shape-vector oligomer-shape)
          :for rotamer-index = (when (slot-boundp monomer-shape 'topology:rotamer-index)
                                 (topology:rotamer-index monomer-shape))
          :for monomer-shape-info :across (topology:monomer-shape-info-vector oligomer-shape)
          :for monomer := (topology:monomer monomer-shape-info)
          :for monomer-shape-kind := (topology:monomer-shape-kind monomer-shape-info)
          :for monomer-context := (topology:monomer-context monomer-shape-info)
          :do (cond
                ((and (eq monomer-shape-kind :sidechain)
                      (eq monomer desired-monomer))
                 (let ((permissible-rotamer
                         (make-permissible-rotamer-for-monomer*
                          oligomer-shape
                          index
                          monomer
                          monomer-context
                          monomer-shape
                          monomer-shape-info
                          rotamer-index
                          rotamers-db)))
                   (return-from make-permissible-sidechain-rotamer-for-monomer permissible-rotamer)))))
    (error "Could not find monomer ~s in ~s" desired-monomer oligomer-shape)))

(defun make-permissible-sidechain-rotamers (oligomer-shape &key monomers)
  "Constructs a permissible-sidechain-rotamers object for the given oligomer-shape.
Optionally filters the process based on a provided list of monomers.
- oligomer-shape: The oligomer-shape to be processed.
- monomers: An optional list of monomers to restrict the processing, nil means all sidechain monomers.
Returns an instance of permissible-sidechain-rotamers."
  (let ((rotamers-db (topology:rotamers-database oligomer-shape))
        (required-backbone-rotamer-indexes (make-array 16 :adjustable t :fill-pointer 0))
        (permissible-rotamer-vector (make-array 16 :adjustable t :fill-pointer 0)))
    ;; Loop over all monomer-shape in the oligomer-shape and downselect the ones that
    ;; are sidechains and if monomers is not nil then only the ones in monomers
    (loop :for index :from 0
          :for monomer-shape :across (topology:monomer-shape-vector oligomer-shape)
          :for rotamer-index = (when (slot-boundp monomer-shape 'topology:rotamer-index)
                                 (topology:rotamer-index monomer-shape))
          :for monomer-shape-info :across (topology:monomer-shape-info-vector oligomer-shape)
          :for monomer := (topology:monomer monomer-shape-info)
          :for monomer-shape-kind := (topology:monomer-shape-kind monomer-shape-info)
          :for monomer-context := (topology:monomer-context monomer-shape-info)
          :do (cond
                ((and (eq monomer-shape-kind :sidechain)
                      (or (null monomers)
                          (member monomer monomers))) ; if monomers defined then use it as a subset
                 (let ((permissible-rotamer
                         (make-permissible-rotamer-for-monomer*
                          oligomer-shape
                          index
                          monomer
                          monomer-context
                          monomer-shape
                          monomer-shape-info
                          rotamer-index
                          rotamers-db)))
                     (vector-push-extend permissible-rotamer permissible-rotamer-vector)))
                ((eq monomer-shape-kind :backbone)
                 (vector-push-extend (cons index (rotamer-index monomer-shape))
                                     required-backbone-rotamer-indexes))
                ))
    (make-instance 'permissible-sidechain-rotamers
                   :permissible-rotamer-vector permissible-rotamer-vector
                   :required-backbone-rotamer-indexes (copy-seq required-backbone-rotamer-indexes))))


(defun monomer-set-from-permissible-rotamers (&rest args)
  (loop for cur = args then (cddr cur)
        for oligomer-shape = (car cur)
        for permissible = (cadr cur)
        until (null cur)
        do (break "Check stuff")))

(defun max-rotamers (permissible-rotamers)
  "Read the maximum rotamer indexes from the stepper and return them in a vector"
  (let* ((len (length (permissible-rotamer-vector permissible-rotamers)))
         (vec (make-array len :element-type 'ext:byte32)))
    (loop for index below len
          for permissible-rotamer = (aref (permissible-rotamer-vector permissible-rotamers) index)
          for max = (length (allowed-rotamer-indexes permissible-rotamer))
          do (setf (aref vec index) max))
    vec))

(defmethod read-rotamers (oligomer-shape (permissible-rotamers permissible-rotamers))
  "Read the rotamer shapes from the stepper and return them in a vector"
  (ensure-oligomer-shape-is-consistent-with-permissible-rotamers oligomer-shape permissible-rotamers)
  (let* ((len (length (permissible-rotamer-vector permissible-rotamers)))
         (vec (make-array len :element-type 'ext:byte32)))
    (loop for index below len
          for permissible-rotamer = (aref (permissible-rotamer-vector permissible-rotamers) index)
          for locus = (monomer-shape-locus permissible-rotamer)
          for monomer-shape = (aref (topology:monomer-shape-vector oligomer-shape) locus)
          for read-rotamer-index = (topology:rotamer-index monomer-shape)
          do (unless (position read-rotamer-index (allowed-rotamer-indexes permissible-rotamer))
               (let ((*print-pretty* nil))
                 (error "In read-rotamers a rotamer-index ~a at index ~a was read that is not one of the expected ones~% ~s - this means that the oligomer-shape got out of sync with the permissible-rotamers" read-rotamer-index index (allowed-rotamer-indexes permissible-rotamer))))
          do (setf (aref vec index) (topology:rotamer-index monomer-shape)))
    vec))

(defmethod write-rotamers (oligomer-shape (permissible-rotamers permissible-rotamers) (vec array))
  "Write a vector of rotamer shapes into the stepper"
  (ensure-oligomer-shape-is-consistent-with-permissible-rotamers oligomer-shape permissible-rotamers)
  (unless (= (length vec)
             (length (permissible-rotamer-vector permissible-rotamers)))
    (error "The number of (permissible-rotamer-vector permissible-rotamers) ~s doesn't match the length of the vec ~s"
           (permissible-rotamer-vector permissible-rotamers)
           vec))
  (let ((len (length (permissible-rotamer-vector permissible-rotamers))))
    (loop for index below len
          for permissible-rotamer = (aref (permissible-rotamer-vector permissible-rotamers) index)
          for locus = (monomer-shape-locus permissible-rotamer)
          for monomer-shape = (aref (topology:monomer-shape-vector oligomer-shape) locus)
          for allowed-rotamer-indexes = (allowed-rotamer-indexes permissible-rotamer)
          for rotamer-index = (aref vec index)
          do (unless (position rotamer-index allowed-rotamer-indexes)
               (error "You are trying to write a rotamer-index ~a into locus ~a but it is not one of the permissible ones: ~s" rotamer-index locus allowed-rotamer-indexes))
          do (setf (topology:rotamer-index monomer-shape) rotamer-index))
    (ensure-oligomer-shape-is-consistent-with-permissible-rotamers oligomer-shape permissible-rotamers)
    oligomer-shape))

(defmethod write-rotamers (oligomer-shape (permissible-backbone-rotamers permissible-backbone-rotamers) (vec array))
  (case (rotamers-state oligomer-shape)
    (:incomplete-no-rotamers
     (set-rotamers-state oligomer-shape :incomplete-backbone-rotamers)
     (call-next-method))
    (:incomplete-backbone-rotamers
     (call-next-method))
    (:complete-sidechain-and-backbone-rotamers
     (set-rotamers-state oligomer-shape :incomplete-backbone-rotamers)
     (call-next-method))
    (t (error "Illegal rotamers-state ~s" (rotamers-state oligomer-shape))))
  )

(defmethod write-rotamers (oligomer-shape (permissible-sidechain-rotamers permissible-sidechain-rotamers) (vec array))
  (case (rotamers-state oligomer-shape)
    (:incomplete-no-rotamers
     (error "You cannot write-rotamers to sidechain-rotamers when rotamers-state is ~s" (rotamers-state oligomer-shape)))
    ((:incomplete-backbone-rotamers :complete-sidechain-and-backbone-rotamers)
     (call-next-method)
     (set-rotamers-state oligomer-shape :complete-sidechain-and-backbone-rotamers))
    (t (error "Illegal rotamers-state ~s" (rotamers-state oligomer-shape)))))

(defmethod write-rotamers (oligomer-shape (permissible-rotamers permissible-rotamers) (val null))
  "Write the first allowed rotamer-index into the rotamer-index of each
monomer-shape in the oligomer-shape. Return the oligomer-shape."
  (let ((len (length (permissible-rotamer-vector permissible-rotamers))))
    (loop for index below len
          for permissible-rotamer = (aref (permissible-rotamer-vector permissible-rotamers) index)
          for locus = (monomer-shape-locus permissible-rotamer)
          for monomer-shape = (aref (topology:monomer-shape-vector oligomer-shape) locus)
          for allowed-rotamer-indexes = (allowed-rotamer-indexes permissible-rotamer)
          for rotamer-index = (progn
                                (unless (> (length allowed-rotamer-indexes) 0)
                                  (error "There must be at least one allowed-rotamer-indexes"))
                                (aref allowed-rotamer-indexes 0))
          do (setf (topology:rotamer-index monomer-shape) rotamer-index)))
  (ensure-oligomer-shape-is-consistent-with-permissible-rotamers oligomer-shape permissible-rotamers)
  oligomer-shape)

(defmethod write-rotamers (oligomer-shape (permissible-backbone-rotamers permissible-backbone-rotamers) (val null))
  (set-rotamers-state oligomer-shape :incomplete-no-rotamers)
  (call-next-method))

(defmethod write-rotamers (oligomer-shape (permissible-sidechain-rotamers permissible-sidechain-rotamers) (val null))
  (case (rotamers-state oligomer-shape)
    ((:incomplete-no-rotamers :incomplete-backbone-rotamers)
     (call-next-method))
    (:complete-sidechain-and-backbone-rotamers
     (set-rotamers-state oligomer-shape :incomplete-backbone-rotamers)
     (call-next-method))
    (t (error "Illegal rotamers-state ~s" (rotamers-state oligomer-shape)))))

;;; Switched to null specializer
#+(or)
(defmethod write-rotamers (oligomer-shape (permissible-rotamers permissible-rotamers) (val integer))
  "Write an integer into the rotamer-index of each monomer-shape in the oligomer-shape. Return the oligomer-shape."
  (error "This function should never be used - it's unsafe, the val would never satisfy the permissible-rotamers.
I suggest changing val to NIL and that forces the first allowed rotamer index to be written")
  (ensure-oligomer-shape-is-consistent-with-permissible-rotamers oligomer-shape permissible-rotamers)
  (let ((len (length (monomer-shape-locus permissible-rotamers))))
    (loop for index below len
          for permissible-rotamer = (aref (permissible-rotamer-vector permissible-rotamers) index)
          for locus = (monomer-locus permissible-rotamer)
          for monomer-shape = (aref (topology:monomer-shape-vector oligomer-shape) locus)
          for allowed-rotamer-indexes = (allowed-rotamer-indexes permissible-rotamer)
          for rotamer-index = val
          do (unless (member rotamer-index allowed-rotamer-indexes)
               (error "You are trying to write a rotamer-index that is not one of the permissible ones: ~s" rotamer-index allowed-rotamer-indexes))
          do (setf (topology:rotamer-index monomer-shape) rotamer-index)))
  oligomer-shape)

#+(or)
(defun number-of-rotamers (permissible-rotamers)
  "Return the total number of rotamer combinations"
  (let ((num 1))
    (loop for permissible-rotamer across (allowed-rotamer-indexes-vector permissible-rotamers)
          for allowed = (allowed-rotamer-indexes permissible-rotamer)
          do (setf num (* num (length allowed))))
    num))

#+(or)
(defun enumerated-rotamer-state (permissible-rotamers index)
  "Return the rotamer-state for a particular enumerated (by index) rotamers.
This is an vector of rotamer-index values enumerated from 0...(number-of-rotamers permissible-rotamers)"
  (error "Fix with permissible-rotamers")
  #+(or)
  (let ((num-rotamers (number-of-rotamers permissible-rotamers)))
    (unless (< index num-rotamers)
      (error "index ~s is out of bounds ~s" index num-rotamers))
    (let* ((num-rotamers (number-of-rotamers permissible-rotamers))
           (max-rotamers (coerce (max-rotamers permissible-rotamers) 'list))
           (rotamer-indices (core:positive-integer-to-mixed-base-digits index max-rotamers))
           (rotamer-vec (coerce (mapcar (lambda (indexes ee) (aref indexes ee))
                                        (coerce (allowed-rotamer-indexes-vector permissible-rotamers) 'list)
                                        rotamer-indices)
                                'vector)))
      rotamer-vec)))

#+(or)
(defun goto-rotamer (permissible-rotamers index)
  "Goto a particular rotamer using an integer index from 0...(number-of-rotamers stepper)"
  (let ((rotamer-state (rotamer-state permissible-rotamers index)))
      (write-rotamers permissible-rotamers rotamer-state)))

(defun zero-rotamers (permissible-rotamers)
  "Return a vector of zero for rotamer index - this shouldn't be written into a permissible-rotamers if you want a correct permissible-rotamers"
  (make-array (length (permissible-rotamer-vector permissible-rotamers)) :initial-element 0))

(defun first-rotamers (oligomer-shape permissible-rotamers)
  "Return a vector of the first rotamer index for each position"
  (declare (ignore oligomer-shape))
  (let ((len (length (permissible-rotamer-vector permissible-rotamers))))
    (make-array len :initial-contents
                (loop :for idx below len
                      :for permissible-rotamer = (aref (permissible-rotamer-vector permissible-rotamers) idx)
                      :for rnd = 0
                      :for index = (aref (allowed-rotamer-indexes permissible-rotamer) rnd)
                      :collect index))))

(defun random-rotamers (oligomer-shape permissible-rotamers)
  "Return a vector of random rotamer indices"
  (declare (ignore oligomer-shape))
  (let ((len (length (permissible-rotamer-vector permissible-rotamers))))
    (make-array len :initial-contents
                (loop :for idx below len
                      :for permissible-rotamer = (aref (permissible-rotamer-vector permissible-rotamers) idx)
                      :for rnd = (random (length (allowed-rotamer-indexes permissible-rotamer)))
                      :for index = (aref (allowed-rotamer-indexes permissible-rotamer) rnd)
                      :collect index))))

(defun apply-rotamers (oligomer-shape permissible-rotamers vec)
  "Return a vector of the first rotamer index for each position"
  (declare (ignore oligomer-shape))
  (let ((len (length (permissible-rotamer-vector permissible-rotamers))))
    (make-array len :initial-contents
                (loop :for idx below len
                      :for permissible-rotamer = (aref (permissible-rotamer-vector permissible-rotamers) idx)
                      :for val = (aref vec idx)
                      :for index = (aref (allowed-rotamer-indexes permissible-rotamer) val)
                      :collect index))))


#+(or)
(defun build-all-sidechains (assembler coords sidechain-stepper)
  "Increment through the sidechain-stepper and build the molecules and accumulate them into an aggregate.
Return that aggregate."
  (let (mols)
  (loop :do (progn
              (apply-permissible-rotamers sidechain-stepper)
              (topology:fill-internals-from-oligomer-shape-and-adjust assembler (oligomer-shape sidechain-stepper))
              (topology:build-all-atom-tree-external-coordinates-and-adjust assembler coords)
              (topology::copy-all-joint-positions-into-atoms assembler coords)
              (push (chem:matter-copy (chem:content-at (topology:aggregate assembler) 0))
                    mols)
              (unless (increment-permissible-rotamers sidechain-stepper)
                (return-from build-all-sidechains (chem:make-aggregate :agg mols)))))))


