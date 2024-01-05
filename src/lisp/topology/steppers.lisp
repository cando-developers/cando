(in-package :topology)

#+(or)
(eval-when (:compile-toplevel :load-toplevel :execute)
  (pushnew :debug-mover *features*))


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

(defun calculate-allowed-rotamers (monomer-shape-index sidechain-monomer-shape oligomer-shape shape-key-to-allowed-rotamers)
  (with-slots (topology:rotamers-map
               topology:oligomer
               ) oligomer-shape
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
        (error 'no-rotamers :message (format nil "Found no allowed-rotamers - phi-deg: ~s  psi-deg: ~s - shape-key ~s -  ~s in ~s - this shouldn't be possible if a backbone-stepper was created and applied to the oligomer-shape" phi-deg psi-deg shape-key sidechain-monomer-shape oligomer-shape)))
      #+debug-mover
      (let ((*print-pretty* nil))
        (format t "        allowed-rotamers -> ~s~%" result))
      result
      )))

(defclass permissible-rotamers ()
  ((monomer-shape-loci :initarg :monomer-shape-loci :reader monomer-shape-loci)
   (allowed-rotamer-indexes-vector :initarg :allowed-rotamer-indexes-vector :reader allowed-rotamer-indexes-vector)
   ))

(defmethod print-object ((obj permissible-rotamers) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (mapcar (lambda (mi ar) (format nil "#~s/~s" mi ar))
                                (coerce (monomer-shape-loci obj) 'list)
                                (coerce (allowed-rotamer-indexes-vector obj) 'list)))))


(defclass permissible-sidechain-rotamers (permissible-rotamers)
  ())

(defclass permissible-backbone-rotamers (permissible-rotamers)
  ())

(defun make-permissible-backbone-rotamers (oligomer-shape)
  "Limit the backbone rotamer indices to the ones that have sidechain rotamers with matching shape-keys"
  (let ((rotamers-db (topology:rotamers-map oligomer-shape))
        (monomer-shape-loci (make-array 16 :adjustable t :fill-pointer 0))
        (allowed-rotamer-indexes-vector (make-array 16 :adjustable t :fill-pointer 0))
        )
    (with-slots (topology:monomer-shape-vector
                 topology:oligomer
                 topology:rotamers-map
                 topology:monomer-shape-map) oligomer-shape
      (let ((special-sidechain-monomer-shapes (make-hash-table)))
        #+(or)(loop :for index below (length topology:monomer-shape-vector)
                    :for monomer-shape = (aref topology:monomer-shape-vector index)
                    :for monomer = (topology:monomer monomer-shape)
                    :for monomer-shape-kind = (topology:monomer-shape-kind monomer-shape)
                    :for monomer-context = (topology:monomer-context monomer-shape)
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
                              (topology:lookup-dihedral-cache oligomer-shape monomer-shape :psi)
                            #+debug-mover
                            (let ((backbone1-monomer-shape-index (position backbone1-monomer-shape topology:monomer-shape-vector)))
                              (format t "     sidechain-monomer-shape-index ~a   backbone1-monomer-shape-index ~a~%" index backbone1-monomer-shape-index))
                            (setf (gethash backbone1-monomer-shape special-sidechain-monomer-shapes) (cons monomer-shape monomer-context))))
        #+debug-mover(format t "-------- Next stage~%")
        (loop :for index below (length topology:monomer-shape-vector)
              :for monomer-shape = (aref topology:monomer-shape-vector index)
              :for monomer-shape-kind = (topology:monomer-shape-kind monomer-shape)
              :for monomer-context = (topology:monomer-context monomer-shape)
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
                    #+(or)((and (eq :backbone monomer-shape-kind)
                                maybe-sidechain)
                           ;; We have a backbone that influences a sidechain
                           ;; identify what backbone rotamers are allowed given the sidechain that it effects
                           (let ((allowed-rotamer-indexes (make-array 16 :adjustable t :fill-pointer 0)))
                             (loop :with backbone-rotamers = rotamers
                                   :for ii :below (length backbone-rotamers)
                                   :for backbone-rotamer = (aref backbone-rotamers ii)
                                   :for backbone-dihedral-cache-deg = (topology:backbone-dihedral-cache-deg backbone-rotamer)
                                   :for phi-1 = (getf backbone-dihedral-cache-deg :phi-1)
                                   :for psi-1 = (getf backbone-dihedral-cache-deg :psi-1)
                                   :for binned-phi-1 = (topology:bin-dihedral-deg phi-1)
                                   :for binned-psi-1 = (topology:bin-dihedral-deg psi-1)
                                   :for sidechain-rotamers = (gethash maybe-sidechain-monomer-context (topology:context-to-rotamers rotamers-db))
                                   :for shape-key-to-index = (topology:shape-key-to-index sidechain-rotamers)
                                   ;; Get the shape-key for the backbone
                                   :for shape-key = (make-phi-psi binned-phi-1 binned-psi-1)
                                   :for has-allowed-sidechain-rotamer-indices = (gethash shape-key shape-key-to-index)
                                   ;; If the sidechain has a shape-key that matches this backbone then
                                   ;;    add the backbone rotamer index to monomer-allowed-rotamer-indexes
                                   :do (progn
                                         #+debug-mover
                                         (let ((*print-pretty* nil))
                                           (format t "For backbone-index ~a rotamer-index-index ~a   shape-key ~a matches-sidechain -> ~a~%" index ii shape-key has-allowed-sidechain-rotamer-indices)))
                                   :when has-allowed-sidechain-rotamer-indices
                                     :do (let ((seen (position ii allowed-rotamer-indexes)))
                                           (unless seen (vector-push-extend ii allowed-rotamer-indexes))
                                           ))
                             (vector-push-extend allowed-rotamer-indexes allowed-rotamer-indexes-vector)
                             (vector-push-extend index monomer-shape-loci)))
                    ((eq :backbone monomer-shape-kind)
                     ;; We have a backbone that doesn't influence any sidechain
                     ;; - so all backbone rotamers are accessible
                     (let ((allowed-rotamer-indexes (make-array 16 :adjustable t :fill-pointer 0)))
                       (loop :with backbone-rotamers = rotamers
                             :for ii :below (length backbone-rotamers)
                             :do (vector-push-extend ii allowed-rotamer-indexes))
                       (vector-push-extend allowed-rotamer-indexes allowed-rotamer-indexes-vector)
                       (vector-push-extend index monomer-shape-loci)))
                    ((eq :sidechain monomer-shape-kind)
                     ;; Do nothing
                     )
                    (t (error "How did we get here?")))))
      (make-instance 'permissible-backbone-rotamers
                     :monomer-shape-loci (copy-seq monomer-shape-loci)
                     :allowed-rotamer-indexes-vector (copy-seq allowed-rotamer-indexes-vector)))))

(defun make-permissible-sidechain-rotamers (oligomer-shape &key monomers)
"Constructs a permissible-sidechain-rotamers object for the given oligomer-shape.
Optionally filters the process based on a provided list of monomers.
- oligomer-shape: The oligomer-shape to be processed.
- monomers: An optional list of monomers to restrict the processing, nil means all sidechain monomers.
Returns an instance of permissible-sidechain-rotamers."
  (let ((rotamers-db (topology:rotamers-map oligomer-shape))
        (monomer-shape-loci (make-array 16 :adjustable t :fill-pointer 0))
        (allowed-rotamer-indexes-vector (make-array 16 :adjustable t :fill-pointer 0))
        )
    ;; Loop over all monomer-shape in the oligomer-shape and downselect the ones that
    ;; are sidechains and if monomers is not nil then only the ones in monomers
    (loop :for index :from 0
          :for monomer-shape :across (topology:monomer-shape-vector oligomer-shape)
          :for monomer := (topology:monomer monomer-shape)
          :for monomer-shape-kind := (topology:monomer-shape-kind monomer-shape)
          :for monomer-context := (topology:monomer-context monomer-shape)
          :when (and (eq monomer-shape-kind :sidechain)
                     (or (null monomers)
                         (member monomer monomers))) ; if monomers defined then use it as a subset
            :do (let* ((rotamers (gethash monomer-context (topology:context-to-rotamers rotamers-db)))
                       (allowed-rotamers (calculate-allowed-rotamers index monomer-shape oligomer-shape (topology:shape-key-to-index rotamers))))
                  (unless allowed-rotamers
                    (error "Could not find any allowed-rotamers for ~s" monomer-shape))
                  (vector-push-extend allowed-rotamers allowed-rotamer-indexes-vector)
                  (vector-push-extend index monomer-shape-loci)))
    (make-instance 'permissible-sidechain-rotamers
                   :monomer-shape-loci monomer-shape-loci
                   :allowed-rotamer-indexes-vector allowed-rotamer-indexes-vector)))


(defun monomer-set-from-permissible-rotamers (&rest args)
  (loop for cur = args then (cddr cur)
        for oligomer-shape = (car cur)
        for permissible = (cadr cur)
        until (null cur)
        do (break "Check stuff")))

(defun max-rotamers (permissible-rotamers)
  "Read the maximum rotamer indexes from the stepper and return them in a vector"
  (let* ((len (length (monomer-shape-loci permissible-rotamers)))
         (vec (make-array len :element-type 'ext:byte32)))
    (loop for index below len
          for max = (length (aref (allowed-rotamer-indexes-vector permissible-rotamers) index))
          do (setf (aref vec index) max))
    vec))

(defmethod read-rotamers (oligomer-shape (permissible-rotamers permissible-rotamers))
  "Read the rotamer shapes from the stepper and return them in a vector"
  (let* ((len (length (monomer-shape-loci permissible-rotamers)))
         (vec (make-array len :element-type 'ext:byte32)))
    (loop for index below len
          for monomer-shape = (aref (topology:monomer-shape-vector oligomer-shape)
                                    (aref (monomer-shape-loci permissible-rotamers) index))
          do (setf (aref vec index) (topology:rotamer-index monomer-shape)))
    vec))

(defmethod write-rotamers (oligomer-shape (permissible-rotamers permissible-rotamers) (vec array))
  "Write a vector of rotamer shapes into the stepper"
  (let ((len (length (monomer-shape-loci permissible-rotamers))))
    (loop for index below len
          for monomer-shape = (aref (topology:monomer-shape-vector oligomer-shape)
                                    (aref (monomer-shape-loci permissible-rotamers) index))
          do (setf (topology:rotamer-index monomer-shape) (aref vec index)))
    oligomer-shape))

(defmethod write-rotamers (oligomer-shape (permissible-rotamers permissible-rotamers) (val integer))
  "Write a vector of rotamer shapes into the stepper. Return the oligomer-shape."
  (let ((len (length (monomer-shape-loci permissible-rotamers))))
    (loop for index below len
          for monomer-shape = (aref (topology:monomer-shape-vector oligomer-shape)
                                    (aref (monomer-shape-loci permissible-rotamers) index))
          do (setf (topology:rotamer-index monomer-shape) val))
    oligomer-shape))

(defun number-of-rotamers (permissible-rotamers)
  "Return the total number of rotamer combinations"
  (let ((num 1))
    (loop for allowed across (allowed-rotamer-indexes-vector permissible-rotamers)
          do (setf num (* num (length allowed))))
    num))

(defun goto-rotamer (permissible-rotamers index)
  "Goto a particular rotamer using an integer index from 0...(number-of-rotamers stepper)"
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
      (write-rotamers permissible-rotamers rotamer-vec))))

(defun zero-rotamers (permissible-rotamers)
  "Return a vector of zero for rotamer index - this shouldn't be written into a permissible-rotamers if you want a correct permissible-rotamers"
  (make-array (length (monomer-shape-loci permissible-rotamers)) :initial-element 0))

(defun first-rotamers (permissible-rotamers)
  "Return a vector of the first rotamer index for each position"
  (map 'vector (lambda (ee) (elt ee 0)) (allowed-rotamer-indexes-vector permissible-rotamers)))

(defun random-rotamers (permissible-rotamers)
  "Return a vector of random rotamer indices"
  (let ((len (length (monomer-shape-loci permissible-rotamers))))
    (make-array len :initial-contents
                (loop :for idx below len
                      :for rnd = (random (length (aref (allowed-rotamer-indexes-vector permissible-rotamers) idx)))
                      :for index = (aref (aref (allowed-rotamer-indexes-vector permissible-rotamers) idx) rnd)
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

