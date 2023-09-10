(in-package :topology)

(defparameter *2pi-s* (* 2.0s0 (float PI 1.0s0)))
(defparameter *pi-s* (float PI 1.0s0))
(defparameter *-pi-s* (- (float PI 1.0s0)))

(defun deg-to-rad (deg)
  (* deg 0.0173433s0))

(defun rad-to-deg (rad)
  (/ rad 0.0173433s0))

(defun radians-difference (b1 b2)
  (let ((diff (float (mod (- b2 b1) *2pi-s*) 1.0s0)))
    (if (< diff *-pi-s*)
        (incf diff *2pi-s*)
        (if (> diff *pi-s*)
            (decf diff *2pi-s*)
            diff))))

(defun radians-add (r1 r2)
  (let ((sum (float (mod (+ r2 r1) *2pi-s*) 1.0s0)))
    (if (< sum *-pi-s*)
        (incf sum *2pi-s*)
        (if (> sum *pi-s*)
            (decf sum *2pi-s*)
            sum))))

(defmacro radians-incf (place rad)
  `(setf ,place (radians-add ,place ,rad)))

(defun degree-difference (b1 b2)
  (let ((diff (float (mod (- b2 b1) 360.0s0) 1.0s0)))
    (if (< diff -180.0s0)
        (incf diff 360.0s0)
        (if (> diff 180.0s0)
            (decf diff 360.0s0)
            diff))))

(defun degrees-add (b1 b2)
  (let ((sum (float (mod (+ b2 b1) 360.0s0) 1.0s0)))
    (if (< sum -180.0s0)
        (incf sum 360.0s0)
        (if (> sum 180.0s0)
            (decf sum 360.0s0)
            sum))))

(defun angle-sub (b2 b1)
  (degree-difference b1 b2))

(defclass internal (serial:serializable)
  ((name :initarg :name :accessor name)
   ))

(defgeneric copy-internal (internal))

(cando.serialize:make-class-save-load
 internal
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (name obj)))))

(defclass jump-internal (internal)
  ())

(cando.serialize:make-class-save-load
 jump-internal
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (name obj)))))

(defclass bonded-internal (internal)
  ((bond :initarg :bond :accessor bond)
   (angle-rad :initarg :angle-rad :accessor angle-rad)
   (dihedral-rad :initarg :dihedral-rad :accessor dihedral-rad)))

(defmethod copy-internal ((internal bonded-internal))
  (when (ext:float-nan-p (dihedral-rad internal))
    (error "Found a NAN dihedral"))
  (make-instance 'bonded-internal
                 :name (name internal)
                 :bond (bond internal)
                 :angle-rad (angle-rad internal)
                 :dihedral-rad (dihedral-rad internal)))

(cando.serialize:make-class-save-load
 bonded-internal
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a :d ~5,2f" (name obj) (rad-to-deg (dihedral-rad obj))))))

(defclass complex-bonded-internal (bonded-internal)
  ())

(cando.serialize:make-class-save-load
 complex-bonded-internal
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (name obj)))))

(defclass fragment-internals (serial:serializable)
  ((trainer-index :initarg :trainer-index :accessor trainer-index)
   (probability :initform 1.0s0 :initarg :probability :accessor probability)
   (energy :initform 0.0s0 :initarg :energy :accessor energy)
   (internals :initarg :internals :accessor internals)
   (shape-key-values :initform nil :initarg :shape-key-values :accessor shape-key-values)
   (coordinates :initarg :coordinates :accessor coordinates)
   (out-of-focus-internals :initarg :out-of-focus-internals :accessor out-of-focus-internals)
   ))

(cando.serialize:make-class-save-load
 fragment-internals
 :print-unreadably
 (lambda (obj stream)
   (let ((*print-pretty* nil))
   (print-unreadable-object (obj stream :type t)
     (let ((internals (coerce (subseq (internals obj) 0 (min (length (internals obj)) 4) ) 'list)))
       (format stream ":trainer-index ~a :energy ~5,2f"
               (trainer-index obj)
               (energy obj)
               )
       #+(or)(format stream ":trainer-index ~a :energy ~10,3f first internals: ~{(~{~s ~6,2f~}) ~}" (trainer-index obj) (energy obj) (mapcar (lambda (x) (list (name x) (rad-to-deg (dihedral-rad x)))) internals )))))))

(defun copy-fragment-internals (fragment-internals)
  (make-instance 'fragment-internals
                 :trainer-index (trainer-index fragment-internals)
                 :probability (probability fragment-internals)
                 :energy (energy fragment-internals)
                 :internals (copy-seq (internals fragment-internals))
                 :shape-key-values (copy-seq (shape-key-values fragment-internals))
                 :coordinates (copy-seq (coordinates fragment-internals))
                 :out-of-focus-internals (let ((ht (make-hash-table)))
                                           (maphash (lambda (key value)
                                                     (setf (gethash key ht) (copy-seq value)))
                                                    (out-of-focus-internals fragment-internals))
                                           ht)))

(defun find-named-fragment-internals (fragment-internals name-or-plug-names)
  (cond
    ((symbolp name-or-plug-names)
     (loop for internal across (internals fragment-internals)
           when (eq name-or-plug-names (name internal))
             do (cond
                  ((typep internal 'bonded-internal)
                   (return (list (topology:dihedral-rad internal))))
                  (t (return (list 0.0))))
           finally (error "Could not find internal with atom-name ~a" name-or-plug-names)))
    ((consp name-or-plug-names)
     (let* ((plug-name (first name-or-plug-names))
            (out-of-focus-internals (gethash plug-name (out-of-focus-internals fragment-internals))))
       (loop for name in (cdr name-or-plug-names)
             collect (loop for internal across out-of-focus-internals
                           when (eq (name internal) name)
                             do (return (topology:dihedral-rad internal))))))
    (t (error "Handle arg ~a" name-or-plug-names))))

(defun cluster-dihedral-vector (fragment-internals names)
  (coerce (loop for name-or-plug-names in names
                append (find-named-fragment-internals fragment-internals name-or-plug-names))
          'vector))

(defun convert-dihedral-rad-vector-to-xy-vector (rad-vector)
  (make-array (* (length rad-vector) 2) :element-type 'single-float
                                        :initial-contents (loop for dihedral in rad-vector
                                                                collect (cos dihedral)
                                                                collect (sin dihedral))))

(defun cluster-dihedral-rad-vector (fragment-internals names)
  (let ((dihedrals (loop for name-or-plug-names in names
                         append (find-named-fragment-internals fragment-internals name-or-plug-names))))
    dihedrals))

(defun calculate-cluster-dihedral-names (focus-monomer oligomer)
  "Calculate the atom-names for dihedrals that are used to cluster fragment-internals"
  (let* ((focus-topology (topology:monomer-topology focus-monomer oligomer))
         (focus-constitution (topology:constitution focus-topology))
         (focus-residue-properties (topology:residue-properties focus-constitution))
         (dihedral-info (getf focus-residue-properties :dihedrals))
         (focus-names (loop for dihedral-name in dihedral-info
                            for dihedral-atom-or-plug-name = (car dihedral-name)
                            collect dihedral-atom-or-plug-name)))
    focus-names
    ))


(defclass context-rotamers (serial:serializable)
  ((focus-monomer-name :initarg :focus-monomer-name :accessor focus-monomer-name)
   (monomer-context :initarg :monomer-context :accessor monomer-context)
   (next-index :initform 0 :initarg :next-index :accessor next-index)
   (fragments :initform (make-array 16 :adjustable t :fill-pointer 0) :initarg :fragments :accessor fragments)))

(cando.serialize:make-class-save-load
 context-rotamers
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~s number-of-fragments ~a" (monomer-context obj) (length (fragments obj))))))

(defmethod context-rotamers ((obj context-rotamers))
  obj)

(defun merge-context-rotamers (context-rotamers-list)
  (unless (> (length context-rotamers-list) 0)
    (error "There must be at least one context-rotamers in list"))
  (let ((unsorted-fragment-internals (make-array 1024 :adjustable t :fill-pointer 0))
        (max-next-index 0))
    (loop for one in context-rotamers-list
          do (setf max-next-index (max max-next-index (next-index one)))
          do (loop for frag-int across (fragments one)
                   do (vector-push-extend frag-int unsorted-fragment-internals)))
    (let ((new-fragment-internals (sort unsorted-fragment-internals #'< :key #'trainer-index)))
      (make-instance 'context-rotamers
                     :focus-monomer-name (focus-monomer-name (first context-rotamers-list))
                     :monomer-context (monomer-context (first context-rotamers-list))
                     :next-index max-next-index
                     :fragments new-fragment-internals))))

(defclass clusterable-context-rotamers (serial:serializable)
  ((foldamer-name :initform :spiroligomer
                   :initarg :foldamer-name :accessor foldamer-name)
   (cluster-dihedral-names :initarg :cluster-dihedral-names :accessor cluster-dihedral-names)
   (context-rotamers :initarg :context-rotamers :accessor context-rotamers)))

(cando.serialize:make-class-save-load clusterable-context-rotamers)

(defun merge-clusterable-context-rotamers (clusterable-context-rotamers-list)
  "Merge a list of clusterable-context-rotamers.
No checking is done to make sure that the list of clusterable-context-rotamers are all compatible"
  (let* ((context-rotamers-list (loop for one in clusterable-context-rotamers-list
                                           collect (context-rotamers one)))
        (context-rotamers (merge-context-rotamers context-rotamers-list))
         (dihedral-names (cluster-dihedral-names (first clusterable-context-rotamers-list)))
         (foldamer-name (foldamer-name (first clusterable-context-rotamers-list)))
         )
    (make-instance 'clusterable-context-rotamers
                   :foldamer-name foldamer-name
                   :cluster-dihedral-names dihedral-names
                   :context-rotamers context-rotamers)))

(defclass context-rotamers-map (serial:serializable)
  ((monomer-context-to-context-rotamers :initform (make-hash-table :test 'equal)
                                              :initarg :monomer-context-to-context-rotamers
                                              :accessor monomer-context-to-context-rotamers)))

(cando.serialize:make-class-save-load context-rotamers-map
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t))))

(defclass rotamer-shape-connections ()
  ())

(defclass sidechain-rotamer-shape-connections ()
  ((fmap :initform (make-hash-table :test 'equal)
         :initarg :fmap
         :accessor fmap)))

(cando.serialize:make-class-save-load
 sidechain-rotamer-shape-connections)

(defun make-sidechain-rotamer-shape-connections ()
  (make-instance 'sidechain-rotamer-shape-connections))

(defclass backbone-rotamer-shape-connections (rotamer-shape-connections)
  ((rotamer-indices :initform (make-array 16 :element-type 'ext:byte32 :fill-pointer 0 :adjustable t)
         :initarg :rotamer-indices
         :accessor rotamer-indices)))

(cando.serialize:make-class-save-load
 backbone-rotamer-shape-connections)

(defun make-backbone-rotamer-shape-connections ()
  (make-instance 'backbone-rotamer-shape-connections))

(defmethod lookup-rotamer-shape-connections ((fsc sidechain-rotamer-shape-connections) key)
  (gethash key (fmap fsc)))

(defmethod lookup-rotamer-shape-connections ((fsc backbone-rotamer-shape-connections) key)
  (rotamer-indices fsc))

(defun append-rotamer-shape-connections (rsc key index)
  (let ((allowed-rotamer-vector (gethash key (fmap rsc))))
    (unless allowed-rotamer-vector
      (setf allowed-rotamer-vector (make-array 16 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
      (setf (gethash key (fmap rsc)) allowed-rotamer-vector))
    (vector-push-extend index allowed-rotamer-vector)))

(defclass rotamer-context-connections ()
  ((fmap :initform (make-hash-table :test 'eq)
         :initarg :fmap
         :accessor fmap)))

(cando.serialize:make-class-save-load
 rotamer-context-connections)

(defun make-rotamer-context-connections ()
  (make-instance 'rotamer-context-connections))

(defun map-rotamer-context-connections (op fcc)
  (maphash (lambda (key value)
             (let ((from (car key))
                   (to (cdr key)))
               (funcall op from to value)))
           (fmap fcc)))

(defun rotamer-context-connections-count (fcc)
  (hash-table-count (fmap fcc)))

(defun set-rotamer-context-connections (fcc from to value)
  (check-type value rotamer-shape-connections)
  (setf (gethash (cons from to) (fmap fcc)) value))

(defun lookup-rotamer-context-connections (fcc key)
  (gethash key (fmap fcc)))

(defun set-rotamer-context-connections (rcc key value)
  (setf (gethash key (fmap rcc)) value))

(defclass connected-rotamers-map (context-rotamers-map)
  ((rotamer-context-connections :initform (make-rotamer-context-connections)
                                :initarg :rotamer-context-connections
                                :accessor rotamer-context-connections)))

(cando.serialize:make-class-save-load
 connected-rotamers-map
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream " rotamer-context-connections-count ~d length-fragments ~d" (rotamer-context-connections-count (rotamer-context-connections obj))
             (hash-table-count (monomer-context-to-context-rotamers obj))))))

(defun matched-context-rotamers-summary (connected-rotamers-map)
  (let ((total-context-rotamers 0)
        (matching-context-rotamers 0)
        (missing-context-rotamers 0))
    (maphash (lambda (key value)
               (declare (ignore key))
               (incf total-context-rotamers (length (fragments value))))
             (monomer-context-to-context-rotamers connected-rotamers-map))
    (map-rotamer-context-connections
     (lambda (from-monomer-context to-monomer-context value)
       (declare (ignore from-monomer-context to-monomer-context))
       (incf matching-context-rotamers (length value)))
     (rotamer-context-connections connected-rotamers-map))
    (map-rotamer-context-connections
     (lambda (from-monomer-context to-monomer-context value)
       (declare (ignore from-monomer-context to-monomer-context))
       (loop for val across value
             when (= (length val) 0)
               do (incf missing-context-rotamers)))
     (rotamer-context-connections connected-rotamers-map))
    (let ((missing-monomer-contexts nil))
      (maphash (lambda (monomer-context context-rotamers)
                 (declare (ignore context-rotamers))
                 (block inner-search
                   (map-rotamer-context-connections
                    (lambda (from-monomer-context to-monomer-context allowed-fragment-indices)
                      (declare (ignore allowed-fragment-indices))
                      (when (or (string= from-monomer-context monomer-context)
                                (string= to-monomer-context monomer-context))
                        (return-from inner-search nil)))
                    (rotamer-context-connections connected-rotamers-map))
                   (push monomer-context missing-monomer-contexts)))
               (monomer-context-to-context-rotamers connected-rotamers-map))
      (values total-context-rotamers matching-context-rotamers missing-context-rotamers missing-monomer-contexts))))

(defparameter *dihedral-threshold* 20.0)
(defparameter *dihedral-rms-threshold* 15.0)

(defun internals-are-different-p (frag1 frag2 &optional )
  (loop for frag1-int across (internals frag1)
        for frag2-int across (internals frag2)
        do (when (and (typep frag1-int 'bonded-internal)
                      (typep frag2-int 'bonded-internal))
             (let* ((dihedral-frag1 (rad-to-deg (dihedral-rad frag1-int)))
                    (dihedral-frag2 (rad-to-deg (dihedral-rad frag2-int)))
                    (aa (degree-difference dihedral-frag1 dihedral-frag2))
                    (different-p (> (abs aa) *dihedral-threshold*)))
               (when different-p
                 (return-from internals-are-different-p t)))))
  nil)

(defun dihedrals-rms (frag1 frag2)
  (let ((sum-of-squares 0.0)
        (num-squares 0))
  (loop for frag1-int across (internals frag1)
        for frag2-int across (internals frag2)
        do (when (and (typep frag1-int 'bonded-internal)
                      (typep frag2-int 'bonded-internal))
             (let* ((dihedral-frag1 (rad-to-deg (dihedral-rad frag1-int)))
                    (dihedral-frag2 (rad-to-deg (dihedral-rad frag2-int)))
                    (aa (degree-difference dihedral-frag1 dihedral-frag2)))
               (incf sum-of-squares (* aa aa))
               (incf num-squares))))
    (let ((rms (sqrt (/ sum-of-squares (float num-squares)))))
      rms)))

(defun rms-internals-are-different-p (frag1 frag2)
  (> (dihedrals-rms frag1 frag2) *dihedral-rms-threshold*))

(defun seen-fragment-internals (context-rotamers fragment-internals)
  (loop for seen-frag across (fragments context-rotamers)
        unless (rms-internals-are-different-p seen-frag fragment-internals)
          do (progn
               (return-from seen-fragment-internals (trainer-index seen-frag))))
  nil)


(defun bad-fragment-internals (fragment-internals)
  (loop for internal across (internals fragment-internals)
        with previous-internal
        do (cond
             ((typep internal 'jump-internal) nil)
             ((> (bond internal) 3.0)
              (return-from bad-fragment-internals
                (if previous-internal
                    (format nil "For atom ~a to ~a bad (bond internal) 3.0 .lt. ~7,2f" (name internal) (name previous-internal) (bond internal))
                    (format nil "For atom ~a bad (bond internal) 3.0 .lt. ~7,2f" (name internal) (bond internal))))))
        do (setf previous-internal internal))
  nil)

(defun save-clusterable-context-rotamers (clusterable-context-rotamers filename)
  (cando.serialize:save-cando clusterable-context-rotamers filename))

(defun load-clusterable-context-rotamers (filename)
  (cando.serialize:load-cando filename))

(defun dump-fragment-internals (fragment-internals finternals)
  (format finternals "begin-conformation ~a~%" (trainer-index fragment-internals))
  (flet ((to-deg (rad)
           (/ rad 0.0174533)))
    (loop for internal across (topology:internals fragment-internals)
          do (cond
               ((typep internal 'topology:jump-internal)
                (format finternals "jump-joint ~a~%" (topology:name internal)))
               ((typep internal 'topology:complex-bonded-internal)
                (format finternals "complex-bonded-joint ~a ~8,3f ~8,3f ~8,3f~%"
                        (topology:name internal)
                        (topology:bond internal)
                        (to-deg (topology:angle-rad internal))
                        (to-deg (topology:dihedral-rad internal))))
               ((typep internal 'topology:bonded-internal)
                (format finternals "bonded-joint ~a ~8,3f ~8,3f ~8,3f~%"
                        (topology:name internal)
                        (topology:bond internal)
                        (to-deg (topology:angle-rad internal))
                        (to-deg (topology:dihedral-rad internal))))
               ))
    (format finternals "end-conformation~%")
    (let ((unique-internals (let (unique-internals)
                              (maphash (lambda (key internals)
                                         (declare (ignore key))
                                         (pushnew internals unique-internals))
                                       (topology:out-of-focus-internals fragment-internals))
                              unique-internals)))
      (loop for internals in unique-internals
            do (when internals
                 (loop for internal across internals
                       do (cond
                            ((typep internal 'topology:jump-internal)
                             (format finternals "jump-joint ~a~%" (topology:name internal)))
                            ((typep internal 'topology:complex-bonded-internal)
                             (format finternals "complex-bonded-joint ~a ~8,3f ~8,3f ~8,3f~%"
                                     (topology:name internal)
                                     (topology:bond internal)
                                     (to-deg (topology:angle-rad internal))
                                     (to-deg (topology:dihedral-rad internal))))
                            ((typep internal 'topology:bonded-internal)
                             (format finternals "bonded-joint ~a ~8,3f ~8,3f ~8,3f~%"
                                     (topology:name internal)
                                     (topology:bond internal)
                                     (to-deg (topology:angle-rad internal))
                                     (to-deg (topology:dihedral-rad internal))))
                            )))))))


(defgeneric fill-joint-internals (joint bond angle-rad dihedral-rad))

(defmethod fill-joint-internals ((joint kin:jump-joint) bond angle-rad dihedral-rad)
  )

(defmethod fill-joint-internals ((joint kin:bonded-joint) bond angle-rad dihedral-rad)
  (kin:set-distance joint bond)
  (kin:set-theta joint angle-rad)
  (kin:set-phi joint dihedral-rad)
  )

(defgeneric extract-bond-angle-rad-dihedral-rad (internal))

(defmethod extract-bond-angle-rad-dihedral-rad ((internal jump-internal))
  (values 0.0 0.0 0.0))

(defmethod extract-bond-angle-rad-dihedral-rad ((internal bonded-internal))
  (break "Check that internal angles are in radians")
  (values (bond internal) (angle-rad internal) (dihedral-rad internal)))

(defgeneric apply-fragment-internals-to-atresidue (fragment-internals atresidue))

(defmethod apply-fragment-internals-to-atresidue ((fragment-internals fragment-internals) atresidue)
 (loop for joint across (joints atresidue)
       for internal across (internals fragment-internals)
       do (multiple-value-bind (bond angle-rad dihedral-rad)
              (extract-bond-angle-rad-dihedral-rad internal)
            (fill-joint-internals joint bond angle-rad dihedral-rad))))

