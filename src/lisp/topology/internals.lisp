(in-package :topology)

(defparameter *2pi-s* (* 2.0s0 (float PI 1.0s0)))
(defparameter *pi-s* (float PI 1.0s0))
(defparameter *-pi-s* (- (float PI 1.0s0)))

(defun deg-to-rad (deg)
  (* deg 0.0173433s0))

(defun rad-to-deg (rad)
  (/ rad 0.0173433s0))

(defun radians-limit (rad)
  "Limit the range of the angle in radians to greater than -PI and less or equal to PI"
  (if (<= rad *-pi-s*)
      (incf rad *2pi-s*)
      (if (> rad *pi-s*)
          (decf rad *2pi-s*)
          rad)))

(defun degrees-limit (deg)
  "Limit the range of the angle in degrees to greater than -180.0 and less or equal to 180.0"
  (if (< deg -180.0s0)
      (incf deg 360.0s0)
      (if (>= deg 180.0s0)
          (decf deg 360.0s0)
          deg)))

(defun radians-sub (b1 b2)
  (let ((diff (float (- b1 b2) 1.0s0)))
    (radians-limit diff)))

(defun radians-add (r1 r2)
  (let ((sum (float (+ r2 r1) 1.0s0)))
    (radians-limit sum)))

(defmacro radians-incf (place rad)
  `(setf ,place (radians-add ,place ,rad)))

(defun degrees-sub (b1 b2)
  (let ((diff (float (- b1 b2) 1.0s0)))
    (degrees-limit diff)))

(defun degrees-add (b1 b2)
  (let ((sum (float (+ b2 b1) 1.0s0)))
    (degrees-limit sum)))

(defun angle-sub (b1 b2)
  (degrees-sub b1 b2))

(defclass dihedral-info ()
  ((name :initarg :name :accessor name)))

(defclass dihedral-info-atom (dihedral-info)
  ((atom-name :initarg :atom-name :accessor atom-name)))

(cando.serialize:make-class-save-load
 dihedral-info-atom
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a" (name obj) (atom-name obj)))))

(defclass dihedral-info-external (dihedral-info)
  ((plug-path :initarg :plug-path :accessor plug-path)
   (external-dihedral-name :initarg :external-dihedral-name :accessor external-dihedral-name)))

(cando.serialize:make-class-save-load
 dihedral-info-external
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a ~a" (name obj) (plug-path obj) (external-dihedral-name obj)))))

(defun create-dihedral-info-from-constitution (constitution)
  ;; Build dihedrals using heavy atoms with more than one bond
  (loop for ca across (constitution-atoms constitution)
        with index = 0
        for element = (element ca)
        when (not (eq element :H))
          collect (make-instance 'dihedral-info-atom
                                 :name (intern (format nil "BB~a" (incf index)) :keyword)
                                 :atom-name (atom-name ca))))

(defun parse-dihedral-info (dihedrals)
  (loop for one in dihedrals
        for name = (second one)
        for info = (first one)
        collect (cond
                  ((symbolp info)
                   (let ((atom-name (intern (string info) :keyword)))
                     (make-instance 'dihedral-info-atom :name name :atom-name atom-name)))
                  ((consp info) (make-instance 'dihedral-info-external
                                               :name name
                                               :plug-path (first info)
                                               :external-dihedral-name (second info)))
                  (t (error "Could not parse ~s as dihedral-info" one)))))

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
  ()
  #+(or)
  ((bond :initarg :bond :accessor bond)
   (angle-rad :initarg :angle-rad :accessor angle-rad)
   (dihedral-rad :initarg :dihedral-rad :accessor dihedral-rad)))

(defmethod copy-internal ((internal bonded-internal))
  (make-instance 'bonded-internal
                 :name (name internal)))

(cando.serialize:make-class-save-load
 bonded-internal
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (name obj)))))

(defclass complex-bonded-internal (bonded-internal)
  ())

(cando.serialize:make-class-save-load
 complex-bonded-internal
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (name obj)))))

(defclass rotamer (cando.serialize:serializable)
  ((internals-values :initarg :internals-values :accessor internals-values)
   (delta-energy :initform 1.0s0 :initarg :delta-energy :accessor delta-energy)
   (probability :initform 1.0s0 :initarg :probability :accessor probability)))

(defclass fragment-internals (rotamer)
  ((trainer-index :initarg :trainer-index :accessor trainer-index)
   (monomer-context :initarg :monomer-context :accessor monomer-context)
   (energy :initform 0.0s0 :initarg :energy :accessor energy)
   (coordinates :initarg :coordinates :accessor coordinates)
   (shape-key-values :initarg :shape-key-values :accessor shape-key-values :documentation "Not needed")
   ))

(defmethod print-object ((obj fragment-internals) stream)
  (if *print-readably*
      (call-next-method)
      (let ((*print-pretty* nil))
        (print-unreadable-object (obj stream :type t)
          (format stream "~a :trainer-index ~a :delta-energy ~5,2f"
                  (monomer-context obj)
                  (trainer-index obj)
                  (delta-energy obj)
                  )
          #+(or)(format stream ":trainer-index ~a :energy ~10,3f first internals: ~{(~{~s ~6,2f~}) ~}" (trainer-index obj) (energy obj) (mapcar (lambda (x) (list (name x) (rad-to-deg (dihedral-rad x)))) internals ))))))

(defun dihedral-rad (internals-values joint-template)
  "Return the dihedral angle corresponding to the joint-template"
  (let ((dihedral-index (+ (* 3 (constitution-atoms-index joint-template)) 2)))
    (aref internals-values dihedral-index)))

(defun (setf dihedral-rad) (value internals-values joint-template)
  "Return the dihedral angle corresponding to the joint-template"
  (let ((dihedral-index (+ (* 3 (constitution-atoms-index joint-template)) 2)))
    (setf (aref internals-values dihedral-index) value)))


(defclass sidechain-rotamer (rotamer)
  (()))

(defclass backbone-rotamer (rotamer)
  ((backbone-dihedral-cache-deg :initarg :backbone-dihedral-cache-deg :accessor backbone-dihedral-cache-deg)))


(defun copy-fragment-internals (fragment-internals)
  "Shallow copy of fragment-internals"
  (make-instance 'fragment-internals
                 :monomer-context (monomer-context fragment-internals)
                 :trainer-index (trainer-index fragment-internals)
                 :probability (probability fragment-internals)
                 :energy (energy fragment-internals)
                 :delta-energy (delta-energy fragment-internals)
                 :internals-values (copy-seq (internals-values fragment-internals))
                 :coordinates (copy-seq (coordinates fragment-internals))
                 ))

(defclass rotamers (cando.serialize:serializable)
  ((rotamers :initarg :rotamers :initform (make-array 16 :adjustable t :fill-pointer 0) :accessor rotamers)))

(defclass sidechain-rotamers (rotamers)
  ((shape-key-to-index :initarg :shape-key-to-index :initform (make-hash-table :test 'equal) :accessor shape-key-to-index
                       :documentation "This is a cons of (phi.psi)")))


(defclass backbone-rotamers (rotamers)
  ())

(defclass rotamers-database (cando.serialize:serializable)
  ((context-to-rotamers :initarg :context-to-rotamers :initform (make-hash-table) :reader context-to-rotamers)
   (foldamer-name :initarg :foldamer-name :reader foldamer-name)))

(defmethod print-object ((obj rotamers-database) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t )
        (format stream "~a contexts" (hash-table-count (context-to-rotamers obj))))))

(defgeneric apply-fragment-internals-to-atresidue (fragment-internals rotamer-index atresidue))

(defmethod apply-fragment-internals-to-atresidue ((obj rotamer) rotamer-index atresidue)
  (setf (rotamer-index atresidue) rotamer-index)
  (loop for joint across (topology:joints atresidue)
        for index from 0
        do (multiple-value-bind (bond angle-rad dihedral-rad)
               (topology:extract-bond-angle-rad-dihedral-rad obj index)
             (when (typep joint 'kin:bonded-joint)
               (topology:fill-joint-internals joint bond angle-rad dihedral-rad)))))


(defmethod monomer-context-to-context-rotamers ((obj rotamers-database))
  (context-to-rotamers obj))

(defgeneric lookup-rotamers-for-context (fcc context &optional errorp))

(defmethod lookup-rotamers-for-context ((db rotamers-database) monomer-context &optional (errorp t))
  (let ((rot (gethash monomer-context (context-to-rotamers db))))
    (if rot
        rot
        (if errorp
            (error "Could not find rotamers for ~s" monomer-context)
            nil))))

(defun extract-dihedral-rad (assembler monomer atom-name)
  (let* ((atres (assembler-atresidue assembler monomer))
         (joint (joint-with-name atres atom-name)))
    (cond
      ((typep joint 'kin:complex-bonded-joint)
       (if (kin:complex-bonded-joint/phi-defined-p joint)
           (kin:bonded-joint/get-phi joint)
           0.0))
      ((typep joint 'kin:bonded-joint)
       (kin:bonded-joint/get-phi joint))
      ;; Otherwise return zero because its
      ;; a jump-joint or a complex-bonded-joint
      ;; and the user said it was an important dihedral
      ;; which it shouldn't be
      (t 0.0))))

(defun find-named-fragment-internals-rad (focused-assembler fragment-internals dihedral-names)
  (let* ((focus-monomer (focus-monomer focused-assembler))
         (oligomer (oligomer-containing-monomer focused-assembler focus-monomer))
         (focus-topology (chem:find-topology (oligomer-monomer-name-for-monomer oligomer focus-monomer)))
         (focus-constitution (constitution focus-topology))
         (focus-dihedrals (getf (residue-properties focus-constitution) :dihedrals))
         (focus-dihedral-info (loop for name in dihedral-names
                                    collect (find name focus-dihedrals :key #'name)))
         )
    (loop for info in focus-dihedral-info
          collect (cond
                    ((typep info 'dihedral-info-external)
                     (with-slots (plug-path external-dihedral-name) info
                       (let ((cur-monomer focus-monomer))
                         (loop for plug-name in plug-path
                               do (setf cur-monomer (monomer-on-other-side cur-monomer plug-name)))
                         (let* ((cur-topology (chem:find-topology (oligomer-monomer-name-for-monomer oligomer cur-monomer)))
                                (cur-constitution (constitution cur-topology))
                                (cur-dihedrals (getf (residue-properties cur-constitution) :dihedrals))
                                (cur-dihedral-info (find external-dihedral-name cur-dihedrals :key #'name)))
                           (with-slots (name atom-name) cur-dihedral-info
                             (extract-dihedral-rad focused-assembler cur-monomer atom-name))))))
                    ((typep info 'dihedral-info-atom)
                     (with-slots (name atom-name) info
                       (extract-dihedral-rad focused-assembler focus-monomer atom-name)))))))


#+(or)
(defun cluster-dihedral-vector (fragment-internals names)
  (error "Where am I called from")
  (coerce (loop for dihedral-name in names
                append (find-named-fragment-internals fragment-internals names))
          'vector))

(defun convert-dihedral-rad-vector-to-xy-vector (rad-vector)
  (make-array (* (length rad-vector) 2) :element-type 'single-float
                                        :initial-contents (loop for dihedral in rad-vector
                                                                collect (cos dihedral)
                                                                collect (sin dihedral))))

(defun cluster-dihedral-rad-vector (focused-assembler fragment-internals names)
(let ((coords (coordinates fragment-internals)))
  (update-ataggregate-joint-tree-internal-coordinates focused-assembler coords)
  (let ((dihedrals (find-named-fragment-internals-rad focused-assembler fragment-internals names)))
    dihedrals)))

(defun calculate-cluster-dihedral-names (focus-monomer oligomer)
  "Calculate the dihedral names for dihedrals that are used to cluster fragment-internals"
  (let* ((focus-topology (topology:monomer-topology focus-monomer oligomer))
         (focus-constitution (topology:constitution focus-topology))
         (focus-residue-properties (topology:residue-properties focus-constitution))
         (cluster-dihedrals-names (getf focus-residue-properties :cluster-dihedrals))
         (dihedral-infos (getf focus-residue-properties :dihedrals))
         ;; If the property cluster-dihedrals exists then use the dihedral names provided
         ;; otherwise use the names in the :dihedrals property
         (dihedral-names (if cluster-dihedrals-names
                             cluster-dihedrals-names
                             (loop for one-dihedral-info in dihedral-infos
                                   collect (name one-dihedral-info)))))
    dihedral-names
    ))


(defclass context-rotamers (serial:serializable)
  ((focus-monomer-name :initarg :focus-monomer-name :accessor focus-monomer-name)
   (monomer-context :initarg :monomer-context :accessor monomer-context)
   (next-index :initform 0 :initarg :next-index :accessor next-index)
   (internals :initarg :internals :accessor :internals)
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
   (print-unreadable-object (obj stream :type t)
     (format stream "~s members" (hash-table-count (monomer-context-to-context-rotamers obj))))))

(defclass rotamer-shape-connections ()
  ())

(defclass sidechain-rotamer-shape-connections ()
  ((phi-psi-map :initform (make-hash-table :test 'equal)
         :initarg :phi-psi-map
         :accessor phi-psi-map)))

(cando.serialize:make-class-save-load
 sidechain-rotamer-shape-connections
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (hash-table-count (phi-psi-map obj))))))

(defun make-sidechain-rotamer-shape-connections ()
  (make-instance 'sidechain-rotamer-shape-connections))

(defclass backbone-rotamer-shape-connections (rotamer-shape-connections)
  ((rotamer-indices :initform (make-array 16 :element-type 'ext:byte32 :fill-pointer 0 :adjustable t)
         :initarg :rotamer-indices
         :accessor rotamer-indices)))

(cando.serialize:make-class-save-load
 backbone-rotamer-shape-connections
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (length (rotamer-indices obj))))))

(defun make-backbone-rotamer-shape-connections (&optional (rotamer-indices nil rotamer-indices-p))
  (if rotamer-indices-p
      (make-instance 'backbone-rotamer-shape-connections :rotamer-indices rotamer-indices)
      (make-instance 'backbone-rotamer-shape-connections)))

(defmethod lookup-rotamer-shape-connections ((fsc sidechain-rotamer-shape-connections) key)
  (gethash key (phi-psi-map fsc)))

(defmethod lookup-rotamer-shape-connections ((fsc backbone-rotamer-shape-connections) key)
  (rotamer-indices fsc))

(defun append-rotamer-shape-connections (rsc key index)
  (let ((allowed-rotamer-vector (gethash key (phi-psi-map rsc))))
    (unless allowed-rotamer-vector
      (setf allowed-rotamer-vector (make-array 16 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
      (setf (gethash key (phi-psi-map rsc)) allowed-rotamer-vector))
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


(defun seen-fragment-internals (context-rotamers fragment-internals)
  (loop for seen-frag across (fragments context-rotamers)
        unless (rms-internals-are-different-p seen-frag fragment-internals)
          do (progn
               (return-from seen-fragment-internals (trainer-index seen-frag))))
  nil)


(defun bad-fragment-internals (internals-vector fragment-internals)
  (loop for internal across internals-vector
        for index from 0
        with previous-internal
        do (multiple-value-bind (bond angle-rad dihedral-rad)
               (topology:extract-bond-angle-rad-dihedral-rad fragment-internals index)
             (declare (ignore angle-rad dihedral-rad))
             (cond
               ((typep internal 'jump-internal) nil)
               ((> bond 3.0)
                (return-from bad-fragment-internals
                  (if previous-internal
                      (format nil "For atom ~a to ~a bad (bond internal) 3.0 .lt. ~7,2f" (name internal) (name previous-internal) bond)
                      (format nil "For atom ~a bad (bond internal) 3.0 .lt. ~7,2f" (name internal) bond))))))
        do (setf previous-internal internal))
  nil)

(defun save-clusterable-context-rotamers (clusterable-context-rotamers filename)
  (cando.serialize:save-cando clusterable-context-rotamers filename))

(defun load-clusterable-context-rotamers (filename)
  (cando.serialize:load-cando filename))

(defun dump-fragment-internals (internals fragment-internals finternals)
  (format finternals "begin-conformation ~a~%" (trainer-index fragment-internals))
  (flet ((to-deg (rad)
           (/ rad 0.0174533)))
    (loop for internal across internals
          for index from 0
          do (multiple-value-bind (bond angle-rad dihedral-rad)
                 (extract-bond-angle-rad-dihedral-rad fragment-internals index)
               (cond
                 ((typep internal 'topology:jump-internal)
                  (format finternals "jump-joint ~a~%" (topology:name internal)))
                 ((typep internal 'topology:complex-bonded-internal)
                  (format finternals "complex-bonded-joint ~a ~8,3f ~8,3f ~8,3f~%"
                          (topology:name internal)
                          bond (to-deg angle-rad) (to-deg dihedral-rad)))
                 ((typep internal 'topology:bonded-internal)
                  (format finternals "bonded-joint ~a ~8,3f ~8,3f ~8,3f~%"
                          (topology:name internal)
                          bond (to-deg angle-rad) (to-deg dihedral-rad)))
                 )))
    (format finternals "end-conformation~%")
    ))


(defgeneric fill-joint-internals (joint bond angle-rad dihedral-rad))

(defmethod fill-joint-internals ((joint kin:jump-joint) bond angle-rad dihedral-rad)
  )

(defmethod fill-joint-internals ((joint kin:xyz-joint) bond angle-rad dihedral-rad)
  "Do nothing with xyz-joints"
  )

(defmethod fill-joint-internals ((joint kin:bonded-joint) bond angle-rad dihedral-rad)
  (kin:set-distance joint bond)
  (kin:set-theta joint angle-rad)
  (kin:set-phi joint dihedral-rad)
  )

(defun extract-bond-angle-rad-dihedral-rad (rotamer index)
  (let ((index3 (* 3 index)))
    (values (aref (internals-values rotamer) index3)
            (aref (internals-values rotamer) (+ 1 index3))
            (aref (internals-values rotamer) (+ 2 index3)))))


(defmethod apply-fragment-internals-to-atresidue ((fragment-internals fragment-internals) rotamer-index atresidue)
  (setf (rotamer-index atresidue) rotamer-index)
  (loop for joint across (joints atresidue)
        for index from 0
        do (multiple-value-bind (bond angle-rad dihedral-rad)
               (extract-bond-angle-rad-dihedral-rad fragment-internals index)
             (fill-joint-internals joint bond angle-rad dihedral-rad))))

(defun analyze-atresidue (atresidue)
  (let ((internals-count 0)
        (internals-defined 0))
    (loop for joint across (joints atresidue)
          do (incf internals-count)
          when (kin:joint/definedp joint)
            do (incf internals-defined))
    (format t "  atresidue internals-count ~d~%" internals-count)
    (format t "  atresidue internals-defined ~d~%" internals-defined)
    ))

#|
(defun defined-dihedrals (assembler atresidue)
  (let ((residue (block find-residue
                  (do-atresidue-residue (search-atresidue search-residue assembler)
                    (when (eq atresidue search-atresidue)
                      (return-from find-residue residue)))
                  (error "Could not find atresidue ~s in ~s" atresidue assembler))))
  (let ((internals-count 0)
        (internals-defined 0))
    (loop for joint across (joints atresidue)
          do (incf internals-count)
          when (kin:joint/definedp joint)
            do (incf internals-defined))
    (format t "  atresidue internals-count ~d~%" internals-count)
    (format t "  atresidue internals-defined ~d~%" internals-defined)
    ))
|#

(defun write-internals (atresidue internals)
  (loop for joint in (joints atresidue)
        for index3 from 0 by 3
        when (typep joint 'kin:bonded-joint)
          do (progn
               (kin:bonded-joint/set-distance joint (aref internals index3))
               (kin:bonded-joint/set-theta joint (aref internals (+ 1 index3)))
               (kin:bonded-joint/set-phi joint (aref internals (+ 2 index3)))
               )))

(defgeneric write-internals-to-vector (joint index3 internals joint-mask))

(defmethod write-internals-to-vector ((joint kin:bonded-joint) index3 internals joint-mask)
  "xyz-joint get their distance,theta,phi coordinates"
  (setf (aref internals index3) (kin:bonded-joint/get-distance joint)
        (aref internals (+ 1 index3)) (kin:bonded-joint/get-theta joint)
        (aref internals (+ 2 index3)) (kin:bonded-joint/get-phi joint)
        (aref joint-mask (/ index3 3)) 1))

(defmethod write-internals-to-vector ((joint kin:xyz-joint) index3 internals joint-mask)
  "xyz-joint get write their x,y,z coordinates"
  (let ((pos (kin:xyz-joint/get-pos joint)))
    (setf (aref internals index3) (geom:vx pos)
          (aref internals (+ 1 index3)) (geom:vy pos)
          (aref internals (+ 2 index3)) (geom:vz pos)
          (aref joint-mask (/ index3 3)) 0)))

(defun extract-internals (atresidue)
  "Extract the internals into a vector of internals and a mask of what internals were defined"
  (let* ((num-joints (length (joints atresidue)))
         (internals (make-array (* 3 num-joints) :element-type (geom:vecreal-type)))
         (joint-mask (make-array num-joints :element-type 'bit :initial-element 0)))
    (loop for joint across (joints atresidue)
          for index from 0
          for index3 from 0 by 3
          when (kin:definedp joint)
            do (write-internals-to-vector joint index3 internals joint-mask))
    (values internals joint-mask)))


(defun indexes-of-internals (atresidue joints)
  (mapcar (lambda (joint)
            (position joint (topology:joints atresidue)))
          joints))

(defun write-merged-internals (target-internals joint-mask children-indexes old-internals new-internals)
  (multiple-value-bind (set-child-indexes unset-child-indexes)
      (loop for child-index in children-indexes
            if (= (aref joint-mask child-index) 1)
              collect child-index into set-child-indexes
            else
              collect child-index into unset-child-indexes
            finally (return (values set-child-indexes unset-child-indexes)))
    (if (null set-child-indexes)
        (loop for child-index in children-indexes
              for index3 = (* 3 child-index)
              do (setf
                  (aref target-internals index3) (aref new-internals index3)
                  (aref target-internals (+ 1 index3)) (aref new-internals (+ 1 index3))
                  (aref target-internals (+ 2 index3)) (aref new-internals (+ 2 index3)))
              )
        (progn
          (when (> (length set-child-indexes) 2)
            (error "set-child-indexes is ~s - it should never be more than 2 elements long" set-child-indexes))
          (let ((delta-dihedrals (loop for child-index in set-child-indexes
                                       for ref-index = (first set-child-indexes) ; index into joint-mask
                                       for ref-index3 = (* 3 ref-index) ; get the index into the internals
                                       for ref-old-dihedral = (aref old-internals (+ 2 ref-index3)) ; The old dihedral we keep
                                       for ref-new-dihedral = (aref new-internals (+ 2 ref-index3)) ; The new dihedral we rotate to old
                                       for delta-dihedral = (topology:radians-sub ref-new-dihedral ref-old-dihedral) ; the delta to add to every other dihedral
                                       collect delta-dihedral)))
            (when (= (length delta-dihedrals) 2)
              (when (> (abs (radians-sub (first delta-dihedrals) (second delta-dihedrals))) 0.1)
                (error "In write-merged-internals the delta-dihedrals are not closely matched: ~s - this will happen if there is a mismatch in the stereochemistry of terminal methyl groups like in LEU and VAL" delta-dihedrals)))
            ;; If delta-dihedrals aren't all the same value then the stereochemistry might be different
            ;; from the rotamer database vs the structure we loaded.
            ;; For example, ILE, VAL or LEU may have the terminal methyl groups swapped.
            (let* ((delta-dihedral (first delta-dihedrals)))
              (loop for child-index in unset-child-indexes
                    for index3 = (* 3 child-index)
                    do (setf
                        (aref target-internals index3) (aref new-internals index3) ; distance
                        (aref target-internals (+ 1 index3)) (aref new-internals (+ 1 index3)) ; angle
                        (aref target-internals (+ 2 index3)) (radians-sub (aref new-internals (+ 2 index3)) delta-dihedral) ; dihedral
                        )))
            )))))

(defun merge-internals (atresidue old-internals joint-mask new-internals)
  (let ((joint-index 0)
        joint
        children
        children-indexes
        (temp-internals (copy-seq old-internals)))
    (tagbody
     top
       (setf joint (aref (topology:joints atresidue) joint-index))
       (setf children (topology:children joint)
             children-indexes (indexes-of-internals atresidue children)
             )
       (when (null children) (go next-iteration))
       (write-merged-internals temp-internals joint-mask children-indexes old-internals new-internals)
       (progn
         (format t "dihedral3 atom ~s in atresidue: ~s~%" joint (name atresidue))
         (format t "  children: ~s~%" children)
         (format t "  updated dihedrals: ~%")
         (loop for joint across (joints atresidue)
               for index from 0
               for index3 from 0 by 3
               for name = (kin:joint/name joint)
               for dihdeg = (rad-to-deg (aref temp-internals (+ 2 index3)))
               for old-dihdeg = (rad-to-deg (aref old-internals (+ 2 index3)))
               do (format t "  ~5s  ~10,3f  ~10,3f~%" name dihdeg old-dihdeg)
               ))
     next-iteration
       (incf joint-index)
       (when (< joint-index (length (joints atresidue))) (go top))
       )
    temp-internals
    ))

