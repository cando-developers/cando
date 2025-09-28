(in-package :topology)

(defparameter *2pi-s* (* 2.0s0 (float PI 1.0s0)))
(defparameter *pi-s* (float PI 1.0s0))
(defparameter *-pi-s* (- (float PI 1.0s0)))
(defconstant +bin-size+ 10 )


(defclass shape-key (cando.serialize:serializable)
  ((parts :initarg :parts :reader parts)))

(defmethod print-object ((shape-key shape-key) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (shape-key stream :type t)
        (format stream "~s" (parts shape-key)))))


(defun make-shape-key (&rest pparts)
  (make-instance 'shape-key :parts pparts))

(defun ensure-shape-key (shape-key)
  (unless (typep shape-key 'shape-key)
    (error "The shape-key ~s must be a shape-key" shape-key))
  shape-key)


(defclass shape-key-map ()
  ((ht :initarg :ht :reader ht)))

(defun shape-key-map-hash (value)
  (sxhash (parts value)))

(defun shape-key-map-test (v1 v2)
  (equal (parts v1) (parts v2)))

(defmethod make-shape-key-map ()
  (make-hash-table :hash-function #'shape-key-map-hash :test #'shape-key-map-test))



(defclass dihedral-range (cando.serialize:serializable)
  ((min-phi :initarg :min-phi :reader :min-phi)
   (max-phi :initarg :max-phi :reader :max-phi)
   (min-psi :initarg :min-psi :reader :min-psi)
   (max-psi :initarg :max-psi :reader :max-psi)))

(defun bin-dihedral-range (dihedral-range)
  "Convert DIGITAL-RANGE into a decimal range and step size"
  (with-slots (min-phi max-phi min-psi max-psi) dihedral-range
    (let ((digital-min-phi (* +bin-size+ (floor min-phi +bin-size+)))
          (digital-max-phi (* +bin-size+ (floor max-phi +bin-size+)))
          (digital-min-psi (* +bin-size+ (ceiling min-psi +bin-size+)))
          (digital-max-psi (* +bin-size+ (ceiling max-psi +bin-size+))))
      (values digital-min-phi digital-max-phi +bin-size+
              digital-min-psi digital-max-psi +bin-size+))))

(defun bin-dihedral-deg (dih-deg)
  "Take a dihedral in degrees and round it to the nearest +bin-size+ degrees"
  (when (and (floatp dih-deg) (ext:float-nan-p dih-deg))
    (error "dih-deg is NAN"))
  (floor (degrees-limit (floor (* +bin-size+ (round dih-deg +bin-size+))))))


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

(defclass one-ring-info (cando.serialize:serializable)
  ((name-to-atom-names :initarg :name-to-atom-names :reader name-to-atom-names)))

(defclass complete-ring-info (one-ring-info) ())
(defclass start-ring-info (one-ring-info) ())
(defclass end-ring-info (one-ring-info) ())



(defclass ring-info (cando.serialize:serializable)
  ((complete :initform (make-hash-table) :initarg :complete :reader complete)
   (start :initform (make-hash-table) :initarg :start :reader :start)
   (end :initform (make-hash-table) :initarg :end :reader end)))

(defun check-ring-args (rings)
  (loop with complete-count = 0
        with start-count = 0
        with end-count = 0
        for (k v) on rings by #'cddr
        unless (member k '(:complete :start :end))
          do (error "Illegal ring argument ~s" k)
        when (eq k :complete)
          do (progn
               (incf complete-count)
               (unless (consp v)
                 (error "The argument of :complete must be a list - not ~s" v)))
        when (eq k :start)
          do (progn
               (incf start-count)
               (unless (consp v)
                 (error "The argument of :start must be a list - not ~s" v)))
        when (eq k :end)
          do (progn
               (incf end-count)
               (unless (consp v)
                 (error "The argument of :end must be a list - not ~s" v)))
        finally (progn
                  (when (> complete-count 1)
                    (error "There should be only one :complete entry"))
                  (when (> start-count 1)
                    (error "There should be only one :start entry"))
                  (when (> end-count 1)
                    (error "There should be only one :end entry")))))

(defun parse-ring-info (rings)
  (check-ring-args rings)
  (let ((complete-plist (getf rings :complete))
        (start-plist (getf rings :start))
        (end-plist (getf rings :end)))
    (make-instance 'ring-info
                   :complete (alexandria:plist-hash-table complete-plist)
                   :start (alexandria:plist-hash-table start-plist)
                   :end (alexandria:plist-hash-table end-plist))))

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
   (free-energy :initform 1.0s0 :initarg :free-energy :accessor free-energy)
   (probability :initform 1.0s0 :initarg :probability :accessor probability)))

(defclass fragment-internals (rotamer)
  ((trainer-index :initarg :trainer-index :accessor trainer-index)
   (monomer-context :initarg :monomer-context :accessor monomer-context)
   (energy :initform 0.0s0 :initarg :energy :accessor energy)
   (coordinates :initarg :coordinates :accessor coordinates)
   ))

(defclass shape-key-cache (cando.serialize:serializable)
  ((cache :initform nil :type list :initarg :cache :accessor cache))
  (:documentation "Maintain cache as a map of shape-key names to shape-keys"))

(defmethod print-object ((obj shape-key-cache) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (format stream "~s" (cache obj)))))

(defgeneric ignore-shape-key-cache (topology foldamer)
  (:documentation "Ignore the shape-key-cache for the topology in the foldamer"))

(defmethod ignore-shape-key-cache (topology foldamer)
  "Generally don't ignore the shape-key-cache"
  (declare (ignore topology foldamer))
  nil)


(defun make-shape-key-cache ()
  "Create an empty SHAPE-KEY-CACHE"
  (make-instance 'shape-key-cache))

(defun add-to-cache (shape-key-cache name shape-key)
  "Add a NAME SHAPE-KEY pair to the SHAPE-KEY-CACHE"
  (unless (or (integerp shape-key)
              (and (consp shape-key)
                   (every (lambda (vv) (integerp vv)) shape-key)))
    (error "shape-key ~s must be an integer or a list of integers" shape-key))
  (push (cons name shape-key) (cache shape-key-cache)))

(defun find-in-cache (shape-key-cache name)
  "Return the dihedral-angle associated with the NAME in SHAPE-KEY-CACHE"
  (let ((entry (assoc name (cache shape-key-cache))))
    (cdr entry)))

(defun cache-alist-copy (shape-key-cache)
  "Return an alist version of the SHAPE-KEY-CACHE that can be altered without
changing the SHAPE-KEY-CACHE."
  (copy-seq (cache shape-key-cache)))

(defclass fragment-internals-with-shape-key-cache (fragment-internals)
  ((shape-key-cache-deg :initarg :shape-key-cache-deg
                                :type shape-key-cache
                                :accessor shape-key-cache-deg)))

(defun make-fragment-internals-with-shape-key-cache-from-fragment-internals (fragment-internals shape-key-cache-deg)
  "Shallow copy of fragment-internals"
  (make-instance 'fragment-internals-with-shape-key-cache
                 :monomer-context (monomer-context fragment-internals)
                 :trainer-index (trainer-index fragment-internals)
                 :probability (probability fragment-internals)
                 :energy (energy fragment-internals)
                 :free-energy (free-energy fragment-internals)
                 :internals-values (copy-seq (internals-values fragment-internals))
                 :coordinates (copy-seq (coordinates fragment-internals))
                 :shape-key-cache-deg shape-key-cache-deg))

(defmethod print-object ((obj fragment-internals) stream)
  (if *print-readably*
      (call-next-method)
      (let ((*print-pretty* nil))
        (print-unreadable-object (obj stream :type t)
          (format stream "~a :trainer-index ~a :free-energy ~5,2f"
                  (monomer-context obj)
                  (trainer-index obj)
                  (free-energy obj)
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

(defmethod print-object ((obj rotamer) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (format stream ":probability ~5,3f" (probability obj)))))

(defun make-sidechain-rotamer (&key internals-values free-energy probability)
  (make-instance 'sidechain-rotamer
                 :internals-values internals-values
                 :free-energy free-energy
                 :probability probability))

(defun make-sidechain-rotamer-from-fragment-internals (fragment-internals)
  (check-type fragment-internals fragment-internals)
  (make-sidechain-rotamer
   :internals-values (internals-values fragment-internals)
   :free-energy (free-energy fragment-internals)
   :probability (probability fragment-internals)))


(defclass backbone-rotamer-base (rotamer)
  ())

(defclass backbone-with-sidechain-rotamer (backbone-rotamer-base)
  ((shape-key-cache-deg :initform (make-shape-key-cache)
                                :initarg :shape-key-cache-deg
                                :type shape-key-cache
                                :accessor shape-key-cache-deg)))

(defmethod print-object ((obj backbone-with-sidechain-rotamer) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (format stream ":shape-key-cache-deg ~s" (shape-key-cache-deg obj)))))

(defun make-backbone-with-sidechain-rotamer (&key internals-values free-energy probability shape-key-cache-deg)
  (make-instance 'backbone-with-sidechain-rotamer
                 :internals-values internals-values
                 :free-energy free-energy
                 :probability probability
                 :shape-key-cache-deg shape-key-cache-deg))

(defun make-backbone-with-sidechain-rotamer-from-fragment-internals (fragment-internals)
  (check-type fragment-internals fragment-internals-with-shape-key-cache)
  (make-backbone-with-sidechain-rotamer
   :internals-values (internals-values fragment-internals)
   :free-energy (free-energy fragment-internals)
   :probability (probability fragment-internals)
   :shape-key-cache-deg (shape-key-cache-deg fragment-internals)))


(defclass backbone-without-sidechain-rotamer (backbone-rotamer-base)
  ())

(defun make-backbone-without-sidechain-rotamer (&key internals-values free-energy probability)
  (make-instance 'backbone-without-sidechain-rotamer
                 :internals-values internals-values
                 :free-energy free-energy
                 :probability probability))

(defun make-backbone-without-sidechain-rotamer-from-fragment-internals (fragment-internals)
  (check-type fragment-internals fragment-internals)
  (make-backbone-without-sidechain-rotamer
   :internals-values (internals-values fragment-internals)
   :free-energy (free-energy fragment-internals)
   :probability (probability fragment-internals)))


(defun copy-fragment-internals (fragment-internals)
  "Shallow copy of fragment-internals"
  (make-instance 'fragment-internals
                 :monomer-context (monomer-context fragment-internals)
                 :trainer-index (trainer-index fragment-internals)
                 :probability (probability fragment-internals)
                 :energy (energy fragment-internals)
                 :free-energy (free-energy fragment-internals)
                 :internals-values (copy-seq (internals-values fragment-internals))
                 :coordinates (copy-seq (coordinates fragment-internals))
                 ))


(defclass rotamers (cando.serialize:serializable)
  ((monomer-context :initarg :monomer-context :reader monomer-context)
   (rotamer-vector :initarg :rotamer-vector
                   :initform (make-array 16 :adjustable t :fill-pointer 0)
                   :accessor rotamer-vector
                   :type (vector rotamer *))))

(defmethod print-object ((object rotamers) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (object stream :type t)
        (format stream "~s" (if (slot-boundp object 'monomer-context)
                                (monomer-context object)
                                "no monomer-context")))))


(defclass sidechain-rotamers (rotamers)
  ((shape-key-to-index :initarg :shape-key-to-index :initform (make-hash-table :test 'equal) :accessor shape-key-to-index
                       :documentation "This is a cons of (phi.psi)")))

(defclass backbone-rotamers-base (rotamers) ())

(defclass backbone-without-sidechain-rotamers (backbone-rotamers-base)
  ())

(defclass backbone-with-sidechain-rotamers (backbone-rotamers-base)
  ())


(defclass rotamers-conformers-database (cando.serialize:serializable)
  ((context-to-rotamers :initarg :context-to-rotamers :initform (make-hash-table) :reader context-to-rotamers)
   (foldamer-name :initarg :foldamer-name :reader foldamer-name)))

(defclass rotamers-database (rotamers-conformers-database) ())

(defclass conformers-database (rotamers-conformers-database) ())

(defun get-conformers (monomer-context rotamers-conformers-database)
  "Return the conformers for the MONOMER-CONTEXT in the ROTAMERS-CONFORMERS-DATABASE"
  (gethash monomer-context (context-to-rotamers rotamers-conformers-database)))

(defmethod print-object ((obj rotamers-conformers-database) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t )
        (format stream "~a contexts" (hash-table-count (context-to-rotamers obj))))))

(defun validate-indexes-for-shape-keys (rotamers-database)
  "Summarize the ROTAMERS-DATABASE. Print monomer-contexts and shape-keys that have empty rotamer-index vectors."
  (let ((empty-indexes 0)
        (filled-indexes 0))
  (maphash (lambda (monomer-context rotamers)
             (when (typep rotamers 'sidechain-rotamers)
               (let ((shape-key-to-index (shape-key-to-index rotamers)))
                 (maphash (lambda (shape-key indexes)
                            (if (= (length indexes) 0)
                                (progn
                                  (incf empty-indexes)
                                  (format t "Empty indexes for ~s ~s~%" monomer-context shape-key))
                                (incf filled-indexes)))
                          shape-key-to-index))))
           (context-to-rotamers rotamers-database))
    (format t "There were ~d empty indexes and ~d filled ones~%" empty-indexes filled-indexes)))


(defgeneric apply-fragment-internals-to-atresidue (assembler assembler-internals fragment-internals rotamer-index atresidue))

(defmethod apply-fragment-internals-to-atresidue (assembler assembler-internals (rotamer rotamer) rotamer-index atresidue)
  (setf (rotamer-index atresidue) rotamer-index)
  (loop for joint across (topology:joints atresidue)
        for index3 from 0 by 3
        do (fill-joint-internals-from-vector assembler assembler-internals joint (internals-values rotamer) index3)))

(defmethod monomer-context-to-context-rotamers ((obj rotamers-database))
  (context-to-rotamers obj))

(defmethod monomer-context-to-context-rotamers ((obj conformers-database))
  (context-to-rotamers obj))

(defgeneric lookup-rotamers-for-context (fcc context &key errorp ignore-sidechains))

(defmethod lookup-rotamers-for-context ((db rotamers-database) monomer-context &key (errorp t) ignore-sidechains)
  (when ignore-sidechains
    (error "ignore-sidechains is only allowed for a conformers-database"))
  (let ((rot (gethash monomer-context (context-to-rotamers db))))
    (if rot
        rot
        (if errorp
            (error "Could not find rotamers for ~s" monomer-context)
            nil))))

(defmethod lookup-rotamers-for-context ((db conformers-database) monomer-context &key (errorp t) (error-value :no-sidechain-rotamer) ignore-sidechains)
  (unless ignore-sidechains
    (error "ignore-sidechains must be set for conformers-database"))
  (let ((rot (gethash monomer-context (context-to-rotamers db))))
    (if rot
        rot
        (if errorp
            (error "Could not find rotamers for ~s" monomer-context)
            error-value))))

(defun extract-dihedral-joint (assembler monomer atom-name)
  (let* ((atres (assembler-atresidue assembler monomer))
         (joint (joint-with-name atres atom-name)))
    joint))

(defun find-named-joints (focused-assembler focus-monomer dihedral-names)
  "Gather up the joints that are named in dihedral-names."
  (let* ((oligomer (oligomer-containing-monomer focused-assembler focus-monomer))
         (focus-topology (chem:find-topology (oligomer-monomer-name-for-monomer oligomer focus-monomer)))
         (focus-constitution (constitution focus-topology))
         (focus-dihedrals (getf (residue-properties focus-constitution) :dihedrals))
         (focus-dihedral-info (loop for name in dihedral-names
                                    for info = (find name focus-dihedrals :key #'name)
                                    if info
                                      collect info
                                    else
                                      do (error "Could not find focus-dihedrals name ~s in ~s" name focus-dihedrals))))
    (loop for info in focus-dihedral-info
          for joint = (cond
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
                                 (extract-dihedral-joint focused-assembler cur-monomer atom-name))))))
                        ((typep info 'dihedral-info-atom)
                         (with-slots (name atom-name) info
                           (extract-dihedral-joint focused-assembler focus-monomer atom-name))))
          if joint
            collect joint
          else
            do (error "About to return a NIL joint for ~s" info))))

(defun extract-dihedral-rad-from-joint (joint internals)
  "Return the internal dihedral value and if not available return (values 0.0 nil).
    If it is available return (values dihedral-rad joint). "
    (cond
      ((typep joint 'kin:complex-bonded-joint)
       (if (kin:complex-bonded-joint/phi-defined-p joint internals)
           (values (kin:bonded-joint/get-phi joint) joint internals)
           (values 0.0 nil)))
      ((typep joint 'kin:bonded-joint)
       (values (kin:bonded-joint/get-phi joint internals) joint))
      ;; Otherwise return zero because its
      ;; a jump-joint or a complex-bonded-joint
      ;; and the user said it was an important dihedral
      ;; which it shouldn't be
      (t (values 0.0 nil))))

(defun extract-dihedral-line-segment-from-joint (joint)
  (let* ((parent-joint (if (kin:joint/parent-bound-p joint)
                           (kin:joint/parent joint)
                           nil))
         (grand-parent-joint (if (and parent-joint (kin:joint/parent-bound-p parent-joint))
                                 (kin:joint/parent parent-joint)
                                 nil))
         (great-grand-parent-joint (if (and grand-parent-joint (kin:joint/parent-bound-p grand-parent-joint))
                                       (kin:joint/parent grand-parent-joint)
                                       nil)))
      (let* ((joint-index-x3 (and joint (kin:joint/position-index-x3 joint)))
             (parent-joint-index-x3 (and parent-joint (kin:joint/position-index-x3 parent-joint)))
             (grand-parent-joint-index-x3 (and grand-parent-joint (kin:joint/position-index-x3 grand-parent-joint)))
             (great-grand-parent-joint-index-x3 (and great-grand-parent-joint (kin:joint/position-index-x3 great-grand-parent-joint))))
        (list joint-index-x3 parent-joint-index-x3 grand-parent-joint-index-x3 great-grand-parent-joint-index-x3))))

(defun find-named-fragment-internals-rad (focused-assembler dihedral-names internals)
  (let ((joints (find-named-joints focused-assembler (focus-monomer focused-assembler) dihedral-names)))
    (mapcar (lambda (joint)
              (extract-dihedral-rad-from-joint joint internals))
            joints)))

(defun joint-to-root-p (joint ht)
  "Return true if a path to the root can be found without hitting any joint in ht."
  (let ((cur joint))
    (loop
      if (kin:joint/parent-bound-p cur)
        do (progn
             (setf cur (kin:joint/parent cur))
             (when (gethash cur ht)
               (return nil)))
      else
        do (return t))))


(defun find-joint-closest-to-root (joints)
  (let ((jht (loop with ht = (make-hash-table)
                   for joint in joints
                   do (setf (gethash joint ht) t)
                   finally (return ht))))
    (loop for joint in joints
          when (joint-to-root-p joint jht)
            do (return-from find-joint-closest-to-root joint)))
  (error "Could not find a joint closest to root joints: ~s" joints))

(defun find-named-fragment-canonical-transform (focused-assembler focus-monomer dihedral-names coords)
  (let* ((joints (find-named-joints focused-assembler focus-monomer dihedral-names))
         (joint-closest-to-root (find-joint-closest-to-root joints))
         (maybe-transform nil))
    (when (kin:joint/parent-bound-p joint-closest-to-root)
      (let ((parent (kin:joint/parent joint-closest-to-root)))
        (when (kin:joint/parent-bound-p parent)
          (let ((grand-parent (kin:joint/parent parent)))
            (let* ((jindex3 (kin:joint/position-index-x3 joint-closest-to-root))
                   (pindex3 (kin:joint/position-index-x3 parent))
                   (gpindex3 (kin:joint/position-index-x3 grand-parent))
                   (jx (aref coords jindex3))
                   (jy (aref coords (+ 1 jindex3)))
                   (jz (aref coords (+ 2 jindex3)))
                   (px (aref coords pindex3))
                   (py (aref coords (+ 1 pindex3)))
                   (pz (aref coords (+ 2 pindex3)))
                   (gpx (aref coords gpindex3))
                   (gpy (aref coords (+ 1 gpindex3)))
                   (gpz (aref coords (+ 2 gpindex3)))
                   (origin (geom:vec jx jy jz))
                   (xpoint (geom:vec px py pz))
                   (xypoint (geom:vec gpx gpy gpz))
                   (xvec (geom:vnormalized (geom:v- xpoint origin)))
                   (xyvec (geom:vnormalized (geom:v- xypoint origin)))
                   (zvec (geom:vnormalized (geom:vcross xvec xyvec)))
                   (yvec (geom:vnormalized (geom:vcross zvec xvec)))
                   (toorigin (geom:make-m4-translate (geom:v* origin -1.0)))
                   (rotate (geom:make-m4-rotation-rows xvec yvec zvec))
                   (transform (geom:m*m rotate toorigin))
                   )
              (setf maybe-transform transform))))))
    (if maybe-transform
        maybe-transform
        (geom:make-matrix t))))

(defun find-named-fragment-line-segments (focused-assembler dihedral-names)
  (let* ((joints (find-named-joints focused-assembler (focus-monomer focused-assembler) dihedral-names)))
    (mapcar (lambda (joint) (extract-dihedral-line-segment-from-joint joint))
            joints)))

#+(or)
(defun cluster-dihedral-vector (fragment-internals names)
  (error "Where am I called from")
  (coerce (loop for dihedral-name in names
                append (find-named-fragment-internals fragment-internals names))
          'vector))

(defun convert-dihedral-rad-vector-to-xy-vector (rad-vector)
  (make-array (* (length rad-vector) 2) :element-type (geom:vecreal-type)
                                        :initial-contents (loop for dihedral in rad-vector
                                                                collect (cos dihedral)
                                                                collect (sin dihedral))))

(defun cluster-dihedral-rad-vector (focused-assembler internals names &optional coords)
  (when coords
    (update-ataggregate-joint-tree-internal-coordinates focused-assembler coords))
  (find-named-fragment-internals-rad focused-assembler names internals))

(defun cluster-dihedral-line-segments (focused-assembler fragment-internals names)
  "For debugging purposes, return line-segments for the dihedrals using in clustering"
  (let ((dihedrals (find-named-fragment-line-segments focused-assembler names)))
    dihedrals))

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
   (fragments :initform (make-array 16 :adjustable t :fill-pointer 0) :initarg :fragments :accessor fragments)))

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
  ((rotamer-indexes :initform (make-array 16 :element-type 'ext:byte32 :fill-pointer 0 :adjustable t)
         :initarg :rotamer-indexes
         :accessor rotamer-indexes)))

(cando.serialize:make-class-save-load
 backbone-rotamer-shape-connections
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (length (rotamer-indexes obj))))))

(defun make-backbone-rotamer-shape-connections (&optional (rotamer-indexes nil rotamer-indexes-p))
  (if rotamer-indexes-p
      (make-instance 'backbone-rotamer-shape-connections :rotamer-indexes rotamer-indexes)
      (make-instance 'backbone-rotamer-shape-connections)))

(defmethod lookup-rotamer-shape-connections ((fsc sidechain-rotamer-shape-connections) key)
  (gethash key (phi-psi-map fsc)))

(defmethod lookup-rotamer-shape-connections ((fsc backbone-rotamer-shape-connections) key)
  (rotamer-indexes fsc))

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
                    (lambda (from-monomer-context to-monomer-context allowed-fragment-indexes)
                      (declare (ignore allowed-fragment-indexes))
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

(defgeneric fill-joint-internals-from-vector (assembler assembler-internals joint vector index3))

(defmethod fill-joint-internals-from-vector (assembler assembler-internals (joint kin:jump-joint) vector index3) #|do nothing|#)
(defmethod fill-joint-internals-from-vector (assembler assembler-internals (joint kin:xyz-joint) vector index3) #|do nothing|#)

(defmethod fill-joint-internals-from-vector (assembler assembler-internals (joint kin:bonded-joint) vector index3)
  (multiple-value-bind (bond angle-rad dihedral-rad)
      (values (aref vector index3)
              (aref vector (+ 1 index3))
              (aref vector (+ 2 index3)))
    (joint-set-distance assembler assembler-internals joint bond)
    (joint-set-theta assembler assembler-internals joint angle-rad)
    (joint-set-phi assembler assembler-internals joint dihedral-rad)))

(defmethod fill-joint-internals-from-vector (assembler assembler-internals (joint kin:bonded-joint) (vector sys:simple-vector-double) index3)
  (kin:fill-internals-from-simple-vector-double joint assembler-internals vector index3))

(defmethod fill-joint-internals-from-vector (assembler assembler-internals (joint kin:bonded-joint) (vector sys:simple-vector-float) index3)
  (kin:fill-internals-from-simple-vector-single joint assembler-internals vector index3))


(defun fill-joint-phi (joint phi assembler internals)
  "Write a dihedral angle in radians in PHI into the JOINT.
I have this to trap writing specific phi values for debugging."
  (kin:bonded-joint/set-phi joint internals phi))

#+(or)
(defmethod fill-joint-internals ((joint kin:bonded-joint) bond angle-rad dihedral-rad)
  (kin:set-distance joint bond)
  (kin:set-theta joint angle-rad)
  (fill-joint-phi joint dihedral-rad)
  )

(defun extract-bond-angle-rad-dihedral-rad (rotamer index)
"Extract the bond, angle-rad, and dihedral-rad from the ROTAMER at INDEX.
INDEX counts by 1 and _not_ by 3."
  (let ((index3 (* 3 index)))
    (values (aref (internals-values rotamer) index3)
            (aref (internals-values rotamer) (+ 1 index3))
            (aref (internals-values rotamer) (+ 2 index3)))))


(defmethod apply-fragment-internals-to-atresidue (assembler assembler-internals (fragment-internals fragment-internals) rotamer-index atresidue)
  (setf (rotamer-index atresidue) rotamer-index)
  (loop for joint across (joints atresidue)
        for index3 from 0 by 3
        do (fill-joint-internals-from-vector assembler assembler-internals joint (internals-values fragment-internals) index3)))

(defun analyze-atresidue (assembler assembler-internals atresidue)
  (let ((internals-count 0)
        (internals-defined 0))
    (loop for joint across (joints atresidue)
          do (incf internals-count)
          when (kin:joint/definedp joint assembler-internals)
            do (incf internals-defined))
    (format t "  atresidue internals-count ~d~%" internals-count)
    (format t "  atresidue internals-defined ~d~%" internals-defined)
    ))

(defun write-internals (assembler assembler-internals atresidue internals-values updated-internals-mask &key verbose)
  "Write the distance,angle,dihedral triplets in INTERNALS-VALUES into ASSEMBLER-INTERNALS"
  (loop for joint across (joints atresidue)
        for index from 0
        for index3 from 0 by 3
        if (or (null updated-internals-mask)
               (= (aref updated-internals-mask index) 1))
          do (typecase joint
               (kin:bonded-joint
                (when verbose (format t "Updating bonded-joint ~s~%" joint))
                (fill-joint-internals-from-vector assembler assembler-internals joint internals-values index3))
               (kin:xyz-joint
                (if (kin:xyz-joint/definedp joint assembler-internals)
                    (warn "Skipping defining xyz-joint")
                    (break "What do we do with an undefined xyz-joint?")))
               (t (error "Add support for ~s~%" joint)))
        else
          do (when verbose (format t "NOT updating joint ~s~%" joint))))

(defgeneric write-internals-to-vector (assembler assembler-internals joint index3 internals joint-mask))

(defmethod write-internals-to-vector (assembler assembler-internals (joint kin:bonded-joint) index3 internals joint-mask)
  "xyz-joint get their distance,theta,phi coordinates"
  (setf (aref internals index3) (kin:bonded-joint/get-distance joint assembler-internals)
        (aref internals (+ 1 index3)) (kin:bonded-joint/get-theta joint assembler-internals)
        (aref internals (+ 2 index3)) (kin:bonded-joint/get-phi joint assembler-internals)
        (aref joint-mask (/ index3 3)) 1))

(defmethod write-internals-to-vector (assembler assembler-internals (joint kin:xyz-joint) index3 internals joint-mask)
  "xyz-joint get write their x,y,z coordinates"
  (let ((pos (kin:xyz-joint/transformed-pos joint)))
    (setf (aref internals index3) (geom:vx pos)
          (aref internals (+ 1 index3)) (geom:vy pos)
          (aref internals (+ 2 index3)) (geom:vz pos)
          (aref joint-mask (/ index3 3)) 1)))

(defun heavy-atom-p (joint constitution-atoms)
  (let* ((name (kin:joint/name joint))
         (ca (let ((ca (find name constitution-atoms :key #'atom-name)))
               (unless ca (error "Could not find constitution atom with name ~s in ~s" name constitution-atoms))
               ca))
         (ca-element (element ca)))
    (not (chem:element-is-hydrogen ca-element))))

(defun extract-internals (assembler assembler-internals atresidue)
  "Extract the internals into a vector of internals and a mask of what internals were defined.
Return the internals, a joint-mask of 1 for each joint that has defined internals and a heavy-atom-mask
that is 1 for each heavy atom."
  (let* ((num-joints (length (joints atresidue)))
         (internals (make-array (* 3 num-joints) :element-type (geom:vecreal-type)))
         (joint-mask (make-array num-joints :element-type 'bit :initial-element 0))
         (heavy-atom-mask (make-array num-joints :element-type 'bit :initial-element 0))
         (topology (chem:find-topology (name atresidue)))
         (constitution-atoms (constitution-atoms (constitution topology))))
    (loop for joint across (joints atresidue)
          for index from 0
          for index3 from 0 by 3
          when (heavy-atom-p joint constitution-atoms)
            do (setf (aref heavy-atom-mask index) 1)
          when (kin:definedp joint assembler-internals)
            do (write-internals-to-vector assembler joint index3 internals joint-mask))
    (values internals joint-mask heavy-atom-mask)))


(defun indexes-of-internals (atresidue joints)
  (mapcar (lambda (joint)
            (position joint (topology:joints atresidue)))
          joints))

(defun maybe-merge-internals (joint target-internals joint-mask children-indexes joints to-internals new-internals heavy-atom-mask updated-internals-mask verbose)
  (multiple-value-bind (set-child-indexes unset-child-indexes)
      (loop for child-index in children-indexes
            if (= (aref joint-mask child-index) 1)
              collect child-index into set-child-indexes
            else
              collect child-index into unset-child-indexes
            finally (return (values set-child-indexes unset-child-indexes)))
    (cond
      ((and (null set-child-indexes)
            (loop for index in unset-child-indexes
                  when (= (aref heavy-atom-mask index) 1) ; a child is a heavy atom and not set
                    do (return nil)
                  finally (return t)))
       ;; all children are unset and hydrogen as in methyl groups - all three hydrogens are undefined
       (loop for child-index in children-indexes
             for index3 = (* 3 child-index)
             do (when verbose (format t "terminal H/methyl case - setting updated-internals-mask to 1~%"))
             do (setf
                 (aref target-internals index3) (aref new-internals index3)
                 (aref target-internals (+ 1 index3)) (aref new-internals (+ 1 index3))
                 (aref target-internals (+ 2 index3)) (aref new-internals (+ 2 index3))
                 (aref updated-internals-mask child-index) 1))
       (when verbose
         (format t "Will set internals for children of ~s~%  children ~s at ~s~%"
                 joint
                 (mapcar (lambda (idx) (elt joints idx)) children-indexes)
                 children-indexes)))
      ((null set-child-indexes)
       ;; do nothing when there are unset heavy atoms
       (when verbose
         (format t "Will NOT set internals for children of joint ~s~%  children ~s at ~s~%"
                 joint
                 (mapcar (lambda (idx) (elt joints idx)) children-indexes)
                 children-indexes))
       )
      ((null unset-child-indexes)  ; all children are set - do nothing
       #| Do nothing |#
       (when verbose
         (format t "For joint ~s~%  All children ~s at ~s are already set~%"
                 joint
                 (mapcar (lambda (idx) (elt joints idx)) children-indexes)
                 children-indexes)))
      ((> (length set-child-indexes) 2)
       (error "set-child-indexes is ~s - it should never be more than 2 elements long: children -> ~s~% Perhaps you need to rebuild the assembler because you may have already setup all the internals for this one" set-child-indexes (topology:children joint)))
      (t (let ((delta-dihedrals (loop for child-index in set-child-indexes
                                      for ref-index = (first set-child-indexes) ; index into joint-mask
                                      for ref-index3 = (* 3 ref-index) ; get the index into the internals
                                      for ref-old-dihedral = (aref to-internals (+ 2 ref-index3)) ; The old dihedral we keep
                                      for ref-new-dihedral = (aref new-internals (+ 2 ref-index3)) ; The new dihedral we rotate to old
                                      for delta-dihedral = (topology:radians-sub ref-new-dihedral ref-old-dihedral) ; the delta to add to every other dihedral
                                      collect delta-dihedral)))
           (when (= (length delta-dihedrals) 2)
             (when (> (abs (radians-sub (first delta-dihedrals) (second delta-dihedrals))) 0.1)
               (error "In maybe-write-merged-internals the delta-dihedrals are not closely matched: ~s - this will happen if there is a mismatch in the stereochemistry of terminal methyl groups like in LEU and VAL" delta-dihedrals)))
           ;; If delta-dihedrals aren't all the same value then the stereochemistry might be different
           ;; from the rotamer database vs the structure we loaded.
           ;; For example, ILE, VAL or LEU may have the terminal methyl groups swapped.
           (let* ((delta-dihedral (first delta-dihedrals)))
             (loop for child-index in unset-child-indexes
                   for index3 = (* 3 child-index)
                   do (when verbose (format t "dependent setting updated-internals-mask to 1~%"))
                   do (setf
                       (aref target-internals index3) (aref new-internals index3) ; distance
                       (aref target-internals (+ 1 index3)) (aref new-internals (+ 1 index3)) ; angle
                       (aref target-internals (+ 2 index3)) (radians-sub (aref new-internals (+ 2 index3)) delta-dihedral) ; dihedral
                       (aref updated-internals-mask child-index) 1
                       ))
             (when verbose
               (format t "Setting children of ~s~%  children ~s at ~s using as reference ~s~%"
                       joint
                       (mapcar (lambda (idx) (elt joints idx)) unset-child-indexes)
                       unset-child-indexes (elt joints (first set-child-indexes))))))))))

(defun merge-internals-relative-to-heavy-atoms (atresidue to-internals rotamer-internals undefined-internals-mask heavy-atom-mask verbose)
  (loop named internals
        with temp-internals = (copy-seq to-internals)
        with updated-internals-mask = (make-array (length undefined-internals-mask) :element-type 'bit :initial-element 0)
        for joint-index below (length (joints atresidue))
        do (when (and (= (aref undefined-internals-mask joint-index) 1)
                      (= (aref heavy-atom-mask joint-index) 1) ; skips hydrogen
                      )
             (let* ((joint (aref (topology:joints atresidue) joint-index))
                    (children (topology:children joint))
                    (children-indexes (indexes-of-internals atresidue children)))
               (unless (null children)
                 (maybe-merge-internals joint temp-internals undefined-internals-mask children-indexes (topology:joints atresidue) to-internals rotamer-internals heavy-atom-mask updated-internals-mask verbose))
               #+(or)
               (progn
                 (format t "dihedral3 atom ~s in atresidue: ~s~%" joint (name atresidue))
                 (format t "  children: ~s~%" children)
                 (format t "  updated dihedrals: ~%")
                 (loop for joint across (joints atresidue)
                       for index from 0
                       for index3 from 0 by 3
                       for name = (kin:joint/name joint)
                       for dihdeg = (rad-to-deg (aref temp-internals (+ 2 index3)))
                       for old-dihdeg = (rad-to-deg (aref to-internals (+ 2 index3)))
                       do (format t "  ~5s  ~10,3f  ~10,3f~%" name dihdeg old-dihdeg)
                       ))))
        finally (return-from internals (values temp-internals updated-internals-mask))))

#+(or)
(defun overwrite-missing-internals (atresidue to-internals rotamer-internals undefined-internals-mask verbose)
  (break "Check arguments")
  (loop named internals
        with temp-internals = (copy-seq to-internals)
        with updated-internals-mask = (make-array (length undefined-internals-mask) :element-type 'bit :initial-element 0)
        for joint-index below (length (joints atresidue))
        do (when (and (= (aref undefined-internals-mask joint-index) 1)
                      (= (aref heavy-atom-mask joint-index) 1) ; skips hydrogen
                      )
             (let* ((joint (aref (topology:joints atresidue) joint-index))
                    (children (topology:children joint))
                    (children-indexes (indexes-of-internals atresidue children)))
               (unless (null children)
                 (maybe-merge-internals joint temp-internals undefined-internals-mask children-indexes (topology:joints atresidue) to-internals rotamer-internals heavy-atom-mask updated-internals-mask verbose))
               #+(or)
               (progn
                 (format t "dihedral3 atom ~s in atresidue: ~s~%" joint (name atresidue))
                 (format t "  children: ~s~%" children)
                 (format t "  updated dihedrals: ~%")
                 (loop for joint across (joints atresidue)
                       for index from 0
                       for index3 from 0 by 3
                       for name = (kin:joint/name joint)
                       for dihdeg = (rad-to-deg (aref temp-internals (+ 2 index3)))
                       for old-dihdeg = (rad-to-deg (aref to-internals (+ 2 index3)))
                       do (format t "  ~5s  ~10,3f  ~10,3f~%" name dihdeg old-dihdeg)
                       ))))
        finally (return-from internals (values temp-internals updated-internals-mask))))

(defun internals-twist-dihedral (assembler dihedral-joint angle-rad internals )
  "Twist around all dihedral angles shared by the DIHEDRAL-JOINT by setting the DIHEDRAL-JOINT dihedral angle to ANGLE-RAD"
  (let* ((parent-joint (kin:joint/parent dihedral-joint))
         (children (kin:joint/joint-children parent-joint))
         (current-dihedral (kin:bonded-joint/get-phi dihedral-joint internals))
         (adjustments (loop for child in children
                            for dih = (kin:bonded-joint/get-phi child internals)
                            for delta = (topology:radians-sub dih current-dihedral)
                            collect (cons child delta)))
         )
    (loop for (joint . dih) in adjustments
          do (kin:bonded-joint/set-phi joint internals (topology:radians-add dih angle-rad)))
    #+(or)
    (progn
      (format t "Setting ~s to ~f~%" dihedral-joint (topology:rad-to-deg angle-rad))
      (loop for child in children
            do (format t "Adjusted child ~s -> ~f~%" child (topology:rad-to-deg (kin:bonded-joint/get-phi child internals)))))))
