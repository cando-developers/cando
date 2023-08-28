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
   (angle :initarg :angle :accessor angle)
   (dihedral :initarg :dihedral :accessor dihedral)))

(defmethod copy-internal ((internal bonded-internal))
  (when (ext:float-nan-p (dihedral internal))
    (error "Found a NAN dihedral"))
  (make-instance 'bonded-internal
                 :name (name internal)
                 :bond (bond internal)
                 :angle (angle internal)
                 :dihedral (dihedral internal)))

(cando.serialize:make-class-save-load
 bonded-internal
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a :d ~5,2f" (name obj) (rad-to-deg (dihedral obj))))))

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
   (internals :initarg :internals :accessor internals)
   (energy :initform 0.0s0 :initarg :energy :accessor energy)
   (out-of-focus-internals :initarg :out-of-focus-internals :accessor out-of-focus-internals)
   (coordinates :initarg :coordinates :accessor coordinates)
   ))

(cando.serialize:make-class-save-load
 fragment-internals
 :print-unreadably
 (lambda (obj stream)
   (let ((*print-pretty* nil))
   (print-unreadable-object (obj stream :type t)
     (let ((internals (coerce (subseq (internals obj) 0 (min (length (internals obj)) 4) ) 'list)))
       (format stream ":trainer-index ~a :energy ~10,3f first internals: ~{(~{~s ~6,2f~}) ~}" (trainer-index obj) (energy obj) (mapcar (lambda (x) (list (name x) (rad-to-deg (dihedral x)))) internals )))))))

(defun copy-fragment-internals (fragment-internals)
  (make-instance 'fragment-internals
                 :trainer-index (trainer-index fragment-internals)
                 :internals (copy-seq (internals fragment-internals))
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
                   (return (list (topology:dihedral internal))))
                  (t (return (list 0.0))))
           finally (error "Could not find internal with atom-name ~a" name-or-plug-names)))
    ((consp name-or-plug-names)
     (let* ((plug-name (first name-or-plug-names))
            (out-of-focus-internals (gethash plug-name (out-of-focus-internals fragment-internals))))
       (loop for name in (cdr name-or-plug-names)
             collect (loop for internal across out-of-focus-internals
                           when (eq (name internal) name)
                             do (return (topology:dihedral internal))))))
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


(defclass fragment-conformations (serial:serializable)
  ((focus-monomer-name :initarg :focus-monomer-name :accessor focus-monomer-name)
   (monomer-context :initarg :monomer-context :accessor monomer-context)
   (next-index :initform 0 :initarg :next-index :accessor next-index)
   (fragments :initform (make-array 16 :adjustable t :fill-pointer 0) :initarg :fragments :accessor fragments)))

(cando.serialize:make-class-save-load
 fragment-conformations
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~s number-of-fragments ~a" (monomer-context obj) (length (fragments obj))))))

(defmethod fragment-conformations ((obj fragment-conformations))
  obj)

(defun merge-fragment-conformations (fragment-conformations-list)
  (unless (> (length fragment-conformations-list) 0)
    (error "There must be at least one fragment-conformations in list"))
  (let ((unsorted-fragment-internals (make-array 1024 :adjustable t :fill-pointer 0))
        (max-next-index 0))
    (loop for one in fragment-conformations-list
          do (setf max-next-index (max max-next-index (next-index one)))
          do (loop for frag-int across (fragments one)
                   do (vector-push-extend frag-int unsorted-fragment-internals)))
    (let ((new-fragment-internals (sort unsorted-fragment-internals #'< :key #'trainer-index)))
      (make-instance 'fragment-conformations
                     :focus-monomer-name (focus-monomer-name (first fragment-conformations-list))
                     :monomer-context (monomer-context (first fragment-conformations-list))
                     :next-index max-next-index
                     :fragments new-fragment-internals))))

(defclass clusterable-fragment-conformations (serial:serializable)
  ((cluster-dihedral-names :initarg :cluster-dihedral-names :accessor cluster-dihedral-names)
   (fragment-conformations :initarg :fragment-conformations :accessor fragment-conformations)))

(cando.serialize:make-class-save-load clusterable-fragment-conformations)

(defclass fragment-conformations-map (serial:serializable)
  ((monomer-context-to-fragment-conformations :initform (make-hash-table :test 'equal)
                                              :initarg :monomer-context-to-fragment-conformations
                                              :accessor monomer-context-to-fragment-conformations)))

(cando.serialize:make-class-save-load fragment-conformations-map
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t))))

(defclass fragment-shape-connections ()
  ((fragment-shape-connections-map :initform (make-hash-table :test 'equal)
        :accessor fragment-shape-connections-map)))

(defun make-fragment-shape-connections ()
  (make-instance 'fragment-shape-connections))

(defclass fragment-context-connections ()
  ((fmap :initform (make-hash-table :test 'eq)
        :accessor fmap)))

(defun make-fragment-context-connections (&optional monomer-context-to-fragment-conformations)
  (make-instance 'fragment-context-connections))

(defun map-fragment-context-connections (op fcc)
  (maphash (lambda (key value)
             (let ((from (car key))
                   (to (cdr key)))
               (funcall op from to value)))
           (fmap fcc)))

(defun fragment-context-connections-count (fcc)
  (hash-table-count (fmap fcc)))

(defun set-fragment-context-connections (fcc from to value)
  (check-type value fragment-shape-connections)
  (setf (gethash (cons from to) (fmap fcc)) value))

(defun get-fragment-context-connections (fcc from to)
  (let* ((key (cons from to))
         (res (gethash (cons from to) (fmap fcc))))
    (if res
        res
        (let ((fsc (make-fragment-shape-connections)))
          (setf (gethash key (fmap fcc)) fsc)))))

(defclass matched-fragment-conformations-map (fragment-conformations-map)
  ((fragment-context-connections :initform (make-fragment-context-connections)
                     :accessor fragment-context-connections)))

(cando.serialize:make-class-save-load
 matched-fragment-conformations-map
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream " fragment-context-connections-count ~d length-fragments ~d" (fragment-context-connections-count (fragment-context-connections obj))
             (hash-table-count (monomer-context-to-fragment-conformations obj))))))

(defun matched-fragment-conformations-summary (matched-fragment-conformations-map)
  (let ((total-fragment-conformations 0)
        (matching-fragment-conformations 0)
        (missing-fragment-conformations 0))
    (maphash (lambda (key value)
               (declare (ignore key))
               (incf total-fragment-conformations (length (fragments value))))
             (monomer-context-to-fragment-conformations matched-fragment-conformations-map))
    (map-fragment-context-connections
     (lambda (from-monomer-context to-monomer-context value)
       (declare (ignore from-monomer-context to-monomer-context))
       (incf matching-fragment-conformations (length value)))
     (fragment-context-connections matched-fragment-conformations-map))
    (map-fragment-context-connections
     (lambda (from-monomer-context to-monomer-context value)
       (declare (ignore from-monomer-context to-monomer-context))
       (loop for val across value
             when (= (length val) 0)
               do (incf missing-fragment-conformations)))
     (fragment-context-connections matched-fragment-conformations-map))
    (let ((missing-monomer-contexts nil))
      (maphash (lambda (monomer-context fragment-conformations)
                 (declare (ignore fragment-conformations))
                 (block inner-search
                   (map-fragment-context-connections
                    (lambda (from-monomer-context to-monomer-context allowed-fragment-indices)
                      (declare (ignore allowed-fragment-indices))
                      (when (or (string= from-monomer-context monomer-context)
                                (string= to-monomer-context monomer-context))
                        (return-from inner-search nil)))
                    (fragment-context-connections matched-fragment-conformations-map))
                   (push monomer-context missing-monomer-contexts)))
               (monomer-context-to-fragment-conformations matched-fragment-conformations-map))
      (values total-fragment-conformations matching-fragment-conformations missing-fragment-conformations missing-monomer-contexts))))

(defparameter *dihedral-threshold* 20.0)
(defparameter *dihedral-rms-threshold* 15.0)

(defun internals-are-different-p (frag1 frag2 &optional )
  (loop for frag1-int across (internals frag1)
        for frag2-int across (internals frag2)
        do (when (and (typep frag1-int 'bonded-internal)
                      (typep frag2-int 'bonded-internal))
             (let* ((dihedral-frag1 (rad-to-deg (dihedral frag1-int)))
                    (dihedral-frag2 (rad-to-deg (dihedral frag2-int)))
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
             (let* ((dihedral-frag1 (rad-to-deg (dihedral frag1-int)))
                    (dihedral-frag2 (rad-to-deg (dihedral frag2-int)))
                    (aa (degree-difference dihedral-frag1 dihedral-frag2)))
               (incf sum-of-squares (* aa aa))
               (incf num-squares))))
    (let ((rms (sqrt (/ sum-of-squares (float num-squares)))))
      rms)))

(defun rms-internals-are-different-p (frag1 frag2)
  (> (dihedrals-rms frag1 frag2) *dihedral-rms-threshold*))

(defun seen-fragment-internals (fragment-conformations fragment-internals)
  (loop for seen-frag across (fragments fragment-conformations)
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

(defun save-clusterable-fragment-conformations (clusterable-fragment-conformations filename)
  (cando.serialize:save-cando clusterable-fragment-conformations filename))

(defun load-clusterable-fragment-conformations (filename)
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
                        (to-deg (topology:angle internal))
                        (to-deg (topology:dihedral internal))))
               ((typep internal 'topology:bonded-internal)
                (format finternals "bonded-joint ~a ~8,3f ~8,3f ~8,3f~%"
                        (topology:name internal)
                        (topology:bond internal)
                        (to-deg (topology:angle internal))
                        (to-deg (topology:dihedral internal))))
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
                                     (to-deg (topology:angle internal))
                                     (to-deg (topology:dihedral internal))))
                            ((typep internal 'topology:bonded-internal)
                             (format finternals "bonded-joint ~a ~8,3f ~8,3f ~8,3f~%"
                                     (topology:name internal)
                                     (topology:bond internal)
                                     (to-deg (topology:angle internal))
                                     (to-deg (topology:dihedral internal))))
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
  (values (bond internal) (angle internal) (dihedral internal)))

(defgeneric apply-fragment-internals-to-atresidue (fragment-internals atresidue))

(defmethod apply-fragment-internals-to-atresidue ((fragment-internals fragment-internals) atresidue)
 (loop for joint across (joints atresidue)
       for internal across (internals fragment-internals)
       do (multiple-value-bind (bond angle-rad dihedral-rad)
              (extract-bond-angle-rad-dihedral-rad internal)
            (fill-joint-internals joint bond angle-rad dihedral-rad))))

