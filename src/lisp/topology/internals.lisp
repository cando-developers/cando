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

(cando.serialize:make-class-save-load jump-internal)

(defclass bonded-internal (internal)
  ((bond :initarg :bond :accessor bond)
   (angle :initarg :angle :accessor angle)
   (dihedral :initarg :dihedral :accessor dihedral)))

(defmethod copy-internal ((internal bonded-internal))
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

(cando.serialize:make-class-save-load complex-bonded-internal)

(defclass fragment-internals (serial:serializable)
  ((index :initarg :index :accessor index)
   (probability :initarg :probability :accessor probability)
   (internals :initarg :internals :accessor internals)
   (out-of-focus-internals :initarg :out-of-focus-internals :accessor out-of-focus-internals)))

(cando.serialize:make-class-save-load fragment-internals)

(defun copy-fragment-internals (fragment-internals)
  (make-instance 'fragment-internals
                 :index (index fragment-internals)
                 :internals (copy-seq (internals fragment-internals))
                 :out-of-focus-internals (let ((ht (make-hash-table)))
                                           (maphash (lambda (key value)
                                                      (setf (gethash key ht) (copy-seq value)))
                                                    (out-of-focus-internals fragment-internals))
                                           ht)))

(defclass fragment-conformations (serial:serializable)
  ((focus-monomer-name :initarg :focus-monomer-name :accessor focus-monomer-name)
   (monomer-context :initarg :monomer-context :accessor monomer-context)
   (total-count :initform 0 :initarg :total-count :accessor total-count)
   (fragments :initform (make-array 16 :adjustable t :fill-pointer 0) :initarg :fragments :accessor fragments)))

(cando.serialize:make-class-save-load fragment-conformations)

(defclass fragment-conformations-map (serial:serializable)
  ((monomer-context-to-fragment-conformations :initform (make-hash-table :test 'equal)
                                              :initarg :monomer-context-to-fragment-conformations
                                              :accessor monomer-context-to-fragment-conformations)))

(cando.serialize:make-class-save-load fragment-conformations-map
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t))))

(defclass matched-fragment-conformations-map (fragment-conformations-map)
  ((fragment-matches :initform (make-hash-table :test 'equal)
                     :initarg :fragment-matches
                     :accessor fragment-matches)))

(cando.serialize:make-class-save-load matched-fragment-conformations-map)

(defun matched-fragment-conformations-summary (matched-fragment-conformations-map)
  (let ((total-fragment-conformations 0)
        (matching-fragment-conformations 0)
        (missing-fragment-conformations 0))
    (maphash (lambda (key value)
               (declare (ignore key))
               (incf total-fragment-conformations (length (fragments value))))
             (monomer-context-to-fragment-conformations matched-fragment-conformations-map))
    (maphash (lambda (key value)
               (declare (ignore key))
               (incf matching-fragment-conformations (length value)))
             (fragment-matches matched-fragment-conformations-map))
    (maphash (lambda (key value)
               (declare (ignorable key))
               (loop for val across value
                     when (= (length val) 0)
                       do (incf missing-fragment-conformations)))
             (fragment-matches matched-fragment-conformations-map))
    (let ((missing-monomer-contexts nil))
      (maphash (lambda (monomer-context fragment-conformations)
                 (block inner-search
                   (maphash (lambda (monomer-context-pair allowed-fragment-indices)
                              (when (or (string= (car monomer-context-pair) monomer-context)
                                        (string= (cdr monomer-context-pair) monomer-context))
                                (return-from inner-search nil)))
                            (fragment-matches matched-fragment-conformations-map))
                   (push monomer-context missing-monomer-contexts)))
               (monomer-context-to-fragment-conformations matched-fragment-conformations-map))
      (values total-fragment-conformations matching-fragment-conformations missing-fragment-conformations missing-monomer-contexts))))

(defconstant +dihedral-threshold+ (* 10.0 0.0174533))

(defun similar-internals-p (frag1 frag2 &optional )
  (loop for frag1-int across (internals frag1)
        for frag2-int across (internals frag2)
        do (when (and (typep frag1-int 'bonded-internal)
                      (typep frag2-int 'bonded-internal))
             (let* ((aa (- (dihedral frag1-int) (dihedral frag2-int)))
                    (aamod (- (mod (+ aa 180) 360) 180)))
               (when (> (abs aamod) +dihedral-threshold+)
                 (return-from similar-internals-p nil)))))
  t)

(defun seen-fragment-internals (fragment-conformations fragment-internals)
  (loop for seen-frag across (fragments fragment-conformations)
        when (similar-internals-p seen-frag fragment-internals)
          do (return-from seen-fragment-internals (index seen-frag)))
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

(defun save-fragment-conformations (fragment-conformations filename)
  (cando.serialize:save-cando fragment-conformations filename))

(defun load-fragment-conformations (filename)
  (cando.serialize:load-cando filename))

(defun dump-fragment-internals (fragment-internals finternals)
  (format finternals "begin-conformation ~a~%" (index fragment-internals))
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

(defgeneric fill-joint-internals (joint internal))

(defmethod fill-joint-internals ((joint kin:jump-joint) (internal jump-internal))
  )


(defmethod fill-joint-internals ((joint kin:bonded-joint) (internal bonded-internal))
  (kin:set-distance joint (bond internal))
  (kin:set-theta joint (angle internal))
  (kin:set-phi joint (dihedral internal))
  )
