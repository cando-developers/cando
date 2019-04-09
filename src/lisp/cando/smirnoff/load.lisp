(in-package :smirnoff)

(defvar *type-index* (ext:make-atomic 0))

(defun next-smirnoff-type-symbol ()
  (let* ((index (loop named get-index
                   for index = (ext:atomic-get *type-index*)
                   for next-index = (1+ index)
                   for changed = (ext:atomic-compare-and-swap-weak *type-index* index next-index)
                   when changed
                     do (return-from get-index index)))
         (index-string (format nil "$~36r" index)))
    (unless (<= (length index-string) 4)
      (error "We have exhausted the numbe of Smirnoff atom types that are available"))
    (intern index-string :keyword)))

(defvar *smirnoff-types* (make-hash-table :test #'equal :thread-safe t))

(defvar *smirnoff*)


(defclass force-term ()
  ((smirks :initarg :smirks :accessor smirks)
   (compiled-smirks :initarg :compiled-smirks :accessor compiled-smirks)
   (id :initarg :id :accessor id)))

(defmethod print-object ((object force-term) stream)
  (print-unreadable-object (object stream)
    (format stream "~a" (class-name (class-of object)))
    (format stream " ~s ~s" :smirks (smirks object))
    (format stream " ~s ~s" :id (id object))))

(defclass bond-term (force-term)
  ((k :initarg :k :accessor k)
   (len :initarg :len :accessor len)))

(defclass angle-term (force-term)
  ((k :initarg :k :accessor k)
   (angle-rad :initarg :angle-rad :accessor angle-rad)))

(defclass harmonic-bond-force ()
  ((k-unit :initarg :k-unit :accessor k-unit)
   (length-unit :initarg :length-unit :accessor length-unit)
   (terms :initform (make-array 16 :fill-pointer 0 :adjustable t) :initarg :terms :accessor terms)))

(defclass harmonic-angle-force ()
  ((k-unit :initarg :k-unit :accessor k-unit)
   (angle-unit :initarg :angle-unit :accessor angle-unit)
   (terms :initform (make-array 16 :fill-pointer 0 :adjustable t) :initarg :terms :accessor terms)))

(defclass torsion-part ()
  ((idivf :initarg :idivf :accessor idivf)
   (k :initarg :k :accessor k)
   (periodicity :initarg :periodicity :accessor periodicity)
   (phase-rad :initarg :phase-rad :accessor phase-rad)))

(defclass proper-term (force-term)
  ((parts :initform nil :initarg :parts :accessor parts)))

(defclass improper-term (force-term)
  ((parts :initform nil :initarg :parts :accessor parts)))

(defclass periodic-torsion-force ()
  ((k-unit :initarg :k-unit :accessor k-unit)
   (phase-unit :initarg :phase-unit :accessor phase-unit)
   (terms :initform (make-array 16 :fill-pointer 0 :adjustable t) :initarg :terms :accessor terms)))
  
(defclass atom-term (force-term)
  ((epsilon :initarg :epsilon :accessor epsilon)
   (rmin-half :initarg :rmin-half :accessor rmin-half)))

(defclass nonbond-term (force-term)
  ((type :initarg :type :accessor type)
   (epsilon :initarg :epsilon :accessor epsilon)
   (rmin-half :initarg :rmin-half :accessor rmin-half)))
  
(defclass nonbonded-force ()
  ((coulomb14scale :initarg :coulomb14scale :accessor coulomb14scale)
   (lj14scale :initarg :lj14scale :accessor lj14scale)
   (sigma-unit :initarg :sigma-unit :accessor sigma-unit)
   (epsilon-unit :initarg :epsilon-unit :accessor epsilon-unit)
   (terms :initform (make-array 16 :fill-pointer 0 :adjustable t) :initarg :terms :accessor terms)))

(defmethod print-object ((object nonbond-term) stream)
  (print-unreadable-object (object stream)
    (format stream "~a ~s ~s" (class-name (class-of object)) :type (type object))
    (format stream " ~s ~s" :smirks (smirks object))
    (format stream " ~s ~s" :id (id object))))

(defclass smirnoff-force-field ()
  ((harmonic-bond-force :initarg :harmonic-bond-force :accessor harmonic-bond-force)
   (harmonic-angle-force :initarg :harmonic-angle-force :accessor harmonic-angle-force)
   (periodic-torsion-force :initarg :periodic-torsion-force :accessor periodic-torsion-force)
   (nonbonded-force :initarg :nonbonded-force :accessor nonbonded-force)))

(defclass combined-force-field ()
  ((parts :initform nil :initarg :parts :accessor parts)))

(defmethod chem:force-fields-as-list ((combined-force-field combined-force-field))
  (parts combined-force-field))


(defun ignore-handler (node)
  (declare (ignore node)))

(defmacro with-force-parser ((node root child-tag) &body body )
  (let ((tag (gensym)))
    `(loop for ,node across (plump:children ,root)
           when (typep ,node 'plump:element)
             do (let ((,tag (plump:tag-name node)))
                  (when (string= ,tag ,child-tag)
                    ,@body)))))
  
(defun parse-harmonic-bond-force (root)
  (let* ((attrs (plump:attributes root))
         (length-unit (gethash "length_unit" attrs))
         (k-unit (gethash "k_unit" attrs)))
    (setf (harmonic-bond-force *smirnoff*)
          (make-instance 'harmonic-bond-force
                         :length-unit length-unit
                         :k-unit k-unit)))
  (with-force-parser (node root "Bond")
    (let* ((attrs (plump:attributes node))
           (smirks (gethash "smirks" attrs))
           (id (gethash "id" attrs))
           (k (fortran:parse-double-float (gethash "k" attrs)))
           (len (fortran:parse-double-float (gethash "length" attrs)))
           (term (make-instance 'bond-term
                                :smirks smirks
                                :compiled-smirks (chem:compile-smarts smirks)
                                :id id
                                :k k
                                :len len)))
      (vector-push-extend term (terms (harmonic-bond-force *smirnoff*))))))

(defun parse-harmonic-angle-force (root)
  (let* ((attrs (plump:attributes root))
         (angle-unit (gethash "angle_unit" attrs))
         (angle-to-rad-multiplier (cond
                                    ((string= "degrees" angle-unit)
                                     0.0174533)
                                    (t (error "Add support for angle-unit ~a" angle-unit))))
         (k-unit (gethash "k_unit" attrs)))
    (setf (harmonic-angle-force *smirnoff*)
          (make-instance 'harmonic-angle-force
                         :angle-unit angle-unit
                         :k-unit k-unit))
    (with-force-parser (node root "Angle")
      (let* ((attrs (plump:attributes node))
             (smirks (gethash "smirks" attrs))
             (id (gethash "id" attrs))
             (k (fortran:parse-double-float (gethash "k" attrs)))
             (angle (fortran:parse-double-float (gethash "angle" attrs)))
             (angle-radians (* angle angle-to-rad-multiplier))
             (term (make-instance 'angle-term
                                  :smirks smirks
                                  :compiled-smirks (chem:compile-smarts smirks)
                                  :id id
                                  :k k
                                  :angle-rad angle-radians)))
        (vector-push-extend term (terms (harmonic-angle-force *smirnoff*)))))))

(defun parse-torsion (node torsion-class-name attributes phase-to-rad-multiplier)
  (let* ((smirks (gethash "smirks" attributes))
         (id (gethash "id" attributes))
         (term (make-instance torsion-class-name
                              :smirks smirks
                              :compiled-smirks (chem:compile-smarts smirks)
                              :id id)))
    (block terms
      (loop for index from 1
            for idivf-name = (format nil "idivf~a" index)
            for idivf-str = (gethash idivf-name attributes)
            if idivf-str
              do (let ((idivf (fortran:parse-double-float idivf-str))
                       (k (fortran:parse-double-float (gethash (format nil "k~a" index) attributes)))
                       (periodicity (parse-integer (gethash (format nil "periodicity~a" index) attributes)))
                       (phase (fortran:parse-double-float (gethash (format nil "phase~a" index) attributes))))
                   (push (make-instance 'torsion-part
                                        :idivf idivf
                                        :k k
                                        :periodicity periodicity
                                        :phase-rad (* phase phase-to-rad-multiplier))
                         (parts term)))
            else
              do (return-from terms)))
    term))

(defun parse-periodic-torsion-force (root)
  (let* ((attrs (plump:attributes root))
         (phase-unit (gethash "phase_unit" attrs))
         (phase-to-rad-multiplier (cond
                                    ((string= "degrees" phase-unit)
                                     0.0174533)
                                    (t (error "Add support for phase-unit ~a" phase-unit))))
         (k-unit (gethash "k_unit" attrs)))
    (setf (periodic-torsion-force *smirnoff*)
          (make-instance 'periodic-torsion-force
                         :phase-unit phase-unit
                         :k-unit k-unit))
    (loop for node across (plump:children root)
          when (typep node 'plump:element)
            do (let ((tag (plump:tag-name node)))
                 (cond
                   ((string= tag "Proper")
                    (vector-push-extend
                     (parse-torsion node 'proper-term (plump:attributes node) phase-to-rad-multiplier)
                     (terms (periodic-torsion-force *smirnoff*))))
                   ((string= tag "Improper")
                    (vector-push-extend
                     (parse-torsion node 'improper-term (plump:attributes node) phase-to-rad-multiplier)
                     (terms (periodic-torsion-force *smirnoff*))))
                   (t (error "Illegal tag ~a" tag)))))))

(defun parse-nonbonded-force (root)
  (let* ((attrs (plump:attributes root))
         (coulomb14scale (fortran:parse-double-float (gethash "coulomb14scale" attrs)))
         (lj14scale (fortran:parse-double-float (gethash "lj14scale" attrs)))
         (sigma-unit (gethash "sigma_unit" attrs))
         (epsilon-unit (gethash "epsilon_unit" attrs)))
    (setf (nonbonded-force *smirnoff*)
          (make-instance 'nonbonded-force
                         :coulomb14scale coulomb14scale
                         :lj14scale lj14scale
                         :sigma-unit sigma-unit
                         :epsilon-unit epsilon-unit))
    (with-force-parser (node root "Atom")
      (let* ((attrs (plump:attributes node))
             (smirks (gethash "smirks" attrs))
             (id (gethash "id" attrs))
             (epsilon (fortran:parse-double-float (gethash "epsilon" attrs)))
             (rmin-half (fortran:parse-double-float (let (num)
                                                  (cond
                                                    ((setf num (gethash "rmin_half" attrs))
                                                     num)
                                                    ((setf num (gethash "sigma" attrs))
                                                     ;; handle sigma
                                                     (error "Handle sigma"))
                                                    (t (error "Neither rmin_half or sigma were provided"))))))
             (type (let ((exists (gethash smirks *smirnoff-types*)))
                     (if exists
                         exists
                         (let ((type (next-smirnoff-type-symbol)))
                           (setf (gethash smirks *smirnoff-types*) type)
                           type))))
             (term (make-instance 'nonbond-term
                                  :type type
                                  :smirks smirks
                                  :compiled-smirks (chem:compile-smarts smirks)
                                  :id id
                                  :epsilon epsilon
                                  :rmin-half rmin-half)))
        (vector-push-extend term (terms (nonbonded-force *smirnoff*)))))))

(defparameter *smirnoff-handlers*
  '(("HarmonicBondForce" . parse-harmonic-bond-force)
    ("PeriodicTorsionForce" . parse-periodic-torsion-force)
    ("HarmonicAngleForce" . parse-harmonic-angle-force)
    ("NonbondedForce" . parse-nonbonded-force)
    ("Date" . ignore-handler)
    ("Author" . ignore-handler)
    ("Smirnoff" . parse-dispatcher )))

(defun parse-dispatcher (root)
  (loop for node across (plump:children root)
        when (typep node 'plump:element)
          do (let* ((tag (plump:tag-name node))
                    (found (assoc tag *smirnoff-handlers* :test #'string=)))
               (when found
                 (let ((handler (cdr found)))
                   (funcall handler node))))))

(defun parse-root (root)
  (loop for node across (plump:children root)
        when (typep node 'plump:element)
          do (parse-dispatcher node)))

(defun load-smirnoff (path)
  (unless (typep path 'pathname)
    (error "Expected pathname - got ~s" path))
  (let ((root (plump:parse path))
        (*smirnoff* (make-instance 'smirnoff-force-field)))
    #+(or)(plump:traverse root #'smirnoff-node)
    (parse-root root)
    *smirnoff*))

(defmethod chem:nonbond-component ((force-field smirnoff-force-field))
  (nonbonded-force force-field))

(defmethod chem:clear-combined-force-field ((combined-force-field combined-force-field))
  (setf (parts combined-force-field) nil))

(defmethod chem:add-shadowing-force-field ((combined-force-field combined-force-field) (smirnoff-force-field smirnoff-force-field) pathname)
  (push smirnoff-force-field (parts combined-force-field)))
