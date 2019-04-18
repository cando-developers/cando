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

(defclass bonds-force ()
  ((k-unit :initarg :k-unit :accessor k-unit)
   (length-unit :initarg :length-unit :accessor length-unit)
   (terms :initform (make-array 16 :fill-pointer 0 :adjustable t) :initarg :terms :accessor terms)))

(defclass angles-force ()
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

(defclass proper-torsion-force ()
  ((k-unit :initarg :k-unit :accessor k-unit)
   (phase-unit :initarg :phase-unit :accessor phase-unit)
   (terms :initform (make-array 16 :fill-pointer 0 :adjustable t) :initarg :terms :accessor terms)))

(defclass improper-torsion-force ()
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
  
(defclass vdw-force ()
  ((potential :initarg :potential :accessor potential)
   (combining-rules :initarg :combining-rules :accessor combining-rules)
   (scale12 :initarg :scale12 :accessor scale12)
   (scale13 :initarg :scale13 :accessor scale13)
   (scale14 :initarg :scale14 :accessor scale14)
   (scale15 :initarg :scale15 :accessor scale15)
   (rmin-half-unit :initarg :rmin-half-unit :accessor rmin-half-unit)
   (epsilon-unit :initarg :epsilon-unit :accessor epsilon-unit)
   (switch-width :initarg :switch-width :accessor switch-width)
   (switch-width-unit :initarg :switch-width-unit :accessor switch-width-unit)
   (cutoff :initarg :cutoff :accessor cutoff)
   (cutoff-unit :initarg :cutoff-unit :accessor cutoff-unit)
   (method :initarg :method :accessor method)
   (terms :initform (make-array 16 :fill-pointer 0 :adjustable t) :initarg :terms :accessor terms)))

(defmethod print-object ((object nonbond-term) stream)
  (print-unreadable-object (object stream)
    (format stream "~a ~s ~s" (class-name (class-of object)) :type (type object))
    (format stream " ~s ~s" :smirks (smirks object))
    (format stream " ~s ~s" :id (id object))))

(defclass smirnoff-force-field ()
  ((bonds-force :initarg :bonds-force :accessor bonds-force)
   (angles-force :initarg :angles-force :accessor angles-force)
   (proper-torsion-force :initarg :proper-torsion-force :accessor proper-torsion-force)
   (improper-torsion-force :initarg :improper-torsion-force :accessor improper-torsion-force)
   (vdw-force :initarg :vdw-force :accessor vdw-force)))

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
  
(defun parse-bonds-force (root)
  (let* ((attrs (plump:attributes root))
         (length-unit (gethash "length_unit" attrs))
         (k-unit (gethash "k_unit" attrs)))
    (setf (bonds-force *smirnoff*)
          (make-instance 'bonds-force
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
      (vector-push-extend term (terms (bonds-force *smirnoff*))))))

(defun parse-angles-force (root)
  (let* ((attrs (plump:attributes root))
         (angle-unit (gethash "angle_unit" attrs))
         (angle-to-rad-multiplier (cond
                                    ((string= "degrees" angle-unit)
                                     0.0174533)
                                    (t (error "Add support for angle-unit ~a" angle-unit))))
         (k-unit (gethash "k_unit" attrs)))
    (setf (angles-force *smirnoff*)
          (make-instance 'angles-force
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
        (vector-push-extend term (terms (angles-force *smirnoff*)))))))

(defun parse-torsion (node torsion-class-name attributes phase-to-rad-multiplier)
  (let* ((smirks (gethash "smirks" attributes))
         (id (gethash "id" attributes))
         (term (make-instance torsion-class-name
                              :smirks smirks
                              :compiled-smirks (chem:compile-smarts smirks)
                              :id id)))
    (block terms
      (loop for index from 1
            for k-name = (format nil "k~a" index)
            for k-str = (gethash k-name attributes)
            if k-str
              do (let ((idivf (let ((idivf-str (gethash (format nil "idivf~a" index) attributes)))
                                (if idivf-str
                                    (fortran:parse-double-float idivf-str)
                                    1.0d0)))
                       (k (fortran:parse-double-float k-str))
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

(defun parse-torsion-force (root name term-name force-name)
  (let* ((attrs (plump:attributes root))
         (phase-unit (gethash "phase_unit" attrs))
         (phase-to-rad-multiplier (cond
                                    ((string= "degrees" phase-unit)
                                     0.0174533)
                                    (t (error "Add support for phase-unit ~a" phase-unit))))
         (k-unit (gethash "k_unit" attrs)))
    (let ((force (make-instance force-name
                                :phase-unit phase-unit
                                :k-unit k-unit)))
      (loop for node across (plump:children root)
            when (typep node 'plump:element)
              do (let ((tag (plump:tag-name node)))
                   (cond
                     ((string= tag name)
                      (vector-push-extend
                       (parse-torsion node term-name (plump:attributes node) phase-to-rad-multiplier)
                       (terms force )))
                     (t (error "Illegal tag ~a" tag)))))
      force)))


(defun parse-proper-torsion-force (root)
  (setf (proper-torsion-force *smirnoff*)
        (parse-torsion-force root "Proper" 'proper-term 'proper-torsion-force)))

(defun parse-improper-torsion-force (root)
  (setf (improper-torsion-force *smirnoff*)
        (parse-torsion-force root "Improper" 'improper-term 'improper-torsion-force)))
  
(defun parse-vdw-force (root)
  (let* ((attrs (plump:attributes root))
         (potential (gethash "potential" attrs))
         (combining-rules (gethash "combining_rules" attrs))
         (scale12 (fortran:parse-double-float (gethash "scale12" attrs)))
         (scale13 (fortran:parse-double-float (gethash "scale13" attrs)))
         (scale14 (fortran:parse-double-float (gethash "scale14" attrs)))
         (scale15 (fortran:parse-double-float (gethash "scale15" attrs)))
         (rmin-half-unit (gethash "rmin_half_unit" attrs))
         (epsilon-unit (gethash "epsilon_unit" attrs))
         (switch-width (fortran:parse-double-float (gethash "switch_width" attrs)))
         (switch-width-unit (gethash "switch_width_unit" attrs))
         (cutoff (fortran:parse-double-float (gethash "cutoff" attrs)))
         (cutoff-unit (gethash "cutoff_unit" attrs))
         (method (gethash "method" attrs))
         )
    (setf (vdw-force *smirnoff*)
          (make-instance 'vdw-force
                         :potential potential
                         :combining-rules combining-rules
                         :scale12 scale12
                         :scale13 scale13
                         :scale14 scale14
                         :scale15 scale15
                         :rmin-half-unit rmin-half-unit
                         :epsilon-unit epsilon-unit
                         :switch-width switch-width
                         :switch-width-unit switch-width-unit
                         :cutoff cutoff
                         :cutoff-unit cutoff-unit
                         :method method))
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
        (vector-push-extend term (terms (vdw-force *smirnoff*)))))))

(defparameter *smirnoff-handlers*
  '(("Bonds" . parse-bonds-force)
    ("ProperTorsions" . parse-proper-torsion-force)
    ("ImproperTorsions" . parse-improper-torsion-force)
    ("Angles" . parse-angles-force)
    ("vdW" . parse-vdw-force)
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
  (vdw-force force-field))

(defmethod chem:clear-combined-force-field ((combined-force-field combined-force-field))
  (setf (parts combined-force-field) nil))

(defmethod chem:add-shadowing-force-field ((combined-force-field combined-force-field) (smirnoff-force-field smirnoff-force-field) pathname)
  (push smirnoff-force-field (parts combined-force-field)))
