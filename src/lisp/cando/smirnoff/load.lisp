(in-package :smirnoff)

(defvar *type-index* 0)

(defun next-smirnoff-type-symbol ()
  (let* ((index (1- (mp:atomic-incf *type-index*)))
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
  (;; Force constant, in kilocalories/(mole*angstrom*angstrom)
   (k :initarg :k :accessor k)
   ;; Equilibrium bond length, in angstroms
   (len :initarg :len :accessor len)))

(defclass angle-term (force-term)
  (;; Force constant in kilocalories/(mole*radian*radian)
   (k :initarg :k :accessor k)
   ;; Equilibrium angle in radians
   (angle :initarg :angle :accessor angle)))

(defclass bonds-force ()
  ((terms :initform (make-array 16 :fill-pointer 0 :adjustable t) :initarg :terms :accessor terms)))

(defclass angles-force ()
  ((terms :initform (make-array 16 :fill-pointer 0 :adjustable t) :initarg :terms :accessor terms)))

(defclass torsion-part ()
  ((idivf :initarg :idivf :accessor idivf)
   ;; Force constant in kilocalories/mol
   (k :initarg :k :accessor k)
   ;; Natural number
   (periodicity :initarg :periodicity :accessor periodicity)
   ;; In radians
   ;; NOTE: PHASE is a CL function, so we can't use that name
   (phase :initarg :phase :accessor phase-angle)))

(defclass proper-term (force-term)
  ((parts :initform nil :initarg :parts :accessor parts)))

(defclass improper-term (force-term)
  ((parts :initform nil :initarg :parts :accessor parts)))

(defclass proper-torsion-force ()
  ((terms :initform (make-array 16 :fill-pointer 0 :adjustable t) :initarg :terms :accessor terms)))

(defclass improper-torsion-force ()
  ((terms :initform (make-array 16 :fill-pointer 0 :adjustable t) :initarg :terms :accessor terms)))

(defclass nonbond-term (force-term)
  ((type :initarg :type :accessor type)
   ;; in kj/mol
   (epsilon :initarg :epsilon :accessor epsilon)
   ;; in angstroms
   (rmin-half :initarg :rmin-half :accessor rmin-half)))
  
(defclass vdw-force ()
  ((potential :initarg :potential :accessor potential)
   (combining-rules :initarg :combining-rules :accessor combining-rules)
   (scale12 :initarg :scale12 :accessor scale12)
   (scale13 :initarg :scale13 :accessor scale13)
   (scale14 :initarg :scale14 :accessor scale14)
   (scale15 :initarg :scale15 :accessor scale15)
   ;; in angstroms
   (switch-width :initarg :switch-width :accessor switch-width)
   ;; in radians
   (cutoff :initarg :cutoff :accessor cutoff)
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

(defclass combined-smirnoff-force-field ()
  ((parts :initform nil :initarg :parts :accessor parts)))

(defmethod chem:force-fields-as-list ((combined-smirnoff-force-field combined-smirnoff-force-field))
  (parts combined-smirnoff-force-field))


;; Given a units:quantity, return a number in canonical units
;; (kcal/(mol*angstrom*angstrom))
(defun canonical-force (quantity)
  (units:value-in-unit quantity
                       (load-time-value
                        (units:/ units:kilocalories-per-mole
                                 (units:* units:angstroms units:angstroms)))
                       1))

;; Ditto for lengths
(defun canonical-length (quantity)
  (units:value-in-unit quantity units:angstroms 1))

;; For rmin-half
(defun nanometer-length (quantity)
  (units:value-in-unit quantity units:nanometers 1))

(defun canonical-angle (quantity)
  (units:value-in-unit quantity units:radians 1))

(defun canonical-angle-force (quantity)
  (units:value-in-unit quantity
                       (load-time-value
                        (units:/ units:kilocalories-per-mole
                                 (units:* units:radians units:radians)))
                       1))

(defun canonical-kj/mol (quantity)
  (units:value-in-unit quantity units:kilojoules-per-mole 1))

(defun canonical-kcal/mol (quantity)
  (units:value-in-unit quantity units:kilocalories-per-mole 1))

(defun ignore-handler (node)
  (declare (ignore node)))

(defmacro with-force-parser ((node root child-tag) &body body )
  `(loop for ,node across (plump:children ,root)
         when (and (typep ,node 'plump:element)
                   (string= (plump:tag-name ,node) ,child-tag))
           do (progn ,@body)))

(defun parse-bonds-force (root)
  (let* (#+(or)(attrs (plump:attributes root))
           ;; We'll use the attrs later.
           )
    (setf (bonds-force *smirnoff*)
          (make-instance 'bonds-force)))
  (with-force-parser (node root "Bond")
    (let* ((attrs (plump:attributes node))
           (smirks (gethash "smirks" attrs))
           (id (gethash "id" attrs))
           (k (canonical-force (parse-quantity (gethash "k" attrs))))
           (len (canonical-length (parse-quantity (gethash "length" attrs))))
           (term (make-instance 'bond-term
                                :smirks smirks
                                :compiled-smirks (chem:compile-smarts smirks)
                                :id id :k k :len len)))
      (vector-push-extend term (terms (bonds-force *smirnoff*))))))

(defun parse-angles-force (root)
  (let* (#+(or)(attrs (plump:attributes root)))
    (setf (angles-force *smirnoff*) (make-instance 'angles-force))
    (with-force-parser (node root "Angle")
      (let* ((attrs (plump:attributes node))
             (smirks (gethash "smirks" attrs))
             (id (gethash "id" attrs))
             (k (canonical-angle-force (parse-quantity (gethash "k" attrs))))
             (angle (canonical-angle (parse-quantity (gethash "angle" attrs))))
             (term (make-instance 'angle-term
                                  :smirks smirks
                                  :compiled-smirks (chem:compile-smarts smirks)
                                  :id id
                                  :k k
                                  :angle angle)))
        (vector-push-extend term (terms (angles-force *smirnoff*)))))))

(defun parse-torsion (node torsion-class-name attributes)
  (let* ((smirks (gethash "smirks" attributes))
         (id (gethash "id" attributes))
         (term (make-instance torsion-class-name
                              :smirks smirks :id id
                              :compiled-smirks (chem:compile-smarts smirks))))
    (block terms
      (loop for index from 1
            for k-name = (format nil "k~a" index)
            for k-str = (gethash k-name attributes)
            if k-str
              do (let ((idivf (let ((idivf-str (gethash (format nil "idivf~a" index) attributes)))
                                (if idivf-str
                                    (fortran:parse-vecreal idivf-str)
                                    1.0d0)))
                       (k (canonical-kcal/mol (parse-quantity k-str)))
                       (periodicity
                         (parse-integer (gethash (format nil "periodicity~a" index) attributes)))
                       (phase
                         (canonical-angle (parse-quantity
                                           (gethash (format nil "phase~a" index) attributes)))))
                   (push (make-instance 'torsion-part
                                        :idivf idivf
                                        :k k
                                        :periodicity periodicity
                                        :phase phase)
                         (parts term)))
            else
              do (return-from terms)))
    term))

(defun parse-torsion-force (root name term-name force-name)
  (let* (#+(or)(attrs (plump:attributes root))
           (force (make-instance force-name)))
      (loop for node across (plump:children root)
            when (typep node 'plump:element)
              do (let ((tag (plump:tag-name node)))
                   (cond
                     ((string= tag name)
                      (vector-push-extend
                       (parse-torsion node term-name (plump:attributes node))
                       (terms force)))
                     (t (error "Illegal tag ~a" tag)))))
      force))


(defun parse-proper-torsion-force (root)
  (setf (proper-torsion-force *smirnoff*)
        (parse-torsion-force root "Proper" 'proper-term 'proper-torsion-force)))

(defun parse-improper-torsion-force (root)
  (setf (improper-torsion-force *smirnoff*)
        (parse-torsion-force root "Improper" 'improper-term 'improper-torsion-force)))

(defun parse-rmin_half (attrs)
  (let ((rmin_half (gethash "rmin_half" attrs))
        (sigma (gethash "sigma" attrs)))
    (cond
      (rmin_half
       (nanometer-length (parse-quantity rmin_half)))
      (sigma
       (let* ((sigma-nm (nanometer-length (parse-quantity sigma)))
              (rmin_half-nm (* (/ sigma-nm 2.0) (expt 2.0 (/ 1.0 6.0)))))
         rmin_half-nm))
      (t (error "Need either rmin_half or sigma for vdw parameter")))))

(defun parse-vdw-force (root)
  (flet ((safe-gethash (key ht)
           (let ((val (gethash key ht)))
             (unless val (error "The offxml file is missing ~a in ~a" key root))
             val)))
    (let* ((attrs (plump:attributes root))
           (potential (safe-gethash "potential" attrs))
           (combining-rules (safe-gethash "combining_rules" attrs))
           (scale12 (fortran:parse-vecreal (safe-gethash "scale12" attrs)))
           (scale13 (fortran:parse-vecreal (safe-gethash "scale13" attrs)))
           (scale14 (fortran:parse-vecreal (safe-gethash "scale14" attrs)))
           (scale15 (fortran:parse-vecreal (safe-gethash "scale15" attrs)))
           (switch-width (canonical-length (parse-quantity (safe-gethash "switch_width" attrs))))
           (cutoff (canonical-length (parse-quantity (safe-gethash "cutoff" attrs))))
           (method (safe-gethash "method" attrs))
           )
      (setf (vdw-force *smirnoff*)
            (make-instance 'vdw-force
                           :potential potential
                           :combining-rules combining-rules
                           :scale12 scale12
                           :scale13 scale13
                           :scale14 scale14
                           :scale15 scale15
                           :switch-width switch-width
                           :cutoff cutoff
                           :method method))
      (with-force-parser (node root "Atom")
        (let* ((attrs (plump:attributes node))
               (smirks (safe-gethash "smirks" attrs))
               (id (safe-gethash "id" attrs))
               (epsilon (canonical-kj/mol (parse-quantity (safe-gethash "epsilon" attrs))))
               (rmin-half (parse-rmin_half attrs))
               (type (or (gethash smirks *smirnoff-types*)
                         (let ((type (next-smirnoff-type-symbol)))
                           (setf (gethash smirks *smirnoff-types*) type)
                           type)))
               (term (make-instance 'nonbond-term
                                    :type type
                                    :smirks smirks
                                    :compiled-smirks (chem:compile-smarts smirks)
                                    :id id
                                    :epsilon epsilon
                                    :rmin-half rmin-half)))
          (vector-push-extend term (terms (vdw-force *smirnoff*))))))))

(defparameter *smirnoff-handlers*
  '(("Bonds" . parse-bonds-force)
    ("ProperTorsions" . parse-proper-torsion-force)
    ("ImproperTorsions" . parse-improper-torsion-force)
    ("Angles" . parse-angles-force)
    ("vdW" . parse-vdw-force)
    ("Electrostatics" . ignore-handler)
    ("ToolkitAM1BCC" . ignore-handler)
    ("Date" . ignore-handler)
    ("Author" . ignore-handler)
    ("Smirnoff" . parse-dispatcher )))

(defun parse-dispatcher (root)
  (loop for node across (plump:children root)
        when (typep node 'plump:element)
          do (let* ((tag (plump:tag-name node))
                    (found (assoc tag *smirnoff-handlers* :test #'string=)))
               (if found
                 (let ((handler (cdr found)))
                   (funcall handler node))
                 #+cando-warnings(warn "Ignoring tag: ~a~%" tag)))))

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

(defmethod chem:clear-combined-force-field ((combined-smirnoff-force-field combined-smirnoff-force-field))
  (setf (parts combined-smirnoff-force-field) nil))

(defmethod chem:add-shadowing-force-field ((combined-smirnoff-force-field combined-smirnoff-force-field) (smirnoff-force-field smirnoff-force-field) pathname)
  (push smirnoff-force-field (parts combined-smirnoff-force-field)))
