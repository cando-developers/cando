(in-package :kin)

(defclass internal ()
  ((atom-index :initarg :atom-index :accessor atom-index)
   (atom-name :initarg :atom-name :accessor atom-name)
   ))

(cando:make-class-save-load internal
                            :print-unreadably
                            (lambda (obj stream)
                              (print-unreadable-object (obj stream :type t)
                                (format stream "~s" (atom-name obj)))))

(defclass bonded-internal (internal)
  ((distance :initform nil :initarg :distance :accessor distance)
   (angle :initform nil :initarg :angle :accessor angle)
   (dihedral :initform nil :initarg :dihedral :accessor dihedral)))

(cando:make-class-save-load
 bonded-internal
 :print-unreadably (lambda (obj stream)
                     (print-unreadable-object (obj stream :type t)
                       (format stream "~s ~5,3f ~5,3f ~5,3f"
                               (atom-name obj)
                               (distance obj)
                               (/ (angle obj) 0.0174533)
                               (/ (dihedral obj) 0.0174533)))))

(defclass jump-internal (internal)
  ())
(cando:make-class-save-load jump-internal
                            :print-unreadably
                            (lambda (obj stream)
                              (print-unreadable-object (obj stream :type t)
                                (format stream "~s" (atom-name obj)))))

(defclass complex-bonded-internal (bonded-internal)
  ())
(cando:make-class-save-load complex-bonded-internal
                            :print-unreadably
                            (lambda (obj stream)
                              (print-unreadable-object (obj stream :type t)
                                (format stream "~s" (atom-name obj)))))

(defgeneric internal-match (inta intb))

(defmethod internal-match ((inta bonded-internal) (intb bonded-internal))
  "Only compare dihedrals"
  (< (abs (- (dihedral inta) (dihedral intb))) (* 10.0 0.0174533)))

(defmethod internal-match ((inta jump-internal) (intb jump-internal))
  "jump-internal's always match"
  t)

(defmethod internal-match ((inta complex-bonded-internal) (intb complex-bonded-internal))
  "complex-bonded-internal only test if there is a dihedral"
  (if (dihedral inta)
      (< (abs (- (dihedral inta) (dihedral intb))) (* 30.0 0.0174533))
      t))

(defgeneric print-internal (obj stream))

(defmethod print-internal ((obj bonded-internal) stream)
  (format stream "~a ~a ~s ~5,3f ~8,3f ~8,3f~%"
          (class-name (class-of obj))
          (atom-index obj)
          (atom-name obj)
          (distance obj)
          (/ (angle obj) 0.0174533)
          (/ (dihedral obj) 0.0174533)))      

(defmethod print-internal ((obj complex-bonded-internal) stream)
  (format stream "~a ~a ~s ~5,3f ~{~8,3f ~}~%"
          (class-name (class-of obj))
          (atom-index obj)
          (atom-name obj)
          (distance obj)
          (cond
            ((and (angle obj) (dihedral obj))
             (list 
              (/ (angle obj) 0.0174533)
              (/ (dihedral obj) 0.0174533))
             )
            ((and (angle obj) (null (dihedral obj)))
             (list 
              (/ (angle obj) 0.0174533)
              ))
            (t nil))))

(defmethod print-internal ((obj jump-internal) stream)
  (format stream "~a ~a ~s~%"
          (class-name (class-of obj))
          (atom-index obj)
          (atom-name obj)))


(defclass rotamer ()
  ((cycle :initarg :cycle :accessor cycle)
   (residue-index :initarg :residue-index :accessor residue-index)
   (internals :initarg :internals :accessor internals)))

(cando:make-class-save-load rotamer)

(defun internals-match (internalsa internalsb)
  (loop for inta in (internals internalsa)
        for intb in (internals internalsb)
        when (not (internal-match inta intb))
          do (return-from internals-match nil))
  t)

(defmethod print-internal ((obj rotamer) stream)
  (format stream "BEGIN-ATRESIDUE~%")
  (format stream "CYCLE ~a~%" (cycle obj))
  (format stream "RESIDUE-INDEX ~a~%" (residue-index obj))
  (loop for internal in (internals obj)
        do (print-internal internal stream))
  (format stream "END-ATRESIDUE~%"))
        
(defclass rotamer-list ()
  ((cycles :initform 0 :initarg :cycles :accessor cycles)
   (entries :initform nil :initarg :entries
               :accessor entries)))

(cando:make-class-save-load rotamer-list)

(defmethod print-internal ((obj rotamer-list) stream)
  (format stream "BEGIN-ROTAMER-LIST :cycles ~a~%" (cycles obj))
  (loop for internals in (entries obj)
        do (print-internal internals stream))
  (format stream "END-ROTAMER-LIST~%"))

(defun save-rotamer-list (data filename)
  (with-open-file (fout filename :direction :output)
    (print-internal data fout)))

(defun load-rotamer-list (filename)
  (with-open-file (fin filename :direction :input)
    (parse-rotamer-list fin)))


(defclass rotamer-library ()
  ((context-to-rotamer-list :initform (make-hash-table)
                            :initarg :context-to-rotamer-list
                            :reader context-to-rotamer-list)))

(cando:make-class-save-load rotamer-library)


(defun number-of-rotamers (rotamer-library)
  (let ((count 0))
    (maphash (lambda (key value)
               (incf count (length (entries value))))
             (context-to-rotamer-list rotamer-library))
    count))


(defgeneric extract-joint-internals (joint index))

(defmethod extract-joint-internals ((joint kin:bonded-joint) index)
  (let ((atom-name (kin:name joint))
        (distance (kin:bonded-joint/get-distance joint))
        (angle (kin:bonded-joint/get-angle joint))
        (dihedral (kin:bonded-joint/get-dihedral joint)))
    (make-instance 'bonded-internal
                   :atom-index index
                   :atom-name atom-name
                   :distance distance
                   :angle angle
                   :dihedral dihedral)))

(defmethod extract-joint-internals ((joint kin:complex-bonded-joint) index)
  (let ((atom-name (kin:name joint))
        (distance (kin:bonded-joint/get-distance joint))
        (angle (if (kin:input-stub-joint1bound-p joint)
                      (kin:bonded-joint/get-angle joint)
                      nil))
        (dihedral (if (kin:input-stub-joint2bound-p joint)
                      (kin:bonded-joint/get-dihedral joint)
                      nil)))
    (make-instance 'complex-bonded-internal
                   :atom-index index
                   :atom-name atom-name
                   :distance distance
                   :angle angle
                   :dihedral dihedral)))

(defmethod extract-joint-internals ((joint kin:jump-joint) index)
  (let ((atom-name (kin:name joint)))
    (make-instance 'jump-internal
                   :atom-index index
                   :atom-name atom-name)))


(defun extract-rotamer (trainer &key index residue-index)
  (with-slots (conformation focus-residue-sequence-number)
      trainer
    (copy-atom-positions-into-joints conformation)
    (update-joint-tree-internal-coordinates conformation)
    (with-slots (ataggregate) conformation
      (let* ((atmolecule (atmolecule-aref ataggregate 0))
             (atresidue (atresidue-aref atmolecule focus-residue-sequence-number))
             (internals (loop for index below (joint-count atresidue)
                              for joint = (joint-aref atresidue index)
                              for internal = (extract-joint-internals joint index)
                              collect internal)))
        (make-instance 'rotamer
                       :cycle index
                       :residue-index focus-residue-sequence-number
                       :internals internals)))))



(defgeneric apply-internal-to-joint (joint internal))

(defmethod apply-internal-to-joint ((joint jump-joint) (internal jump-internal))
  ;; Do nothing
  )

(defmethod apply-internal-to-joint ((joint bonded-joint) (internal bonded-internal))
  (let ((distance (distance internal))
        (angle (angle internal))
        (dihedral (dihedral internal)))
    (kin:set-distance joint distance)
    (kin:set-theta joint angle)
    (kin:set-phi joint dihedral)))


(defmethod apply-internal-to-joint ((joint complex-bonded-joint) (internal complex-bonded-internal))
  (let ((distance (distance internal))
        (angle (angle internal))
        (dihedral (dihedral internal)))
    (kin:set-distance joint distance)
    (when (kin:input-stub-joint1bound-p joint)
      (kin:set-theta joint angle)
      (when (kin:input-stub-joint2bound-p joint)
        (kin:set-phi joint dihedral)))))

(defun apply-rotamer-to-atresidue (atresidue rotamer)
  (let ((joints (joints atresidue))
        (internals (internals rotamer)))
    (loop for index below (length joints)
          for joint = (aref joints index)
          for internal = (elt internals index)
          do (apply-internal-to-joint joint internal))))
