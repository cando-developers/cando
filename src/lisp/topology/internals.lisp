(in-package :topology)

(defclass internal ()
  ((name :initarg :name :accessor name)
   ))

(cando:make-class-save-load internal)

(defclass jump-internal (internal)
  ())

(cando:make-class-save-load jump-internal)

(defclass bonded-internal (internal)
  ((bond :initarg :bond :accessor bond)
   (angle :initarg :angle :accessor angle)
   (dihedral :initarg :dihedral :accessor dihedral)))

(cando:make-class-save-load bonded-internal)

(defclass complex-bonded-internal (bonded-internal)
  ())

(cando:make-class-save-load complex-bonded-internal)

(defclass out-of-focus-internal ()
  ((name :initarg :name :accessor name)
   (atres-name :initarg :atres-name :accessor atres-name)
   (p-name :initarg :p-name :accessor p-name)
   (p-atres-name :initarg :p-atres-name :accessor p-atres-name)
   (gp-name :initarg :gp-name :accessor gp-name)
   (gp-atres-name :initarg :gp-atres-name :accessor gp-atres-name)
   (ggp-name :initarg :ggp-name :accessor ggp-name)
   (ggp-atres-name :initarg :ggp-atres-name :accessor ggp-atres-name)
   (dihedral-rad :initarg :dihedral-rad :accessor dihedral-rad)))

(cando:make-class-save-load
 out-of-focus-internal
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a(~a) ~a(~a) ~a(~a) ~a(~a) ~5,2f"
             (name obj)
             (atres-name obj)
             (p-name obj)
             (p-atres-name obj)
             (gp-name obj)
             (gp-atres-name obj)
             (ggp-name obj)
             (ggp-atres-name obj)
             (dihedral-rad obj)))))

(defclass fragment-internals ()
  ((index :initarg :index :accessor index)
   (internals :initarg :internals :accessor internals)
   (out-of-focus-internals :initarg :out-of-focus-internals :accessor out-of-focus-internals)))

(cando:make-class-save-load fragment-internals)

(defclass fragment-conformations ()
  ((monomer-context :initarg :monomer-context :accessor monomer-context)
   (total-count :initform 0 :initarg :total-count :accessor total-count)
   (fragments :initform nil :initarg :fragments :accessor fragments)))

(cando:make-class-save-load fragment-conformations)

(defclass fragment-conformations ()
  ((monomer-context-to-fragment-conformations :initform (make-hash-table :test 'equal)
                                     :initarg :monomer-context-to-fragment-conformations
                                     :accessor monomer-context-to-fragment-conformations)))

(cando:make-class-save-load
 fragment-conformations
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t))))


(defconstant +dihedral-threshold+ (* 10.0 0.0174533))

(defun similar-internals-p (frag1 frag2 &optional )
  (loop for frag1-int in (internals frag1)
        for frag2-int in (internals frag2)
        do (when (and (typep frag1-int 'bonded-internal)
                      (typep frag2-int 'bonded-internal))
             (let* ((aa (- (dihedral frag1-int) (dihedral frag2-int)))
                    (aamod (- (mod (+ aa 180) 360) 180)))
               (when (> (abs aamod) +dihedral-threshold+)
                 (return-from similar-internals-p nil)))))
  t)

(defun seen-fragment-internals (fragment-conformations fragment-internals)
  (loop for seen-frag in (fragments fragment-conformations)
        when (similar-internals-p seen-frag fragment-internals)
          do (return-from seen-fragment-internals (index seen-frag)))
  nil)


(defun good-fragment-internals (fragment-internals)
  (loop for internal in (internals fragment-internals)
        do (cond
             ((typep internal 'jump-internal))
             ((> (bond internal) 1.8)
              (return-from good-fragment-internals nil))
             ))
  t)



(defun save-fragment-conformations (fragment-conformations filename)
  (cando:save-cando fragment-conformations filename))

(defun load-fragment-conformations (filename)
  (cando:load-cando filename))

(defun dump-fragment-internals (fragment-internals finternals)
  (format finternals "begin-conformation ~a~%" (index fragment-internals))
  (flet ((to-deg (rad)
           (/ rad 0.0174533)))
    (loop for internal in (topology:internals fragment-internals)
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
               )))
  (format finternals "end-conformation~%"))


(defgeneric fill-joint-internals (joint internal))

(defmethod fill-joint-internals ((joint kin:jump-joint) (internal jump-internal))
  )


(defmethod fill-joint-internals ((joint kin:bonded-joint) (internal bonded-internal))
  (kin:set-distance joint (bond internal))
  (kin:set-theta joint (angle internal))
  (kin:set-phi joint (dihedral internal))
  )
