(in-package :cando-user)
(progn
  (defparameter c1 (chem:make-atom :c1 :c))
  (defparameter h11 (chem:make-atom :h11 :h))
  (defparameter h12 (chem:make-atom :h12 :h))
  (defparameter h13 (chem:make-atom :h13 :h))
  (defparameter c2 (chem:make-atom :c2 :c))
  (defparameter h21 (chem:make-atom :h21 :h))
  (defparameter h22 (chem:make-atom :h22 :h))
  (defparameter c3 (chem:make-atom :c3 :c))
  (defparameter h31 (chem:make-atom :h31 :h))
  (defparameter h32 (chem:make-atom :h32 :h))
  (defparameter c4 (chem:make-atom :c4 :c))
  (defparameter h41 (chem:make-atom :h41 :h))
  (defparameter h42 (chem:make-atom :h42 :h))
  (defparameter h43 (chem:make-atom :h43 :h))
  (chem:bond-to c1 h11 :single-bond)
  (chem:bond-to c1 h12 :single-bond)
  (chem:bond-to c1 h13 :single-bond)
  (chem:bond-to c1 c2 :single-bond)
  (chem:bond-to c2 h21 :single-bond)
  (chem:bond-to c2 h22 :single-bond)
  (chem:bond-to c2 c3 :single-bond)
  (chem:bond-to c3 h31 :single-bond)
  (chem:bond-to c3 h32 :single-bond)
  (chem:bond-to c3 c4 :single-bond)
  (chem:bond-to c4 h41 :single-bond)
  (chem:bond-to c4 h42 :single-bond)
  (chem:bond-to c4 h43 :single-bond)
  (defparameter res (chem:make-residue :butane))
  (loop for atom in (list c1 h11 h12 h13 c2 h21 h22 c3 h31 h32 c4 h41 h42 h43)
        do (chem:add-matter res atom))
  (defparameter mol (chem:make-molecule :butane))
  (chem:add-matter mol res)
  (cando:jostle mol)
  (leap:easy-gaff2)
  (leap:assign-atom-types mol)
  (leap:minimize mol)
  )




(progn
  
  (defparameter aimj (chem:make-atom-id-map))
  (chem:resize-aggregate aimj 1)
  (chem:resize-molecule aimj 0 1)
  (chem:resize-residue aimj 0 0 14)

  (defparameter aima (chem:make-atom-id-map))
  (chem:resize-aggregate aima 1)
  (chem:resize-molecule aima 0 1)
  (chem:resize-residue aima 0 0 14)

  

  (defun jj (atom-id atom)
    (let ((joint (kin:make-jump-joint atom-id (chem:get-name atom))))
      (chem:atom-id-map-set aimj atom-id joint)
      (chem:atom-id-map-set aima atom-id atom)
      joint))

  (defun cbj (atom-id atom distance-joint &optional angle-joint dihedral-joint)
    (let ((joint (kin:make-complex-bonded-joint atom-id (chem:get-name atom))))
      (kin:joint/add-child distance-joint joint)
      (when dihedral-joint
        (kin:set-input-stub-joint2 joint dihedral-joint))
      (when angle-joint
        (kin:set-input-stub-joint1 joint angle-joint))
      (chem:atom-id-map-set aimj atom-id joint)
      (chem:atom-id-map-set aima atom-id atom)
      joint))

  (defun bj (atom-id atom parent)
    (let ((joint (kin:make-bonded-joint atom-id (chem:get-name atom))))
      (kin:joint/add-child parent joint)
      (chem:atom-id-map-set aimj atom-id joint)
      (chem:atom-id-map-set aima atom-id atom)
      joint))

  (defparameter jh11 (jj  '(0 0 0) h11))
  (defparameter jc1  (cbj '(0 0 1) c1  jh11))
  (defparameter jh12 (cbj '(0 0 2) h12 jc1 jh11))
  (defparameter jh13 (cbj '(0 0 3) h13 jc1 jh12 jh11 ))
  (defparameter jc2  (cbj '(0 0 4) c2  jc1 jh12 jh11 ))
  (defparameter jh21 (bj  '(0 0 5) h21 jc2))
  (defparameter jh22 (bj  '(0 0 6) h22 jc2))
  (defparameter jc3  (bj  '(0 0 7) c3  jc2))
  (defparameter jh31 (bj  '(0 0 8) h31 jc3))
  (defparameter jh32 (bj  '(0 0 9) h32 jc3))
  (defparameter jc4  (bj  '(0 0 10) c4 jc3))
  (defparameter jh41 (bj  '(0 0 11) h41 jc4))
  (defparameter jh42 (bj  '(0 0 12) h42 jc4))
  (defparameter jh43 (bj  '(0 0 13) h43 jc4))
  (defparameter jroot jh11)
  )

(defun copy-externals-to-joints (aima aimj)
  (chem:walk aima (lambda (atom-id atom)
                    (let ((pos (chem:get-position atom))
                          (joint (chem:atom-id-map-get aimj atom-id)))
                      (kin:set-position joint pos)))))


(defun copy-externals-to-atoms (aimj aima)
  (chem:walk aimj (lambda (atom-id joint)
                    (let ((pos (kin:get-position joint))
                          (atom (chem:atom-id-map-get aima atom-id)))
                      (chem:set-position atom pos)))))



(copy-externals-to-joints aima aimj)

(kin:update-internal-coords jroot)

(defgeneric dump-internal (atom-id joint stream))

(defmethod dump-internal (atom-id (joint kin:jump-joint) stream)
  (format stream "jj  ~a ~3a~%" atom-id (kin:name joint)))

(defmethod dump-internal (atom-id (joint kin:complex-bonded-joint) stream)
  (format stream "cbj ~a ~a ~3,2f ~a ~3,2f ~a ~3,2f ~a~%" atom-id (kin:name joint)
          (kin:get-distance joint) (kin:input-stub-joint0 joint)
          (/ (kin:get-theta joint) 0.0174533) (when (kin:input-stub-joint1bound-p joint)
                                  (kin:input-stub-joint1 joint))
          (/ (kin:get-phi joint) 0.0174533)  (when (kin:input-stub-joint2bound-p joint)
                                (kin:input-stub-joint2 joint))
          ))


(defmethod dump-internal (atom-id (joint kin:bonded-joint) stream)
  (format stream "bj  ~a ~a ~3,2f ~a ~3,2f ~a ~3,2f ~a~%" atom-id
          (kin:name joint)
          (kin:get-distance joint)
          (kin:input-stub-joint0 joint)
          (/ (kin:get-theta joint) 0.0174533)
          (kin:input-stub-joint1 joint)
          (/ (kin:get-phi joint) 0.0174533)
          (kin:input-stub-joint1 joint)))

(defun dump-internals (aimj stream)
  (chem:walk aimj (lambda (atom-id joint)
                    (dump-internal atom-id joint stream))))


(defun input-stub (joint)
  (let ((input-stub (kin:get-input-stub joint)))
    input-stub
    ))

(defun input-stub-shapes (joint)
  (let* ((radius 0.1)
         (input-stub (kin:get-input-stub joint))
         (trans (geom:get-translation input-stub))
         (xcol (geom:get-x-column input-stub))
         (xpos (geom:v+ trans xcol))
         (ycol (geom:get-y-column input-stub))
         (ypos (geom:v+ trans ycol))
         (zcol (geom:get-z-column input-stub))
         (zpos (geom:v+ trans zcol))
         )
    `(( :type "cone"
        :position1 ,(geom:vlist trans)
        :position2 ,(geom:vlist xpos)
        :color (1 0 0)
        :radius ,radius)
      ( :type "cone"
        :position1 ,(geom:vlist trans)
        :position2 ,(geom:vlist ypos)
        :color (0 1 0)
        :radius ,radius)
      ( :type "cone"
        :position1 ,(geom:vlist trans)
        :position2 ,(geom:vlist zpos)
        :color (0 0 1)
        :radius ,radius))))


(defclass joint-builder ()
  ((atom-to-joint :initform (make-hash-table) :initarg :atom-to-joint :reader atom-to-joint)
   (atom-to-residue :initform (make-hash-table) :initarg :atom-to-residue :reader atom-to-residue)
   (residue-to-id :initform (make-hash-table) :initarg :residue-to-id :reader residue-to-id))
  )


(defun build-atom-to-residue-map (mol)
  (let ((a2r (make-hash-table)))
    (cando:do-residues (res mol)
      (chem:do-atoms (atm res)
        (setf (gethash atm a2r) res)))
    a2r))

(defun atom-tree-from-spanning-tree (atom)
  (let ((span (chem:make-spanning-loop atom)))
    (loop for atm = (chem:next span (lambda (&rest r) t))
          while atm do (let ((back (chem:get-back-span span atm)))
                         (print (list atm back))))))

#|
(kin:update-xyz-coords jroot)

(copy-externals-to-atoms aimj aima)


jroot
             
|#





(asdf:load-asd "/home/meister/Development/cando-llvm14/extensions/cando/src/lisp/design/design.asd")

(ql:quickload :design)
