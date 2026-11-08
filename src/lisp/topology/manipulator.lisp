(in-package :topology)

(defparameter +fscale-frac-bin-size+ 0.1) ; Bin Cremer-Pople q-values in 0.05 value bin sizes - it was 0.01 but that generates too many shape-key
(defparameter +fscale-deg-bin-size+ 10)     ; Bin angles in 5-degree bin sizes  - it was 1-degree but that generates too many shape-key


(defmacro with-cpring5 ((vec idx &optional suffix) &body body)
  "Inside BODY expose two SETF-able symbol-macros:

      @5q[-SUFFIX]   ≡ (aref VEC IDX)
      @5phi[-SUFFIX] ≡ (aref VEC (1+ IDX))

Supply SUFFIX (a symbol or string) when you need distinct names
while nesting several WITH-CPRING5 forms."
  (let* ((pkg *package*)
         (tail  (if suffix (format nil "-~A" suffix) ""))
         (qsym  (intern (concatenate 'string "@5Q"  tail)  pkg))
         (psym  (intern (concatenate 'string "@5PHI" tail) pkg))
         (v     (gensym "VEC-"))
         (i     (gensym "IDX-")))
    `(let ((,v ,vec)
           (,i ,idx))
       (symbol-macrolet ((,qsym  (aref ,v ,i))
                         (,psym  (aref ,v (1+ ,i))))
         ,@body))))

;;; Test code for with-cpring5
#+(or)
(let ((v (make-array 10 :element-type 'double-float :initial-element 0d0)))
  (with-cpring5 (v 2 outer)
    (with-cpring5 (v 4 inner)
      (setf @5q-outer 1d0   @5phi-outer 2d0   ; v[2], v[3]
            @5q-inner 3d0   @5phi-inner 4d0)  ; v[4], v[5]
      (+ @5q-outer @5phi-outer @5q-inner @5phi-inner))))
;; ⇒ 10.0d0


(defmacro with-cpring6-old ((cpinternals cpindex) &body body)
  "Expose three SETF-able symbol-macros inside BODY:
     @6q     ⇔ (aref CPINTERNALS CPINDEX)
     @6theta ⇔ (aref CPINTERNALS (+ 1 CPINDEX)
     @6phi   ⇔ (aref CPINTERNALS (+ 2 CPINDEX))"
  (let ((vec (gensym "VEC-"))
        (idx (gensym "IDX-")))
    `(let ((,vec ,cpinternals)
           (,idx ,cpindex))
       (symbol-macrolet ((@6q     (aref ,vec ,idx))
                         (@6theta (aref ,vec (+ 1 ,idx)))
                         (@6phi   (aref ,vec (+ 2 ,idx)))
                         )
                  ,@body))))

(defmacro with-cpring6 ((vec idx &optional suffix) &body body)
  "Inside BODY expose three SETF-able places:

     @6q[-SUFFIX]      ≡ (aref VEC IDX)
     @6theta[-SUFFIX]  ≡ (aref VEC (1+  IDX))
     @6phi[-SUFFIX]    ≡ (aref VEC (+ 2 IDX))"
  (let* ((pkg *package*)
         (tail (if suffix (format nil "-~A" suffix) ""))
         (qsym     (intern (concatenate 'string "@6Q"     tail) pkg))
         (thetasym (intern (concatenate 'string "@6THETA" tail) pkg))
         (phisym   (intern (concatenate 'string "@6PHI"   tail) pkg))
         (v  (gensym "VEC-"))
         (i  (gensym "IDX-")))
    `(let ((,v ,vec)
           (,i ,idx))
       (symbol-macrolet ((,qsym     (aref ,v ,i))
                         (,thetasym (aref ,v (1+ ,i)))
                         (,phisym   (aref ,v (+ 2 ,i))))
         ,@body))))


;;; Test code
#+(or)
(let ((v (make-array 10 :element-type 'double-float :initial-element 0d0)))
  (with-cpring6 (v 1 outer)
    (with-cpring6 (v 4 inner)
      (setf @6q-outer 1d0  @6theta-outer 2d0  @6phi-outer 3d0
            @6q-inner 4d0 @6theta-inner 5d0 @6phi-inner 6d0)
      (+ @6q-outer @6theta-outer @6phi-outer
         @6q-inner @6theta-inner @6phi-inner))))
;; ⇒ 21.0d0



(defmacro with-zmatrix-internals-old ((cpinternals cpindex) &body body)
  "Expose three SETF-able symbol-macros inside BODY:
     @distance    ⇔ (aref CPINTERNALS CPINDEX)
     @angle       ⇔ (aref CPINTERNALS (+ 1 CPINDEX)
     @dihedral    ⇔ (aref CPINTERNALS (+ 2 CPINDEX))"
  (let ((vec (gensym "VEC-"))
        (idx (gensym "IDX-")))
    `(let ((,vec ,cpinternals)
           (,idx ,cpindex))
       (symbol-macrolet ((@distance (aref ,vec ,idx))
                         (@angle    (aref ,vec (+ 1 ,idx)))
                         (@dihedral (aref ,vec (+ 2 ,idx)))
                         )
         ,@body))))


(defmacro with-zmatrix-internals ((vec idx &optional suffix) &body body)
  "Inside BODY treat the vector entry IDX and the two that follow as
   @distance[-SUFFIX]  ≡ (aref VEC IDX)
   @angle[-SUFFIX]     ≡ (aref VEC (1+  IDX))
   @dihedral[-SUFFIX]  ≡ (aref VEC (+ 2 IDX))"
  (let* ((pkg        *package*)
         (tail       (if suffix (format nil "-~A" suffix) ""))
         (dsym       (intern (concatenate 'string "@DISTANCE" tail) pkg))
         (asym       (intern (concatenate 'string "@ANGLE"    tail) pkg))
         (dhsym      (intern (concatenate 'string "@DIHEDRAL" tail) pkg))
         (v          (gensym "VEC-"))
         (i          (gensym "IDX-")))
    `(let ((,v ,vec)
           (,i ,idx))
       (symbol-macrolet ((,dsym  (aref ,v ,i))
                         (,asym  (aref ,v (1+ ,i)))
                         (,dhsym (aref ,v (+ 2 ,i))))
         ,@body))))

;; Quick sanity-check (no external test library needed)
#+(or)
(let ((v (make-array 7 :element-type 'double-float :initial-element 0d0)))
  (with-zmatrix-internals (v 1)            ; @distance @angle @dihedral
    (with-zmatrix-internals (v 4 sub)      ; @distance-sub @angle-sub @dihedral-sub
      (setf @distance      1d0
            @angle         2d0
            @dihedral      3d0
            @distance-sub  4d0
            @angle-sub     5d0
            @dihedral-sub  6d0)
      (assert (= (+ @distance @angle @dihedral
                    @distance-sub @angle-sub @dihedral-sub)
                 21d0)))))
;; ⟹ returns NIL if all assertions pass


;;; The old way of extracting the parameters
(defun extract-ring6-cps-old (cpinternals cpindex)
  (values (aref cpinternals cpindex)         ; q
          (aref cpinternals (+ 1 cpindex))   ; theta
          (aref cpinternals (+ 2 cpindex)))) ; phi

;;; The more readable and more future proof
(defun extract-ring6-cps (cpinternals cpindex)
  (with-cpring6 (cpinternals cpindex)
    (values @6q @6theta @6phi)))

;;; The old way of extracting the parameters
(defun extract-ring5-cps-old (cpinternals cpindex)
  (values (aref cpinternals cpindex)         ; q
          (aref cpinternals (+ 1 cpindex)))) ; phi

;;; The more readable and more future proof
(defun extract-ring5-cps (cpinternals cpindex)
  (with-cpring5 (cpinternals cpindex)
    (values @5q @5phi)))

(defclass manipulator ()
  ((assembler :initarg :assembler :reader assembler)
   (monomer-to-drivers-map :initform (make-hash-table) :reader monomer-to-drivers-map)
   (internal-index :initform 0 :initarg :internal-index :accessor internal-index)
   (other-drivers :initform nil :initarg :other-drivers :accessor other-drivers)
   (static-drivers :initform nil :initarg :static-drivers :accessor static-drivers)))

(defgeneric add-to-drivers (manipulator driver monomer))
(defgeneric sort-drivers (manipulator))

(defmethod add-to-drivers ((manipulator manipulator) driver monomer)
  (push driver (other-drivers manipulator))
  (push driver (gethash monomer (monomer-to-drivers-map manipulator)))
  )

(defmethod sort-drivers ((manipulator manipulator))
  (let ((sorted-other-drivers (sort (copy-seq (other-drivers manipulator)) #'< :key #'driver-index3)))
    (setf (other-drivers manipulator) sorted-other-drivers)))

(defclass focused-manipulator (manipulator)
  ((focus-monomer :initarg :focus-monomer :reader focus-monomer)
   (focus-drivers :initform nil :initarg :focus-drivers :accessor focus-drivers)
   ))

(defmethod add-to-drivers ((manipulator focused-manipulator) driver monomer)
  (if (eq monomer (focus-monomer manipulator))
      (push driver (focus-drivers manipulator))
      (push driver (other-drivers manipulator)))
  (push driver (gethash monomer (monomer-to-drivers-map manipulator)))
  )

(defmethod sort-drivers ((manipulator focused-manipulator))
  (call-next-method)
  (let ((sorted-focus-drivers (sort (copy-seq (focus-drivers manipulator)) #'< :key #'driver-index3)))
    (setf (focus-drivers manipulator) sorted-focus-drivers)))

(defclass range ()
  ((begin :initarg :begin :reader begin)
   (end :initarg :end :reader end)))

(defun new-range (manipulator size)
  "Advance the internal-index of the DRIVER by SIZE entries and return the range"
  (let ((begin (internal-index manipulator))
        (end (incf (internal-index manipulator) size)))
    (make-instance 'range :begin begin :end end)))

(defun range-size (range)
  "Return the size of the range"
  (- (end range) (begin range)))

(defclass ring-driver ()
  ((monomer :initarg :monomer :reader monomer)
   (cremer-pople-range :initarg :cremer-pople-range :reader cremer-pople-range)
   (lengths-range :initarg :lengths-range :reader lengths-range)
   (angles-range :initarg :angles-range :reader angles-range)
   (ring-joints :initarg :ring-joints :reader ring-joints)
   (exocyclic-joints :initarg :exocyclic-joints :reader exocyclic-joints)
   (build-joints :initarg :build-joints :initform (make-hash-table) :reader build-joints)
   ))

(defmethod print-object ((driver ring-driver) stream)
  (print-unreadable-object (driver stream :type t)
    (format stream "~s ~s" (first (ring-joints driver)) (monomer driver))))

(defun ring-driver-distance-angle-arrays (driver cpinternals)
  "Return displaced arrays into the distances and angles for the RING-DRIVER"
    (let ((distances (make-array (range-size (lengths-range driver)) :element-type 'double-float :displaced-to cpinternals :displaced-index-offset (begin (lengths-range driver))))
          (angles (make-array (range-size (angles-range driver)) :element-type 'double-float :displaced-to cpinternals :displaced-index-offset (begin (angles-range driver)))))
      (values distances angles)))

(defclass ring5-driver (ring-driver) ())

(defclass ring6-driver (ring-driver) ())

(defclass fused-ring ()
  ((ring1 :initarg :ring1 :reader ring1)
   (ring2 :initarg :ring2 :reader ring2)))

(defclass fused-ring-driver ()
  ((rings :initarg :rings :reader rings)))


(defmethod monomer ((driver fused-ring-driver))
  (monomer (ring2 (rings driver))))


(defmethod print-object ((driver fused-ring-driver) stream)
  (print-unreadable-object (driver stream :type t)
    (format stream "~s" (ring2 (rings driver)))))


(defgeneric driver-index3 (obj)
  (:documentation "Return the lowest index3 for the object"))

(defmethod driver-index3 ((ring-driver ring-driver))
  (kin:joint/position-index-x3 (first (ring-joints ring-driver))))

(defmethod driver-index3 ((fused-ring-driver fused-ring-driver))
  (kin:joint/position-index-x3 (first (ring-joints (ring2 (rings fused-ring-driver))))))

(defclass dihedral-driver ()
  ((name :initarg :name :initform nil :reader name)
   (monomer :initarg :monomer :reader monomer)
   (internal-range :initarg :internal-range :reader internal-range)
   (dihedral-joint :initarg :dihedral-joint :reader dihedral-joint)))

(defmethod print-object ((driver dihedral-driver) stream)
  (print-unreadable-object (driver stream :type t)
    (format stream "~s" (dihedral-joint driver))))

(defmethod driver-index3 ((dihedral-driver dihedral-driver))
  (kin:joint/position-index-x3 (dihedral-joint dihedral-driver)))

(defclass exocyclic-joint ()
  ((exo-range :initarg :exo-range :reader exo-range)
   (joints :initarg :joints :reader joints)))

(defclass ring-atom-and-exocyclic-atoms ()
  ((atm :initarg :atm :reader atm)
   (prev-atom :initarg :prev-atom :reader prev-atom)
   (next-atom :initarg :next-atom :reader next-atom)
   (exo-atoms :initarg :exo-atoms :reader exo-atoms)))

(defmethod print-object ((obj ring-atom-and-exocyclic-atoms) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~s" (atm obj))))


(defun identify-and-order-rings (assembler joint-to-atom)
  "Identify the rings in the ASSEMBLER and sort the atoms in each ring in order of how atoms are built"
  (let* ((aggregate (topology:aggregate assembler))
         (root-map (topology:root-map (topology:joint-tree assembler)))
         (root-joint (gethash (first (topology:always-oligomers assembler)) root-map))
         (stack (list (first root-joint)))
         (atoms (make-hash-table))
         (atom-index 0))
    (loop while
          (let* ((node (pop stack))
                 (atm (gethash node joint-to-atom)))
            (setf (gethash atm atoms) atom-index)
            (incf atom-index)
            (loop for child in (reverse (kin:joint/joint-children node))
                  do (push child stack))
            stack))
    ;; Order the atoms in rings
    (let ((mol (chem:content-at aggregate 0)))
      (chem-info:with-rings-do (rings (chem:identify-rings mol))
        (chem:with-aromaticity-information (aromaticity-info mol :am1bcc)
          (multiple-value-bind (non-aromatic-rings aromatic-rings)
              (loop for ring in rings
                    if (every (lambda (atm) (gethash atm aromaticity-info)) ring)
                      collect ring into aromatic-rings
                    else
                      collect ring into non-aromatic-rings
                    finally (return (values non-aromatic-rings aromatic-rings)))
            (flet ((order-atoms (xx yy)
                     (< (gethash xx atoms) (gethash yy atoms))))
              (let ((sorted-rings (loop for ring in non-aromatic-rings
                                        for sorted = (sort (copy-seq ring) #'order-atoms)
                                        collect sorted)))
                (values sorted-rings aromatic-rings)))))))))


(defparameter *manip-atom-joint-map* nil)
(defparameter *manip-rings* nil)
(defparameter *manip-aromatic-rings* nil)

(defmacro with-identified-and-ordered-rings ((rings aromatic-rings atom-joint-map assembler) &body body)
  "Find rings and aromatic-rings and run code within their scope"
  (let ((ht (gensym "HT-"))
        (joint (gensym "JOINT-"))
        (atm (gensym "ATM-"))
        (run (gensym "RUN-"))
        )
    `(flet ((,run (,rings ,aromatic-rings ,atom-joint-map)
              ,@body))
       (if *manip-atom-joint-map*
           (let ((,atom-joint-map *manip-atom-joint-map*)
                 (,rings *manip-rings*)
                 (,aromatic-rings *manip-aromatic-rings*))
             (,run ,rings ,aromatic-rings ,atom-joint-map))
           (let ((,atom-joint-map (let ((,ht (make-hash-table)))
                                    (topology:do-joint-atom (,joint ,atm ,assembler)
                                      (setf (gethash ,joint ,ht) ,atm)
                                      (setf (gethash ,atm ,ht) ,joint))
                                    ,ht)))
             (multiple-value-bind (,rings ,aromatic-rings)
                 (identify-and-order-rings ,assembler ,atom-joint-map)
               (let ((*manip-atom-joint-map* ,atom-joint-map)
                     (*manip-rings* ,rings)
                     (*manip-aromatic-rings* ,aromatic-rings))
                 (,run ,rings ,aromatic-rings ,atom-joint-map))))))))

(defun rings-and-exocyclic-atoms (rings)
  (flet ((find-exo-atoms (atm prev-atom next-atom)
           (loop for bonded in (chem:atom/bonded-atoms-as-list atm)
                 unless (or (eq prev-atom bonded) (eq next-atom bonded))
                   collect bonded)))
      (loop for ring in rings
            for ring-atm-ht = (let ((ht (make-hash-table)))
                                (loop for atm in ring
                                      do (setf (gethash atm ht) t))
                                ht)
            collect (loop with prev-atom = (car (last ring))
                          for cur on ring
                          for atm = (car cur)
                          for next-atom = (or (cadr cur)
                                              (car ring))
                          for exo-atoms = (find-exo-atoms atm prev-atom next-atom)
                          until (null cur)
                          collect (make-instance 'ring-atom-and-exocyclic-atoms
                                                 :atm atm
                                                 :prev-atom prev-atom
                                                 :next-atom next-atom
                                                 :exo-atoms exo-atoms)
                          do (setf prev-atom atm)))))

(defun build-one-ring-driver (manipulator ring-and-exo-atoms atom-joint-map focus-monomer joint-to-monomer)
  (declare (ignore focus-monomer))
  (let* ((cremer-pople-range (new-range manipulator (cond
                                                 ((= (length ring-and-exo-atoms) 5) 2)
                                                 ((= (length ring-and-exo-atoms) 6) 3)
                                                 (t (error "Illegal ring-and-exo-atoms length ~d" ring-and-exo-atoms)))))
         (lengths-range (new-range manipulator (length ring-and-exo-atoms)))
         (angles-range (new-range manipulator (length ring-and-exo-atoms)))
         (all-joints (make-hash-table))
         (ring-joints (loop for ring-atom-and-exo-atoms in ring-and-exo-atoms
                            for ring-atom = (atm ring-atom-and-exo-atoms)
                            for ring-joint = (gethash ring-atom atom-joint-map)
                            do (setf (gethash ring-joint all-joints) t)
                            collect ring-joint))
         (exocyclic-joints (loop for ring-atom-and-exo-atoms in ring-and-exo-atoms
                                 for prev-atom = (prev-atom ring-atom-and-exo-atoms)
                                 for next-atom = (next-atom ring-atom-and-exo-atoms)
                                 for cur-atom = (atm ring-atom-and-exo-atoms)
                                 append (loop for exo-atom in (exo-atoms ring-atom-and-exo-atoms)
                                              for exo-range = (new-range manipulator 3) ; store bond, angle, improper dihedral
                                              for exo-joint = (gethash exo-atom atom-joint-map)
                                              do (setf (gethash exo-joint all-joints) t)
                                              collect (make-instance 'exocyclic-joint
                                                                     :exo-range exo-range
                                                                     :joints (list (gethash prev-atom atom-joint-map)
                                                                                   (gethash next-atom atom-joint-map)
                                                                                   (gethash cur-atom atom-joint-map)
                                                                                   exo-joint)))))
         (build-joints (let (bj)
                         (flet ((in-ring-system-p (jnt)
                                  (etypecase jnt
                                    ;; Cartesian joints store their reconstructed position
                                    ;; directly; they need no internal-coordinate ancestry.
                                    ((or kin:jump-joint kin:xyz-joint) t)
                                    (kin:complex-bonded-joint
                                     (let* ((jp (kin:joint/parent jnt))
                                            (jgp (when (typep jp 'kin:complex-bonded-joint) (kin:joint/parent jp)))
                                            (jggp (when (typep jgp 'kin:complex-bonded-joint)  (kin:joint/parent jgp))))
                                       (cond
                                         (jggp (return-from in-ring-system-p (and (gethash jggp all-joints)
                                                                                  (gethash jgp all-joints)
                                                                                  (gethash jp all-joints))))
                                         (jgp (return-from in-ring-system-p (and (gethash jp all-joints)
                                                                                 (gethash jgp all-joints))))
                                         (jp (return-from in-ring-system-p (gethash jp all-joints)))
                                         (t (error "This should never happen")))))
                                    (kin:bonded-joint
                                     (let* ((jp (kin:joint/parent jnt))
                                            (jgp (kin:joint/parent jp))
                                            (jggp (kin:joint/parent jgp)))
                                       (and
                                        (gethash jp all-joints)
                                        (gethash jgp all-joints)
                                        (gethash jggp all-joints)))))))
                           (maphash (lambda (jnt dummy)
                                      (declare (ignore dummy))
                                      (when (in-ring-system-p jnt)
                                        (push jnt bj)))
                                    all-joints))
                         bj))
         (first-joint (first ring-joints))
         (monomer (gethash first-joint joint-to-monomer))
         (ring-driver (make-instance (cond
                                       ((= (length ring-joints) 5) 'ring5-driver)
                                       ((= (length ring-joints) 6) 'ring6-driver)
                                       (t (error "Wrong ring size - must be 5 or 6 - was ~d" (length ring-joints))))
                                     :monomer monomer
                                     :build-joints build-joints
                                     :cremer-pople-range cremer-pople-range
                                     :lengths-range lengths-range
                                     :angles-range angles-range
                                     :ring-joints ring-joints
                                     :exocyclic-joints exocyclic-joints)))
    (values ring-driver monomer)))


(defun count-shared-atoms (ring1 ring2)
  (let ((ring1-atoms (make-hash-table))
        (count 0))
    (loop for ring-atom-and-exocyclic-atoms in ring1
          do (setf (gethash (atm ring-atom-and-exocyclic-atoms) ring1-atoms) t))
    (loop for ring-atom-and-exocyclic-atoms in ring2
          when (gethash (atm ring-atom-and-exocyclic-atoms) ring1-atoms)
            do (incf count))
    count))

(defun build-all-ring-drivers (manipulator rings-and-exo-atoms atom-joint-map focus-monomer joints-to-monomers)
  (let ((rings-to-fused-rings (make-hash-table))
        (fused-rings nil))
    (loop for (ring1 . rest) on rings-and-exo-atoms while rest
          do (loop for ring2 in rest
                   for shared-atoms = (count-shared-atoms ring1 ring2)
                   when (= shared-atoms 2)
                     do (let* ((ring1-atom1 (atm (first ring1)))
                               (ring2-atom1 (atm (first ring2)))
                               (ring1-joint (gethash ring1-atom1 atom-joint-map))
                               (ring2-joint (gethash ring2-atom1 atom-joint-map))
                               (ring1-index3 (kin:joint/position-index-x3 ring1-joint))
                               (ring2-index3 (kin:joint/position-index-x3 ring2-joint)))
                          (if (< ring1-index3 ring2-index3)
                              (let ((fused-ring-pair (cons ring1 ring2)))
                                (push fused-ring-pair fused-rings)
                                (setf (gethash ring1 rings-to-fused-rings) fused-ring-pair
                                      (gethash ring2 rings-to-fused-rings) fused-ring-pair))
                              (let ((fused-ring-pair (cons ring2 ring1)))
                                (push fused-ring-pair fused-rings)
                                (setf (gethash ring1 rings-to-fused-rings) fused-ring-pair
                                      (gethash ring2 rings-to-fused-rings) fused-ring-pair))))))
    ;; Add the sington rings to the manipulator
    (loop for ring in rings-and-exo-atoms
          when (and (not (gethash ring rings-to-fused-rings))
                    (>= (length ring) 5))  ; ignore rings less than 5 atoms
            do (multiple-value-bind (ring-driver monomer)
                   (build-one-ring-driver manipulator ring atom-joint-map focus-monomer joints-to-monomers)
                 (unless monomer
                   (error "No monomer for joint ~s" ring-driver))
                 (add-to-drivers manipulator ring-driver monomer)
                 ))
    ;; Add the fused rings to the manipulator
    (loop for (ring1 . ring2 ) in fused-rings
          do (multiple-value-bind (ring1-driver ring1-monomer)
                 (build-one-ring-driver manipulator ring1 atom-joint-map focus-monomer joints-to-monomers)
               (declare (ignore ring1-monomer))
               (multiple-value-bind (ring2-driver ring2-monomer)
                   (build-one-ring-driver manipulator ring2 atom-joint-map focus-monomer joints-to-monomers)
                 (unless ring2-monomer
                   (error "No monomer for driver~s" ring2-driver))
                 (let ((fused-ring-driver (make-instance 'fused-ring-driver
                                                         :rings (make-instance 'fused-ring :ring1 ring1-driver :ring2 ring2-driver))))
                   (add-to-drivers manipulator fused-ring-driver ring2-monomer)
                   ))))
    ))


#+(or)
(defun build-dihedral-driver (static-p driver rotatable-joint focus-monomer joints-to-monomers)
  (let ((internal-range (new-range driver 3)))
    (let* ((monomer (gethash rotatable-joint joints-to-monomers))
           (dihedral-driver (make-instance 'dihedral-driver
                                           :monomer monomer
                                           :internal-range internal-range
                                           :dihedral-joint rotatable-joint)))
      (unless monomer
        (error "No monomer for joint ~s" rotatable-joint))
      (cond
        (static-p
         (push dihedral-driver (static-drivers driver)))
        ((eq monomer focus-monomer)
         (push dihedral-driver (focus-drivers driver)))
        (t (push dihedral-driver (other-drivers driver))))
      dihedral-defun)))

(defun dihedral-component-name (manipulator monomer joint)
  "Prefer the constitution's dihedral label, falling back to the atom name."
  (let* ((atom-name (kin:joint/name joint))
         (oligomer (oligomer-containing-monomer (assembler manipulator) monomer t))
         (topology (chem:find-topology (oligomer-monomer-name oligomer monomer)))
         (dihedrals (getf (residue-properties (constitution topology)) :dihedrals))
         (info (find-if (lambda (info)
                          (and (typep info 'dihedral-info-atom)
                               (eq (atom-name info) atom-name)))
                        dihedrals)))
    (if info (name info) atom-name)))

(defun make-dihedral-driver (manipulator dihedral-joint joints-to-monomers)
  (let* ((internal-range (new-range manipulator 3))
         (monomer (gethash dihedral-joint joints-to-monomers)))
    (values (make-instance 'dihedral-driver
                           :monomer monomer
                           :name (dihedral-component-name manipulator monomer dihedral-joint)
                           :internal-range internal-range
                           :dihedral-joint dihedral-joint)
            monomer)))

(defun build-rest-of-manipulator (manipulator rings-and-exo-atoms rotatable-dihedrals static-dihedrals joints-to-monomers atom-joint-map &key focus-monomer)
  (build-all-ring-drivers manipulator rings-and-exo-atoms atom-joint-map focus-monomer joints-to-monomers)
  (loop for rotatable-dihedral in rotatable-dihedrals
        do (multiple-value-bind (dihedral-driver monomer)
               (make-dihedral-driver manipulator rotatable-dihedral joints-to-monomers)
             (add-to-drivers manipulator dihedral-driver monomer)
             ))
  (sort-drivers manipulator)
  (loop for static-dihedral in static-dihedrals
        do (let ((dihedral-driver (make-dihedral-driver manipulator static-dihedral joints-to-monomers)))
             (push dihedral-driver (static-drivers manipulator))))
  manipulator)

(defun identify-build-joints (rings-and-exocyclic-atoms rotatable-dihedrals atom-joint-map)
  (let ((build-joints (make-hash-table)))
    (loop for ring-and-exocyclic-atoms in rings-and-exocyclic-atoms
          do (loop for ring-atom-and-exocyclic-atoms in ring-and-exocyclic-atoms
                   for atm = (atm ring-atom-and-exocyclic-atoms)
                   for atm-joint = (gethash atm atom-joint-map)
                   do (setf (gethash atm-joint build-joints) t)
                   do (loop for exo-atom in (exo-atoms ring-atom-and-exocyclic-atoms)
                            for exo-joint = (gethash exo-atom atom-joint-map)
                            do (setf (gethash exo-joint build-joints) t))))
    (loop for rotatable-dihedral in rotatable-dihedrals
          do (setf (gethash rotatable-dihedral build-joints) t))
    build-joints
    ))

(defun static-dihedral-joint-p (joint)
  "A static dihedral needs a non-complex bonded joint and three bound ancestors."
  (and (typep joint 'kin:bonded-joint)
       (not (typep joint 'kin:complex-bonded-joint))
       (loop for ancestor = joint then (kin:joint/parent ancestor)
             repeat 3
             always (kin:joint/parent-bound-p ancestor))))

(defun inter-ring-dihedral-allowed-p (atmb atmc rings-ht)
  "Allow a candidate torsion unless its central atoms share a ring."
  ;; Being in rings is not enough to exclude a bond: the connecting bond of
  ;; two separate rings must remain a rotation driver. Check every membership
  ;; so bonds within fused rings are still excluded. Ring objects in this table
  ;; are shared by identity, as populated by CALCULATE-ROTATABLE-DIHEDRALS.
  (not (intersection (gethash atmb rings-ht) (gethash atmc rings-ht) :test #'eq)))

(defun calculate-rotatable-dihedrals (assembler)
  (let ((rotatable-dihedrals nil)
        (aggregate (topology:aggregate assembler)))
    (with-identified-and-ordered-rings (rings aromatic-rings atom-joint-map assembler)
      (let ((rings-ht (let ((ht (make-hash-table)))
                        (loop for ring in rings
                              do (loop for atm in ring
                                       do (push ring (gethash atm ht))))
                        (loop for ring in aromatic-rings
                              do (loop for atm in ring
                                       do (push ring (gethash atm ht))))
                        ht)))
        (loop for atmol across (topology:atmolecules (topology:ataggregate assembler))
              do (loop for atres across (topology:atresidues atmol)
                       do (loop for joint across (topology:joints atres)
                                for id = (kin:joint/id joint)
                                for mol = (chem:content-at aggregate (first id))
                                for res = (chem:content-at mol (second id))
                                for atm = (chem:content-at res (third id))
                                when (and (eq (class-of joint) (load-time-value (find-class 'kin:bonded-joint)))
                                          (chem:heavy-element (chem:get-element atm)))
                                  do (let* ((parent (kin:joint/parent joint))
                                            (grand-parent (kin:joint/parent parent))
                                            #+(or)(great-grand-parent (kin:joint/parent grand-parent)))
                                       (flet ((rotatablep (jointb jointc)
                                                (let* ((atmb (gethash jointb atom-joint-map))
                                                       (atmc (gethash jointc atom-joint-map)))
                                                  (inter-ring-dihedral-allowed-p atmb atmc rings-ht))))
                                         (when (rotatablep parent grand-parent)
                                           #+(or)(format t "rotatable joint = ~s  parent = ~s~%" joint parent)
                                           (push joint rotatable-dihedrals)))))))
        rotatable-dihedrals))))


(defgeneric make-manipulator (assembler))


(defun initialize-manipulator-drivers (manipulator oligomer &key focus-monomer)
  "Populate a fresh MANIPULATOR with drivers for OLIGOMER and return it.
The caller chooses the manipulator class and supplies its assembler. Reuse
ring discovery bindings when available, and classify static joints here so
all callers apply the same joint-type and ancestry requirements."
  (let* ((assembler (assembler manipulator))
         (atmolecule (aref (atmolecules (ataggregate assembler)) 0))
         (joints-to-monomers (make-hash-table)))
    (loop for monomer across (monomers (oligomer-space oligomer))
          for monomer-pos = (gethash monomer (monomer-positions assembler))
          for atmol = (aref (atmolecules (ataggregate assembler))
                           (molecule-index monomer-pos))
          for atres = (aref (atresidues atmol) (residue-index monomer-pos))
          do (loop for joint across (joints atres)
                   do (setf (gethash joint joints-to-monomers) monomer)))
    (with-identified-and-ordered-rings (rings aromatic-rings atom-joint-map assembler)
      (declare (ignore aromatic-rings))
      (let* ((rotatable-dihedrals (calculate-rotatable-dihedrals assembler))
             (rings-and-exo-atoms (rings-and-exocyclic-atoms rings))
             (build-joints (identify-build-joints rings-and-exo-atoms
                                                 rotatable-dihedrals atom-joint-map))
             (static-dihedrals nil))
        (walk-atmolecule-joints
         atmolecule
         (lambda (joint pos)
           (declare (ignore pos))
           (when (and (static-dihedral-joint-p joint)
                      (null (gethash joint build-joints)))
             (push joint static-dihedrals))))
        (build-rest-of-manipulator manipulator rings-and-exo-atoms
                                   rotatable-dihedrals static-dihedrals
                                   joints-to-monomers atom-joint-map
                                   :focus-monomer focus-monomer)))))

(defmethod make-manipulator ((assembler topology:assembler))
  (let* ((oligomer-shape (first (oligomer-shapes assembler)))
         (oligomer (oligomer oligomer-shape)))
    (initialize-manipulator-drivers
     (make-instance 'manipulator :assembler assembler) oligomer)))

(defmethod make-manipulator ((assembler topology:training-assembler))
  (let* ((oligomer (first (oligomers assembler)))
         (focus-monomer (focus-monomer assembler))
         (manipulator (make-instance 'focused-manipulator
                                     :focus-monomer focus-monomer
                                     :assembler assembler)))
    (initialize-manipulator-drivers manipulator oligomer
                                    :focus-monomer focus-monomer)
    (values manipulator oligomer focus-monomer)))

(defun get-indexes3 (ring-joints)
  (let ((indexes3 (loop with indexes3 = (make-array (length ring-joints) :element-type 'ext:byte32)
                        for idx below (length ring-joints)
                        for joint in ring-joints
                        do (setf (aref indexes3 idx) (kin:joint/position-index-x3 joint))
                        finally (return indexes3))))
    indexes3))


(defun extract-exocyclic-joint (exocyclic-joint driver cpinternals externals)
  (declare (ignore driver))
  (let* ((joints (joints exocyclic-joint))
         (ja (first joints))
         (jb (second joints))
         (jc (third joints))
         (jd (fourth joints))
         (iva (kin:joint/position-index-x3 ja))
         (ivb (kin:joint/position-index-x3 jb))
         (ivc (kin:joint/position-index-x3 jc))
         (ivd (kin:joint/position-index-x3 jd))
         (distance (geom:calculate-distance-array ivc ivd externals))
         (angle (geom:calculate-angle-array ivb ivc ivd externals))
         (dihedral (geom:calculate-dihedral-array iva ivb ivc ivd externals))
         (exo-index (begin (exo-range exocyclic-joint))))
    (with-zmatrix-internals (cpinternals exo-index)
      (setf @distance distance
            @angle angle
            @dihedral dihedral))))

(defmethod fill-cpinternals-with-cremer-pople ((driver ring5-driver) cpinternals externals indexes3)
  (let ((cremer-pople-index (begin (cremer-pople-range driver))))
    (multiple-value-bind (q2 phi2)
        (chem:cpring-5-coordinates-to-cremer-pople externals indexes3)
      #+(or)
      (progn
        (format t "debug fill-cpinternals-with-cremer-pople --------------------------~%")
        (format t "coordinates = ~a~%" (chem:cpring-mathematica externals indexes3))
        (format t "cremer-pople = ~s~%" (list q2 phi2)))
      (with-cpring5 (cpinternals cremer-pople-index)
        (setf @5q q2
              @5phi phi2)))))

(defmethod fill-cpinternals-with-cremer-pople ((driver ring6-driver) cpinternals externals indexes3)
  (let ((cremer-pople-index (begin (cremer-pople-range driver))))
    (with-cpring6 (cpinternals cremer-pople-index)
      (multiple-value-setq (@6q @6theta @6phi)
        (chem:cpring-6-coordinates-to-cremer-pople externals indexes3)))))

(defmethod driver-extract-cpinternals-from-coordinates ((driver ring-driver) cpinternals externals)
  (let* ((ring-joints (ring-joints driver))
         (indexes3 (get-indexes3 ring-joints)))
    ;; Get the cremer-pople parameters
    (fill-cpinternals-with-cremer-pople driver cpinternals externals indexes3)
    ;; Get the bond lengths and angles
    (multiple-value-bind (distances angles)
        (ring-driver-distance-angle-arrays driver cpinternals)
      (chem:cpring-coordinates-to-distances-and-angles distances angles externals indexes3)
      (loop for exocyclic-joint in (exocyclic-joints driver)
            do (extract-exocyclic-joint exocyclic-joint driver cpinternals externals))
      )))

(defmethod driver-extract-cpinternals-from-coordinates ((driver fused-ring-driver) cpinternals externals)
  (driver-extract-cpinternals-from-coordinates (ring1 (rings driver)) cpinternals externals)
  (driver-extract-cpinternals-from-coordinates (ring2 (rings driver)) cpinternals externals))

(defmethod generate-coordinates-from-cpinternals ((driver ring-driver) externals cpinternals)
  (let* ((ring-joints (ring-joints driver))
         (indexes3 (get-indexes3 ring-joints))
         (cremer-pople-index (begin (cremer-pople-range driver))))
    ;; Get the cremer-pople parameters
    (cond
      ((= (length ring-joints) 5)
       (with-cpring5 (cpinternals cremer-pople-index)
         (multiple-value-setq (@5q @5phi)
           (chem:cpring-5-coordinates-to-cremer-pople externals indexes3))))
      ((= (length ring-joints) 6)
       (with-cpring6 (cpinternals cremer-pople-index)
         (multiple-value-setq (@6q @6theta @6phi)
           (chem:cpring-6-coordinates-to-cremer-pople externals indexes3))))
      (t (error "Bad ring size")))
    ;; Get the bond lengths and angles
    (multiple-value-bind (distances angles)
        (ring-driver-distance-angle-arrays driver cpinternals)
      (chem:cpring-coordinates-to-distances-and-angles distances angles externals indexes3)
      (loop for exocyclic-joint in (exocyclic-joints driver)
            do (extract-exocyclic-joint exocyclic-joint driver cpinternals externals))
      )))

(defmethod driver-extract-cpinternals-from-coordinates ((driver dihedral-driver) cpinternals externals)
  (let* ((joint (dihedral-joint driver))
         (jd joint)
         (jc (progn
               (assert (kin:joint/parent-bound-p jd))
               (kin:joint/parent jd)))
         (jb (progn
               (assert (kin:joint/parent-bound-p jc))
               (kin:joint/parent jc)))
         (ja (progn
               (assert (kin:joint/parent-bound-p jb))
               (kin:joint/parent jb)))
         (ivd (kin:joint/position-index-x3 jd))
         (ivc (kin:joint/position-index-x3 jc))
         (ivb (kin:joint/position-index-x3 jb))
         (iva (kin:joint/position-index-x3 ja))
         (distance (geom:calculate-distance-array ivc ivd externals))
         (angle (geom:calculate-angle-array ivb ivc ivd externals))
         (dihedral (geom:calculate-dihedral-array iva ivb ivc ivd externals))
         (internal-index (begin (internal-range driver))))
    (with-zmatrix-internals (cpinternals internal-index)
      (setf @distance distance
            @angle angle
            @dihedral dihedral))))

(defun make-cpinternals-for-manipulator (manipulator)
  (make-array (internal-index manipulator) :element-type 'double-float))

(defgeneric extract-cpinternals-from-coordinates (manipulator cpinternals externals)
  (:documentation "Extract into CPINTERNALS internal coordinates derived from EXTERNALS"))

(defmethod extract-cpinternals-from-coordinates ((manipulator manipulator) cpinternals externals)
  (loop for driver in (other-drivers manipulator)
        do (driver-extract-cpinternals-from-coordinates driver cpinternals externals))
  (loop for driver in (static-drivers manipulator)
        do (driver-extract-cpinternals-from-coordinates driver cpinternals externals)))

(defmethod extract-cpinternals-from-coordinates ((manipulator focused-manipulator) cpinternals externals)
  "In addition to the cpinternals for manipulator instances focused-manipulators have focus-drivers"
  (loop for driver in (focus-drivers manipulator)
        do (driver-extract-cpinternals-from-coordinates driver cpinternals externals))
  (call-next-method))


(defun build-exocyclic-joints (exocyclic-joints temp-externals cpinternals)
  (loop for exocyclic-joint in exocyclic-joints
        for exo-index = (begin (exo-range exocyclic-joint))
        for ja = (first (joints exocyclic-joint))
        for jb = (second (joints exocyclic-joint))
        for jc = (third (joints exocyclic-joint))
        for jd = (fourth (joints exocyclic-joint))
        for ia = (kin:joint/position-index-x3 ja)
        for ib = (kin:joint/position-index-x3 jb)
        for ic = (kin:joint/position-index-x3 jc)
        for id = (kin:joint/position-index-x3 jd)
        do (with-zmatrix-internals (cpinternals exo-index)
             (geom:in-place-build-using-bond-angle-dihedral-array temp-externals id
                                                                  @distance ic
                                                                  @angle ib
                                                                  @dihedral ia))

        ))

#+(or)(defparameter *rings* nil)

(defun finish-ring-internals (assembler internals zjs lengths angles temp-externals indexes3 cpinternals driver &key debug-cremer-pople)
  (declare (ignore assembler debug-cremer-pople))
  #+(or)(setf temp-externals (chem:make-nvector (length temp-externals)))
  (chem:cpring-generate-coordinates zjs lengths angles temp-externals indexes3)
  #+(or)
  (when debug-cremer-pople
    (format t "---------------------------------------~%")
    (format t "    zjs = ~s~%" zjs)
    (let* ((new-zjs (chem:cpring-coordinates-to-zjs temp-externals indexes3))
           (new-cpvalues (cond
                           ((= (length zjs) 5)
                            (multiple-value-list (chem:cpring-5-coordinates-to-cremer-pople temp-externals indexes3)))
                           ((= (length zjs) 6)
                            (multiple-value-list (chem:cpring-6-coordinates-to-cremer-pople temp-externals indexes3)))
                           (t (error "Bad ring size ~d" (length zjs))))))
      (format t "new-zjs = ~s~%" new-zjs)
      (format t "~a~%" (chem:cpring-mathematica temp-externals indexes3))
      (format t "original cpvalues = ~s~%" debug-cremer-pople)
      (format t "     new cpvalues = ~s~%" new-cpvalues)))
  (build-exocyclic-joints (exocyclic-joints driver) temp-externals cpinternals)
  #+(or)(push (cons assembler (copy-seq temp-externals)) *rings*)
  (loop for build-joint in (build-joints driver)
        do (kin:joint/update-internal-coord build-joint internals temp-externals))
  )

(defgeneric build-only-ring-externals-from-cpinternals (assembler driver internals temp-externals cpinternals))

(defmethod build-only-ring-externals-from-cpinternals (assembler (driver ring5-driver) internals temp-externals cpinternals)
  (let* ((ring-joints (ring-joints driver))
         (indexes3 (get-indexes3 ring-joints))
         (cremer-pople-index (begin (cremer-pople-range driver)))
         (zjs (with-cpring5 (cpinternals cremer-pople-index)
                (chem:cpring-5-cremer-pople-to-zjs @5q @5phi))))
    (multiple-value-bind (distances angles)
        (ring-driver-distance-angle-arrays driver cpinternals)
      (chem:cpring-generate-coordinates zjs distances angles temp-externals indexes3))))

(defmethod build-only-ring-externals-from-cpinternals (assembler (driver ring6-driver) internals temp-externals cpinternals)
  (let* ((ring-joints (ring-joints driver))
         (indexes3 (get-indexes3 ring-joints))
         (cremer-pople-index (begin (cremer-pople-range driver)))
         (zjs (with-cpring6 (cpinternals cremer-pople-index)
                (chem:cpring-6-cremer-pople-to-zjs @6q @6theta @6phi))))
    (multiple-value-bind (distances angles)
        (ring-driver-distance-angle-arrays driver cpinternals)
      (chem:cpring-generate-coordinates zjs distances angles temp-externals indexes3))))


(defmethod build-internals-from-cpinternals (assembler (driver ring5-driver) internals temp-externals cpinternals)
  (let* ((ring-joints (ring-joints driver))
         (indexes3 (get-indexes3 ring-joints))
         (cremer-pople-index (begin (cremer-pople-range driver)))
         (zjs (with-cpring5 (cpinternals cremer-pople-index)
                (chem:cpring-5-cremer-pople-to-zjs @5q @5phi))))
    (multiple-value-bind (distances angles)
        (ring-driver-distance-angle-arrays driver cpinternals)
    (finish-ring-internals assembler internals zjs distances angles temp-externals indexes3 cpinternals driver))))

(defmethod build-internals-from-cpinternals (assembler (driver ring6-driver) internals temp-externals cpinternals)
  (let* ((ring-joints (ring-joints driver))
         (indexes3 (get-indexes3 ring-joints))
         (cremer-pople-index (begin (cremer-pople-range driver)))
         (zjs (with-cpring6 (cpinternals cremer-pople-index)
                (chem:cpring-6-cremer-pople-to-zjs @6q @6theta @6phi))))
    (multiple-value-bind (distances angles)
        (ring-driver-distance-angle-arrays driver cpinternals)
      (finish-ring-internals assembler internals zjs distances angles temp-externals indexes3 cpinternals driver))))

(defmethod build-internals-from-cpinternals (assembler (driver fused-ring-driver) internals temp-externals cpinternals)
  (build-internals-from-cpinternals assembler (ring1 (rings driver)) internals temp-externals cpinternals)
  (build-internals-from-cpinternals assembler (ring2 (rings driver)) internals temp-externals cpinternals))

(defmethod build-internals-from-cpinternals (assembler (driver dihedral-driver) internals temp-externals cpinternals)
  (declare (ignore assembler temp-externals))
  (let ((cpinternals-index3 (begin (internal-range driver)))
        (internals-index3 (kin:joint/position-index-x3 (dihedral-joint driver))))
    (setf (aref internals internals-index3) (aref cpinternals cpinternals-index3)
          (aref internals (+ 1 internals-index3)) (aref cpinternals (+ 1 cpinternals-index3))
          (aref internals (+ 2 internals-index3)) (aref cpinternals (+ 2 cpinternals-index3)))))

(defgeneric manipulator-build-internals-from-cpinternals (manipulator temp-externals cpinternals internals))

(defmethod manipulator-build-internals-from-cpinternals ((manipulator manipulator) temp-externals cpinternals internals)
  "Build the INTERALS from CPINTERNALS and use TEMP-EXTERNALS as a scratch pad of coordinates."
  (let ((assembler (assembler manipulator)))
    (loop for driver in (other-drivers manipulator)
          do (build-internals-from-cpinternals (assembler manipulator) driver internals temp-externals cpinternals))
    (loop for driver in (static-drivers manipulator)
          do (build-internals-from-cpinternals (assembler manipulator) driver internals temp-externals cpinternals))
    ))

(defmethod manipulator-build-internals-from-cpinternals ((manipulator focused-manipulator) temp-externals cpinternals internals)
  "Build the INTERALS from CPINTERNALS and use TEMP-EXTERNALS as a scratch pad of coordinates."
  (let ((assembler (assembler manipulator)))
    (loop for driver in (focus-drivers manipulator)
          do (build-internals-from-cpinternals (assembler manipulator) driver internals temp-externals cpinternals)))
  (call-next-method))

(defun fscale-frac (num &optional (bin-size +fscale-frac-bin-size+))
  (floor (fround (/ num bin-size))))

(defun fscale-frac-reverse (int &optional (bin-size +fscale-frac-bin-size+))
  (* int bin-size))

(defun fscale-deg (num &optional (bin-size +fscale-deg-bin-size+))
  (floor (fround (/ (topology:rad-to-deg num) bin-size))))

(defun fscale-deg-reverse (int &optional (bin-size +fscale-deg-bin-size+))
  (topology:deg-to-rad (* int bin-size)))


;;; SHAPE-KEY-LEGACY: moved to shape-key.lisp; regression: tests/shape-key.lisp.
#+(or)
(defclass binned-shape-key-cache (cando.serialize:serializable)
  ((entries :initarg :entries :initform nil :reader binned-shape-key-entries))
  (:documentation "An alist mapping symbolic names to nonempty integer-bin lists.
A dihedral has one bin; five- and six-membered rings have two and three bins.
Treat the entries as read-only. Ring labels count rings; dihedrals retain their names."))

;;; SHAPE-KEY-LEGACY: moved to shape-key.lisp.
#+(or)
(defclass binned-shape-key (binned-shape-key-cache) ()
  (:documentation "A selected shape key, distinct from the complete driver cache.
Shares the cache's labeled components, validation, and lookup protocol."))

;;; SHAPE-KEY-LEGACY: replaced by canonical construction in shape-key.lisp.
#+(or)
(defun make-binned-shape-key (entries)
  "Construct a selected key from an alist of names and bin lists."
  (make-instance 'binned-shape-key :entries entries))

;;; SHAPE-KEY-LEGACY: replaced by canonical construction in shape-key.lisp.
#+(or)
(defun make-binned-shape-key-cache (entries)
  "Construct a cache from an alist of names and bin lists."
  (make-instance 'binned-shape-key-cache :entries entries))

;;; SHAPE-KEY-LEGACY: moved to shape-key.lisp.
#+(or)
(defun binned-shape-key-proper-list-p (value)
  (and (listp value)
       (handler-case (integerp (list-length value))
         (type-error () nil))))

;;; SHAPE-KEY-LEGACY: validation and duplicate rejection now in shape-key.lisp.
#+(or)
(defmethod initialize-instance :after ((key binned-shape-key-cache) &key)
  (let ((entries (binned-shape-key-entries key)))
    (unless (and (binned-shape-key-proper-list-p entries)
                 (every (lambda (entry)
                          (and (consp entry)
                               (car entry) (symbolp (car entry))
                               (binned-shape-key-proper-list-p (cdr entry))
                               (<= 1 (length (cdr entry)) 3)
                               (every #'integerp (cdr entry))))
                        entries))
      (let ((*print-circle* t))
        (error "Invalid binned shape key entries: ~s" entries)))))

;;; SHAPE-KEY-LEGACY: moved to shape-key.lisp.
#+(or)
(defmethod binned-shape-key-names ((key binned-shape-key-cache))
  (mapcar #'car (binned-shape-key-entries key)))

;;; SHAPE-KEY-LEGACY: defensive reader now in shape-key.lisp.
#+(or)
(defmethod binned-shape-key-components ((key binned-shape-key-cache))
  (mapcar #'cdr (binned-shape-key-entries key)))

;;; SHAPE-KEY-LEGACY: moved to shape-key.lisp.
#+(or)
(defmethod print-object ((key binned-shape-key-cache) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (key stream :type t)
        (write (binned-shape-key-entries key) :stream stream))))

;;; SHAPE-KEY-LEGACY: duplicates now rejected at construction/insertion.
#+(or)
(defun binned-shape-key-component (key name)
  "Return the bin list labeled NAME; signal missing or ambiguous labels."
  (let ((entries (binned-shape-key-entries key)))
    (unless (= 1 (count name entries :key #'car :test #'eq))
      (error "Expected one component named ~s in binned shape key; available names: ~s"
             name (binned-shape-key-names key)))
    (cdr (assoc name entries :test #'eq))))

;;; SHAPE-KEY-LEGACY: positional ring labels; use declared geometry and validated cache insertion.
#+(or)
(defun combine-binned-shape-key-caches (keys)
  "Concatenate complete components, numbering rings independently of dihedrals."
  (error "Do I want to hard code ring component names")
  (let ((entries nil) (ring-index 0))
    (dolist (key keys)
      (loop for (name . component) in (binned-shape-key-entries key)
            do (push (cons (if (> (length component) 1)
                               (case (incf ring-index)
                                 (1 :first-ring)
                                 (2 :second-ring)
                                 (otherwise (intern (format nil "RING-~d" ring-index) :keyword)))
                               name)
                           component)
                     entries)))
    (make-binned-shape-key-cache (nreverse entries))))

;;; SHAPE-KEY-LEGACY: use measure-binned-shape-key-cache with a foldamer policy.
#+(or)
(defgeneric binned-shape-key-cache (driver cpinternals)
  (:documentation "Return a BINNED-SHAPE-KEY-CACHE containing labeled bin lists."))
;;; SHAPE-KEY-LEGACY: use write-ring-cache-component with a foldamer policy.
#+(or)
(defgeneric binned-shape-key-write (driver cpinternals shape-key))

;;; SHAPE-KEY-LEGACY: declared names and foldamer binning replace this method.
#+(or)
(defmethod binned-shape-key-cache ((driver ring5-driver) cpinternals)
  (let* ((cpindex (begin (cremer-pople-range driver))))
    (with-cpring5 (cpinternals cpindex)
      (make-binned-shape-key-cache
       (list (list :first-ring (fscale-frac @5q)
                   (fscale-deg @5phi)))))))

;;; SHAPE-KEY-LEGACY: replaced by write-ring-cache-component; tests/cache-geometry.
#+(or)
(defmethod binned-shape-key-write ((driver ring5-driver) cpinternals shape-key-part)
  (let* ((cpindex (begin (cremer-pople-range driver))))
    (with-cpring5 (cpinternals cpindex)
      (let* ((5q (fscale-frac-reverse (first shape-key-part)))
             (5phi (fscale-deg-reverse (second shape-key-part))))
        (setf @5q 5q
              @5phi 5phi)))))

;;; SHAPE-KEY-LEGACY: declared names and foldamer binning replace this method.
#+(or)
(defmethod binned-shape-key-cache ((driver ring6-driver) cpinternals)
  (let* ((cpindex (begin (cremer-pople-range driver))))
    (with-cpring6 (cpinternals cpindex)
      (make-binned-shape-key-cache
       (list (list :first-ring (fscale-frac @6q)
                   (fscale-deg @6theta)
                   (fscale-deg @6phi)))))))

;;; SHAPE-KEY-LEGACY: replaced by write-ring-cache-component; tests/cache-geometry.
#+(or)
(defmethod binned-shape-key-write ((driver ring6-driver) cpinternals shape-key-part)
  (let* ((cpindex (begin (cremer-pople-range driver))))
    (with-cpring6 (cpinternals cpindex)
      (let* ((6q (fscale-frac-reverse (first shape-key-part)))
             (6theta (fscale-deg-reverse (second shape-key-part)))
             (6phi (fscale-deg-reverse (third shape-key-part))))
        (setf @6q 6q
              @6theta 6theta
              @6phi 6phi)))))

;;; SHAPE-KEY-LEGACY: cache declarations select individual rings; coupling is deferred.
#+(or)
(defmethod binned-shape-key-cache ((driver fused-ring-driver) cpinternals)
  (let ((rings (rings driver)))
    (binned-shape-key-cache (list (ring1 rings) (ring2 rings)) cpinternals)))

;;; SHAPE-KEY-LEGACY: declared names and foldamer binning replace this method.
#+(or)
(defmethod binned-shape-key-cache ((driver dihedral-driver) cpinternals)
  (let ((index (begin (internal-range driver))))
    (with-zmatrix-internals (cpinternals index)
      (make-binned-shape-key-cache
       (list (list (or (name driver)
                       (kin:joint/name (dihedral-joint driver)))
                   (fscale-deg @dihedral)))))))

;;; SHAPE-KEY-LEGACY: replaced by measure-binned-shape-key-cache; tests/cache-geometry.
#+(or)
(defmethod binned-shape-key-cache ((drivers list) cpinternals)
  (combine-binned-shape-key-caches
   (mapcar (lambda (driver) (binned-shape-key-cache driver cpinternals)) drivers)))

;;; SHAPE-KEY-LEGACY: use binned-shape-key-pathname-part for complete selected keys.
#+(or)
(defmethod binned-shape-key-string ((driver ring5-driver) cpinternals)
  (destructuring-bind (q2i phi2i)
      (binned-shape-key-component (binned-shape-key-cache driver cpinternals) :first-ring)
    (format nil "q~dp~d" q2i phi2i)))

;;; SHAPE-KEY-LEGACY: use binned-shape-key-pathname-part for complete selected keys.
#+(or)
(defmethod binned-shape-key-string ((driver ring6-driver) cpinternals)
  (destructuring-bind (qqi thetai phi2i)
      (binned-shape-key-component (binned-shape-key-cache driver cpinternals) :first-ring)
    (format nil "Q~dt~dp~d" qqi thetai phi2i)))

(defmethod dump-cpinternals* ((driver ring5-driver) cpinternals)
  (let ((cpindex (begin (cremer-pople-range driver))))
    (with-cpring5 (cpinternals cpindex)
      (format nil "~8,5f ~8,1fdeg"
              @5q
              (topology:rad-to-deg @5phi)))))

(defmethod dump-cpinternals* ((driver ring6-driver) cpinternals)
  (let ((cpindex (begin (cremer-pople-range driver))))
    (with-cpring6 (cpinternals cpindex)
      (format nil "~8,5f ~8,1fdeg ~8,1fdeg"
              @6q
              (topology:rad-to-deg @6theta)
              (topology:rad-to-deg @6phi)))))

(defmethod dump-cpinternals* ((driver fused-ring-driver) cpinternals)
  (format nil "~a ~a"
          (dump-cpinternals* (ring1 (rings driver)) cpinternals)
          (dump-cpinternals* (ring2 (rings driver)) cpinternals)))


(defmethod dump-cpinternals* ((driver dihedral-driver) cpinternals)
  (declare (ignore driver cpinternals stream))
  nil
)

(defgeneric dump-manipulator-cpinternals (manipulator cpinternals &optional stream))

(defmethod dump-manipulator-cpinternals ((manipulator manipulator) cpinternals &optional (stream t))
  (format stream "other-drivers---~%")
  (loop for driver in (other-drivers manipulator)
        for repr = (dump-cpinternals* driver cpinternals)
        when repr
             do (format stream "~a~%" repr)))

(defmethod dump-manipulator-cpinternals ((manipulator focused-manipulator) cpinternals &optional (stream t))
  (format stream "focus-drivers---~%")
  (loop for driver in (focus-drivers manipulator)
        for repr = (dump-cpinternals* driver cpinternals)
        when repr
             do (format stream "~a : ~s~%" repr driver))
  (call-next-method))

(defparameter *manipulate-quoted-keywords* '(:ring :dih))
(defun process-clause (monomer
                       &rest clause
                       &key &allow-other-keys)
  (list* `',monomer
         (loop for (key value) on clause by #'cddr
               collect key
               when (member key *manipulate-quoted-keywords*)
                 collect `',value
               else
                 collect value)))


(defgeneric do-manipulation (driver cpinternals target target-value &key &allow-other-keys))

(defmethod do-manipulation ((driver fused-ring-driver) cpinternals (target (eql :ring)) target-value &key qq theta phi theta-deg phi-deg)
  (when (and theta theta-deg)
    (error "Only one of :theta or :theta-deg are allowed"))
  (when (and phi phi-deg)
    (error "Only one of :phi or :phi-deg are allowed"))
  (cond
    ((string-equal "dkp" (string target-value))
     (with-cpring6 (cpinternals (begin (cremer-pople-range (ring1 (rings driver)))))
       (when qq (setf @6q qq))
       (when theta (setf @6theta theta))
       (when theta-deg (setf @6theta (deg-to-rad theta-deg)))
       (when phi (setf @6phi phi))
       (when phi-deg (setf @6phi (deg-to-rad phi-deg)))
       ))
    ((string-equal "bb" (string target-value))
     (let ((ring (ring2 (rings driver))))
       (etypecase ring
         (ring5-driver
          (with-cpring5 (cpinternals (begin (cremer-pople-range ring)))
            (when qq (setf @5q qq))
            (when phi (setf @5phi phi))
            (when phi-deg (setf @5phi (deg-to-rad phi-deg)))
            ))
         (ring6-driver
          (with-cpring6 (cpinternals (begin (cremer-pople-range ring)))
            (when qq (setf @6q qq))
            (when theta (setf @6theta theta))
            (when theta-deg (setf @6theta (deg-to-rad theta-deg)))
            (when phi (setf @6phi phi))
            (when phi-deg (setf @6phi (deg-to-rad phi-deg)))
            )))))
    (t (error ":ring argument must be one of dkp or bb - saw ~s" target-value)))
  )

(defmethod do-manipulation (driver cpinternals (target (eql :dih)) target-value &rest args)
  (break "Check driver and clause for dih"))


(defun %manipulate (manipulator cpinternals monomer &rest clause &key &allow-other-keys)
  (let* ((assembler (topology:assembler manipulator))
         (oligomer-shape (first (topology:oligomer-shapes assembler)))
         (oligomer (topology:oligomer oligomer-shape))
         (monomer-string (symbol-name monomer)))
    (loop for driver in (other-drivers manipulator)
          for driver-monomer = (monomer driver)
          for driver-monomer-symbol = (topology:oligomer-monomer-name oligomer driver-monomer)
          for driver-monomer-string = (symbol-name driver-monomer-symbol)
          when (string= driver-monomer-string monomer-string)
            do (apply 'do-manipulation driver cpinternals clause)
          )))


#|
(manipulate manipulator
  ( pro4ss :ring dkp :qq 0.32 :phi (deg 25.0) :theta (deg 90.0))
  ( pro4rr :ring dkp :phi (deg (+ -180.0 25.0)))
  ( 4pr :dih c1 :angle (deg 180))
  ( foo :dih c2 :angle (deg 180))
)
|#

(defmacro manipulate (manipulator cpinternals &body clauses)
  (let ((man (gensym "MANIPULATOR"))
        (cpi (gensym "CPINTERNALS")))
    `(let ((,man ,manipulator)
           (,cpi ,cpinternals))
       ,@(loop for clause in clauses
               collect `(%manipulate ,man ,cpi ,@(apply 'process-clause clause)))))
  )



(defun externals-from-cpinternals (manipulator cpinternals assembler-internals &key 
                                                             (externals (make-coordinates-for-assembler (assembler manipulator)))
                                                             (temp-externals externals))
  (manipulator-build-internals-from-cpinternals manipulator temp-externals cpinternals assembler-internals)
  (update-externals (assembler manipulator) assembler-internals :coords externals :internals assembler-internals)
  externals)

(defun aggregate-from-cpinternals (manipulator cpinternals &key (externals (make-coordinates-for-assembler (assembler manipulator))))
  (let* ((assembler (assembler manipulator))
         (aggregate (chem:matter-copy (topology:aggregate assembler))))
    (externals-from-cpinternals manipulator cpinternals :externals externals)
    (chem:matter/apply-coordinates aggregate externals)
    aggregate))
