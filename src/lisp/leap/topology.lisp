(in-package :leap.topology)

(defclass amber-topology ()
  ((aggregate :initarg :aggregate :accessor aggregate)
   ))


(defun find-unique-name (name-string-set prefix name-index)
  "Assemble a name from prefix and name-index and if it is not in name-string-set
then return it and the incremented name-index.  If it is in the name-string-set
keep incrementing name-index until a name is found that is NOT in name-string-set."
  (let ((*print-base* 36))
    (loop for tries below 1000
          do (let ((new-name (format nil "~a~s" prefix (incf name-index))))
               (unless (gethash new-name name-string-set)
                 (return-from find-unique-name (values new-name name-index)))))
    (error "More than 1000 names were generated and no unique names were found")))
    
(defun compress-atom-names (atom-names)
  "Construct a hash-table that maps current atom-names to 4 character atom names.
Return (values compressed-atom-name-map max-atom-name-length). "
  (let ((name-map (make-hash-table))
        (name-string-set (make-hash-table :test #'equal))
        (name-index-map (make-hash-table :test #'equal))
        (max-name-length 0)
        (pathological-case-counter 0))
    (loop for name across atom-names
          do (setf (gethash name name-map) nil))
    (loop for name being the hash-keys in name-map
          do (setf (gethash (string name) name-string-set) T))
    (loop for name being the hash-keys in name-map
          for name-string = (string name)
          for non-alpha = (position-if-not #'alpha-char-p name-string)
          for prefix = (subseq name-string 0 non-alpha)
          for name-index = (gethash prefix name-index-map 0)
          do (when (> (length (string name)) max-name-length)
               (setf max-name-length (length (string name))))
             (if (> (length name-string) 4)
                 (multiple-value-bind (unique-name new-name-index)
                     (find-unique-name name-string-set prefix name-index)
                   (when (> (length unique-name) 4)
                     ;; We have a pathological case where a name is still >4 chars
                     (multiple-value-setq (unique-name new-name-index)
                       (find-unique-name name-string-set "%" pathological-case-counter))
                     (setf pathological-case-counter new-name-index)
                     (when (> (length unique-name) 4)
                       (error "A hyper-pathological event has taken place - we could not identify a unique atom name with < 4 characters for ~a - we got ~a" name unique-name)))
                   (setf (gethash name name-map) unique-name)
                   (setf (gethash prefix name-index-map) new-name-index))
                 (setf (gethash name name-map) (string name))))
    (values name-map max-name-length)))
    
(defun collapse-stretch-parameters (kb-vec r0-vec atom1-vec atom2-vec)
  (let ((j-vec (make-array 256 :fill-pointer 0 :adjustable t))
        (jnext 0)
        jtemp
        (uniques (make-hash-table :test #'eq))
        (kbj-vec (make-array 256 :element-type 'double-float :fill-pointer 0 :adjustable t))
        (r0j-vec (make-array 256 :element-type 'double-float :fill-pointer 0 :adjustable t)))
    (loop for i from 0 below (length kb-vec)
       for atom1 = (aref atom1-vec i)
       for atom2 = (aref atom2-vec i)
       for key = (chem:canonical-stretch-key (chem:get-type atom1) (chem:get-type atom2))
       for kbi = (aref kb-vec i)
       for r0i = (aref r0-vec i)
       do (if (setf jtemp (gethash key uniques))
              (vector-push-extend jtemp j-vec)
              (progn
                (setf (gethash key uniques) jnext)
                (vector-push-extend jnext j-vec)
                (vector-push-extend kbi kbj-vec)
                (vector-push-extend r0i r0j-vec)
                (incf jnext))))
    (values j-vec kbj-vec r0j-vec)))

(defun prepare-amber-energy-stretch (energy-function)
  (let* ((energy-stretch (chem:get-stretch-component energy-function))
         (stretch-vectors (chem:extract-vectors-as-alist energy-stretch))
         (kb-vector (cdr (assoc :kb stretch-vectors)))
         (r0-vector (cdr (assoc :r0 stretch-vectors)))
         (i1-vector (cdr (assoc :i1 stretch-vectors)))
         (i2-vector (cdr (assoc :i2 stretch-vectors)))
         (atom1-vector (cdr (assoc :atom1 stretch-vectors)))
         (atom2-vector (cdr (assoc :atom2 stretch-vectors)))
         )
#||    (Format t "Number of stretch terms -> ~a~%" (length kb-vector))
    (format t "kb-vector ~a~%" kb-vector)
    (format t "r0-vector ~a~%" r0-vector)
    (format t "i1-vector ~a~%" i1-vector)
    (format t "i2-vector ~a~%" i2-vector)
    (format t "atom1-vector ~a~%" atom1-vector)
    (format t "atom2-vector ~a~%" atom2-vector)
||#
    (multiple-value-bind (j-vec kbj-vec r0j-vec)
        (collapse-stretch-parameters kb-vector r0-vector atom1-vector atom2-vector)
      (Let ((with-h 0)
            (without-h 0))
        (loop for i from 0 below (length kb-vector)
           for atom1 = (aref atom1-vector i)
           for atom2 = (aref atom2-vector i)
           do (if (or (eq (chem:get-element atom1) :h)
                      (eq (chem:get-element atom2) :h))
                  (incf with-h)
                  (incf without-h)))
        (let ((ibh (make-array with-h :element-type (array-element-type i1-vector)))
              (jbh (make-array with-h :element-type (array-element-type i2-vector)))
              (icbh (make-array with-h :element-type '(signed-byte 32)))
              (ib (make-array without-h :element-type (array-element-type i1-vector)))
              (jb (make-array without-h :element-type (array-element-type i2-vector)))
              (icb (make-array without-h :element-type '(signed-byte 32)))
              (curh 0)
              (cur 0))
          (loop for i from 0 below (length kb-vector)
             for kb = (aref kb-vector i)
             for r0 = (aref r0-vector i)
             for i1 = (aref i1-vector i)
             for i2 = (aref i2-vector i)
             for atom1 = (aref atom1-vector i)
             for atom2 = (aref atom2-vector i)
;;;             do (format t "atom1: ~a atom2: ~a i1: ~a i2: ~a~% " atom1 atom2 i1 i2)
             do (if (or (eq (chem:get-element atom1) :h)
                        (eq (chem:get-element atom2) :h))
                    (progn
                      (setf (aref ibh curh) i1)
                      (setf (aref jbh curh) i2)
                      (setf (aref icbh curh) (+ (aref j-vec i) 1))
                      (incf curh))
                    (progn
                      (setf (aref ib cur) i1)
                      (setf (aref jb cur) i2)
                      (setf (aref icb cur) (+ (aref j-vec i) 1))
                      (incf cur))))
          (values with-h without-h ibh jbh icbh ib jb icb kbj-vec r0j-vec))))))

(defun canonical-angle-key (type1 type2 type3)
  (declare (symbol type1 type2 type3))
  (intern
   (if (string< (string type1) (string type3))
       (format nil "~a-~a-~a" (string type1) (string type2) (string type3))
       (format nil "~a-~a-~a" (string type3) (string type2) (string type1)))
   :keyword))

(defun collapse-angle-parameters (kt-vec t0-vec atom1-vec atom2-vec atom3-vec)
  (let ((j-vec (make-array 256 :fill-pointer 0 :adjustable t))
        (jnext 0)
        jtemp
        (uniques (make-hash-table :test #'eq))
        (ktj-vec (make-array 256 :element-type 'double-float :fill-pointer 0 :adjustable t))
        (t0j-vec (make-array 256 :element-type 'double-float :fill-pointer 0 :adjustable t)))
    (loop for i from 0 below (length kt-vec)
       for atom1 = (aref atom1-vec i)
       for atom2 = (aref atom2-vec i)
       for atom3 = (aref atom3-vec i)
       for key = (canonical-angle-key (chem:get-type atom1) (chem:get-type atom2) (chem:get-type atom3))
       for kti = (aref kt-vec i)
       for t0i = (aref t0-vec i)
;;;       do (format t "atom1: ~a atom2: ~a atom3: ~a  key: ~a~% " atom1 atom2 atom3 key)
       do (if (setf jtemp (gethash key uniques))
              (vector-push-extend jtemp j-vec)
              (progn
                (setf (gethash key uniques) jnext)
                (vector-push-extend jnext j-vec)
                (vector-push-extend kti ktj-vec)
                (vector-push-extend t0i t0j-vec)
                (incf jnext))))
    (values j-vec ktj-vec t0j-vec)))

(defun prepare-amber-energy-angle (energy-function)
  (let* ((energy-angle (chem:get-angle-component energy-function))
         (angle-vectors (chem:extract-vectors-as-alist energy-angle))
         (kt-vector (cdr (assoc :kt angle-vectors)))
         (t0-vector (cdr (assoc :t0 angle-vectors)))
         (i1-vector (cdr (assoc :i1 angle-vectors)))
         (i2-vector (cdr (assoc :i2 angle-vectors)))
         (i3-vector ( cdr (assoc :i3 angle-vectors)))
         (atom1-vector (cdr (assoc :atom1 angle-vectors)))
         (atom2-vector (cdr (assoc :atom2 angle-vectors)))
         (atom3-vector (cdr (assoc :atom3 angle-vectors)))
         )
    (multiple-value-bind (j-vec ktj-vec t0j-vec)
        (collapse-angle-parameters kt-vector t0-vector atom1-vector atom2-vector atom3-vector)
    (Let ((with-h 0)
          (without-h 0))
      (loop for i from 0 below (length kt-vector)
         for atom1 = (aref atom1-vector i)
         for atom2 = (aref atom2-vector i)
         for atom3 = (aref atom3-vector i)
         for has-h = (or (eq (chem:get-element atom1) :h)
                         (eq (chem:get-element atom3) :h))
;;;         do (format t "Angle ~a ~a ~a  has-h -> ~a~%" atom1 atom2 atom3 has-h)
         do (if has-h
                (incf with-h)
                (incf without-h)))
      (let ((ith (make-array with-h :element-type (array-element-type i1-vector)))
            (jth (make-array with-h :element-type (array-element-type i2-vector)))
            (kth (make-array with-h :element-type (array-element-type i3-vector)))
            (icth (make-array with-h :element-type '(signed-byte 32)))
            (it (make-array without-h :element-type (array-element-type i1-vector)))
            (jt (make-array without-h :element-type (array-element-type i2-vector)))
            (kt1 (make-array without-h :element-type (array-element-type i3-vector)))
            (ict (make-array without-h :element-type '(signed-byte 32)))
            (curh 0)
            (cur 0))
        (loop for i from 0 below (length kt-vector)
           for kt = (aref kt-vector i)
           for t0 = (aref t0-vector i)
           for i1 = (aref i1-vector i)
           for i2 = (aref i2-vector i)
           for i3 = (aref i3-vector i)
           for atom1 = (aref atom1-vector i)
           for atom2 = (aref atom2-vector i)
           for atom3 = (aref atom3-vector i)
           do (if (or (eq (chem:get-element atom1) :h)
                      (eq (chem:get-element atom3) :h))
                  (progn
                    (setf (aref ith curh) i1)
                    (setf (aref jth curh) i2)
                    (setf (aref kth curh) i3)
                    (setf (aref icth curh) (+ (aref j-vec i) 1))
                    (incf curh))
                  (progn
                    (setf (aref it cur) i1)
                    (setf (aref jt cur) i2)
                    (setf (aref kt1 cur) i3)
                    (setf (aref ict cur) (+ (aref j-vec i) 1))
                    (incf cur))))
        (values with-h without-h ith jth kth icth it jt kt1 ict ktj-vec t0j-vec))))))

(defun canonical-dihedral-key (type1 type2 type3 type4 n)
  (declare (symbol type1 type2 type3 type4))
  (intern
   (if (string< (string type1) (string type4))
       (format nil "p-~a-~a-~a-~a-~d" (string type1) (string type2) (string type3) (string type4) n)
       (format nil "p-~a-~a-~a-~a-~d" (string type4) (string type3) (string type2) (string type1) n))
   :keyword))

(defun canonical-improper-key (atom1 atom2 atom3 atom4)
  (multiple-value-bind (satom1 satom2 satom4)
      (chem:improper-atom-sort atom1 atom2 atom4)
    (intern
     (format nil "i-~a-~a-~a-~a"
             (string (chem:get-type satom1))
             (string (chem:get-type satom2))
             (string (chem:get-type atom3))
             (string (chem:get-type satom4)))
     :keyword)))

(defun collapse-dihedral-parameters (v-vec in-vec phase-vec proper-vec atom1-vec atom2-vec atom3-vec atom4-vec)
  (let ((j-vec (make-array 256 :fill-pointer 0 :adjustable t))
        (jnext 0)
        jtemp
        (uniques (make-hash-table :test #'eq))
        (vj-vec (make-array 256 :element-type 'double-float :fill-pointer 0 :adjustable t))
        (inj-vec (make-array 256 :element-type '(integer 1 6) :fill-pointer 0 :adjustable t))
        (phasej-vec (make-array 256 :element-type 'double-float :fill-pointer 0 :adjustable t))
        (properj-vec (make-array 256 :fill-pointer 0 :adjustable t)))
    (loop for i from 0 below (length v-vec)
       for atom1 = (aref atom1-vec i)
       for atom2 = (aref atom2-vec i)
       for atom3 = (aref atom3-vec i)
       for atom4 = (aref atom4-vec i)
       for type1 = (chem:get-type atom1)
       for type2 = (chem:get-type atom2)
       for type3 = (chem:get-type atom3)
       for type4 = (chem:get-type atom4)
;       for in = (aref in-vec i)
       for proper = (aref proper-vec i)
;       for key = (if proper
;                     (canonical-dihedral-key type1 type2 type3 type4 in)
;                     (canonical-improper-key atom1 atom2 atom3 atom4))
       for vi = (aref v-vec i)
       for ini = (aref in-vec i)
       for phasei = (aref phase-vec i)
       for key = (intern (format nil "~15f-~d-~f" vi ini phasei) :keyword) 
;;;       do (format t "atom1: ~a atom2: ~a atom3: ~a atom4 ~a key: ~a vi: ~a ini ~a phase ~a~% " atom1 atom2 atom3 atom4 key vi ini phasei)
       do (if (setf jtemp (gethash key uniques))
              (vector-push-extend jtemp j-vec)
              (progn
                (setf (gethash key uniques) jnext)
                (vector-push-extend jnext j-vec)
                (vector-push-extend vi vj-vec)
                (vector-push-extend ini inj-vec)
                (vector-push-extend phasei phasej-vec)
                (vector-push-extend proper properj-vec)
                (incf jnext))))
    (values j-vec vj-vec inj-vec phasej-vec properj-vec)))

(defun dihedral-sort-order (x y pv iv jv kv lv nv)
  (let ((px (aref pv x))
        (py (aref pv y)))
    (cond
      ;; If x and y are both proper dihedrals then we have to apply further tests to order them
      ((and px py)
       (let ((ix (aref iv x))
             (iy (aref iv y)))
         (if (<= ix iy)
             (if (< ix iy)
                 t
                 (let ((jx (aref jv x))
                       (jy (aref jv y)))
                   (if (<= jx jy)
                       (if (< jx jy)
                           t
                           (let ((kx (aref kv x))
                                 (ky (aref kv y)))
                             (if (<= kx ky)
                                 (if (< kx ky)
                                     t
                                     (let ((lx (aref lv x))
                                           (ly (aref lv y)))
                                       (if (<= lx ly)
                                           (if (< lx ly)
                                               t
                                               (let ((nx (aref nv x))
                                                     (ny (aref nv y)))
                                                 (if (<= nx ny)
                                                     (if (< nx ny)
                                                         t
                                                         (error "Identical terms"))
                                                     nil)))
                                           nil)))
                                 nil)))
                       nil)))
             nil)))
      ;; If x is a proper dihedral and y is an improper then they are in order
      ((and px (not py)) t)
      ;; If x is an improper and y is a proper then they are not in order
      ((and (not px) py) nil)
      ;; If x and y are impropers then they are in order
      ;;   - we dont care about ordering improper dihedrals
      (t t))))

(defun check-if-dihedrals-are-ordered (p iv jv kv lv nv)
  (loop for x from 0 below (1- (length iv))
     for y = (1+ x)
     always (dihedral-sort-order x y p iv jv kv lv nv)))

(defun make-indices (num)
  (coerce (loop for i below num collect i) '(vector fixnum)))

(defun sort-dihedrals (pv iv jv kv lv nv)
  (let ((indices (make-indices (length iv))))
    (flet ((order (x y)
             (dihedral-sort-order x y pv iv jv kv lv nv)))
      (stable-sort indices #'order))))

(Defun prepare-amber-energy-dihedral (energy-function ib jb it kt)
  "Prepare the amber energy dihedral.  
Use the ib,jb (the atoms of bonds without hydrogen) and
it,kt (the terminal atoms of angles without hydrogen) and
if the terminal atoms of a proper dihedral are in one of those two lists
then don't calculate 1,4 interactions"
  (let* ((energy-dihedral (chem:get-dihedral-component energy-function))
         (dihedral-vectors (chem:extract-vectors-as-alist energy-dihedral))
         (v-vector (cdr (assoc :v dihedral-vectors)))
         (in-vector (cdr (assoc :in dihedral-vectors)))
         (phase-vector (cdr (assoc :phase dihedral-vectors)))
         (i1-vector (cdr (assoc :i1 dihedral-vectors)))
         (i2-vector (cdr (assoc :i2 dihedral-vectors)))
         (i3-vector (cdr (assoc :i3 dihedral-vectors)))
         (i4-vector (cdr (assoc :i4 dihedral-vectors)))
         (proper-vector (cdr (assoc :proper dihedral-vectors)))         
         (atom1-vector (cdr (assoc :atom1 dihedral-vectors)))
         (atom2-vector (cdr (assoc :atom2 dihedral-vectors)))
         (atom3-vector (cdr (assoc :atom3 dihedral-vectors)))
         (atom4-vector (cdr (assoc :atom4 dihedral-vectors)))
         )
    #+(or)
    (progn
      (format t "v-vector ~s~%" v-vector)
      (format t "in-vector ~s~%" in-vector)
      (format t "phase-vector ~s~%" phase-vector)
      (format t "i1-vector ~s~%" i1-vector)
      (format t "i2-vector ~s~%" i2-vector)
      (format t "i3-vector ~s~%" i3-vector)
      (format t "i4-vector ~s~%" i4-vector)
      (format t "proper-vector ~s~%" proper-vector))
    (let ((indices (if (check-if-dihedrals-are-ordered proper-vector i1-vector i2-vector i3-vector i4-vector in-vector)
                       (make-indices (length v-vector))
                       (sort-dihedrals proper-vector i1-vector i2-vector i3-vector i4-vector in-vector)))))
    (let (proper-prev i1prev i2prev i3prev i4prev)
      ;; If i3x or i4x are zero - then they can not be negated
      ;; we have to flip the dihedral order to ensure this does not happen
      (loop for x below (length v-vector)
            for properx = (aref proper-vector x)
            for i1x = (aref i1-vector x)
            for i2x = (aref i2-vector x)
            for i3x = (aref i3-vector x)
            for i4x = (aref i4-vector x)
            do (when (or (= i3x 0) (= i4x 0))
                 (psetf (aref i1-vector x) i4x
                        (aref i2-vector x) i3x
                        (aref i3-vector x) i2x
                        (aref i4-vector x) i1x)))
      ;; Modify i1-vector, i2-vector, i3-vector, i4-vector to indicate
      ;; whether 1-4 interaction should be suppressed and  (i3-vector value is negated)
      ;; whether the dihedral is improper or not  (i4-vector value is negated)
      (flet ((in-same-bond-or-angle (i1x i4x ib jb it kt)
               ;; Return T if i1x,i4x are in a bond
               (loop for bond-index below (length ib)
                     for ib-index = (aref ib bond-index)
                     for jb-index = (aref jb bond-index)
                     when (or (and (= ib-index i1x) (= jb-index i4x))
                              (and (= ib-index i4x) (= jb-index i1x)))
                       do (return-from in-same-bond-or-angle t))
               ;; Return T if i1x,i4x are in an angle
               (loop for angle-index below (length it)
                     for it-index = (aref it angle-index)
                     for kt-index = (aref kt angle-index)
                     when (or (and (= it-index i1x) (= kt-index i4x))
                              (and (= it-index i4x) (= kt-index i1x)))
                       do (return-from in-same-bond-or-angle t))
               nil))
        (loop for x below (length v-vector)
              for properx = (aref proper-vector x)
              for i1x = (aref i1-vector x)
              for i2x = (aref i2-vector x)
              for i3x = (aref i3-vector x)
              for i4x = (aref i4-vector x)
              do (cond
                   ((and properx
                         proper-prev
                         (eq i1x i1prev)
                         (eq i2x i2prev)
                         (eq i3x i3prev)
                         (eq i4x i4prev))
                    (setf (aref i3-vector x) (- i3x)))
                   (properx
                    (when (in-same-bond-or-angle i1x i4x ib jb it kt)
                      (setf (aref i3-vector x) (- i3x))))
                   (t                   ; It's an improper
                    (setf (aref i3-vector x) (- i3x)
                          (aref i4-vector x) (- i4x))))
                 (setf proper-prev properx
                       i1prev i1x
                       i2prev i2x
                       i3prev i3x
                       i4prev i4x)))
      (multiple-value-bind (j-vec vj-vec inj-vec phasej-vec properj-vec)
          (collapse-dihedral-parameters v-vector in-vector phase-vector proper-vector atom1-vector atom2-vector atom3-vector atom4-vector)
        (let ((with-h 0)
              (without-h 0))
          (loop for i from 0 below (length v-vector)
                for atom1 = (aref atom1-vector i)
                for atom2 = (aref atom2-vector i)
                for atom3 = (aref atom3-vector i)
                for atom4 = (aref atom4-vector i)
                do (if (or (eq (chem:get-element atom1) :h)
                           (eq (chem:get-element atom2) :h)
                           (eq (chem:get-element atom3) :h)
                           (eq (chem:get-element atom4) :h))
                       (incf with-h)
                       (incf without-h)))
          (let ((iph (make-array with-h :element-type (array-element-type i1-vector)))
                (jph (make-array with-h :element-type (array-element-type i2-vector)))
                (kph (make-array with-h :element-type (array-element-type i3-vector)))
                (lph (make-array with-h :element-type (array-element-type i4-vector)))
                (icph (make-array with-h :element-type '(signed-byte 32)))
                (ip (make-array without-h :element-type (array-element-type i1-vector)))
                (jp (make-array without-h :element-type (array-element-type i2-vector)))
                (kp (make-array without-h :element-type (array-element-type i3-vector)))
                (lp (make-array without-h :element-type (array-element-type i4-vector)))
                (icp (make-array without-h :element-type '(signed-byte 32)))
                (curh 0)
                (cur 0))
            (loop for i from 0 below (length v-vector)
                  for v = (aref v-vector i)
                  for in = (aref in-vector i)
                  for i1 = (aref i1-vector i)
                  for i2 = (aref i2-vector i)
                  for i3 = (aref i3-vector i)
                  for i4 = (aref i4-vector i)
                  for atom1 = (aref atom1-vector i)
                  for atom2 = (aref atom2-vector i)
                  for atom3 = (aref atom3-vector i)
                  for atom4 = (aref atom4-vector i)
                  do (if (or (eq (chem:get-element atom1) :h)
                             (eq (chem:get-element atom2) :h)
                             (eq (chem:get-element atom3) :h)
                             (eq (chem:get-element atom4) :h))
                         (progn
                           (setf (aref iph curh) i1)
                           (setf (aref jph curh) i2)
                           (setf (aref kph curh) i3)
                           (setf (aref lph curh) i4)
                           (setf (aref icph curh) (+ (aref j-vec i) 1))
                           (incf curh))
                         (progn
                           (setf (aref ip cur) i1)
                           (setf (aref jp cur) i2)
                           (setf (aref kp cur) i3)
                           (setf (aref lp cur) i4)
                           (setf (aref icp cur) (+ (aref j-vec i) 1))
                           (incf cur))))
            (values with-h without-h iph jph kph lph icph ip jp kp lp icp vj-vec inj-vec phasej-vec properj-vec)))))))

(defun canonical-nonbond-key (type1 type2)
  (declare (symbol type1 type2))
  (intern (format nil "~a-~a" (string type1) (string type2))
          :keyword))

(defun generate-nonbond-parameters (ffnonbond-db type-index-vector)
  (let ((iac-vec (make-array 256 :element-type '(signed-byte 32) :fill-pointer 0 :adjustable t))
        (jnext 1)
        jtemp
        (uniques (make-hash-table :test #'eq))
        (type-indexj-vec (make-array 256 :element-type '(signed-byte 32) :fill-pointer 0 :adjustable t))
        (local-typej-vec (make-array 256 :element-type '(signed-byte 32) :fill-pointer 0 :adjustable t))
        (ntypes 0))
    (loop for i from 0 below (length type-index-vector)
       for index = (aref type-index-vector i)
       for ffnonbond = (chem:get-ffnonbond-using-type-index ffnonbond-db index)
       for epsilon = (chem:get-epsilon-k-cal ffnonbond)
       for rm =  (chem:get-radius-angstroms ffnonbond)
       for key = (intern (format nil "~f-~f" epsilon rm) :keyword)
       do (if (setf jtemp (gethash key uniques))        
;       do (if (setf jtemp (gethash index uniques))
              (vector-push-extend jtemp iac-vec)
              (progn
;                (setf (gethash index uniques) jnext)
                (setf (gethash key uniques) jnext)
                (vector-push-extend jnext iac-vec)
                (vector-push-extend index type-indexj-vec)
                (vector-push-extend jnext local-typej-vec)
                (incf jnext))))
    (setf ntypes (length local-typej-vec))
    (let ((ico-vec (make-array (* ntypes ntypes) :element-type '(signed-byte 32)))
;          (ffnonbond1 (chem:get-ffnonbond-using-type-index ffnonbond-db type-index1))
          (cn1-vec (make-array (/ (* ntypes (+ ntypes 1)) 2) :element-type 'double-float)) 
          (cn2-vec (make-array (/ (* ntypes (+ ntypes 1)) 2) :element-type 'double-float))
          (type1 0)
          (type2 0)
          (rstar 0)
          (epsilonij 0))
;;      (format *debug-io* "In generate-nonbond-parameters type-indexj-vec -> ~s~%" type-indexj-vec)
      (loop for i from 0 below ntypes
         for type1 = (aref type-indexj-vec i)
         do (loop for j from 0 below ntypes
               for type2 = (aref type-indexj-vec j)
               do (when (<= i j)
;;;                    (format *debug-io* "In generate-nonbond-parameters i->~a j->~a type1->~a type2->~a~%" i j type1 type2)
                    (let ((ffnonbond1 (chem:get-ffnonbond-using-type-index ffnonbond-db type1))
                          (ffnonbond2 (chem:get-ffnonbond-using-type-index ffnonbond-db type2)))
                      (setf rstar (+ (chem:get-radius-angstroms ffnonbond1) (chem:get-radius-angstroms ffnonbond2))
                            epsilonij (sqrt (* (chem:get-epsilon-k-cal ffnonbond1) (chem:get-epsilon-k-cal ffnonbond2)))
                            (aref cn1-vec (- (+ (/ (* (+ j 1) j) 2) (+ i 1)) 1)) (* epsilonij (expt rstar 12.0))
                            (aref cn2-vec (- (+ (/ (* (+ j 1) j) 2) (+ i 1)) 1)) (* 2.0 epsilonij (expt rstar 6.0))
                            (aref ico-vec (+ (* ntypes i) j)) (+ (/ (* (+ j 1) j) 2) (+ i 1)))
                      (if (< i j)
                          (setf (aref ico-vec (+ (* ntypes j) i)) (+ (/ (* (+ j 1) j) 2) (+ i 1))))))))
;;;                      (format t "type1 ~a type2 ~a~%" type1 type2)    
      (values ntypes ico-vec iac-vec local-typej-vec cn1-vec cn2-vec))))

(defun chem:prepare-amber-energy-nonbond (energy-function ffnonbond-db)
  (let* ((atom-table (chem:atom-table energy-function))
         (natom (chem:get-number-of-atoms atom-table))
         (atom-name-vector (make-array natom))
         (atom-type-vector (make-array natom))
         (charge-vector (make-array natom :element-type 'double-float))
         (mass-vector (make-array natom :element-type 'double-float))
         (type-index-vector (make-array natom :element-type '(signed-byte 32)))
         (atomic-number-vector (make-array natom :element-type '(signed-byte 32)))
         (atom-radius-vector (make-array natom :element-type 'double-float))
         (energy-nonbond (chem:get-nonbond-component energy-function)))
;;;    (format t "energy-nonbond -> ~a~%" energy-nonbond)
;;;    (format t "ffnonbond-db -> ~a~%" ffnonbond-db)
    (loop for i from 0 below natom
       for atom-name = (chem:elt-atom-name atom-table i)
       for atom-type = (chem:elt-atom-type atom-table i)
       for charge = (chem:elt-charge atom-table i)
       ;;for mass = (chem:elt-mass atom-table i)
       for type-index = (chem:elt-type-index atom-table i)
       for atomic-number = (chem:elt-atomic-number atom-table i)
       for ffnonbondi = (chem:get-ffnonbond-using-type-index ffnonbond-db type-index)
       for atom-radius = (chem:get-radius-angstroms ffnonbondi)
       for mass = (chem:get-mass ffnonbondi)
       do (setf (aref atom-name-vector i) atom-name
                (aref atom-type-vector i) atom-type
                (aref charge-vector i) charge
                (aref mass-vector i) mass
                (aref type-index-vector i) type-index
                (aref atomic-number-vector i) atomic-number
                (aref atom-radius-vector i) atom-radius))
;    (format t "atom-name-vector -> ~s~%" atom-name-vector)
;    (format t "atom-type-vector -> ~s~%" atom-type-vector)
;    (format t "charge-vector -> ~s~%" charge-vector)
;    (format t "type-index-vector -> ~s~%" type-index-vector)
                                        ;    (format t "atom-radius-vector -> ~s~%" atom-radius-vector)
    (multiple-value-bind (ntypes ico-vec iac-vec local-typej-vec cn1-vec cn2-vec)
        (generate-nonbond-parameters ffnonbond-db type-index-vector)
      (list (cons :ntypes ntypes)
            (cons :atom-name-vector (copy-seq atom-name-vector))
            (cons :atom-type-vector (copy-seq atom-type-vector))
            (cons :charge-vector (copy-seq charge-vector))
            (cons :mass-vector (copy-seq mass-vector))
            (cons :atomic-number-vector (copy-seq atomic-number-vector))
            (cons :atom-radius-vector (copy-seq atom-radius-vector))
            (cons :ico-vec (copy-seq ico-vec))
            (cons :iac-vec (copy-seq iac-vec))
            (cons :local-typej-vec (copy-seq local-typej-vec))
            (cons :cn1-vec (copy-seq cn1-vec))
            (cons :cn2-vec (copy-seq cn2-vec))))))


(defun prepare-residue (energy-function)
  (let* ((atom-table (chem:atom-table energy-function))
         (residue-vector (chem:atom-table-residues atom-table))
         (residue-name-vector (chem:atom-table-residue-names atom-table))
         (molecule-vector (chem:atom-table-atoms-per-molecule atom-table))
         (residue-pointer-vector (make-array 256 :element-type '(signed-byte 32) :fill-pointer 0 :adjustable t))
         (atoms-per-molecule (make-array 256 :element-type '(signed-byte 32) :fill-pointer 0 :adjustable t))
         (nmolecule 0)
         (nresidue 0)
         (nmxrs 0))
    (setf nmolecule (- (length molecule-vector) 1))
    (setf nresidue (length residue-name-vector))
    (loop for i from 0 below (length residue-vector)
       do (if (= i 0)
              (vector-push-extend (+ (aref residue-vector 0) 1) residue-pointer-vector)
              (if (/= (aref residue-vector (- i 1)) (aref residue-vector i))
                  (vector-push-extend (+ (aref residue-vector i) 1) residue-pointer-vector))))
    ;;      for fresidue = (+ (aref residue-vector i) 1)
    ;;      do (setf (aref residue-vector i) fresidue))
    (loop for i from 0 below nresidue
       for inmxrs = (- (aref residue-vector (+ i 1)) (aref residue-vector i))
       do (if (> inmxrs nmxrs)
              (setf nmxrs inmxrs)))
    (loop for i from 0 below nmolecule
       for natom = (- (aref molecule-vector (+ i 1)) (aref molecule-vector i))
       do (vector-push-extend natom atoms-per-molecule))
    (values nresidue nmxrs residue-pointer-vector residue-name-vector atoms-per-molecule)))

;; for now, hardwire the Bondi radii
(defun prepare-generalized-born (atomic-number-vector)
  (let ((generalized-born-radius (make-array (length atomic-number-vector)))
        (generalized-born-screen (make-array (length atomic-number-vector)))
        (atomic-number 0)
        (radius 0)
        (screen 0))
    (loop for i from 0 below (length atomic-number-vector)
       for atomic-number = (aref atomic-number-vector i)
       do (case atomic-number (1 (setf radius 1.3))
                (6 (setf radius 1.7))
                (7 (setf radius 1.55))
                (8 (setf radius 1.5))
                (9 (setf radius 1.5))
                (14 (setf radius 2.1))
                (15 (setf radius 1.85))
                (16 (setf radius 1.8))
                (17 (setf radius 1.7))
                (otherwise  (setf radius 1.5)))
       do (setf (aref generalized-born-radius i) radius)
       do (case atomic-number (1 (setf screen 0.85))
                (6 (setf screen 0.72))
                (7 (setf screen 0.79))
                (8 (setf screen 0.85))
                (9 (setf screen 0.88))
                (15 (setf screen 0.86))
                (16 (setf screen 0.96))
                (otherwise  (setf screen 0.8)))
       do (setf (aref generalized-born-screen i) screen))
    (values generalized-born-radius generalized-born-screen)))

(defun solvent-pointers (aggregate)
  (let ((residue-count 0)
        (molecule-count 0))
    (chem:map-residues
     nil
     (lambda (r)
       (unless (string= (chem:get-name r) "WAT")
                (incf residue-count)))
     aggregate)
    (chem:map-molecules
     nil
     (lambda (r)
       (incf molecule-count))
     aggregate)
    (values residue-count molecule-count)
    ))

(defun save-amber-parm-format (aggregate topology-pathname coordinate-pathname &key system assign-types)
  (let* ((energy-function (chem:make-energy-function aggregate system
                                                     :use-excluded-atoms t
                                                     :assign-types assign-types))
         (nonbonds (chem:get-nonbond-component energy-function))
         (number-excluded-atoms (chem:number-excluded-atoms nonbonds))
         (excluded-atom-list (chem:excluded-atom-list nonbonds))
         (topology-pathname (merge-pathnames topology-pathname))
         (coordinate-pathname (merge-pathnames coordinate-pathname))
         (natom (chem:get-number-of-atoms (chem:atom-table energy-function))))
    ;; Skip assigning MarkMainChainAtoms and MarkSideChain atoms for now
    ;; see (unitio.c:4889).  This won't mean anything for spiroligomers.
    (fortran:with-fortran-output-file (ftop topology-pathname :direction :output)
      (fortran:debug "-1-")             ;
      (fortran:fformat 1 "%-80s")
      ;;                  (fortran:fwrite (core:strftime 81 "%%VERSION  VERSION_STAMP = V0002.000  DATE = %m/%d/%y  %H:%M:%S"))
;;; temporary!!!
      (multiple-value-bind (second minute hour date month year)
          (get-decoded-time)
        (fortran:fwrite (format
                         nil
                         "%VERSION  VERSION_STAMP = V0001.000  DATE = ~2,'0d/~2,'0d/~2,'0d  ~2,'0d:~2,'0d:~2,'0d"
                         month date (- year 2000) hour minute second)))
      (fortran:fwrite "%FLAG TITLE")
      (fortran:fwrite "%FORMAT(20a4)")
      (fortran:fwrite (string (chem:get-name aggregate)))
      ;; This function will be very, very long (for Common Lisp)
      ;; To avoid lots of nested scopes we will declare one large scope
      ;;   and declare all of the variables in that scope here at the top
      (let (atom-vectors
            ntypes atom-name atom-type charge mass atomic-number atom-radius ico iac local-typej-vec cn1-vec cn2-vec #|nonbonds|#
            nbonh mbona ibh jbh icbh ib jb icb kbj-vec r0j-vec #|stretches|#
            ntheth mtheta ith jth kth icth it jt kt1 ict ktj-vec t0j-vec #|angles|#
            nphih mphia iph jph kph lph icph ip jp kp lp icp vj-vec inj-vec phasej-vec properj-vec #|dihedrals|#
            nhparm nparm  nnb nres residue-vec residue-name-vec atoms-per-molecule
            generalized-born-radius generalized-born-screen
            residue-count molecule-count
            nbona    ntheta nphia  NUMBND NUMANG NPTRA
            NATYP    NPHB   IFPERT NBPER  NGPER  NDPER
            MBPER    MGPER  MDPER  IFBOX  nmxrs  IFCAP
            NUMEXTRA NCOPY
            )
        ;; Here we need to calculate all of the values for %FLAG POINTERS
        (multiple-value-setq (nbonh mbona ibh jbh icbh ib jb icb kbj-vec r0j-vec)
          (prepare-amber-energy-stretch energy-function))
        (multiple-value-setq (ntheth mtheta ith jth kth icth it jt kt1 ict ktj-vec t0j-vec)
          (prepare-amber-energy-angle energy-function))
        (multiple-value-setq (nphih mphia iph jph kph lph icph ip jp kp lp icp vj-vec inj-vec phasej-vec properj-vec)
          (prepare-amber-energy-dihedral energy-function ib jb it kt1))
                                        ;        (multiple-value-setq (ntypes atom-name charge mass atomic-number ico iac local-typej-vec cn1-vec cn2-vec)
                                        ;          (chem:prepare-amber-energy-nonbond energy-function))
        (multiple-value-setq (nres nmxrs residue-vec residue-name-vec atoms-per-molecule)
          (prepare-residue energy-function))
        (multiple-value-setq (residue-count molecule-count)
          (solvent-pointers aggregate))
        (setf atom-vectors (chem:prepare-amber-energy-nonbond
                            energy-function
                            (chem:lookup-nonbond-force-field-for-aggregate aggregate system)))
        (setf ntypes (cdr (assoc :ntypes atom-vectors)))
        (setf atom-name (cdr (assoc :atom-name-vector atom-vectors)))
        (setf atom-type (cdr (assoc :atom-type-vector atom-vectors)))
        (setf charge (cdr (assoc :charge-vector atom-vectors)))
        (setf mass (cdr (assoc :mass-vector atom-vectors)))
        (setf atomic-number (cdr (assoc :atomic-number-vector atom-vectors)))
        (setf atom-radius (cdr (assoc :atom-radius-vector atom-vectors)))
        (setf ico (cdr (assoc :ico-vec atom-vectors)))
        (setf iac (cdr (assoc :iac-vec atom-vectors)))
        (setf local-typej-vec (cdr (assoc :local-typej-vec atom-vectors)))
        (setf cn1-vec (cdr (assoc :cn1-vec atom-vectors)))
        (setf cn2-vec (cdr (assoc :cn2-vec atom-vectors)))
        (multiple-value-setq (generalized-born-radius generalized-born-screen)
          (prepare-generalized-born atomic-number)) 
                                        ;        (setf generalized-born-screen (prepare-generalized-born atomic-number))
        (setf nhparm 0)
        (setf nparm 0)
        (setf nnb (length excluded-atom-list))
        (setf nbona mbona)
        (setf ntheta mtheta)
        (setf nphia mphia)
        (setf numbnd (length kbj-vec))
        (setf numang (length ktj-vec))
        (setf nptra (length vj-vec))
        (setf natyp (length (remove-duplicates atom-type)))
        (setf nphb 0)
        (setf ifpert 0)    
        (setf nbper 0)     
        (setf ngper 0)     
        (setf ndper 0)     
        (setf mbper 0)     
        (setf mgper 0)    
        (setf mdper 0)
        (if (chem:has-property aggregate :bounding-box)
            (setf ifbox 1)
            (setf ifbox 0))
        (setf ifcap 0)     
        (setf numextra 0)  
        (setf ncopy 0)
                                        ;       (setf atom-type-name (chem:get-type atoms))
 
        #| dihedrals, nonbonds, others??? |#
        ;; --- Done calculating all of the values
        ;; --- Now write out all of the values
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG POINTERS")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-2-")
        (fortran:fformat 10 "%8d")
        ;; NATOM
        (fortran:fwrite natom)
        (fortran:fwrite ntypes)
        (fortran:fwrite nbonh)    ;number of bonds containing hydrogen
        (fortran:fwrite mbona) ;number of bonds not containing hydrogen
        (fortran:fwrite ntheth) ; number of angles containing hydrogen
        (fortran:fwrite mtheta) ; number of angles not containing hydrogen
        (fortran:fwrite nphih) ; number of dihedrals containing hydrogen
        (fortran:fwrite mphia) ; number of dihedrals not containing hydrogen
        (fortran:fwrite nhparm)         ; currently not used
        (fortran:fwrite nparm) ; used to determine if addles created prmtop
        (fortran:fwrite nnb)   ; number of excluded atoms
        (fortran:fwrite nres)  ; number of residues
        (fortran:fwrite nbona) ; MBONA + number of constraint bonds
        (fortran:fwrite ntheta) ; MTHETA + number of constraint angles
        (fortran:fwrite nphia) ; MPHIA + number of constraint dihedrals
        (fortran:fwrite numbnd)      ; number of unique bond types
        (fortran:fwrite numang)      ; number of unique angle types
        (fortran:fwrite nptra)       ; number of unique dihedral types
        (fortran:fwrite natyp) ; number of atom types in parameter file, see SOLTY below
        (fortran:fwrite nphb) ; number of distinct 10-12 hydrogen bond pair types
        (fortran:fwrite ifpert) ; set to 1 if perturbation info is to be read in
        (fortran:fwrite nbper)  ; number of bonds to be perturbed
        (fortran:fwrite ngper)  ; number of angles to be perturbed
        (fortran:fwrite ndper)  ; number of dihedrals to be perturbed
        (fortran:fwrite mbper) ; number of bonds with atoms completely in perturbed group
        (fortran:fwrite mgper) ; number of angles with atoms completely in perturbed group
        (fortran:fwrite mdper) ; number of dihedrals with atoms completely in perturbed groups
        (fortran:fwrite ifbox) ; set to 1 if standard periodic box, 2 when truncated octahedral
        (fortran:fwrite nmxrs) ; number of atoms in the largest residue
        (fortran:fwrite ifcap) ; set to 1 if the CAP option from edit was specified
        (fortran:fwrite numextra) ; number of extra points found in topology
        (fortran:fwrite ncopy) ; number of PIMD slices / number of beads
        (fortran:end-line)

        ;; Next) 
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG ATOM_NAME")
        (fortran:fwrite "%FORMAT(20a4)")
        (fortran:debug "-3-")
        (fortran:fformat 20 "%-4s")
        (multiple-value-bind (compressed-atom-names max-name-length)
            (compress-atom-names atom-name)
          (loop for name across atom-name
                for compressed-name = (gethash name compressed-atom-names)
                when (> (length compressed-name) 4)
                  do (error "There is an atom name ~a that cannot be made unique and less than 4 characters we generated ~a !!!" name compressed-name)
                do (fortran:fwrite compressed-name)))
        (fortran:end-line)
        ;; write the atom names

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG CHARGE")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-4-")
        (fortran:fformat 5 "%16.8e")
        (loop for ch across charge
              do (fortran:fwrite (* ch 18.2223)))
        (fortran:end-line)
        ;; write the atom charges

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG ATOMIC_NUMBER")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-5-")
        (fortran:fformat 10 "%8d")
        (loop for number across atomic-number
              do (fortran:fwrite number))
        (fortran:end-line)
        ;; write the atomic number of each atom

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG MASS")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-6-")
        (fortran:fformat 5 "%16.8e")
        (loop for ma across mass
              do (fortran:fwrite ma))
        (fortran:end-line)
        ;; write the atom masses

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG ATOM_TYPE_INDEX")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-7-")
        (fortran:fformat 10 "%8d")
        (loop for ia across iac
              do (fortran:fwrite ia))
        (fortran:end-line)
        ;; write the index fot the atom types

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG NUMBER_EXCLUDED_ATOMS")
        (fortran:end-line)
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-8-")
        (fortran:fformat 10 "%8d")
        (loop for na across number-excluded-atoms
              do (fortran:fwrite na))
        (fortran:end-line)
        ;; write the total number of excluded atoms for atom "i"

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG NONBONDED_PARM_INDEX")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-9-")
        (fortran:fformat 10 "%8d")
        (loop for ic across ico
              do (fortran:fwrite ic))
        (fortran:end-line)
        ;; provides the index to the nobon parameter arrays CN1, CN2 and ASOL, BSOL.

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG RESIDUE_LABEL")
        (fortran:fwrite "%FORMAT(20A4)")
        (fortran:debug "-10-")
        (fortran:fformat 20 "%-4s")
        (loop for ren across residue-name-vec
              do (fortran:fwrite (string ren)))
        (fortran:end-line)
        ;; write the name of each of the residues

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG RESIDUE_POINTER")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-11-")
        (fortran:fformat 10 "%8d")
        (loop for re from 0 below (- (length residue-vec) 1)
              do (fortran:fwrite (aref residue-vec re)))
        (fortran:end-line)
        ;; write the atoms in each residue are listed for atom "1" in IPRES(i) to IPRES(i+1)-1

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG BOND_FORCE_CONSTANT")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-12-")
        (fortran:fformat 5 "%16.8e")
        (loop for kb across kbj-vec
              do (fortran:fwrite kb))
        (fortran:end-line)
        
        ;; write the force constant for the bonds of each type

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG BOND_EQUIL_VALUE")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-13-")
        (fortran:fformat 5 "%16.8e")
        (loop for r0 across r0j-vec
              do (fortran:fwrite r0))
        (fortran:end-line)
        ;; write the equilibrium bond length for the bonds of each type

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG ANGLE_FORCE_CONSTANT")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-14-")
        (fortran:fformat 5 "%16.8e")
        (loop for kt across ktj-vec
              do (fortran:fwrite kt))
        (fortran:end-line)
        ;; write the force constant for the angles of each type

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG ANGLE_EQUIL_VALUE")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-15-")
        (fortran:fformat 5 "%16.8e")
        (loop for t0 across t0j-vec
              do (fortran:fwrite t0))
        (fortran:end-line)
        ;; write the equilibrium angle for the angles of each type

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG DIHEDRAL_FORCE_CONSTANT")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-16-")
        (fortran:fformat 5 "%16.8e")
        (loop for v0 across vj-vec
              do (fortran:fwrite v0))
        (fortran:end-line)
        ;; write the force constant for the dihedral of a given type

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG DIHEDRAL_PERIODICITY")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-17-")
        (fortran:fformat 5 "%16.8e")
        (loop for in0 across inj-vec
              do (fortran:fwrite (float in0)))
        (fortran:end-line)
        ;; write the periodicity of the dihedral of a given type

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG DIHEDRAL_PHASE")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-18-")
        (fortran:fformat 5 "%16.8e")
        (loop for p0 across phasej-vec
              do (fortran:fwrite p0))
        (fortran:end-line)
        ;; write the phase of the dihedral of a given type

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG SCEE_SCALE_FACTOR")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-19-")
        (fortran:fformat 5 "%16.8e")
        (loop for pr0 across properj-vec
           do (if pr0
                  (fortran:fwrite 1.2)
                  (fortran:fwrite 0.0)))
        (fortran:end-line)
        ;; write the 1-4 electrostatic scaling constant

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG SCNB_SCALE_FACTOR")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-20-")
        (fortran:fformat 5 "%16.8e")
        (loop for pr0 across properj-vec
           do (if pr0
                  (fortran:fwrite 2.0)
                  (fortran:fwrite 0.0)))
        (fortran:end-line)
        ;; write the 1-4 vdw scaling constant

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG SOLTY")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-21-")
        (fortran:fformat 5 "%16.8e")
        (loop repeat natyp
              do (fortran:fwrite 0.0))
        (fortran:end-line)
        ;; currently unused

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG LENNARD_JONES_ACOEF")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-22-")
        (fortran:fformat 5 "%16.8e")
        (loop for cn1 across cn1-vec
              do (fortran:fwrite cn1))
        (fortran:end-line)
        ;; write the Lennard Jones r**12 terms for all possible atom type interactions

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG LENNARD_JONES_BCOEF")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-23-")
        (fortran:fformat 5 "%16.8e")
        (loop for cn2 across cn2-vec
              do (fortran:fwrite cn2))
        (fortran:end-line)
        ;; write the Lennard Jones r**6 terms for all possible atom type interactions

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG BONDS_INC_HYDROGEN")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-24-")
        (fortran:fformat 10 "%8d")
        (loop for i below (length ibh)
              for ibhi = (aref ibh i)
              for jbhi = (aref jbh i)
              for icbhi = (aref icbh i)
              do (fortran:fwrite ibhi)
              do (fortran:fwrite jbhi)
              do (fortran:fwrite icbhi))
        (fortran:end-line)
        
        ;; write IBH, JBH, ICBH

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG BONDS_WITHOUT_HYDROGEN")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-25-")
        (fortran:fformat 10 "%8d")
        (loop for i below (length ib)
              for ibi = (aref ib i)
              for jbi = (aref jb i)
              for icbi = (aref icb i)
              do (fortran:fwrite ibi)
              do (fortran:fwrite jbi)
              do (fortran:fwrite icbi))
        (fortran:end-line)
        ;; write IB, JB, ICB

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG ANGLES_INC_HYDROGEN")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-26-")
        (fortran:fformat 10 "%8d")
        (loop for i below (length ith)
              for ithi = (aref ith i)
              for jthi = (aref jth i)
              for kthi = (aref kth i)
              for icthi = (aref icth i)
              do (fortran:fwrite ithi)
              do (fortran:fwrite jthi)
              do (fortran:fwrite kthi)
              do (fortran:fwrite icthi))
        (fortran:end-line)
        ;; write ITH, JTH, KTH, ICTH

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG ANGLES_WITHOUT_HYDROGEN")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-27-")
        (fortran:fformat 10 "%8d")
        (loop for i below (length it)
              for iti = (aref it i)
              for jti = (aref jt i)
              for kti = (aref kt1 i)
              for icti = (aref ict i)
              do (fortran:fwrite iti)
              do (fortran:fwrite jti)
              do (fortran:fwrite kti)
              do (fortran:fwrite icti))
        (fortran:end-line)
        ;; write IT, JT, KT, ICT

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG DIHEDRALS_INC_HYDROGEN")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-28-")
        (fortran:fformat 10 "%8d")
        (loop for i below (length iph)
              for iphi = (aref iph i)
              for jphi = (aref jph i)
              for kphi = (aref kph i)
              for lphi = (aref lph i)
              for icphi = (aref icph i)
              do (fortran:fwrite iphi)
              do (fortran:fwrite jphi)
              do (fortran:fwrite kphi)
              do (fortran:fwrite lphi)
              do (fortran:fwrite icphi))
        (fortran:end-line)
        ;; write IPH, JPH, KPH, LPH, ICPH

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG DIHEDRALS_WITHOUT_HYDROGEN")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-29-")
        (fortran:fformat 10 "%8d")
        (loop for i below (length ip)
              for ipi = (aref ip i)
              for jpi = (aref jp i)
              for kpi = (aref kp i)
              for lpi = (aref lp i)
              for icpi = (aref icp i)
              do (fortran:fwrite ipi)
              do (fortran:fwrite jpi)
              do (fortran:fwrite kpi)
              do (fortran:fwrite lpi)
              do (fortran:fwrite icpi))
        (fortran:end-line)
        ;; write IP, JP, KP, LP, ICP

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG EXCLUDED_ATOMS_LIST")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-30-")
        (fortran:fformat 10 "%8d")
        (loop for atom across excluded-atom-list
              do (fortran:fwrite (+ atom 1)))
        (fortran:end-line)
        ;; write excluded atoms list

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG HBOND_ACOEF")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-31-")
        (fortran:fformat 5 "%16.8e")
        (fortran:end-line)
        ;;This term has been dropped from most modern force fields.

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG HBOND_BCOEF")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-32-")
        (fortran:fformat 5 "%16.8e")
        (fortran:end-line)
        ;;This term has been dropped from most modern force fields.

        ;; Next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG HBCUT")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-33-")
        (fortran:fformat 5 "%16.8e")
        (fortran:end-line)
        ;;no longer used for anything.

        ;;next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG AMBER_ATOM_TYPE")
        (fortran:fwrite "%FORMAT(20A4)")
        (fortran:debug "-34-")
        (fortran:fformat 20 "%-4s")
        (loop for type across atom-type
              do (fortran:fwrite (string type)))
        (fortran:end-line)

        ;;next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG TREE_CHAIN_CLASSIFICATION")
        (fortran:fwrite "%FORMAT(20A4)")
        (fortran:debug "-35-")
        (fortran:fformat 20 "%-4s")
        (loop repeat natom
              do (fortran:fwrite "M"))
        (fortran:end-line)
        ;; We are not considering protein/DNA thing, so just put "main chain" for all atoms.

        ;;next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG JOIN_ARRAY")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-36-")
        (fortran:fformat 10 "%8d")
        (loop repeat natom
              do (fortran:fwrite "0"))
        (fortran:end-line)
        ;;This section is no longer used and is currently just filled with zeros.

        ;;next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG IROTAT")
        (fortran:fwrite "%FORMAT(10I8)")
        (fortran:debug "-37-")
        (fortran:fformat 10 "%8d")
        (loop repeat natom
              do (fortran:fwrite "0"))
        (fortran:end-line)
        ;;This section is not used and is currently just filled with zeros.

        ;;next
        (if (chem:has-property aggregate :bounding-box)
            (progn
              (fortran:fformat 1 "%-80s")
              (fortran:fwrite "%FLAG SOLVENT_POINTERS")
              (fortran:fwrite "%FORMAT(3I8)")
              (fortran:debug "-38-")
              (fortran:fformat 3 "%8d")
;;;temporary!!!
              (fortran:fwrite residue-count) 
              (fortran:fwrite molecule-count) 
              (fortran:fwrite "2") 
              (fortran:end-line)
              (fortran:fformat 1 "%-80s")
              (fortran:fwrite "%FLAG ATOMS_PER_MOLECULE")
              (fortran:fwrite "%FORMAT(10I8)")
              (fortran:debug "-39-")
              (fortran:fformat 10 "%8d")
              (loop for natom across atoms-per-molecule
                    do (fortran:fwrite natom))
              (fortran:end-line)
              ;; number of atoms per molecule
              
              ;;next
              (fortran:fformat 1 "%-80s")
              (fortran:fwrite "%FLAG BOX_DIMENSIONS")
              (fortran:fwrite "%FORMAT(5E16.8)")
              (fortran:debug "-40-")
              (fortran:fformat 5 "%16.8e")
              (let ((solvent-box (chem:matter-get-property aggregate :bounding-box)))
                (unless (and solvent-box (listp solvent-box) (= (length solvent-box) 3))
                  (error "There must be a solvent-box property in the aggregate properties and it must be a list of length three numbers"))
                (fortran:fwrite "90.0000000") ;box angle
                (fortran:fwrite (float (first solvent-box)))
                (fortran:fwrite (float (second solvent-box)))
                (fortran:fwrite (float (third solvent-box)))) 
              (fortran:end-line)))

        ;;next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG RADIUS_SET")
        (fortran:fwrite "%FORMAT(1a80)")
        (fortran:debug "-41-")
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "modified Bondi radii (mbondi)") ;default in leap
        (fortran:end-line)

        ;;next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG RADII")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-42-")
        (fortran:fformat 5 "%16.8e")
        (loop for radius across generalized-born-radius
              do (fortran:fwrite radius))
        (fortran:end-line)
        ;;Generalized Born intrinsic dielectric radii

        ;;next
        (fortran:fformat 1 "%-80s")
        (fortran:fwrite "%FLAG SCREEN")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-43-")
        (fortran:fformat 5 "%16.8e")
        (loop for screen across generalized-born-screen
              do (fortran:fwrite screen))
        (fortran:end-line)
        ))
;;;    (format *debug-io* "coordinate-pathname -> ~s~%" coordinate-pathname)
    (fortran:with-fortran-output-file (ftop coordinate-pathname :direction :output :if-exists :supersede)
      (fortran:fformat 20 "%-4s")
      (fortran:fwrite (string (chem:get-name aggregate)))
      (fortran:end-line)
                                        ;      (fortran:fformat 1 "%5d")
                                        ;      (fortran:fwrite natom)
                                        ;      (fortran:fformat 5 "%15.7lf")
                                        ;      (fortran:fwrite 0.0)
                                        ;      (fortran:fwrite 0.0)
      (fortran:fwrite (format nil  " ~5d~%" natom))
      (fortran:fformat 6 "%12.7lf")
      (let ((atom-table (chem:atom-table energy-function)))
        (loop for i from 0 below natom
              for atom = (chem:elt-atom atom-table i)
              for atom-coordinate-index-times3 = (chem:elt-atom-coordinate-index-times3 atom-table i)
              for pos = (chem:get-position atom)
;;;           do (format *debug-io* "atom-coordinate-index-times3 -> ~a~%" atom-coordinate-index-times3)
              do (progn
                   (fortran:fwrite (geom:vx pos))
                   (fortran:fwrite (geom:vy pos))
                   (fortran:fwrite (geom:vz pos)))))
      ;; write out the solvent box
      (if (chem:has-property aggregate :bounding-box)
          (let ((solvent-box (chem:matter-get-property aggregate :bounding-box)))
            (unless (and solvent-box (listp solvent-box) (= (length solvent-box) 3))
              (error "There must be a solvent-box property in the aggregate properties and it must be a list of length three numbers"))
            (fortran:fwrite (float (first solvent-box)))
            (fortran:fwrite (float (second solvent-box)))
            (fortran:fwrite (float (third solvent-box)))
            (fortran:fwrite "90.0000000") ;box angle
            (fortran:fwrite "90.0000000") 
            (fortran:fwrite "90.0000000")))
      (fortran:end-line)))
  t)

(defvar %flag-title "%FLAG TITLE")
(defvar %flag-pointers "%FLAG POINTERS")
(defvar %flag-atom-name "%FLAG ATOM_NAME")
(defvar %flag-charge "%FLAG CHARGE")
(defvar %flag-atomic-number "%FLAG ATOMIC_NUMBER")
(defvar %flag-mass "%FLAG MASS")
(defvar %flag-atom-type-index "%FLAG ATOM_TYPE_INDEX")
(defvar %flag-number-excluded-atoms "%FLAG NUMBER_EXCLUDED_ATOMS")
(defvar %flag-nonbonded-parm-index "%FLAG NONBONDED_PARM_INDEX")
(defvar %flag-residue-label "%FLAG RESIDUE_LABEL")
(defvar %flag-residue-pointer "%FLAG RESIDUE_POINTER")
(defvar %flag-bond-force-constant "%FLAG BOND_FORCE_CONSTANT")
(defvar %flag-bond-equil-value "%FLAG BOND_EQUIL_VALUE")
(defvar %flag-angle-force-constant "%FLAG ANGLE_FORCE_CONSTANT")
(defvar %flag-angle-equil-value "%FLAG ANGLE_EQUIL_VALUE")
(defvar %flag-dihedral-force-constant "%FLAG DIHEDRAL_FORCE_CONSTANT")
(defvar %flag-dihedral-periodicity "%FLAG DIHEDRAL_PERIODICITY")
(defvar %flag-dihedral-phase "%FLAG DIHEDRAL_PHASE")
(defvar %flag-scee-scale-factor "%FLAG SCEE_SCALE_FACTOR")
(defvar %flag-scnb-scale-factor "%FLAG SCNB_SCALE_FACTOR")
(defvar %flag-solty "%FLAG SOLTY")
(defvar %flag-lennard-jones-acoef "%FLAG LENNARD_JONES_ACOEF")
(defvar %flag-lennard-jones-bcoef "%FLAG LENNARD_JONES_BCOEF")
(defvar %flag-bonds-inc-hydrogen "%FLAG BONDS_INC_HYDROGEN")
(defvar %flag-bonds-without-hydrogen "%FLAG BONDS_WITHOUT_HYDROGEN")
(defvar %flag-angles-inc-hydrogen "%FLAG ANGLES_INC_HYDROGEN")
(defvar %flag-angles-without-hydrogen "%FLAG ANGLES_WITHOUT_HYDROGEN")
(defvar %flag-dihedrals-inc-hydrogen "%FLAG DIHEDRALS_INC_HYDROGEN")
(defvar %flag-dihedrals-without-hydrogen "%FLAG DIHEDRALS_WITHOUT_HYDROGEN")
(defvar %flag-excluded-atoms-list "%FLAG EXCLUDED_ATOMS_LIST")
(defvar %flag-amber-atom-type "%FLAG AMBER_ATOM_TYPE")

(defun verify-%flag-line (line)
  (unless (string-equal line "%FLAG" :start1 0 :end1 5)
    (error "Expected %FLAG at the start of the line - got: ~s" line)))

(defmacro rlog (fmt &rest args)
  `(progn
     (cl:format *debug-io* ,fmt ,@args)
     (finish-output *debug-io*)))

;(defun read-amber-parm-format (stream)
;  (let ((fif (fortran:make-fortran-input-file :stream stream))
(defun read-amber-parm-format (fif)
  (let (natom ntypes nbonh mbona ntheth mtheta nphih mphia nhparm nparm
        nnb nres nbona ntheta nphia numbnd numang nptra
        natyp nphb ifpert nbper ngper ndper
        mbper mgper mdper ifbox nmxrs ifcap numextra ncopy
        atom-name charge atomic-number mass atom-type-index number-excluded-atoms
        nonbonded-parm-index residue-label residue-pointer
        bond-force-constant bond-equil-value
        angle-force-constant  angle-equil-value
        dihedral-force-constant dihedral-periodicity dihedral-phase
        scee-scale-factor scnb-scale-factor solty
        lennard-jones-acoef lennard-jones-bcoef
        bonds-inc-hydrogen bonds-without-hydrogen
        angles-inc-hydrogen angles-without-hydrogen
        dihedrals-inc-hydrogen dihedrals-without-hydrogen
        excluded-atoms-list amber-atom-type)
    (rlog "Starting read-amber-parm-format~%")
    (fortran:fread-line fif)     ; Skip the version and timestamp line
    (fortran:fread-line fif)     ; read the first %FLAG line
    ;; From here on down - read the input file
    (loop for line = (fortran:fortran-input-file-look-ahead fif)
       while line
       do (rlog "line ~a~%" line)
       do (verify-%flag-line line)
       do (cond
            ((string-equal %flag-title line :end2 (length %flag-title))
             (fortran:fread-line-or-error fif) ; read %FORMAT(20a4)...
             (fortran:fread-line-or-error fif) ; read ITITL : title
             (fortran:fread-line-or-error fif)) ; Read the next %flag line or eof
            ((string-equal %flag-pointers line :end2 (length %flag-pointers)) ; line is %FLAG POINTERS
             (fortran:fread-line-or-error fif) ; read %FORMAT(10i8) NATOM ...
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) ; read first line of data
               (let ((pointers (fortran:fread-vector fif per-line format-char width)))
                 (setf natom (aref pointers 0) ; natom
                       ntypes (aref pointers 1)
                       nbonh (aref pointers 2)
                       mbona (aref pointers 3)
                       ntheth (aref pointers 4)
                       mtheta (aref pointers 5)
                       nphih (aref pointers 6)
                       mphia (aref pointers 7)
                       nhparm (aref pointers 8)
                       nparm  (aref pointers 9)
                       nnb  (aref pointers 10)
                       nres  (aref pointers 11)
                       nbona  (aref pointers 12)
                       ntheta  (aref pointers 13)
                       nphia  (aref pointers 14)
                       numbnd  (aref pointers 15)
                       numang  (aref pointers 16)
                       nptra (aref pointers 17)
                       natyp  (aref pointers 18)
                       nphb  (aref pointers 19)
                       ifpert (aref pointers 20)
                       nbper (aref pointers 21)
                       ngper  (aref pointers 22)
                       ndper  (aref pointers 23)
                       mbper  (aref pointers 24)
                       mgper (aref pointers 25)
                       mdper (aref pointers 26)
                       ifbox (aref pointers 27)
                       nmxrs (aref pointers 28)
                       ifcap (aref pointers 29)
                       numextra (aref pointers 30)
                       ncopy (aref pointers 31))))
             )
            ((string-equal %flag-atom-name line :end2 (length %flag-atom-name))
             (fortran:fread-line-or-error fif) 
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif)
               (setf atom-name
                     (map 'vector
                          (lambda (string)
                            (intern string "KEYWORD"))
                          (fortran:fread-vector fif per-line format-char width)))))
            ((string-equal %flag-charge line :end2 (length %flag-charge))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf charge (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-atomic-number line :end2 (length %flag-atomic-number))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf atomic-number (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-mass line :end2 (length %flag-mass))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf mass (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-atom-type-index line :end2 (length %flag-atom-type-index))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf atom-type-index (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-number-excluded-atoms line :end2 (length %flag-number-excluded-atoms))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf number-excluded-atoms (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-nonbonded-parm-index line :end2 (length %flag-nonbonded-parm-index))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf nonbonded-parm-index (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-residue-label line :end2 (length %flag-residue-label))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf residue-label
                     (map 'vector
                          (lambda (string)
                            (intern string "KEYWORD"))
                          (fortran:fread-vector fif per-line format-char width)))))
            ((string-equal %flag-residue-pointer line :end2 (length %flag-residue-pointer))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf residue-pointer (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-bond-force-constant line :end2 (length %flag-bond-force-constant))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf bond-force-constant (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-bond-equil-value line :end2 (length %flag-bond-equil-value))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf bond-equil-value (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-angle-force-constant line :end2 (length %flag-angle-force-constant))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf angle-force-constant (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-angle-equil-value line :end2 (length %flag-angle-equil-value))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf angle-equil-value (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-dihedral-force-constant line :end2 (length %flag-dihedral-force-constant))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf dihedral-force-constant (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-dihedral-periodicity line :end2 (length %flag-dihedral-periodicity))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf dihedral-periodicity (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-dihedral-phase line :end2 (length %flag-dihedral-phase))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf dihedral-phase (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-scee-scale-factor line :end2 (length %flag-scee-scale-factor))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf scee-scale-factor (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-scnb-scale-factor line :end2 (length %flag-scnb-scale-factor))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf scnb-scale-factor (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-solty line :end2 (length %flag-solty))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf solty (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-lennard-jones-acoef line :end2 (length %flag-lennard-jones-acoef))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf lennard-jones-acoef (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-lennard-jones-bcoef line :end2 (length %flag-lennard-jones-bcoef))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf lennard-jones-bcoef (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-bonds-inc-hydrogen line :end2 (length %flag-bonds-inc-hydrogen))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf bonds-inc-hydrogen (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-bonds-without-hydrogen line :end2 (length %flag-bonds-without-hydrogen))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf bonds-without-hydrogen (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-angles-inc-hydrogen line :end2 (length %flag-angles-inc-hydrogen))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf angles-inc-hydrogen (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-angles-without-hydrogen line :end2 (length %flag-angles-without-hydrogen))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf angles-without-hydrogen (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-dihedrals-inc-hydrogen line :end2 (length %flag-dihedrals-inc-hydrogen))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf dihedrals-inc-hydrogen (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-dihedrals-without-hydrogen line :end2 (length %flag-dihedrals-without-hydrogen))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf dihedrals-without-hydrogen (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-excluded-atoms-list line :end2 (length %flag-excluded-atoms-list))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf excluded-atoms-list (fortran:fread-vector fif per-line format-char width))))
            ((string-equal %flag-amber-atom-type line :end2 (length %flag-amber-atom-type))
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (setf amber-atom-type (fortran:fread-vector fif per-line format-char width))))
            (t
             (cl:format t "Unknown flag ~a~%" line)
             (fortran:fread-line-or-error fif)  
             (multiple-value-bind (per-line format-char width decimal)
                 (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
               (fortran:fread-line-or-error fif) 
               (fortran:fread-vector fif per-line format-char width)))))
    (rlog "natom -> ~s~%" natom)
    (rlog "ntypes -> ~s~%" ntypes)
    (rlog "nbonh -> ~s~%" nbonh)
    (rlog "mbona -> ~s~%" mbona)
    (rlog "ntheth -> ~s~%" ntheth)
    (rlog "mtheta -> ~s~%" mtheta)
    (rlog "nphih -> ~s~%" nphih)
    (rlog "mphia -> ~s~%" mphia)
    (rlog "nhparm -> ~s~%" nhparm)
    (rlog "nparm -> ~s~%" nparm)
    (rlog "nnb -> ~s~%" nnb)
    (rlog "nres -> ~s~%" nres)
    (rlog "nbona -> ~s~%" nbona)
    (rlog "ntheta -> ~s~%" ntheta)
    (rlog "nphia -> ~s~%" nphia)
    (rlog "numbnd -> ~s~%" numbnd)
    (rlog "numang -> ~s~%" numang)
    (rlog "nptra -> ~s~%" nptra)
    (rlog "natyp -> ~s~%" natyp)
    (rlog "nphb -> ~s~%" nphb)
    (rlog "ifpert -> ~s~%" ifpert)
    (rlog "nbper -> ~s~%" nbper)
    (rlog "ngper -> ~s~%" ngper)
    (rlog "ndper -> ~s~%" ndper)
    (rlog "mbper -> ~s~%" mbper)
    (rlog "mgper -> ~s~%" mgper)
    (rlog "mdper -> ~s~%" mdper)
    (rlog "ifbox -> ~s~%" ifbox)
    (rlog "nmxrs -> ~s~%" nmxrs)
    (rlog "ifcap -> ~s~%" ifcap)
    (rlog "numextra -> ~s~%" numextra)
    (rlog "ncopy -> ~s~%" ncopy)
    (rlog "atom-name -> ~s~%" atom-name)
    (rlog "charge -> ~s~%" charge)
    (rlog "atomic-number -> ~s~%" atomic-number)
    (rlog "mass -> ~s~%" mass)
    (rlog "atom-type-index -> ~s~%" atom-type-index)
    (rlog "number-excluded-atoms -> ~s~%" number-excluded-atoms)
    (rlog "nonbonded-parm-index -> ~s~%" nonbonded-parm-index)
    (rlog "residue-label -> ~s~%" residue-label)
    (rlog "residue-pointer -> ~s~%" residue-pointer)
    (rlog "bond-force-constant -> ~s~%" bond-force-constant)
    (rlog "bond-equil-value -> ~s~%" bond-equil-value)
    (rlog "angle-force-constant -> ~s~%" angle-force-constant)
    (rlog "angle-equil-value -> ~s~%" angle-equil-value)
    (rlog "dihedral-force-constant -> ~s~%" dihedral-force-constant)
    (rlog "dihedral-periodicity -> ~s~%" dihedral-periodicity)
    (rlog "dihedral-phase -> ~s~%" dihedral-phase)
    (rlog "scee-scale-factor -> ~s~%" scee-scale-factor)
    (rlog "scnb-scale-factor -> ~s~%" scnb-scale-factor)
    (rlog "solty -> ~s~%" solty)
    (rlog "lennard-jones-acoef -> ~s~%" lennard-jones-acoef)
    (rlog "lennard-jones-bcoef -> ~s~%" lennard-jones-bcoef)
    (rlog "bonds-inc-hydrogen -> ~s~%" bonds-inc-hydrogen)
    (rlog "bonds-without-hydrogen -> ~s~%" bonds-without-hydrogen)
    (rlog "angles-inc-hydrogen -> ~s~%" angles-inc-hydrogen)
    (rlog "angles-without-hydrogen -> ~s~%" angles-without-hydrogen)
    (rlog "dihedrals-inc-hydrogen -> ~s~%" dihedrals-inc-hydrogen)
    (rlog "dihedrals-without-hydrogen -> ~s~%" dihedrals-without-hydrogen)
    (rlog "excluded-atoms-list -> ~s~%" excluded-atoms-list)
    (rlog "amber-atom-type -> ~s~%" amber-atom-type)
    (let ((energy-stretch (core:make-cxx-object 'chem:energy-stretch))
          (energy-angle (core:make-cxx-object 'chem:energy-angle))
          (energy-dihedral (core:make-cxx-object 'chem:energy-dihedral))
          (energy-nonbond (core:make-cxx-object 'chem:energy-nonbond))
          (atom-table (core:make-cxx-object 'chem:atom-table))
          ;; ... more of these
          (kbs-vec (make-array (+ nbonh mbona) :element-type 'double-float))
          (r0s-vec (make-array (+ nbonh mbona) :element-type 'double-float))
          (i1s-vec (make-array (+ nbonh mbona) :element-type '(signed-byte 32)))
          (i2s-vec (make-array (+ nbonh mbona) :element-type '(signed-byte 32)))
          (atom1s-vec (make-array (+ nbonh mbona)))
          (atom2s-vec (make-array (+ nbonh mbona)))
          (kta-vec (make-array (+ ntheth mtheta) :element-type 'double-float))
          (t0a-vec (make-array (+ ntheth mtheta) :element-type 'double-float))
          (i1a-vec (make-array (+ ntheth mtheta) :element-type '(signed-byte 32)))
          (i2a-vec (make-array (+ ntheth mtheta) :element-type '(signed-byte 32)))
          (i3a-vec (make-array (+ ntheth mtheta) :element-type '(signed-byte 32)))
          (atom1a-vec (make-array (+ ntheth mtheta)))
          (atom2a-vec (make-array (+ ntheth mtheta)))
          (atom3a-vec (make-array (+ ntheth mtheta)))
          (vd-vec (make-array (+ nphih mphia) :element-type 'double-float))
          (ind-vec (make-array (+ nphih mphia) :element-type '(signed-byte 32)))
          (dnd-vec (make-array (+ nphih mphia) :element-type 'double-float))
          (phased-vec (make-array (+ nphih mphia):element-type 'double-float))
          (i1d-vec (make-array (+ nphih mphia) :element-type '(signed-byte 32)))
          (i2d-vec (make-array (+ nphih mphia) :element-type '(signed-byte 32)))
          (i3d-vec (make-array (+ nphih mphia) :element-type '(signed-byte 32)))
          (i4d-vec (make-array (+ nphih mphia) :element-type '(signed-byte 32)))
          (properd-vec (make-array (+ nphih mphia)))
          (atom1d-vec (make-array (+ nphih mphia)))
          (atom2d-vec (make-array (+ nphih mphia)))
          (atom3d-vec (make-array (+ nphih mphia)))
          (atom4d-vec (make-array (+ nphih mphia)))
          (counts 0)
          (counta 0)
          (countd 0)                               
            stretch-vectors angle-vectors dihedral-vectors nonbond-vectors atom-table-vectors
          (atoms (make-array natom))
          )
      ;; Create a vector of atoms to pass to the atom-table and to set the atoms for stretch, angle, dihedral etc.
      (loop for i from 0 below natom
         for name = (aref atom-name i)
         for atom-element = (chem:element-for-atomic-number (aref atomic-number i))
         do (setf (aref atoms i)  (chem:make-atom name atom-element))
           )
      (loop for i from 0 below (length residue-label)
         for name = (aref residue-label i)
         for begin-atom-index = (1- (aref residue-pointer i))
         for end-atom-index = (1- (aref residue-pointer (1+ i)))
         do (let ((residue (chem:make-residue name)))
              (loop for atomi from begin-atom-index below end-atom-index
                   do (chem:add-matter residue (aref atoms atomi)))))
      (rlog "atoms -> ~s~%" atoms)
      (loop for i from 0 below numbnd
            do (loop for j from 0 below nbonh
                     for jicbh = (aref bonds-inc-hydrogen (+ (* j 3) 2))
                                        ;               do (rlog "i ~a j ~a jicbh ~a~%" i j jicbh)
                     do (when (= jicbh (+ i 1))
                          (setf (aref kbs-vec counts) (aref bond-force-constant i)
                                (aref r0s-vec counts) (aref bond-equil-value i)
                                (aref i1s-vec counts) (aref bonds-inc-hydrogen (* j 3))
                                (aref i2s-vec counts) (aref bonds-inc-hydrogen (+ (* j 3) 1))
                                (aref atom1s-vec counts) (aref atoms (/ (aref bonds-inc-hydrogen (* j 3)) 3))
                                (aref atom2s-vec counts) (aref atoms (/ (aref bonds-inc-hydrogen (+ (* j 3) 1)) 3)))
                          (incf counts)))
            do (loop for j from 0 below mbona
                     for jicb = (aref bonds-without-hydrogen (+ (* j 3) 2))           
                     do (when (= jicb (+ i 1))
                          (setf (aref kbs-vec counts) (aref bond-force-constant i)
                                (aref r0s-vec counts) (aref bond-equil-value i)
                                (aref i1s-vec counts) (aref bonds-without-hydrogen (* j 3))
                                (aref i2s-vec counts) (aref bonds-without-hydrogen (+ (* j 3) 1))
                                (aref atom1s-vec counts) (aref atoms (/ (aref bonds-without-hydrogen (* j 3)) 3))
                                (aref atom2s-vec counts) (aref atoms (/ (aref bonds-without-hydrogen (+ (* j 3) 1)) 3)))
                          (incf counts))) 
            )
      (setf stretch-vectors (acons :kb kbs-vec stretch-vectors))
      (setf stretch-vectors (acons :r0 r0s-vec stretch-vectors))
      (setf stretch-vectors (acons :i1 i1s-vec stretch-vectors))
      (setf stretch-vectors (acons :i2 i2s-vec stretch-vectors))
      (setf stretch-vectors (acons :atom1 atom1s-vec stretch-vectors))
      (setf stretch-vectors (acons :atom2 atom2s-vec stretch-vectors))
      (rlog "stretch-vectors -> ~s~%" stretch-vectors)
      (chem:fill-from-vectors-in-alist energy-stretch stretch-vectors)

      (loop for i from 0 below numang
            do (loop for j from 0 below ntheth
                     for jicth = (aref angles-inc-hydrogen (+ (* j 4) 3))
                                        ;               do (rlog "i ~a j ~a jicth ~a~%" i j jicth)
                     do (when (= jicth (+ i 1))
                          (setf (aref kta-vec counta) (aref angle-force-constant i)
                                (aref t0a-vec counta) (aref angle-equil-value i)
                                (aref i1a-vec counta) (aref angles-inc-hydrogen (* j 4))
                                (aref i2a-vec counta) (aref angles-inc-hydrogen (+ (* j 4) 1))
                                (aref i3a-vec counta) (aref angles-inc-hydrogen (+ (* j 4) 2))
                                (aref atom1a-vec counta) (aref atoms (/ (aref angles-inc-hydrogen (* j 4)) 3))
                                (aref atom2a-vec counta) (aref atoms (/ (aref angles-inc-hydrogen (+ (* j 4) 1)) 3))
                                (aref atom3a-vec counta) (aref atoms (/ (aref angles-inc-hydrogen (+ (* j 4) 2)) 3)))
                          (incf counta)))
            do (loop for j from 0 below mtheta
                     for jict = (aref angles-without-hydrogen (+ (* j 4) 3))
                     do (when (= jict (+ i 1))
                          (setf (aref kta-vec counta) (aref angle-force-constant i)
                                (aref t0a-vec counta) (aref angle-equil-value i)
                                (aref i1a-vec counta) (aref angles-without-hydrogen (* j 4))
                                (aref i2a-vec counta) (aref angles-without-hydrogen (+ (* j 4) 1))
                                (aref i3a-vec counta) (aref angles-without-hydrogen (+ (* j 4) 2))
                                (aref atom1a-vec counta) (aref atoms (/ (aref angles-without-hydrogen (* j 4)) 3))
                                (aref atom2a-vec counta) (aref atoms (/ (aref angles-without-hydrogen (+ (* j 4) 1)) 3))
                                (aref atom3a-vec counta) (aref atoms (/ (aref angles-without-hydrogen (+ (* j 4) 2)) 3)))
                          (incf counta))))
      (setf angle-vectors (acons :kt kta-vec  angle-vectors))
      (setf angle-vectors (acons :t0 t0a-vec angle-vectors))
      (setf angle-vectors (acons :i1 i1a-vec angle-vectors))
      (setf angle-vectors (acons :i2 i2a-vec angle-vectors))
      (setf angle-vectors (acons :i3 i3a-vec angle-vectors))
      (setf angle-vectors (acons :atom1 atom1a-vec angle-vectors))
      (setf angle-vectors (acons :atom2 atom2a-vec angle-vectors))
      (setf angle-vectors (acons :atom3 atom3a-vec angle-vectors))
      (rlog "angle-vectors -> ~s~%" angle-vectors)
      (chem:fill-from-vectors-in-alist energy-angle angle-vectors)

      (loop for i from 0 below nptra
         do (loop for j from 0 below nphih
               for jicph = (aref dihedrals-inc-hydrogen (+ (* j 5) 4))
                                        ;              do (rlog "i ~a j ~a jicth ~a~%" i j jicph)
               do (when (= jicph (+ i 1))
                    (setf (aref vd-vec countd) (aref dihedral-force-constant i)
                          (aref ind-vec countd) (round (aref dihedral-periodicity i)) 
                          (aref dnd-vec countd) (aref dihedral-periodicity i) 
                          (aref phased-vec countd) (aref dihedral-phase i)
                          (aref i1d-vec countd) (aref dihedrals-inc-hydrogen (* j 5)) 
                          (aref i2d-vec countd) (aref dihedrals-inc-hydrogen (+ (* j 5) 1))
                          (aref atom1d-vec countd) (aref atoms (/ (aref dihedrals-inc-hydrogen (* j 5)) 3))
                          (aref atom2d-vec countd) (aref atoms (/ (aref dihedrals-inc-hydrogen (+ (* j 5) 1)) 3)))
                    (if (< (aref dihedrals-inc-hydrogen (+ (* j 5) 2)) 0)
                        (setf (aref i3d-vec countd) (- (aref dihedrals-inc-hydrogen (+ (* j 5) 2)))
                              (aref atom3d-vec countd) (aref atoms (/ (- (aref dihedrals-inc-hydrogen (+ (* j 5) 2))) 3)))
                        (setf (aref i3d-vec countd) (aref dihedrals-inc-hydrogen (+ (* j 5) 2))
                              (aref atom3d-vec countd) (aref atoms (/ (aref dihedrals-inc-hydrogen (+ (* j 5) 2)) 3))))
                    (if (< (aref dihedrals-inc-hydrogen (+ (* j 5) 3)) 0)
                        (setf (aref i4d-vec countd) (- (aref dihedrals-inc-hydrogen (+ (* j 5) 3)))
                              (aref atom4d-vec countd) (aref atoms (/ (- (aref dihedrals-inc-hydrogen (+ (* j 5) 3))) 3)))
                        (setf (aref i4d-vec countd) (aref dihedrals-inc-hydrogen (+ (* j 5) 3))
                              (aref atom4d-vec countd) (aref atoms (/ (aref dihedrals-inc-hydrogen (+ (* j 5) 3)) 3))
                              (aref properd-vec countd) T))
                    (incf countd)))
         do (loop for j from 0 below mphia
               for jicp = (aref dihedrals-without-hydrogen (+ (* j 5) 4))
               do (when (= jicp (+ i 1))
                    (setf (aref vd-vec countd) (aref dihedral-force-constant i)
                          (aref ind-vec countd) (round (aref dihedral-periodicity i)) 
                          (aref dnd-vec countd) (aref dihedral-periodicity i) 
                          (aref phased-vec countd) (aref dihedral-phase i)
                          (aref i1d-vec countd) (aref dihedrals-without-hydrogen (* j 5)) 
                          (aref i2d-vec countd) (aref dihedrals-without-hydrogen (+ (* j 5) 1))
                          (aref atom1d-vec countd) (aref atoms (/ (aref dihedrals-without-hydrogen (* j 5)) 3))
                          (aref atom2d-vec countd) (aref atoms (/ (aref dihedrals-without-hydrogen (+ (* j 5) 1)) 3)))
                    (if (< (aref dihedrals-without-hydrogen (+ (* j 5) 2)) 0)
                        (setf (aref i3d-vec countd) (- (aref dihedrals-without-hydrogen (+ (* j 5) 2)))
                              (aref atom3d-vec countd) (aref atoms (/ (- (aref dihedrals-without-hydrogen (+ (* j 5) 2))) 3)))
                        (setf (aref i3d-vec countd) (aref dihedrals-without-hydrogen (+ (* j 5) 2))
                              (aref atom3d-vec countd) (aref atoms (/ (aref dihedrals-without-hydrogen (+ (* j 5) 2)) 3))))
                    (if (< (aref dihedrals-without-hydrogen (+ (* j 5) 3)) 0)
                        (setf (aref i4d-vec countd) (- (aref dihedrals-without-hydrogen (+ (* j 5) 3)))
                              (aref atom4d-vec countd) (aref atoms (/ (- (aref dihedrals-without-hydrogen (+ (* j 5) 3))) 3)))
                        (setf (aref i4d-vec countd) (aref dihedrals-without-hydrogen (+ (* j 5) 3))
                              (aref atom4d-vec countd) (aref atoms (/ (aref dihedrals-without-hydrogen (+ (* j 5) 3)) 3))
                              (aref properd-vec countd) T))
                    (incf countd)))
           )
      (setf dihedral-vectors (acons :v vd-vec dihedral-vectors))
      (setf dihedral-vectors (acons :in ind-vec dihedral-vectors))
      (setf dihedral-vectors (acons :dn dnd-vec dihedral-vectors))
      (setf dihedral-vectors (acons :phase phased-vec dihedral-vectors))
      (setf dihedral-vectors (acons :i1 i1d-vec dihedral-vectors))
      (setf dihedral-vectors (acons :i2 i2d-vec dihedral-vectors))
      (setf dihedral-vectors (acons :i3 i3d-vec dihedral-vectors))
      (setf dihedral-vectors (acons :i4 i4d-vec dihedral-vectors))
      (setf dihedral-vectors (acons :proper properd-vec dihedral-vectors))
      (setf dihedral-vectors (acons :atom1 atom1d-vec dihedral-vectors))
      (setf dihedral-vectors (acons :atom2 atom2d-vec dihedral-vectors))
      (setf dihedral-vectors (acons :atom3 atom3d-vec dihedral-vectors))
      (setf dihedral-vectors (acons :atom4 atom4d-vec dihedral-vectors))
      (rlog "dihedral-vectors -> ~s~%" dihedral-vectors)
      (chem::fill-from-vectors-in-alist energy-dihedral dihedral-vectors)
      ;;atom-table
      (setf atom-table-vectors (acons :atom-name-vector atom-name atom-table-vectors))
      (setf atom-table-vectors (acons :atom-type-vector amber-atom-type atom-table-vectors))
      (setf atom-table-vectors (acons :charge-vector charge atom-table-vectors))
      (setf atom-table-vectors (acons :mass-vector mass atom-table-vectors))
      (setf atom-table-vectors (acons :atomic-number-vector atomic-number atom-table-vectors))
 ;     (setf atom-table-vectors (acons :residues residue-pointer atom-table-vectors))
 ;     (setf atom-table-vectors (acons :residue-names residue-labels atom-table-vectors))
          
      (rlog "atom-table-vectors -> ~s~%" atom-table-vectors)
      (dolist (entry atom-table-vectors)
        (format *debug-io* "entry -> ~s  (type-of (cdr entry)) -> ~s~%" entry (type-of (cdr entry))))
      (chem::fill-atom-table-from-vectors atom-table atom-table-vectors)
      ;;nonbond
      (setf nonbond-vectors (acons :ntypes ntypes nonbond-vectors))
      (setf nonbond-vectors (acons :atom-name-vector (copy-seq atom-name) nonbond-vectors))
      (setf nonbond-vectors (acons :charge-vector (copy-seq charge) nonbond-vectors))
      (setf nonbond-vectors (acons :mass-vector (copy-seq mass) nonbond-vectors))
      (setf nonbond-vectors (acons :atomic-number-vector (copy-seq atomic-number) nonbond-vectors))
      (setf nonbond-vectors (acons :ico-vec (copy-seq nonbonded-parm-index) nonbond-vectors))
      (setf nonbond-vectors (acons :iac-vec (copy-seq atom-type-index) nonbond-vectors))
      (setf nonbond-vectors (acons :cn1-vec (copy-seq lennard-jones-acoef) nonbond-vectors))
      (setf nonbond-vectors (acons :cn2-vec (copy-seq lennard-jones-bcoef) nonbond-vectors))
      (rlog "nonbond-vectors -> ~s~%" nonbond-vectors)
      (chem:construct-nonbond-terms-from-aList energy-nonbond nonbond-vectors)
      (loop for i from 0 below (length excluded-atoms-list)
         for atom = (- (aref excluded-atoms-list i) 1)
          do (setf (aref excluded-atoms-list i) atom))
      (chem:set-nonbond-excluded-atom-info energy-nonbond atom-table (copy-seq excluded-atoms-list) (copy-seq number-excluded-atoms))
      ;;
      ;; Now we have energy-stretch, energy-angle, and energy-dihedral
      ;;   we want to put them into an energy-function
      (rlog "atom-table-vectors -> ~s~%" atom-table-vectors)
      ;;
      #+(or) ;; the old code
      (let ((energy-function (core:make-cxx-object 'chem:energy-function
                                                   :stretch energy-stretch
                                                   :angle energy-angle
                                                   :dihedral energy-dihedral
                                                   :nobond energy-nonbond
                                                   :atom-table atom-table)))
        energy-function #| <-- This was missing before |# )
      (let ((alist (list (cons :atom-table atom-table)
                         (cons :stretch energy-stretch)
                         (cons :angle energy-angle)
                         (cons :dihedral energy-dihedral)
                         (cons :nonbond energy-nonbond)))
            (energy-function (core:make-cxx-object 'chem:energy-function)))
        (chem:fill-energy-function-from-alist energy-function alist)
        energy-function)
      )))

(defun read-amber-coordinate-file (fif)
  (fortran:fread-line fif)       ; Skip the version and timestamp line
  (let* ((line (fortran:fread-line fif))
         (results (make-array 3 :element-type 't :adjustable t :fill-pointer 0))
         ;; Read format FORMAT(I5,5E15.7) NATOM,TIME,TEMP
         (natoms (parse-integer line :start 0 :end 5))
         (time (if (>= (length line) 15)
                   (fortran::parse-double-float (subseq line 5 nil) :start 0 :end 15)
                   nil))
         (temp (if (>= (length line) 35)
                   (fortran::parse-double-float (subseq line (+ 5 15) nil) :start 0 :end 15)
                   nil)))
    (fortran:fread-line fif) 
    (copy-seq (fortran:fread-vector fif 6 #\F 12))))


;;; The following code is to generate a human readable representation of an energy-function
;;; with everything sorted so that the terms can be compared side-by-side using something like
;;; 'ediff'

(defstruct stretch-term
  atom1-name atom2-name atom1-type atom2-type kb r0)

#|(defun extract-energy-stretch (energy-function)
  (let* ((energy-stretch (chem:get-stretch-component energy-function))
         (stretch-vectors (chem:extract-vectors-as-alist energy-stretch))
         (kb-vector (cdr (assoc :kb stretch-vectors)))
         (r0-vector (cdr (assoc :r0 stretch-vectors)))
         (i1-vector (cdr (assoc :i1 stretch-vectors)))
         (i2-vector (cdr (assoc :i2 stretch-vectors)))
         (atom1-vector (cdr (assoc :atom1 stretch-vectors)))
         (atom2-vector (cdr (assoc :atom2 stretch-vectors)))
         )
    (rlog "stretch-vectors ~s~%" stretch-vectors)  
    (let ((stretches (loop for atom1 across atom1-vector
                        for atom2 across atom2-vector
                        for kb across kb-vector
                        for r0 across r0-vector
                        for atom1-name = (chem:get-name atom1)
                        for atom2-name = (chem:get-name atom2)
                        for atom1-type = (chem:get-type atom1)
                        for atom2-type = (chem:get-type atom2)
                        collect (make-stretch-term :atom1-name atom1-name
                                                   :atom2-name atom2-name
                                                   :atom1-type atom1-type
                                                   :atom2-type atom2-type
                                                   :kb kb
                                                   :r0 r0))))
      (rlog "stretch-vectors ~s~%" stretch-vectors) 
      (loop for stretch in stretches
         do (when (string> (string (stretch-term-atom1-name stretch))
                           (string (stretch-term-atom2-name stretch)))
              (rotatef (stretch-term-atom1-name stretch)
                       (stretch-term-atom2-name stretch))
              (rotatef (stretch-term-atom1-type stretch)
                       (stretch-term-atom2-type stretch))))
      (flet ((order-stretch (s1 s2)
               (if (string< (string (stretch-term-atom1-name s1))
                            (string (stretch-term-atom1-name s2)))
                   t
                   (if (string> (string (stretch-term-atom1-name s1))
                                (string (stretch-term-atom1-name s2)))
                       nil
                       (string< (string (stretch-term-atom2-name s1))
                                (string (stretch-term-atom2-name s2)))))))
        (sort stretches #'order-stretch)))
(rlog "stretch-vectors ~s~%" stretch-vectors))) |#

(defun extract-energy-stretch (energy-function)
  (let* ((energy-stretch (chem:get-stretch-component energy-function))
         (stretch-vectors (chem:extract-vectors-as-alist energy-stretch))
         (kb-vector (cdr (assoc :kb stretch-vectors)))
         (r0-vector (cdr (assoc :r0 stretch-vectors)))
         (i1-vector (cdr (assoc :i1 stretch-vectors)))
         (i2-vector (cdr (assoc :i2 stretch-vectors)))
         (atom1-vector (cdr (assoc :atom1 stretch-vectors)))
         (atom2-vector (cdr (assoc :atom2 stretch-vectors)))
         )
    (let ((stretches (loop for atom1 across atom1-vector
                        for atom2 across atom2-vector
                        for kb across kb-vector
                        for r0 across r0-vector
                        for i1 across i1-vector
                        for i2 across i2-vector
                        for atom1-name = atom1
                        for atom2-name = atom2
                        for atom1-type = i1
                        for atom2-type = i2
                        collect (make-stretch-term :atom1-name atom1-name
                                                   :atom2-name atom2-name
                                                   :atom1-type atom1-type
                                                   :atom2-type atom2-type
                                                   :kb kb
                                                   :r0 r0))))
      (rlog "stretch-vectors ~s~%" stretch-vectors) 
      (loop for stretch in stretches
         do (when (string> (string (stretch-term-atom1-name stretch))
                           (string (stretch-term-atom2-name stretch)))
              (rotatef (stretch-term-atom1-name stretch)
                       (stretch-term-atom2-name stretch))
              (rotatef (stretch-term-atom1-type stretch)
                       (stretch-term-atom2-type stretch))))
      (flet ((order-stretch (s1 s2)
               (if (string< (string (stretch-term-atom1-name s1))
                            (string (stretch-term-atom1-name s2)))
                   t
                   (if (string> (string (stretch-term-atom1-name s1))
                                (string (stretch-term-atom1-name s2)))
                       nil
                       (string< (string (stretch-term-atom2-name s1))
                                (string (stretch-term-atom2-name s2)))))))
        (sort stretches #'order-stretch)))
    (rlog "stretch-vectors ~s~%" stretch-vectors)))


(defstruct angle-term
  atom1-name atom2-name atom3-name atom1-type atom2-type atom3-type kt t0)

(defun extract-energy-angle (energy-function)
  (let* ((energy-angle (chem:get-angle-component energy-function))
         (angle-vectors (chem:extract-vectors-as-alist energy-angle))
         (kt-vector (cdr (assoc :kt angle-vectors)))
         (t0-vector (cdr (assoc :t0 angle-vectors)))
         (i1-vector (cdr (assoc :i1 angle-vectors)))
         (i2-vector (cdr (assoc :i2 angle-vectors)))
         (i3-vector (cdr (assoc :i3 angle-vectors)))
         (atom1-vector (cdr (assoc :atom1 angle-vectors)))
         (atom2-vector (cdr (assoc :atom2 angle-vectors)))
         (atom3-vector (cdr (assoc :atom3 angle-vectors)))
         )
    (let ((angles (loop for atom1 across atom1-vector
                     for atom2 across atom2-vector
                     for atom3 across atom3-vector
                     for kt across kt-vector
                     for t0 across t0-vector
                     for atom1-name = (chem:get-name atom1)
                     for atom2-name = (chem:get-name atom2)
                     for atom3-name = (chem:get-name atom3)
                     for atom1-type = (chem:get-type atom1)
                     for atom2-type = (chem:get-type atom2)
                     for atom3-type = (chem:get-type atom3)
                     collect (make-angle-term :atom1-name atom1-name
                                              :atom2-name atom2-name
                                              :atom3-name atom3-name
                                              :atom1-type atom1-type
                                              :atom2-type atom2-type
                                              :atom3-type atom3-type
                                              :kt kt
                                              :t0 t0))))
      (rlog "angle-vectors ~s~%" angle-vectors) 
      (loop for angle in angles
         do (when (string> (string (angle-term-atom1-name angle))
                           (string (angle-term-atom3-name angle)))
              (rotatef (angle-term-atom1-name angle)
                       (angle-term-atom3-name angle))
              (rotatef (angle-term-atom1-type angle)
                       (angle-term-atom3-type angle))))
      (flet ((order-angle (a1 a2)
               (if (string< (string (angle-term-atom1-name a1))
                            (string (angle-term-atom1-name a2)))
                   t
                   (if (string< (string (angle-term-atom2-name a1))
                                (string (angle-term-atom2-name a2)))
                       t
                       (if (string> (string (angle-term-atom2-name a1))
                                    (string (angle-term-atom2-name a2)))
                           nil
                           (string< (string (angle-term-atom3-name a1))
                                    (string (angle-term-atom3-name a2))))))))
        (sort angles #'order-angle)))
          (rlog "angle-vectors ~s~%" angle-vectors)))

(defstruct dihedral-term
  atom1-name atom2-name atom3-name atom4-name atom1-type atom2-type atom3-type atom4-type v in phase proper)

(defun extract-energy-dihedral (energy-function)
  (let* ((energy-dihedral (chem:get-dihedral-component energy-function))
         (dihedral-vectors (chem:extract-vectors-as-alist energy-dihedral))
         (v-vector (cdr (assoc :v dihedral-vectors)))
         (in-vector (cdr (assoc :in dihedral-vectors)))
         (phase-vector (cdr (assoc :phase dihedral-vectors)))
         (proper-vector (cdr (assoc :prpoer dihedral-vectors)))
         (i1-vector (cdr (assoc :i1 dihedral-vectors)))
         (i2-vector (cdr (assoc :i2 dihedral-vectors)))
         (i3-vector (cdr (assoc :i3 dihedral-vectors)))
         (i4-vector (cdr (assoc :i4 dihedral-vectors)))
         (atom1-vector (cdr (assoc :i1 dihedral-vectors)))
         (atom2-vector (cdr (assoc :i2 dihedral-vectors)))
         (atom3-vector (cdr (assoc :i3 dihedral-vectors)))
         (atom4-vector (cdr (assoc :i4 dihedral-vectors)))
         )
    (let ((dihedrals (loop for atom1 across atom1-vector
                        for atom2 across atom2-vector
                        for atom3 across atom3-vector
                        for atom4 across atom4-vector
                        for v across v-vector
                        for in across in-vector
                        for phase across phase-vector
                        for proper across proper-vector
                        for atom1-name = (chem:get-name atom1)
                        for atom2-name = (chem:get-name atom2)
                        for atom3-name = (chem:get-name atom3)
                        for atom4-name = (chem:get-name atom4)
                        for atom1-type = (chem:get-type atom1)
                        for atom2-type = (chem:get-type atom2)
                        for atom3-type = (chem:get-type atom3)
                        for atom4-type = (chem:get-type atom4)
                        collect (make-dihedral-term :atom1-name atom1-name
                                                    :atom2-name atom2-name
                                                    :atom3-name atom3-name
                                                    :atom4-name atom4-name
                                                    :atom1-type atom1-type
                                                    :atom2-type atom2-type
                                                    :atom3-type atom3-type
                                                    :atom4-type atom4-type
                                                    :v v
                                                    :in in
                                                    :phase phase
                                                    :propeer proper))))
      (rlog "dihedral-vectors ~s~%" dihedral-vectors) 
      (loop for dihedral in dihedrals
         do (when (string> (string (dihedral-term-atom1-name dihedral))
                           (string (dihedral-term-atom4-name dihedral)))
              (rotatef (dihedral-term-atom1-name dihedral)
                       (dihedral-term-atom4-name dihedral))
              (rotatef (dihedral-term-atom1-type dihedral)
                       (dihedral-term-atom4-type dihedral))
              (rotatef (dihedral-term-atom2-name dihedral)
                       (dihedral-term-atom3-name dihedral))
              (rotatef (dihedral-term-atom2-type dihedral)
                       (dihedral-term-atom3-type dihedral))))
      (flet ((order-dihedral (d1 d2)
               (if (string< (string (dihedral-term-atom1-name d1))
                            (string (dihedral-term-atom1-name d2)))
                   t
                   (if (string< (string (dihedral-term-atom2-name d1))
                                (string (dihedral-term-atom2-name d2)))
                       t
                       (if (string< (string (dihedral-term-atom3-name d1))
                                    (string (dihedral-term-atom3-name d2)))
                           t
                           (if (string> (string (dihedral-term-atom3-name d1))
                                        (string (dihedral-term-atom3-name d2)))
                               nil
                               (string< (string (dihedral-term-atom4-name d1))
                                        (string (dihedral-term-atom4-name d2)))))))))
        (sort dihedrals #'order-dihedral)))
    (rlog "dihedral-vectors ~s~%" dihedral-vectors)))

                              
