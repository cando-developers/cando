(in-package :leap.topology)


(defconstant +amber-charge-conversion-18.2223+ 18.2223)

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

(defun bond-order-int (atom1 atom2)
  "Return the bond order as an integer"
  (let* ((bond (chem:get-bond-to atom1 atom2)))
    (chem:get-order-as-int bond)))
               
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
              (bond-orders (make-array without-h :element-type '(signed-byte 8)))
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
                      (setf (aref bond-orders cur) (bond-order-int atom1 atom2))
                      (incf cur))))
          (values with-h without-h ibh jbh icbh ib jb icb kbj-vec r0j-vec bond-orders))))))

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
  (let ((j-vec (make-array (length v-vec) :fill-pointer 0))
        (jnext 0)
        jtemp
        (uniques (make-hash-table :test #'equal))
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
          for key = (list vi ini phasei)
          ;;          for key = (intern (format nil "~15f-~d-~f" vi ini phasei) :keyword) 
 ;;;         do (format t "atom1: ~a atom2: ~a atom3: ~a atom4 ~a key: ~a vi: ~a ini ~a phase ~a~% " atom1 atom2 atom3 atom4 key vi ini phasei)
          do (if (setf jtemp (gethash key uniques))
                 (vector-push jtemp j-vec)
                 (progn
                   (vector-push jnext j-vec)
                   (setf (gethash key uniques) jnext)
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

(defun prepare-amber-energy-dihedral (energy-function ib jb it kt)
  "Prepare the amber energy dihedral.  
Use the ib,jb (the atoms of bonds without hydrogen) and
it,kt (the terminal atoms of angles without hydrogen) and
if the terminal atoms of a proper dihedral are in one of those two lists
then don't calculate 1,4 interactions"
  (format t "Starting prepare-amber-energy-dihedral~%")
  (finish-output)
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
      (format t "Ordering i1,i2,i3,i4 prepare-amber-energy-dihedral~%")
      (finish-output)
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
      (format t "1-4 interactions prepare-amber-energy-dihedral~%")
      (finish-output)
      (let ((in-bond-or-angle (make-hash-table :test #'equal)))
        (progn
          (loop for bond-index below (length ib)
                for ib-index = (aref ib bond-index)
                for jb-index = (aref jb bond-index)
                do (setf (gethash (cons ib-index jb-index) in-bond-or-angle) t)
                   (setf (gethash (cons jb-index ib-index) in-bond-or-angle) t))
          (loop for angle-index below (length it)
                for it-index = (aref it angle-index)
                for kt-index = (aref kt angle-index)
                do (setf (gethash (cons it-index kt-index) in-bond-or-angle) t)
                   (setf (gethash (cons kt-index it-index) in-bond-or-angle) t)))
        (flet ((in-same-bond-or-angle (i1x i4x ib jb it kt)
                 ;; Return T if i1x,i4x are in a bond
                 ;; Use hash table
                 (let ((fast-result (gethash (cons i1x i4x) in-bond-or-angle))
                       #+(or) (slow-result (block slow-path
                                      (loop for bond-index below (length ib)
                                            for ib-index = (aref ib bond-index)
                                            for jb-index = (aref jb bond-index)
                                            when (or (and (= ib-index i1x) (= jb-index i4x))
                                                     (and (= ib-index i4x) (= jb-index i1x)))
                                              do (return-from slow-path t))
                                      ;; Return T if i1x,i4x are in an angle
                                      (loop for angle-index below (length it)
                                            for it-index = (aref it angle-index)
                                            for kt-index = (aref kt angle-index)
                                            when (or (and (= it-index i1x) (= kt-index i4x))
                                                     (and (= it-index i4x) (= kt-index i1x)))
                                              do (return-from slow-path t))
                                      nil)))
                   #+(or)  (unless (eq fast-result slow-result)
                     (error "The fast result and the slow result don't match"))
                   fast-result)))
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
                     (t                 ; It's an improper
                      (setf (aref i3-vector x) (- i3x)
                            (aref i4-vector x) (- i4x))))
                   (setf proper-prev properx
                         i1prev i1x
                         i2prev i2x
                         i3prev i3x
                         i4prev i4x))))
      (multiple-value-bind (j-vec vj-vec inj-vec phasej-vec properj-vec)
          (collapse-dihedral-parameters v-vector in-vector phase-vector proper-vector atom1-vector atom2-vector atom3-vector atom4-vector)
        (let ((with-h 0)
              (without-h 0))
          (format t "Counting w and w/o water prepare-amber-energy-dihedral~%")
          (finish-output)
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
            (format t "Extracting prepare-amber-energy-dihedral~%")
            (finish-output)
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
                          (ffnonbond2 (chem:get-ffnonbond-using-type-index ffnonbond-db type2))
                          ;; The index is calculated with the same equation in leap
                          ;; iIndex = iY * (iY + 1) / 2 + iX + 1;        /* +1 because they are FORTRAN */
                          ;; See unitio.c line 3639
                          (index (- (+ (/ (* (+ j 1) j) 2) (+ i 1)) 1)) ; (((((j + 1) * j) / 2) + (i + 1)) - 1)
                          (ico-index (+ (* ntypes i) j)))
                      (setf rstar (+ (chem:get-radius-angstroms ffnonbond1) (chem:get-radius-angstroms ffnonbond2))
                            epsilonij (sqrt (* (chem:get-epsilon-k-cal ffnonbond1) (chem:get-epsilon-k-cal ffnonbond2)))
                            (aref cn1-vec index) (* epsilonij (expt rstar 12.0))
                            (aref cn2-vec index) (* 2.0 epsilonij (expt rstar 6.0))
                            (aref ico-vec ico-index) (+ (/ (* (+ j 1) j) 2) (+ i 1)))
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

(defun prepare-molecules (energy-function)
  "For each molecule assign an integer index for the force-field and return (values molecule-ff-index-vec ff-name-vec."
  (let* ((atom-table (chem:atom-table energy-function))
         (molecules (chem:atom-table-molecules atom-table))
         (force-field-names (make-hash-table))
         (force-field-name-index 0))
    (loop for molecule across molecules
          for force-field-name = (chem:force-field-name molecule)
          unless (gethash force-field-name force-field-names)
            do (setf (gethash force-field-name force-field-names)
                     (prog1 force-field-name-index
                       (incf force-field-name-index))))
    (let ((force-field-names-vec (make-array (hash-table-count force-field-names))))
      (maphash (lambda (ff-name ff-index)
                 (setf (elt force-field-names-vec ff-index) ff-name))
               force-field-names)
      (let ((molecule-force-field-names-index-vec (make-array (length molecules))))
        (loop for mol across molecules
              for index from 0
              for ff-name = (chem:force-field-name mol)
              for ff-index = (gethash ff-name force-field-names)
              do (setf (elt molecule-force-field-names-index-vec index) ff-index))
        (values molecule-force-field-names-index-vec force-field-names-vec)))))

(defun prepare-residue (energy-function residue-name-to-pdb-alist)
  (let* ((atom-table (chem:atom-table energy-function))
         (residue-vector (chem:atom-table-residues atom-table))
         (residue-pointer-prepare-vector (chem:atom-table-residue-pointers atom-table))
         (residue-name-vector (chem:atom-table-residue-names atom-table))
         (molecule-vector (chem:atom-table-atoms-per-molecule atom-table))
         (residue-pointer-vector (make-array 256 :element-type '(signed-byte 32) :fill-pointer 0 :adjustable t))
         (atoms-per-molecule (make-array 256 :element-type '(signed-byte 32) :fill-pointer 0 :adjustable t))
         (nmolecule 0)
         (nresidue 0)
         (nmxrs 0))
    ;; Fix residue names that are too long
    (loop for index from 0 below (length residue-name-vector)
          for residue-name = (aref residue-name-vector index)
          when (> (length (string residue-name)) 3)
            do (let ((short-residue-name-pair (assoc residue-name residue-name-to-pdb-alist)))
                 (unless short-residue-name-pair
                   (error "There is must be a short residue name for ~a" residue-name))
                 (let ((short-name (cdr short-residue-name-pair)))
                   (format t "Using short name ~s in place of ~s~%" short-name residue-name)
                   (when (< (length (string short-name)) 3)
                     (warn "The short name ~s should be no more than three characters long - residue-name-to-pdb-alist -> ~s" short-name residue-name-to-pdb-alist)
                     (setf short-name (intern (subseq (string short-name) 0 3) :keyword)))
                   (setf (elt residue-name-vector index) short-name))))
    (setf nmolecule (length molecule-vector))
    (setf nresidue (length residue-name-vector))
    (loop for i from 0 below (length residue-pointer-prepare-vector)
       do (if (= i 0)
              (vector-push-extend (+ (aref residue-pointer-prepare-vector 0) 1) residue-pointer-vector)
              (if (/= (aref residue-pointer-prepare-vector (- i 1)) (aref residue-pointer-prepare-vector i))
                  (vector-push-extend (+ (aref residue-pointer-prepare-vector i) 1) residue-pointer-vector))))
    ;;      for fresidue = (+ (aref residue-pointer-prepare-vector i) 1)
    ;;      do (setf (aref residue-pointer-prepare-vector i) fresidue))
    (loop for i from 0 below nresidue
       for inmxrs = (- (aref residue-pointer-prepare-vector (+ i 1)) (aref residue-pointer-prepare-vector i))
       do (if (> inmxrs nmxrs)
              (setf nmxrs inmxrs)))
    (let ((prev-molv 0))
      (loop for i from 0 below nmolecule
            for molv = (aref molecule-vector i)
            for natom = (- molv prev-molv)
            do (vector-push-extend natom atoms-per-molecule)
            (setf prev-molv molv)))
    (values nresidue nmxrs residue-pointer-vector residue-name-vector atoms-per-molecule residue-vector)))

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

(defmacro outline-progn (&body body)
  "Break some code out into a separate function"
  `(funcall
    (lambda ()
      ,@body)))

     
(defun save-amber-parm-format-using-energy-function (energy-function topology-pathname coordinate-pathname residue-name-to-pdb-alist &key (cando-extensions t))
  "Generate an AMBER topology/coordinate file pair using the energy function. 
Arguments:
energy-function : The energy-function to generate the topology from.
topology-pathname : Where to write the topology file.
coordinate-pathname : Where to write the coordinate file (ascii).
residue-name-to-pdb-alist : An alist of long residue names to short PDB residue names.
cando-extensions               : T if you want cando-extensions written to the topology file."
  (let* ((bar (cando:make-progress-bar :style :bar :message "Saving" :total 41 :width 41 :divisions 41))
         (bar-counter 0)
         (nonbonds (chem:get-nonbond-component energy-function))
         (number-excluded-atoms (chem:number-excluded-atoms nonbonds))
         (excluded-atom-list (chem:excluded-atom-list nonbonds))
         (topology-pathname (merge-pathnames topology-pathname))
         (coordinate-pathname (merge-pathnames coordinate-pathname))
         (atom-table (chem:atom-table energy-function))
         (natom (chem:get-number-of-atoms atom-table))
         residue-vec)
    ;; Skip assigning MarkMainChainAtoms and MarkSideChain atoms for now
    ;; see (unitio.c:4889).  This won't mean anything for spiroligomers.
    (format t "Writing to ~a~%" topology-pathname)
    (finish-output)
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
      (fortran:fwrite (string (chem:aggregate-name atom-table)))
      ;; This function will be very, very long (for Common Lisp)
      ;; To avoid lots of nested scopes we will declare one large scope
      ;;   and declare all of the variables in that scope here at the top
      (let (atom-vectors
            molecule-force-field-name-indices force-field-names-vec
            ntypes atom-name atom-type charge mass atomic-number atom-radius ico iac local-typej-vec cn1-vec cn2-vec #|nonbonds|#
            nbonh mbona ibh jbh icbh ib jb icb kbj-vec r0j-vec #|stretches|#
            ntheth mtheta ith jth kth icth it jt kt1 ict ktj-vec t0j-vec #|angles|#
            nphih mphia iph jph kph lph icph ip jp kp lp icp vj-vec inj-vec phasej-vec properj-vec #|dihedrals|#
            nhparm nparm  nnb nres residue-pointer-vec residue-name-vec atoms-per-molecule
            generalized-born-radius generalized-born-screen
            nbona    ntheta nphia  NUMBND NUMANG NPTRA
            NATYP    NPHB   IFPERT NBPER  NGPER  NDPER
            MBPER    MGPER  MDPER  IFBOX  nmxrs  IFCAP
            NUMEXTRA NCOPY
            non-h-bond-orders
            )
        ;; Here we need to calculate all of the values for %FLAG POINTERS
        (multiple-value-setq (nbonh mbona ibh jbh icbh ib jb icb kbj-vec r0j-vec non-h-bond-orders)
          (prepare-amber-energy-stretch energy-function))
        (multiple-value-setq (ntheth mtheta ith jth kth icth it jt kt1 ict ktj-vec t0j-vec)
          (prepare-amber-energy-angle energy-function))
        (multiple-value-setq (nphih mphia iph jph kph lph icph ip jp kp lp icp vj-vec inj-vec phasej-vec properj-vec)
          (prepare-amber-energy-dihedral energy-function ib jb it kt1))
                                        ;        (multiple-value-setq (ntypes atom-name charge mass atomic-number ico iac local-typej-vec cn1-vec cn2-vec)
                                        ;          (chem:prepare-amber-energy-nonbond energy-function))
        (multiple-value-setq (nres nmxrs residue-pointer-vec residue-name-vec atoms-per-molecule residue-vec)
          (prepare-residue energy-function residue-name-to-pdb-alist))
        (multiple-value-setq (molecule-force-field-name-indices force-field-names-vec)
          (prepare-molecules energy-function))
        (setf atom-vectors (chem:prepare-amber-energy-nonbond energy-function (chem:nonbond-force-field-for-aggregate atom-table)))
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
        (setf ifbox (chem:bounding-box-bound-p atom-table))
        (setf ifcap 0)     
        (setf numextra 0)  
        (setf ncopy 0)
                                        ;       (setf atom-type-name (chem:get-type atoms))
 
        #| dihedrals, nonbonds, others??? |#
        ;; --- Done calculating all of the values
        ;; --- Now write out all of the values
        (fortran:fformat 1 "%-80s")
        (cando:progress-advance bar (incf bar-counter))
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
        (fortran:fwrite (if ifbox 1 0)) ; set to 1 if standard periodic box, 2 when truncated octahedral
        (fortran:fwrite nmxrs) ; number of atoms in the largest residue
        (fortran:fwrite ifcap) ; set to 1 if the CAP option from edit was specified
        (fortran:fwrite numextra) ; number of extra points found in topology
        (fortran:fwrite ncopy) ; number of PIMD slices / number of beads
        (fortran:end-line)

        ;; Next)
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
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
         (fortran:end-line))
        ;; write the atom names

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG CHARGE")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-4-")
         (fortran:fformat 5 "%16.8e")
         (loop for ch across charge
               do (fortran:fwrite (* ch +amber-charge-conversion-18.2223+)))
         (fortran:end-line))
        ;; write the atom charges

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG ATOMIC_NUMBER")
         (fortran:fwrite "%FORMAT(10I8)")
         (fortran:debug "-5-")
         (fortran:fformat 10 "%8d")
         (loop for number across atomic-number
               do (fortran:fwrite number))
         (fortran:end-line))
        ;; write the atomic number of each atom

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG MASS")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-6-")
         (fortran:fformat 5 "%16.8e")
         (loop for ma across mass
               do (fortran:fwrite ma))
         (fortran:end-line))
        ;; write the atom masses

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG ATOM_TYPE_INDEX")
         (fortran:fwrite "%FORMAT(10I8)")
         (fortran:debug "-7-")
         (fortran:fformat 10 "%8d")
         (loop for ia across iac
               do (fortran:fwrite ia))
         (fortran:end-line))
        ;; write the index fot the atom types

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG NUMBER_EXCLUDED_ATOMS")
         (fortran:end-line)
         (fortran:fwrite "%FORMAT(10I8)")
         (fortran:debug "-8-")
         (fortran:fformat 10 "%8d")
         (loop for na across number-excluded-atoms
               do (fortran:fwrite na))
         (fortran:end-line))
        ;; write the total number of excluded atoms for atom "i"

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG NONBONDED_PARM_INDEX")
         (fortran:fwrite "%FORMAT(10I8)")
         (fortran:debug "-9-")
         (fortran:fformat 10 "%8d")
         (loop for ic across ico
               do (fortran:fwrite ic))
         (fortran:end-line))
        ;; provides the index to the nobon parameter arrays CN1, CN2 and ASOL, BSOL.

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG RESIDUE_LABEL")
         (fortran:fwrite "%FORMAT(20A4)")
         (fortran:debug "-10-")
         (fortran:fformat 20 "%-4s")
         (loop for ren across residue-name-vec
               do (fortran:fwrite (string ren)))
         (fortran:end-line))
        ;; write the name of each of the residues

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG RESIDUE_POINTER")
         (fortran:fwrite "%FORMAT(10I8)")
         (fortran:debug "-11-")
         (fortran:fformat 10 "%8d")
         (loop for re from 0 below (- (length residue-pointer-vec) 1)
               do (fortran:fwrite (aref residue-pointer-vec re)))
         (fortran:end-line))
        ;; write the atoms in each residue are listed for atom "1" in IPRES(i) to IPRES(i+1)-1

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG BOND_FORCE_CONSTANT")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-12-")
         (fortran:fformat 5 "%16.8e")
         (loop for kb across kbj-vec
               do (fortran:fwrite kb))
         (fortran:end-line))
        
        ;; write the force constant for the bonds of each type

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG BOND_EQUIL_VALUE")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-13-")
         (fortran:fformat 5 "%16.8e")
         (loop for r0 across r0j-vec
               do (fortran:fwrite r0))
         (fortran:end-line))
        ;; write the equilibrium bond length for the bonds of each type

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG ANGLE_FORCE_CONSTANT")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-14-")
         (fortran:fformat 5 "%16.8e")
         (loop for kt across ktj-vec
               do (fortran:fwrite kt))
         (fortran:end-line))
        ;; write the force constant for the angles of each type

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG ANGLE_EQUIL_VALUE")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-15-")
         (fortran:fformat 5 "%16.8e")
         (loop for t0 across t0j-vec
               do (fortran:fwrite t0))
         (fortran:end-line))
        ;; write the equilibrium angle for the angles of each type

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG DIHEDRAL_FORCE_CONSTANT")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-16-")
         (fortran:fformat 5 "%16.8e")
         (loop for v0 across vj-vec
               do (fortran:fwrite v0))
         (fortran:end-line))
        ;; write the force constant for the dihedral of a given type

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG DIHEDRAL_PERIODICITY")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-17-")
         (fortran:fformat 5 "%16.8e")
         (loop for in0 across inj-vec
               do (fortran:fwrite (float in0)))
         (fortran:end-line))
        ;; write the periodicity of the dihedral of a given type

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG DIHEDRAL_PHASE")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-18-")
         (fortran:fformat 5 "%16.8e")
         (loop for p0 across phasej-vec
               do (fortran:fwrite p0))
         (fortran:end-line))
        ;; write the phase of the dihedral of a given type

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG SCEE_SCALE_FACTOR")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-19-")
         (fortran:fformat 5 "%16.8e")
         (loop for pr0 across properj-vec
               do (if pr0
                      (fortran:fwrite 1.2)
                      (fortran:fwrite 0.0)))
         (fortran:end-line))
        ;; write the 1-4 electrostatic scaling constant

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG SCNB_SCALE_FACTOR")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-20-")
         (fortran:fformat 5 "%16.8e")
         (loop for pr0 across properj-vec
               do (if pr0
                      (fortran:fwrite 2.0)
                      (fortran:fwrite 0.0)))
         (fortran:end-line))
        ;; write the 1-4 vdw scaling constant

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG SOLTY")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-21-")
         (fortran:fformat 5 "%16.8e")
         (loop repeat natyp
               do (fortran:fwrite 0.0))
         (fortran:end-line))
        ;; currently unused

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG LENNARD_JONES_ACOEF")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-22-")
         (fortran:fformat 5 "%16.8e")
         (loop for cn1 across cn1-vec
               do (fortran:fwrite cn1))
         (fortran:end-line))
        ;; write the Lennard Jones r**12 terms for all possible atom type interactions

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG LENNARD_JONES_BCOEF")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-23-")
         (fortran:fformat 5 "%16.8e")
         (loop for cn2 across cn2-vec
               do (fortran:fwrite cn2))
         (fortran:end-line))
        ;; write the Lennard Jones r**6 terms for all possible atom type interactions

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
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
         (fortran:end-line))
        
        ;; write IBH, JBH, ICBH

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
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
         (fortran:end-line))
        ;; write IB, JB, ICB

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
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
         (fortran:end-line))
        ;; write ITH, JTH, KTH, ICTH

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
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
         (fortran:end-line))
        ;; write IT, JT, KT, ICT

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
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
         (fortran:end-line))
        ;; write IPH, JPH, KPH, LPH, ICPH

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
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
         (fortran:end-line))
        ;; write IP, JP, KP, LP, ICP

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG EXCLUDED_ATOMS_LIST")
         (fortran:fwrite "%FORMAT(10I8)")
         (fortran:debug "-30-")
         (fortran:fformat 10 "%8d")
         (loop for atom across excluded-atom-list
               do (fortran:fwrite (+ atom 1)))
         (fortran:end-line))
        ;; write excluded atoms list

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG HBOND_ACOEF")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-31-")
         (fortran:fformat 5 "%16.8e")
         (fortran:fwrite 0.0)
         (fortran:end-line))
        ;;This term has been dropped from most modern force fields.

        ;; Next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG HBOND_BCOEF")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-32-")
         (fortran:fformat 5 "%16.8e")
         (fortran:fwrite 0.0)
         (fortran:end-line))
        ;;This term has been dropped from most modern force fields.

        ;; Next
        (fortran:fformat 1 "%-80s")
        (cando:progress-advance bar (incf bar-counter))
        (fortran:fwrite "%FLAG HBCUT")
        (fortran:fwrite "%FORMAT(5E16.8)")
        (fortran:debug "-33-")
        (fortran:fformat 5 "%16.8e")
        (fortran:fwrite 0.0)
        (fortran:end-line)
        ;;no longer used for anything.

        ;;next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG AMBER_ATOM_TYPE")
         (fortran:fwrite "%FORMAT(20A4)")
         (fortran:debug "-34-")
         (fortran:fformat 20 "%-4s")
         (loop for type across atom-type
               do (fortran:fwrite (string type)))
         (fortran:end-line))

        ;;next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG TREE_CHAIN_CLASSIFICATION")
         (fortran:fwrite "%FORMAT(20A4)")
         (fortran:debug "-35-")
         (fortran:fformat 20 "%-4s")
         (loop repeat natom
               do (fortran:fwrite "M"))
         (fortran:end-line))
        ;; We are not considering protein/DNA thing, so just put "main chain" for all atoms.

        ;;next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG JOIN_ARRAY")
         (fortran:fwrite "%FORMAT(10I8)")
         (fortran:debug "-36-")
         (fortran:fformat 10 "%8d")
         (loop repeat natom
               do (fortran:fwrite "0"))
         (fortran:end-line))
        ;;This section is no longer used and is currently just filled with zeros.

        ;;next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG IROTAT")
         (fortran:fwrite "%FORMAT(10I8)")
         (fortran:debug "-37-")
         (fortran:fformat 10 "%8d")
         (loop repeat natom
               do (fortran:fwrite "0"))
         (fortran:end-line))
        ;;This section is not used and is currently just filled with zeros.

        ;;next
        (outline-progn
         (if (chem:bounding-box-bound-p atom-table)
             (progn
               (fortran:fformat 1 "%-80s")
               (cando:progress-advance bar (incf bar-counter))
               (format t "solvent_pointers info: iptres nspm nspsol -> ~a ~a ~a~%"
                       (chem:final-solute-residue-iptres-bound-p atom-table)
                       (chem:total-number-of-molecules-nspm-bound-p atom-table)
                       (chem:first-solvent-molecule-nspsol-bound-p atom-table))
               (when (and (chem:final-solute-residue-iptres-bound-p atom-table)
                          (chem:total-number-of-molecules-nspm-bound-p atom-table)
                          (chem:first-solvent-molecule-nspsol-bound-p atom-table))
                 (format t "Generating SOLVENT_POINTERS~%")
                 (let ((final-solute-residue-iptres (chem:final-solute-residue-iptres atom-table))
                       (total-number-of-molecules-nspm (chem:total-number-of-molecules-nspm atom-table))
                       (first-solvent-molecule-nspsol (chem:first-solvent-molecule-nspsol atom-table)))
                   (fortran:fwrite "%FLAG SOLVENT_POINTERS")
                   (fortran:fwrite "%FORMAT(3I8)")
                   (fortran:debug "-38-")
                   (fortran:fformat 3 "%8d")
                   (fortran:fwrite final-solute-residue-iptres)
                   (fortran:fwrite total-number-of-molecules-nspm)
                   (fortran:fwrite first-solvent-molecule-nspsol)
                   (fortran:end-line)))
               (fortran:fformat 1 "%-80s")
               (cando:progress-advance bar (incf bar-counter))
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
               (cando:progress-advance bar (incf bar-counter))
               (fortran:fwrite "%FLAG BOX_DIMENSIONS")
               (fortran:fwrite "%FORMAT(5E16.8)")
               (fortran:debug "-40-")
               (fortran:fformat 5 "%16.8e")
               (let ((solvent-box (chem:bounding-box atom-table)))
                 (unless (and solvent-box (listp solvent-box)
                              (or (= (length solvent-box) 3)
                                  (= (length solvent-box) 6)))
                   (error "There must be a solvent-box property in the aggregate properties and it must be a list of length three (x,y,z) or six (x,y,z,angle1,angle2,angle3) numbers"))
                 (if (> (length solvent-box) 3)
                     (fortran:fwrite (float (fourth solvent-box)))
                     (fortran:fwrite "90.0000000") ;box angle
                     )
                 (fortran:fwrite (float (first solvent-box)))
                 (fortran:fwrite (float (second solvent-box)))
                 (fortran:fwrite (float (third solvent-box)))) 
               (fortran:end-line))))

        ;;next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG RADIUS_SET")
         (fortran:fwrite "%FORMAT(1a80)")
         (fortran:debug "-41-")
         (fortran:fformat 1 "%-80s")
         (fortran:fwrite "modified Bondi radii (mbondi)") ;default in leap
         (fortran:end-line))

        ;;next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG RADII")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-42-")
         (fortran:fformat 5 "%16.8e")
         (loop for radius across generalized-born-radius
               do (fortran:fwrite radius))
         (fortran:end-line))
        ;;Generalized Born intrinsic dielectric radii

        ;;next
        (outline-progn
         (fortran:fformat 1 "%-80s")
         (cando:progress-advance bar (incf bar-counter))
         (fortran:fwrite "%FLAG SCREEN")
         (fortran:fwrite "%FORMAT(5E16.8)")
         (fortran:debug "-43-")
         (fortran:fformat 5 "%16.8e")
         (loop for screen across generalized-born-screen
               do (fortran:fwrite screen))
         (fortran:end-line))
        (when cando-extensions
          (outline-progn
           (fortran:fformat 1 "%-80s")
           (cando:progress-advance bar (incf bar-counter))
           (fortran:fwrite "%FLAG BOND_ORDERS")
           (fortran:fwrite "%FORMAT(40I2)")
           (fortran:debug "-44-")
           (fortran:fformat 40 "%2d")
           (loop for orderi across non-h-bond-orders
                 do (fortran:fwrite orderi))
           (fortran:end-line))
          ;; If there is more than one force field or the only force-field is not
          ;; :DEFAULT - then write out the per-molecule force-field name
          (when (or (> (length force-field-names-vec) 1)
                    (not (eq (elt force-field-names-vec 0) :default)))
            (outline-progn
             (let ((max-force-field-name-len 0))
               (loop for name across force-field-names-vec
                     do (setf max-force-field-name-len (max max-force-field-name-len
                                                            (length (string name)))))
               (fortran:fformat 1 (format nil "%-80s"))
               (cando:progress-advance bar (incf bar-counter))
               (fortran:fwrite "%FLAG FORCE_FIELD_NAMES")
               (fortran:fwrite (format nil "%FORMAT(1a~d)" (1+ max-force-field-name-len)))
               (fortran:debug "-45-")
               (fortran:fformat 1 (format nil "%-~ds" (1+ max-force-field-name-len)))
               (loop for name across force-field-names-vec
                     do (fortran:fwrite (string name)))
               (fortran:end-line)))
            (outline-progn
             (fortran:fformat 1 "%-80s")
             (cando:progress-advance bar (incf bar-counter))
             (fortran:fwrite "%FLAG MOLECULE_FORCE_FIELD_INDEX")
             (fortran:fwrite (format nil "%FORMAT(20I3)"))
             (fortran:debug "-46-")
             (fortran:fformat 20 "%3d")
             (loop for index across molecule-force-field-name-indices
                   do (fortran:fwrite (1+ index)))
             (fortran:end-line))))
        ))
;;;    (format *debug-io* "coordinate-pathname -> ~s~%" coordinate-pathname)
      (fortran:with-fortran-output-file (ftop coordinate-pathname :direction :output :if-exists :supersede)
      (fortran:fformat 20 "%-4s")
      (fortran:fwrite (string (chem:aggregate-name atom-table)))
      (fortran:end-line)
                                        ;      (fortran:fformat 1 "%5d")
                                        ;      (fortran:fwrite natom)
                                        ;      (fortran:fformat 5 "%15.7lf")
                                        ;      (fortran:fwrite 0.0)
                                        ;      (fortran:fwrite 0.0)
      (fortran:fwrite (format nil  " ~5d~%" natom))
      (fortran:fformat 6 "%12.7lf")
      (let ((ox 0.0)
            (oy 0.0)
            (oz 0.0))
        (when (chem:bounding-box-bound-p atom-table)
          (let ((solvent-box (chem:bounding-box atom-table)))
            (unless (and solvent-box (listp solvent-box) (= (length solvent-box) 3))
              (error "There must be a solvent-box property in the aggregate properties and it must be a list of length three numbers"))
            (setf ox (/ (float (first solvent-box)) 2.0)
                  oy (/ (float (second solvent-box)) 2.0)
                  oz (/ (float (third solvent-box)) 2.0))))
        (loop for i from 0 below natom
              for atom = (chem:elt-atom atom-table i)
              for atom-coordinate-index-times3 = (chem:elt-atom-coordinate-index-times3 atom-table i)
              for pos = (chem:get-position atom)
;;;           do (format t "atom-coordinate-index-times3 -> ~a~%" atom-coordinate-index-times3)
              do (progn
                   (fortran:fwrite (+ (geom:vx pos) ox))
                   (fortran:fwrite (+ (geom:vy pos) oy))
                   (fortran:fwrite (+ (geom:vz pos) oz))))
        (fortran:end-line))
      ;; write out the solvent box
      (if (chem:bounding-box-bound-p atom-table)
          (let ((solvent-box (chem:bounding-box atom-table)))
            (unless (and solvent-box (listp solvent-box) (= (length solvent-box) 3))
              (error "There must be a solvent-box property in the aggregate properties and it must be a list of length three numbers"))
            (fortran:fwrite (float (first solvent-box)))
            (fortran:fwrite (float (second solvent-box)))
            (fortran:fwrite (float (third solvent-box)))
            (fortran:fwrite "90.0000000") ;box angle
            (fortran:fwrite "90.0000000") 
            (fortran:fwrite "90.0000000")))
      (fortran:end-line))
    (cando:progress-done bar)
    (values energy-function)))

(defun save-amber-parm-format (aggregate topology-pathname coordinate-pathname &key assign-types residue-name-to-pdb-alist (cando-extensions t))
  (format t "Constructing energy function~%")
  (finish-output)
  (let* ((energy-function (chem:make-energy-function aggregate 
                                                     :use-excluded-atoms t
                                                     :assign-types assign-types)))
    ;;; We need to:
    ;;;  (1) make sure energy function copies :bounding-box property from aggregate
    ;;;  (2) Copy the name of the aggregate into the energy function
    ;;;  (3) Separate the solvent molecules from solute molecules and order them in the energy-function
    ;;;  (4) Copy the result of (chem:lookup-nonbond-force-field-for-aggregate aggregate force-field) into the energy-function
    (save-amber-parm-format-using-energy-function energy-function
                                                  topology-pathname
                                                  coordinate-pathname
                                                  residue-name-to-pdb-alist
                                                  :cando-extensions cando-extensions)))

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
(defvar %flag-solvent-pointers "%FLAG SOLVENT_POINTERS")
(defvar %flag-atoms-per-molecule "%FLAG ATOMS_PER_MOLECULE")
(defvar %flag-bond-orders "%FLAG BOND_ORDERS")
(defvar %flag-force-field-names "%FLAG FORCE_FIELD_NAMES")
(defvar %flag-molecule-force-field-index "%FLAG MOLECULE_FORCE_FIELD_INDEX")

(defun verify-%flag-line (line)
  (unless (string-equal line "%FLAG" :start1 0 :end1 5)
    (error "Expected %FLAG at the start of the line - got: ~s" line)))

#+(or)
(defmacro rlog (fmt &rest args)
  `(progn
           (cl:format *debug-io* ,fmt ,@args)
           (finish-output *debug-io*)))
(defmacro rlog (fmt &rest args)
  nil)

(defun generate-aggregate-for-energy-function (energy-function)
  (let* ((aggregate (chem:make-aggregate))
         (atom-table (chem:atom-table energy-function))
         (molecules (chem:atom-table-molecules atom-table))
         (aggregate (chem:make-aggregate)))
    (loop for mol across molecules
          do (chem:add-matter aggregate mol))
    aggregate))

;(defun read-amber-parm-format (stream)
;  (let ((fif (fortran:make-fortran-input-file :stream stream))
(defun read-amber-parm-format (topology-pathname)
  "Return (values energy-function) - use generate-aggregate-for-energy-function to get an aggregate"
  (fortran:with-fortran-input-file (fif topology-pathname :direction :input)
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
          excluded-atoms-list amber-atom-type solvent-pointers atoms-per-molecule
          non-h-bond-orders
          molecules-vec residues-vec
          force-field-names molecule-force-field-index
          )
      (rlog "Starting read-amber-parm-format~%")
      (fortran:fread-line fif)   ; Skip the version and timestamp line
      (fortran:fread-line fif)   ; read the first %FLAG line
      ;; From here on down - read the input file
      (outline-progn
       (loop for line = (fortran:fortran-input-file-look-ahead fif)
             while line
             ;;do (rlog "line ~a~%" line)
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
                                  (intern (string-trim " " string) "KEYWORD"))
                                (fortran:fread-vector fif per-line format-char width)))))
                  ((string-equal %flag-charge line :end2 (length %flag-charge))
                   (fortran:fread-line-or-error fif)  
                   (multiple-value-bind (per-line format-char width decimal)
                       (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
                     (fortran:fread-line-or-error fif) 
                     (setf charge (fortran:fread-vector fif per-line format-char width))
                     (loop for chargei from 0 below (length charge)
                           do (setf (aref charge chargei) (/ (aref charge chargei) +amber-charge-conversion-18.2223+)))))
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
                                  (intern (string-trim " " string) "KEYWORD"))
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
                     (let ((type-strings (fortran:fread-vector fif per-line format-char width)))
                       (setf amber-atom-type (make-array (length type-strings) :adjustable nil))
                       (loop for index from 0 below (length type-strings)
                             for type-string = (string-trim " " (elt type-strings index))
                             do (setf (elt amber-atom-type index) (intern type-string :keyword))))))
                  ((string-equal %flag-solvent-pointers line :end2 (length %flag-solvent-pointers))
                   (fortran:fread-line-or-error fif)  
                   (multiple-value-bind (per-line format-char width decimal)
                       (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
                     (fortran:fread-line-or-error fif) 
                     (setf solvent-pointers (fortran:fread-vector fif per-line format-char width))))
                  ((string-equal %flag-atoms-per-molecule line :end2 (length %flag-atoms-per-molecule))
                   (fortran:fread-line-or-error fif)  
                   (multiple-value-bind (per-line format-char width decimal)
                       (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
                     (fortran:fread-line-or-error fif) 
                     (setf atoms-per-molecule (fortran:fread-vector fif per-line format-char width))))
                  ((string-equal %flag-bond-orders line :end2 (length %flag-bond-orders))
                   (fortran:fread-line-or-error fif)  
                   (multiple-value-bind (per-line format-char width decimal)
                       (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
                     (fortran:fread-line-or-error fif) 
                     (setf non-h-bond-orders (fortran:fread-vector fif per-line format-char width))))
                  ((string-equal %flag-force-field-names line :end2 (length %flag-force-field-names))
                   (fortran:fread-line-or-error fif)  
                   (multiple-value-bind (per-line format-char width decimal)
                       (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
                     (fortran:fread-line-or-error fif) 
                     (setf force-field-names (fortran:fread-vector fif per-line format-char width))
                     (loop for index from 0
                           for ffname across force-field-names
                           do (setf (elt force-field-names index) (intern (string-trim " " ffname) :keyword)))))
                  ((string-equal %flag-molecule-force-field-index line :end2 (length %flag-molecule-force-field-index))
                   (fortran:fread-line-or-error fif)  
                   (multiple-value-bind (per-line format-char width decimal)
                       (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
                     (fortran:fread-line-or-error fif) 
                     (setf molecule-force-field-index (fortran:fread-vector fif per-line format-char width))))
                  (t
                   (if chem:*verbose* (cl:format t "Unknown flag ~a~%" line))
                   (fortran:fread-line-or-error fif)  
                   (multiple-value-bind (per-line format-char width decimal)
                       (fortran:parse-fortran-format-line (fortran:fortran-input-file-look-ahead fif))
                     (fortran:fread-line-or-error fif) 
                     (fortran:fread-vector fif per-line format-char width)))))
       )
      ;;(rlog "natom -> ~s~%" natom)
      ;;(rlog "ntypes -> ~s~%" ntypes)
      ;;(rlog "nbonh -> ~s~%" nbonh)
      ;;(rlog "mbona -> ~s~%" mbona)
      ;;(rlog "ntheth -> ~s~%" ntheth)
      ;;(rlog "mtheta -> ~s~%" mtheta)
      ;;(rlog "nphih -> ~s~%" nphih)
      ;;(rlog "mphia -> ~s~%" mphia)
      ;;(rlog "nhparm -> ~s~%" nhparm)
      ;;(rlog "nparm -> ~s~%" nparm)
      ;;(rlog "nnb -> ~s~%" nnb)
      ;;(rlog "nres -> ~s~%" nres)
      ;;(rlog "nbona -> ~s~%" nbona)
      ;;(rlog "ntheta -> ~s~%" ntheta)
      ;;(rlog "nphia -> ~s~%" nphia)
      ;;(rlog "numbnd -> ~s~%" numbnd)
      ;;(rlog "numang -> ~s~%" numang)
      ;;(rlog "nptra -> ~s~%" nptra)
      ;;(rlog "natyp -> ~s~%" natyp)
      ;;(rlog "nphb -> ~s~%" nphb)
      ;;(rlog "ifpert -> ~s~%" ifpert)
      ;;(rlog "nbper -> ~s~%" nbper)
      ;;(rlog "ngper -> ~s~%" ngper)
      ;;(rlog "ndper -> ~s~%" ndper)
      ;;(rlog "mbper -> ~s~%" mbper)
      ;;(rlog "mgper -> ~s~%" mgper)
      ;;(rlog "mdper -> ~s~%" mdper)
      ;;(rlog "ifbox -> ~s~%" ifbox)
      ;;(rlog "nmxrs -> ~s~%" nmxrs)
      ;;(rlog "ifcap -> ~s~%" ifcap)
      ;;(rlog "numextra -> ~s~%" numextra)
      ;;(rlog "ncopy -> ~s~%" ncopy)
      ;;(rlog "atom-name -> ~s~%" atom-name)
      ;;(rlog "charge -> ~s~%" charge)
      ;;(rlog "atomic-number -> ~s~%" atomic-number)
      ;;(rlog "mass -> ~s~%" mass)
      ;;(rlog "atom-type-index -> ~s~%" atom-type-index)
      ;;(rlog "number-excluded-atoms -> ~s~%" number-excluded-atoms)
      ;;(rlog "nonbonded-parm-index -> ~s~%" nonbonded-parm-index)
      ;;(rlog "residue-label -> ~s~%" residue-label)
      ;;(rlog "residue-pointer -> ~s~%" residue-pointer)
      ;;(rlog "bond-force-constant -> ~s~%" bond-force-constant)
      ;;(rlog "bond-equil-value -> ~s~%" bond-equil-value)
      ;;(rlog "angle-force-constant -> ~s~%" angle-force-constant)
      ;;(rlog "angle-equil-value -> ~s~%" angle-equil-value)
      ;;(rlog "dihedral-force-constant -> ~s~%" dihedral-force-constant)
      ;;(rlog "dihedral-periodicity -> ~s~%" dihedral-periodicity)
      ;;(rlog "dihedral-phase -> ~s~%" dihedral-phase)
      ;;(rlog "scee-scale-factor -> ~s~%" scee-scale-factor)
      ;;(rlog "scnb-scale-factor -> ~s~%" scnb-scale-factor)
      ;;(rlog "solty -> ~s~%" solty)
      ;;(rlog "lennard-jones-acoef -> ~s~%" lennard-jones-acoef)
      ;;(rlog "lennard-jones-bcoef -> ~s~%" lennard-jones-bcoef)
      ;;(rlog "bonds-inc-hydrogen -> ~s~%" bonds-inc-hydrogen)
      ;;(rlog "bonds-without-hydrogen -> ~s~%" bonds-without-hydrogen)
      ;;(rlog "angles-inc-hydrogen -> ~s~%" angles-inc-hydrogen)
      ;;(rlog "angles-without-hydrogen -> ~s~%" angles-without-hydrogen)
      ;;(rlog "dihedrals-inc-hydrogen -> ~s~%" dihedrals-inc-hydrogen)
      ;;(rlog "dihedrals-without-hydrogen -> ~s~%" dihedrals-without-hydrogen)
      ;;(rlog "excluded-atoms-list -> ~s~%" excluded-atoms-list)
      ;;(rlog "amber-atom-type -> ~s~%" amber-atom-type)
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
            (atoms (make-array natom)))
        ;; Create a vector of atoms to pass to the atom-table and to set the atoms for stretch, angle, dihedral etc.
        (rlog "Create atoms vector~%")
        (loop for i from 0 below natom
              for name = (aref atom-name i)
              for atom-element = (chem:element-for-atomic-number (aref atomic-number i))
              do (setf (aref atoms i)  (chem:make-atom name atom-element))
              )
        (setf (aref residue-pointer (length residue-pointer)) (+ 1 natom))
        ;;(format t "residue-pointer ~s~%" residue-pointer)
	(setf residues-vec (make-array (length residue-label) :element-type t :adjustable nil))
        (setf molecules-vec (make-array 256 :element-type t :fill-pointer 0 :adjustable t))
        (let (residue-accumulate
              (atoms-in-molecule 0)
              (molecule-index 0))
          ;; Figure out the residues and molecules for the resulting aggregate
          (loop for i from 0 below (length residue-label)
                for name = (aref residue-label i)
                for begin-atom-index = (1- (aref residue-pointer i))
                for end-atom-index = (1- (aref residue-pointer (1+ i)))
                do (let ((residue (chem:make-residue name)))
                     (loop for atomi from begin-atom-index below end-atom-index
                           for atom = (aref atoms atomi)
                           do (chem:add-matter residue atom))
		     (setf (aref residues-vec i) residue)
                     (push residue residue-accumulate)
                     (incf atoms-in-molecule (- end-atom-index begin-atom-index))
                     (when (= atoms-in-molecule (elt atoms-per-molecule molecule-index))
                       (let ((molecule (chem:make-molecule)))
                         (when solvent-pointers
                           (if (>= (1+ molecule-index) (elt solvent-pointers 2))
                               (chem:setf-molecule-type molecule :solvent)
                               (chem:setf-molecule-type molecule :solute)))
                         (mapc (lambda (res) (chem:add-matter molecule res)) (nreverse residue-accumulate))
                         (vector-push-extend molecule molecules-vec)
                         (setf atoms-in-molecule 0)
                         (setf residue-accumulate nil)
                         (incf molecule-index))))))
        (when molecule-force-field-index
          (loop for moli from 0 below (length molecules-vec)
                for molecule = (elt molecules-vec moli)
                for force-field-index = (1- (elt molecule-force-field-index moli))
                for force-field-name = (elt force-field-names force-field-index)
                do (chem:setf-force-field-name molecule force-field-name)))
        ;; Identify the force-field types
        ;; (warn "This is where I set the force-field name for each molecule")
        #+(or)(let ((force-field-index 0))
          (setf force-field-names (make-hash-table))
          (loop for molecule across molecules-vec
                for force-field-name = (chem:get-force-field molecule)
                do (if (gethash force-field-name force-field-names)
                       nil
                       (setf (gethash force-field-name force-field-names) (prog1 force-field-index
                                                                            (incf force-field-index))))))
        ;;(rlog "atoms -> ~s~%" atoms)
        (rlog "Create stretch vectors~%")
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
                            (incf counts))))
        ;; Now form the bonds with the correct bond orders
        ; (format t "read-amber-parm about to create bonds~%")
        (loop for bondi from 0 below nbonh
              for atom1-index = (/ (aref bonds-inc-hydrogen (* bondi 3)) 3)
              for atom2-index = (/ (aref bonds-inc-hydrogen (+ (* bondi 3) 1)) 3)
              for atom1 = (aref atoms atom1-index)
              for atom2 = (aref atoms atom2-index)
              do (chem:bond-to atom1 atom2 :single-bond)) ; All bonds to hydrogen are single bonds
        (if non-h-bond-orders
            (loop for bondi from 0 below mbona
              for atom1-index = (/ (aref bonds-without-hydrogen (* bondi 3)) 3)
              for atom2-index = (/ (aref bonds-without-hydrogen (+ (* bondi 3) 1)) 3)
              for atom1 = (aref atoms atom1-index)
              for atom2 = (aref atoms atom2-index)
              for bond-order-int = (aref non-h-bond-orders bondi)
                  do (chem:bond-to-order-int atom1 atom2 bond-order-int))
            (loop for bondi from 0 below mbona
              for atom1 = (aref atoms (/ (aref bonds-without-hydrogen (* bondi 3)) 3))
              for atom2 = (aref atoms (/ (aref bonds-without-hydrogen (+ (* bondi 3) 1)) 3))
                  do (chem:bond-to atom1 atom2 :unknown-order-bond)))
        ;; (format t "read-amber-parm done creating bonds~%")
        (setf stretch-vectors (acons :kb kbs-vec stretch-vectors))
        (setf stretch-vectors (acons :r0 r0s-vec stretch-vectors))
        (setf stretch-vectors (acons :i1 i1s-vec stretch-vectors))
        (setf stretch-vectors (acons :i2 i2s-vec stretch-vectors))
        (setf stretch-vectors (acons :atom1 atom1s-vec stretch-vectors))
        (setf stretch-vectors (acons :atom2 atom2s-vec stretch-vectors))
        ;;(rlog "stretch-vectors -> ~s~%" stretch-vectors)
        (chem:fill-from-vectors-in-alist energy-stretch stretch-vectors)
        (rlog "Create angle vectors~%")
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
        ;;(rlog "angle-vectors -> ~s~%" angle-vectors)
        (chem:fill-from-vectors-in-alist energy-angle angle-vectors)

        (rlog "Create dihedral vectors~%")
        (outline-progn
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
               ))
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
        ;;(rlog "dihedral-vectors -> ~s~%" dihedral-vectors)
        (chem::fill-from-vectors-in-alist energy-dihedral dihedral-vectors)
        ;;atom-table
        (rlog "Create atom-table vectors~%")
        (setf atom-table-vectors (acons :atom-vector atoms atom-table-vectors))
        (setf atom-table-vectors (acons :atom-name-vector atom-name atom-table-vectors))
        (setf atom-table-vectors (acons :atom-type-vector amber-atom-type atom-table-vectors))
        (setf atom-table-vectors (acons :charge-vector charge atom-table-vectors))
        (setf atom-table-vectors (acons :mass-vector mass atom-table-vectors))
        (setf atom-table-vectors (acons :atomic-number-vector atomic-number atom-table-vectors))
                                        ;     (setf atom-table-vectors (acons :residues residue-pointer atom-table-vectors))
                                        ;     (setf atom-table-vectors (acons :residue-names residue-labels atom-table-vectors))
          
        ;;(rlog "atom-table-vectors -> ~s~%" atom-table-vectors)
        (chem::fill-atom-table-from-vectors atom-table atom-table-vectors)

        ;;nonbond
        (rlog "Create nonbond vectors~%")
        (setf nonbond-vectors (acons :ntypes ntypes nonbond-vectors))
        (setf nonbond-vectors (acons :atom-name-vector (copy-seq atom-name) nonbond-vectors))
        (setf nonbond-vectors (acons :atom-type-vector amber-atom-type nonbond-vectors))
        (setf nonbond-vectors (acons :charge-vector (copy-seq charge) nonbond-vectors))
        (setf nonbond-vectors (acons :mass-vector (copy-seq mass) nonbond-vectors))
        (setf nonbond-vectors (acons :atomic-number-vector (copy-seq atomic-number) nonbond-vectors))
        (setf nonbond-vectors (acons :ico-vec (copy-seq nonbonded-parm-index) nonbond-vectors))
        (setf nonbond-vectors (acons :iac-vec (copy-seq atom-type-index) nonbond-vectors))
        (setf nonbond-vectors (acons :cn1-vec (copy-seq lennard-jones-acoef) nonbond-vectors))
        (setf nonbond-vectors (acons :cn2-vec (copy-seq lennard-jones-bcoef) nonbond-vectors))
        
        ;;(rlog "nonbond-vectors -> ~s~%" nonbond-vectors)
        (chem:construct-nonbond-terms-from-aList energy-nonbond nonbond-vectors)
        (loop for i from 0 below (length excluded-atoms-list)
              for atom = (- (aref excluded-atoms-list i) 1)
              do (setf (aref excluded-atoms-list i) atom))
        (chem:set-nonbond-excluded-atom-info energy-nonbond atom-table (copy-seq excluded-atoms-list) (copy-seq number-excluded-atoms))

        ;; for energy nonbond test
        #+(or)(chem:expand-excluded-atoms-to-terms energy-nonbond)
        ;;

        ;; Now we have energy-stretch, energy-angle, and energy-dihedral
        ;;   we want to put them into an energy-function
        ;;        (rlog "atom-table-vectors -> ~s~%" atom-table-vectors)
        ;;
        #+(or) ;; the old code
        (let ((energy-function (core:make-cxx-object 'chem:energy-function
                                                     :stretch energy-stretch
                                                     :angle energy-angle
                                                     :dihedral energy-dihedral
                                                     :nobond energy-nonbond
                                                     :atom-table atom-table)))
          energy-function #| <-- This was missing before |# )
	;; fill in atom-table information
	(chem:setf-atom-table-residue-pointers atom-table residue-pointer)
        (chem:setf-atom-table-residue-names atom-table (make-array (length residue-label) :adjustable t :initial-contents residue-label))
	(chem:setf-atom-table-residue-pointers atom-table residue-pointer)
        (if solvent-pointers
            (progn
              (chem:set-first-solvent-molecule-nspsol atom-table (elt solvent-pointers 2))
              (chem:set-final-solute-residue-iptres atom-table (elt solvent-pointers 0))
              (chem:set-total-number-of-molecules-nspm atom-table (elt solvent-pointers 1))))
  	(chem:setf-atom-table-atoms-per-molecule atom-table atoms-per-molecule)
	(chem:setf-atom-table-residues atom-table residues-vec)
        (chem:setf-atom-table-molecules atom-table (copy-seq molecules-vec))
	;; more here
        (let ((alist (list (cons :atom-table atom-table)
                           (cons :stretch energy-stretch)
                           (cons :angle energy-angle)
                           (cons :dihedral energy-dihedral)
                           (cons :nonbond energy-nonbond)))
              (energy-function (core:make-cxx-object 'chem:energy-function)))
          (chem:fill-energy-function-from-alist energy-function alist)
          (values energy-function (generate-aggregate-for-energy-function energy-function )))))))


(defun read-amber-ascii-restart-file (coordinate-filename &key read-velocities (read-bounding-box t))
  "Return (values number-of-atoms coordinates bounding-box-or-nil velocities-or-nil)"
  (fortran:with-fortran-input-file (fif coordinate-filename :direction :input)
    (fortran:fread-line fif)     ; Skip the version and timestamp line
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
      (let* ((coordinates (copy-seq (fortran:fread-double-float-vector fif 6 12 (* natoms 3))))
             (velocities (when read-velocities
                           (copy-seq (fortran:fread-double-float-vector fif 6 12 (* natoms 3)))))
             (bounding-box (when read-bounding-box
                             (copy-seq (fortran:fread-double-float-vector fif 6 12 6)))))
        (values natoms coordinates (coerce bounding-box 'list) velocities)))))

;;; The following code is to generate a human readable representation of an energy-function
;;; with everything sorted so that the terms can be compared side-by-side using something like
;;; 'ediff'

(defstruct stretch-term
  atom1-name atom2-name atom1-type atom2-type kb r0)


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
                           for atom1-name = (chem:get-name atom1)
                           for atom2-name = (chem:get-name atom2)
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
                 (rotatef (stretch-term-atom2-name stretch)
                          (stretch-term-atom1-name stretch))
                 (rotatef (stretch-term-atom2-type stretch)
                          (stretch-term-atom1-type stretch))))
      (flet ((order-stretch (s1 s2)
               (if (string< (string (stretch-term-atom1-name s1))
                            (string (stretch-term-atom1-name s2)))
                   t
                   (if (string> (string (stretch-term-atom1-name s1))
                                (string (stretch-term-atom1-name s2)))
                       nil
                       (string< (string (stretch-term-atom2-name s1))
                                (string (stretch-term-atom2-name s2)))))))
        (sort stretches #'order-stretch))
      (with-open-file (stream "/tmp/stretch.dat" :direction :output)
        (format stream "stretch ~s~%" stretches)))))
;;    (rlog "stretch-vectors2 ~s~%" stretches)


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

(defun float= (x y &optional (delta 0.01))
  (let ((absx (abs x)))
    (if (< absx 1e-6)
        (< (abs y) 1e-6)
        (< (/ (abs (- x y)) absx) delta))))
               
(defun compare-energy-nonbond-parameters (energy-function1 energy-function2)
  (let* ((energy-nonbond1 (chem:get-nonbond-component  energy-function1))
         (energy-nonbond2 (chem:get-nonbond-component energy-function2))
         (nonbond-vectors1 (chem:extract-vectors-as-alist energy-nonbond1))
         (nonbond-vectors2 (chem:extract-vectors-as-alist energy-nonbond2)))
    (loop for i from 0 below (length (cdr (assoc :da nonbond-vectors1)))
          for da1 = (aref (cdr (assoc :da nonbond-vectors1)) i)
          for da2 = (aref (cdr (assoc :da nonbond-vectors2)) i)
          for dc1 = (aref (cdr (assoc :dc nonbond-vectors1)) i)
          for dc2 = (aref (cdr (assoc :dc nonbond-vectors2)) i)
          for charge11 = (aref (cdr (assoc :charge1 nonbond-vectors1)) i)
          for charge12 = (aref (cdr (assoc :charge2 nonbond-vectors1)) i)
          for charge21 = (aref (cdr (assoc :charge1 nonbond-vectors2)) i)
          for charge22 = (aref (cdr (assoc :charge2 nonbond-vectors2)) i)
          for i11 = (aref (cdr (assoc :i1 nonbond-vectors1)) i)
          for i12 = (aref (cdr (assoc :i2 nonbond-vectors1)) i)
          for i21 = (aref (cdr (assoc :i1 nonbond-vectors2)) i)
          for i22 = (aref (cdr (assoc :i2 nonbond-vectors2)) i)
          for atomname11 = (aref (cdr (assoc :atomname1 nonbond-vectors1)) i)
          for atomname12 = (aref (cdr (assoc :atomname2 nonbond-vectors1)) i)
          for atomname21 = (aref (cdr (assoc :atomname1 nonbond-vectors2)) i)
          for atomname22 = (aref (cdr (assoc :atomname2 nonbond-vectors2)) i)
           do (when (or (not (float= da1 da2))
                       (not (float= dc1 dc2))
                       (not (float= charge11 charge21))
                       (not (float= charge12 charge22))
                       (/= i11 i21)
                       (/= i12 i22)
                       )
               (format t "line ~a da1: ~a dc1: ~a charge11: ~a charge12: ~a atomname11: ~a atomname12: ~a~%"
                       i da1 dc1 charge11 charge12 atomname11 atomname12)
               (format t "line ~a da2: ~a dc2: ~a charge21: ~a charge22: ~a atomname21: ~a atomname22: ~a~%"
                       i da2 dc2 charge21 charge22 atomname21 atomname22)))))


(defclass amber-topology-coord-pair ()
 ((topology-filename :initarg :topology-filename :accessor topology-filename)
   (coordinate-filename :initarg :coordinate-filename :accessor coordinate-filename)
   (energy-function :initarg :energy-function :accessor energy-function)
   (aggregate :initarg :aggregate :accessor aggregate)
   (number-of-atoms :initarg :number-of-atoms :accessor number-of-atoms)
   (current-coordinates :initarg :current-coordinates :accessor current-coordinates)))

(defmethod cando:agg ((object amber-topology-coord-pair))
  (aggregate object))

(defclass amber-topology-restart-pair (amber-topology-coord-pair)
  ())

(defmethod cando:agg ((object amber-topology-restart-pair))
  (aggregate object))

(defclass amber-topology-trajectory-pair (amber-topology-coord-pair)
  ((netcdf :initarg :netcdf :accessor netcdf)))

(defun write-coordinates-into-energy-function-atom-table (energy-function coordinates)
  (let ((atom-table (chem:atom-table energy-function)))
    (loop for index from 0 below (/ (length coordinates) 3)
          for coord-index = 0 then (+ coord-index 3)
          for xpos = (elt coordinates coord-index)
          for ypos = (elt coordinates (+ 1 coord-index))
          for zpos = (elt coordinates (+ 2 coord-index))
          for atom = (chem:elt-atom atom-table index)
          do (chem:set-position-xyz atom xpos ypos zpos))))


(defun cell-lengths (netcdf)
  (let ((sv (static-vectors:make-static-vector 3 :element-type 'double-float)))
    (netcdf:get-vara-double netcdf "cell_lengths" (vector 0) (vector 3) sv)
    (copy-seq sv)))

(defun cell-angles (netcdf)
  (let ((sv (static-vectors:make-static-vector 3 :element-type 'double-float)))
    (netcdf:get-vara-double netcdf "cell_angles" (vector 0) (vector 3) sv)
    (copy-seq sv)))


(defun read-bounding-box (netcdf)
  (let* ((lengths (cell-lengths netcdf))
         (angles (cell-angles netcdf))
         (bounding-box (append (coerce lengths 'list) (coerce angles 'list))))
    bounding-box))

(defun read-amber-restart-file (coordinate-filename)
  (unless (probe-file coordinate-filename)
    (error "Could not open restart file ~a" coordinate-filename))
  (let ((fin (open coordinate-filename :direction :input)))
    (let* ((c1 (read-char fin))
           (c2 (read-char fin))
           (c3 (read-char fin))
           (c4 (read-char fin))
           (is-netcdf (and (char= c1 #\C)
                           (char= c2 #\D)
                           (char= c3 #\F)
                           (<= (char-int c4) 10))))
      (close fin)
      (if is-netcdf
          (let* ((netcdf (netcdf:nc-open coordinate-filename :mode netcdf-cffi:+nowrite+))
                 (number-of-atoms (netcdf:len (gethash "atom" (netcdf:dimensions netcdf))))
                 (coords (static-vectors:make-static-vector (* number-of-atoms 3) :element-type 'single-float))
                 (bounding-box (read-bounding-box netcdf))
                 (result (netcdf:get-vara-float netcdf
                                                "coordinates"
                                                (vector 0 0)
                                                (vector number-of-atoms 3)
                                                coords)))
            (unless (= result 0)
              (error "Could not read coordinates result -> ~a" result))
            (netcdf:nc-close netcdf)
            (values number-of-atoms coords bounding-box))
          (read-amber-ascii-restart-file coordinate-filename)))))

(defun load-amber-topology-restart-pair (&key topology-filename coordinate-filename)
  (multiple-value-bind (energy-function aggregate)
      (read-amber-parm-format topology-filename)
    (multiple-value-bind (number-of-atoms coordinates bounding-box)
        (read-amber-restart-file coordinate-filename)
      (write-coordinates-into-energy-function-atom-table energy-function coordinates)
      (when bounding-box
        (chem:set-property aggregate :bounding-box bounding-box))
      (make-instance 'amber-topology-restart-pair
                     :topology-filename topology-filename
                     :coordinate-filename coordinate-filename
                     :energy-function energy-function
                     :aggregate aggregate
                     :number-of-atoms number-of-atoms
                     :current-coordinates coordinates))))

(defun load-amber-topology-trajectory-pair (&key topology-filename coordinate-filename)
  (multiple-value-bind (energy-function aggregate)
      (read-amber-parm-format topology-filename)
    (let* ((crd (netcdf:nc-open coordinate-filename :mode netcdf-cffi:+nowrite+))
           (number-of-atoms (netcdf:len (gethash "atom" (netcdf:dimensions crd))))
           (coords (static-vectors:make-static-vector (* number-of-atoms 3) :element-type 'single-float))
           (pair (make-instance 'amber-topology-trajectory-pair
                                :topology-filename topology-filename
                                :coordinate-filename coordinate-filename
                                :netcdf crd
                                :energy-function energy-function
                                :aggregate aggregate
                                :number-of-atoms number-of-atoms
                                :current-coordinates coords)))
      (gctools:finalize pair (lambda ()
                               (format t "Closing netcdf file~%")
                               (netcdf:nc-close (netcdf pair))))
      pair)))

(defun change-coordinate-file (amber-topology-pair coordinate-filename)
  "Switch to another coordinate file"
  (when (netcdf amber-topology-pair)
    (netcdf:nc-close (netcdf amber-topology-pair)))
  (let ((crd (netcdf:nc-open coordinate-filename :mode netcdf-cffi:+nowrite+)))
    (setf (netcdf amber-topology-pair) crd
          (coordinate-filename amber-topology-pair) coordinate-filename)))
  
(defun number-of-frames (amber-topology-pair)
  "Return the number of frames in the coordinate file.
If it's a restart file then return NIL"
  (let ((dim (gethash "frame" (netcdf:dimensions (netcdf amber-topology-pair)))))
    (if dim
        (netcdf:len dim)
        nil)))

(defun read-frame-into-atoms (amber-topology-pair &optional frame-index)
  "Read a frame of coordinates into the atom positions and return the aggregate with the new atom positions. DO NOT MODIFY THIS AGGREGATE"
  (let ((number-of-frames (number-of-frames amber-topology-pair)))
    (when number-of-frames
      (when (>= frame-index number-of-frames)
        (error "Tried to load a frame beyond the end of the netcdf file - only ~a frames available" number-of-frames)))
    (let ((result (netcdf:get-vara-float (netcdf amber-topology-pair)
                                         "coordinates"
                                         (vector frame-index 0 0)
                                         (vector 1 (number-of-atoms amber-topology-pair) 3)
                                         (current-coordinates amber-topology-pair))))
      (write-coordinates-into-energy-function-atom-table (energy-function amber-topology-pair) (current-coordinates amber-topology-pair))
      (aggregate amber-topology-pair))))

