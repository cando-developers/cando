(in-package :leap.topology)

(defclass amber-topology ()
  ((aggregate :initarg :aggregate :accessor aggregate)
   ))

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
    (format t "Number of stretch terms -> ~a~%" (length kb-vector))
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
	     do (format t "atom1: ~a atom2: ~a i1: ~a i2: ~a~% " atom1 atom2 i1 i2)
	     do (if (or (eq (chem:get-element atom1) :h)
			(eq (chem:get-element atom2) :h))
		    (progn
		      (setf (aref ibh curh) i1)
		      (setf (aref jbh curh) i2)
		      (setf (aref icbh curh) (aref j-vec i))
		      (incf curh))
		    (progn
		      (setf (aref ib cur) i1)
		      (setf (aref jb cur) i2)
		      (setf (aref icb cur) (aref j-vec i))
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
       do (format t "atom1: ~a atom2: ~a atom3: ~a  key: ~a~% " atom1 atom2 atom3 key)
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
                    (setf (aref icth curh) (aref j-vec i))
                    (incf curh))
                  (progn
                    (setf (aref it cur) i1)
                    (setf (aref jt cur) i2)
                    (setf (aref kt1 cur) i3)
                    (setf (aref ict cur) (aref j-vec i))
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
       for in = (aref in-vec i)
       for proper = (aref proper-vec i)
       for key = (if proper
                     (canonical-dihedral-key type1 type2 type3 type4 in)
                     (canonical-improper-key atom1 atom2 atom3 atom4))
       for vi = (aref v-vec i)
       for ini = (aref in-vec i)
       for phasei = (aref phase-vec i)
       do (format t "atom1: ~a atom2: ~a atom3: ~a atom4 ~a key: ~a vi: ~a ini ~a phase ~a~% " atom1 atom2 atom3 atom4 key vi ini phasei)
       do (if (setf jtemp (gethash key uniques))
              (vector-push-extend jtemp j-vec)
              (progn
                (setf (gethash key uniques) jnext)
                (vector-push-extend jnext j-vec)
                (vector-push-extend vi vj-vec)
                (vector-push-extend ini inj-vec)
                (vector-push-extend phasei phasej-vec)
                (vector-push-extend proper properj-vec)
                (format t " ~a ~a~% " vi ini)
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

(Defun prepare-amber-energy-dihedral (energy-function)
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
    (let ((indices (if (check-if-dihedrals-are-ordered proper-vector i1-vector i2-vector i3-vector i4-vector in-vector)
                       (make-indices (length v-vector))
                       (sort-dihedrals proper-vector i1-vector i2-vector i3-vector i4-vector in-vector)))))
      (let (proper-prev i1prev i2prev i3prev i4prev)
        (loop for x below (length v-vector)
           for properx = (aref proper-vector x)
           for i1x = (aref i1-vector x)
           for i2x = (aref i2-vector x)
           for i3x = (aref i3-vector x)
           for i4x = (aref i4-vector x)
           do (when (and (eq properx proper-prev)
                         (eq i1x i1prev)
                         (eq i2x i2prev)
                         (eq i3x i3prev)
                         (eq i4x i4prev))
                (setf (aref i3-vector x) (- i3x)))
           do (when (not properx)
                (setf (aref i4-vector x) (- i4x)))
           do (setf proper-prev properx
                    i1prev i1x
                    i2prev i2x
                    i3prev i3x
                    i4prev i4x))
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
                                (setf (aref icph curh) (aref j-vec i))
                                (incf curh))
                              (progn
                                (setf (aref ip cur) i1)
                                (setf (aref jp cur) i2)
                                (setf (aref kp cur) i3)
                                (setf (aref lp cur) i4)
                                (setf (aref icp cur) (aref j-vec i))
                                (incf cur))))
                    (values with-h without-h iph jph kph lph icph ip jp kp lp icp vj-vec inj-vec phasej-vec)))))))

(defun canonical-nonbond-key (type1 type2)
  (declare (symbol type1 type2))
  (intern (format nil "~a-~a" (string type1) (string type2))
          :keyword))

(defun generate-nonbond-parameters (ffnonbond-db type-index-vector)
;(defun generate-nonbond-parameters (type-index-vector)
  (let ((iac-vec (make-array 256 :fill-pointer 0 :adjustable t))
        (jnext 1)
        jtemp
        (uniques (make-hash-table :test #'eq))
        (type-indexj-vec (make-array 256 :element-type 'fixnum :fill-pointer 0 :adjustable t))
        (local-typej-vec (make-array 256 :element-type 'fixnum :fill-pointer 0 :adjustable t))
        (ntypes 0))
    (loop for i from 0 below (length type-index-vector)
       for index = (aref type-index-vector i)
       do (if (setf jtemp (gethash index uniques))
              (vector-push-extend jtemp iac-vec)
              (progn
                (setf (gethash index uniques) jnext)
                (vector-push-extend jnext iac-vec)
                (vector-push-extend index type-indexj-vec)
                (vector-push-extend jnext local-typej-vec)
                (incf jnext))))
    (setf ntypes (length local-typej-vec))
    (let ((ico-vec (make-array (* ntypes ntypes) :element-type 'fixnum))
;          (ffnonbond1 (chem:get-ffnonbond-using-type-index ffnonbond-db type-index1))
          (cn1-vec (make-array (/ (* ntypes (+ ntypes 1)) 2) :element-type 'double-float)) 
          (cn2-vec (make-array (/ (* ntypes (+ ntypes 1)) 2) :element-type 'double-float))
          (type1 0)
          (type2 0)
          (cur 0)
          (curcn 0)
          (rstar 0)
          (epsilonij 0))
      (loop for i from 0 below ntypes
         for iaci = (- (aref local-typej-vec i) 1)
         do (loop for j from 0 below ntypes
               for iacj = (aref local-typej-vec j)
               for icox = (* ntypes (+ iaci iacj))
               do  (progn
                     (setf (aref ico-vec cur) icox)
                     (incf cur))))
      (loop for i from 0 below (length type-indexj-vec)
         for type1 = (aref type-indexj-vec i)
         do (loop for j from 0 below (length type-indexj-vec)
               for type2 = (aref type-indexj-vec j)
               do (when (<= i j)
                    (let ((ffnonbond1 (chem:get-ffnonbond-using-type-index ffnonbond-db type1))
                          (ffnonbond2 (chem:get-ffnonbond-using-type-index ffnonbond-db type2)))
                      (setf rstar (+ (chem:get-radius-angstroms ffnonbond1) (chem:get-radius-angstroms ffnonbond2))
                            epsilonij (sqrt (* (chem:get-epsilon-k-cal ffnonbond1) (chem:get-epsilon-k-cal ffnonbond2)))
                            (aref cn1-vec curcn) (* epsilonij (expt rstar 12.0))
                            (aref cn2-vec curcn) (* 2.0 epsilonij (expt rstar 6.0)))
                      (format t "type1 ~a type2 ~a~%" type1 type2)
                      (incf curcn)))))
      (values ntypes ico-vec iac-vec type-indexj-vec local-typej-vec cn1-vec cn2-vec))))

(defun prepare-amber-energy-nonbond (energy-function)
  (let* ((atom-table (chem:atom-table energy-function))
         (natom (chem:get-number-of-atoms atom-table))
         (atom-name-vector (make-array natom))
         (charge-vector (make-array natom :element-type 'double-float))
         (mass-vector (make-array natom :element-type 'double-float))
         (type-index-vector (make-array natom :element-type 'fixnum))
         (atomic-number-vector (make-array natom :element-type 'fixnum))
         (energy-nonbond (chem:get-nonbond-component energy-function))
         (ffnonbond-db (chem:get-ffnonbond-db energy-nonbond)))
    (format t "energy-nonbond -> ~a~%" energy-nonbond)
    (format t "ffnonbond-db -> ~a~%" ffnonbond-db)
    (loop for i from 0 below natom
       for atom-name = (chem:elt-atom-name atom-table i)
       for charge = (chem:elt-charge atom-table i)
       for mass = (chem:elt-mass atom-table i)
       for type-index = (chem:elt-type-index atom-table i)
       for atomic-number = (chem:elt-atomic-number atom-table i)
       do (setf (aref atom-name-vector i) atom-name
                (aref charge-vector i) charge
                (aref mass-vector i) mass
                (aref type-index-vector i) type-index
                (aref atomic-number-vector i) atomic-number))
    (multiple-value-bind (ntypes ico-vec iac-vec type-indexj-vec local-typej-vec cn1-vec cn2-vec)
        (generate-nonbond-parameters ffnonbond-db type-index-vector)
      (values ntypes atom-name-vector charge-vector mass-vector atomic-number-vector ico-vec iac-vec local-typej-vec cn1-vec cn2-vec))
    )  
  )
                                                                      

(defun save-amber-parm-format (aggregate topology-pathname coordinate-pathname force-field)
  (let* ((energy-function (chem:make-energy-function aggregate force-field
                                                     :use-excluded-atoms t))
         (nonbonds (chem:get-nonbond-component energy-function))
         (number-excluded-atoms (chem:number-excluded-atoms nonbonds))
         (excluded-atom-list (chem:excluded-atom-list nonbonds))
         (topology-pathname (merge-pathnames topology-pathname))
         (coordinate-pathname (merge-pathnames coordinate-pathname)))
    ;; Skip assigning MarkMainChainAtoms and MarkSideChain atoms for now
    ;; see (unitio.c:4889).  This won't mean anything for spiroligomers.
    (fortran:with-fortran-output-file (ftop topology-pathname :direction :output)
      (fortran:debug-on ftop)
      (fortran:debug "-1-")             ;
      (fortran:format 1 "%-80s")
      ;;      (fortran:write (core:strftime 81 "%%VERSION  VERSION_STAMP = V0002.000  DATE = %m/%d/%y  %H:%M:%S"))
      (fortran:write "-----insert the version and time stamp---")
      (fortran:write "%FLAG TITLE")
      (fortran:write "%FORMAT(20a4)")
;      (fortran:write (chem:get-name aggregate))
      ;; This function will be very, very long (for Common Lisp)
      ;; To avoid lots of nested scopes we will declare one large scope
      ;;   and declare all of the variables in that scope here at the top
      (let (natom
            ntypes atom-name charge mass atomic-number ico iac local-typej-vec cn1-vec cn2-vec #|nonbonds|#
            nbonh mbona ibh jbh icbh ib jb icb kbj-vec r0j-vec #|stretches|#
            ntheth mtheta ith jth kth icth it jt kt1 ict ktj-vec t0j-vec #|angles|#
            nphih mphia iph jph kph lph icph ip jp kp lp icp vj-vec inj-vec phasej-vec #|dihedrals|#
            nhparm
            NPARM  NNB    NRES
            NBONA    NTHETA NPHIA  NUMBND NUMANG NPTRA
            NATYP    NPHB   IFPERT NBPER  NGPER  NDPER
            MBPER    MGPER  MDPER  IFBOX  NMXRS  IFCAP
            NUMEXTRA NCOPY
            )
        ;; Here we need to calculate all of the values for %FLAG POINTERS
        (setf natom (chem:get-number-of-atoms (chem:atom-table energy-function))) ;total number of atoms
        (multiple-value-setq (nbonh mbona ibh jbh icbh ib jb icb kbj-vec r0j-vec)
          (prepare-amber-energy-stretch energy-function))
        (multiple-value-setq (ntheth mtheta ith jth kth icth it jt kt1 ict ktj-vec t0j-vec)
          (prepare-amber-energy-angle energy-function))
        (multiple-value-setq (nphih mphia iph jph kph lph icph ip jp kp lp icp vj-vec inj-vec phasej-vec)
          (prepare-amber-energy-dihedral energy-function))
        (multiple-value-setq (ntypes atom-name charge mass atomic-number ico iac local-typej-vec cn1-vec cn2-vec)
          (prepare-amber-energy-nonbond energy-function))
        (setf nhparm 0)
        (setf numbnd (length kbj-vec))
        (setf numang (length ktj-vec))
        (setf nptra (length vj-vec))
 
        #| dihedrals, nonbonds, others??? |#
        ;; --- Done calculating all of the values
        ;; --- Now write out all of the values
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG POINTERS")
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-2-")
        (fortran:format 10 "%8d")
        ;; NATOM
        (fortran:write natom)
        (fortran:write ntypes)
        (fortran:write nbonh)       ;number of bonds containing hydrogen
        (fortran:write mbona)       ;number of bonds not containing hydrogen
        (fortran:write ntheth)      ; number of angles containing hydrogen
        (fortran:write mtheta)      ; number of angles not containing hydrogen
        (fortran:write nphih)       ; number of dihedrals containing hydrogen
        (fortran:write mphia)       ; number of dihedrals not containing hydrogen
        (fortran:write nhparm)    ; currently not used
        ;;(fortran:write NPARM)     ; used to determine if addles created prmtop
        ;;(fortran:write NNB)       ; number of excluded atoms
        ;;(fortran:write NRES)      ; number of residues
        ;;(fortran:write NBONA)     ; MBONA + number of constraint bonds
        ;;(fortran:write NTHETA)    ; MTHETA + number of constraint angles
        ;;(fortran:write NPHIA)     ; MPHIA + number of constraint dihedrals
        (fortran:write numbnd)      ; number of unique bond types
        (fortran:write numang)      ; number of unique angle types
        (fortran:write nptra)       ; number of unique dihedral types
        ;;(fortran:write NATYP)     ; number of atom types in parameter file, see SOLTY below
        ;;(fortran:write NPHB)      ; number of distinct 10-12 hydrogen bond pair types
        ;;(fortran:write IFPERT)    ; set to 1 if perturbation info is to be read in
        ;;(fortran:write NBPER)     ; number of bonds to be perturbed
        ;;(fortran:write NGPER)     ; number of angles to be perturbed
        ;;(fortran:write NDPER)     ; number of dihedrals to be perturbed
        ;;(fortran:write MBPER)     ; number of bonds with atoms completely in perturbed group
        ;;(fortran:write MGPER)     ; number of angles with atoms completely in perturbed group
        ;;(fortran:write MDPER)     ; number of dihedrals with atoms completely in perturbed groups
        ;;(fortran:write IFBOX)     ; set to 1 if standard periodic box, 2 when truncated octahedral
        ;;(fortran:write NMXRS)     ; number of atoms in the largest residue
        ;;(fortran:write IFCAP)     ; set to 1 if the CAP option from edit was specified
        ;;(fortran:write NUMEXTRA)  ; number of extra points found in topology
        ;;(fortran:write NCOPY)     ; number of PIMD slices / number of beads
        (fortran:end-line)

        ;; Next) 
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG ATOM_NAME")
        (fortran:write "%FORMAT(20a4)")
        (fortran:debug "-3-")
        (fortran:format 20 "%4s")
        (loop for name across atom-name
           do (fortran:write name))
        (fortran:end-line)
        ;; write the atom names

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG CHARGE")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-4-")
        (fortran:format 5 "%16.8f")
        (loop for ch across charge
           do (fortran:write ch))
        (fortran:end-line)
        ;; write the atom charges

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG ATOMIC_NUMBER")
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-5-")
        (fortran:format 10 "%8d")
        (loop for number across atomic-number
           do (fortran:write number))
        (fortran:end-line)
        ;; write the atomic number of each atom

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG MASS")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-6-")
        (fortran:format 5 "%16.8f")
        (loop for ma across mass
           do (fortran:write ma))
        (fortran:end-line)
        ;; write the atom masses

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG ATOMIC_TYPE_INDEX")
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-7-")
        (fortran:format 10 "%8d")
        (loop for ia across iac
           do (fortran:write ia))
        (fortran:end-line)
        ;; write the index fot the atom types

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG NUMBER_EXCLUDED_ATOMS")
        (fortran:end-line)
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-8-")
        (fortran:format 10 "%8d")
        (fortran:end-line)
        ;; write the total number of excluded atoms for atom "i"

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG NONBONDED_PARM_INDEX")
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-9-")
        (fortran:format 10 "%8d")
        (loop for ic across ico
           do (fortran:write ic))
        (fortran:end-line)
        ;; provides the index to the nobon parameter arrays CN1, CN2 and ASOL, BSOL.

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG RESIDUE_LABEL")
        (fortran:write "%FORMAT(20A4)")
        (fortran:debug "-10-")
        (fortran:format 20 "%4s")
        (fortran:end-line)
        ;; write the name of each of the residues

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG RESIDUE_POINTER")
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-11-")
        (fortran:format 10 "%8d")
        (fortran:end-line)
        ;; write the atoms in each residue are listed for atom "1" in IPRES(i) to IPRES(i+1)-1

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG BOND_FORCE_CONSTANT")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-12-")
        (fortran:format 5 "%16.8f")
        (loop for kb across kbj-vec
           do (fortran:write kb))
        (fortran:end-line)
        
        ;; write the force constant for the bonds of each type

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%BOND_EQUIL_VALUE")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-13-")
        (fortran:format 5 "%16.8f")
        (loop for r0 across r0j-vec
           do (fortran:write r0))
        (fortran:end-line)
        ;; write the equilibrium bond length for the bonds of each type

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG ANGLE_FORCE_CONSTANT")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-14-")
        (fortran:format 5 "%16.8f")
        (loop for kt across ktj-vec
           do (fortran:write kt))
        (fortran:end-line)
        ;; write the force constant for the angles of each type

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG ANGLE_EQUIL_VALUE")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-15-")
        (fortran:format 5 "%16.8f")
        (loop for t0 across t0j-vec
           do (fortran:write t0))
        (fortran:end-line)
        ;; write the equilibrium angle for the angles of each type

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG DIHEDRAL_FORCE_CONSTANT")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-16-")
        (fortran:format 5 "%16.8f")
        (loop for v0 across vj-vec
           do (fortran:write v0))
        (fortran:end-line)
        ;; write the force constant for the dihedral of a given type

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG DIHEDRAL_PERIODICITY")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-17-")
        (fortran:format 5 "%16.8f")
        (loop for in0 across inj-vec
           do (fortran:write (float in0)))
        (fortran:end-line)
        ;; write the periodicity of the dihedral of a given type

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG DIHEDRAL_PHASE")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-18-")
        (fortran:format 5 "%16.8f")
        (loop for p0 across phasej-vec
           do (fortran:write p0))
        (fortran:end-line)
        ;; write the phase of the dihedral of a given type

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG SCEE_SCALE_FACTOR")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-19-")
        (fortran:format 5 "%16.8f")
        (fortran:end-line)
        ;; write the 1-4 electrostatic scaling constant

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG SCNB_SCALE_FACTOR")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-20-")
        (fortran:format 5 "%16.8f")
        (fortran:end-line)
        ;; write the 1-4 vdw scaling constant

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG SOLTY")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-21-")
        (fortran:format 5 "%16.8f")
        (fortran:end-line)
        ;; currently unused

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG LENNARD_JONES_ACOEF")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-22-")
        (fortran:format 5 "%16.8f")
        (loop for cn1 across cn1-vec
           do (fortran:write cn1))
        (fortran:end-line)
        ;; write the Lennard Jones r**12 terms for all possible atom type interactions

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG LENNARD_JONES_BCOEF")
        (fortran:write "%FORMAT(5E16.8)")
        (fortran:debug "-23-")
        (fortran:format 5 "%16.8f")
        (loop for cn2 across cn2-vec
           do (fortran:write cn2))
        (fortran:end-line)
       ;; write the Lennard Jones r**6 terms for all possible atom type interactions

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG BONDS_INC_HYDROGEN")
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-24-")
        (fortran:format 10 "%8d")
        (loop for i below (length ibh)
           for ibhi = (aref ibh i)
           for jbhi = (aref jbh i)
           for icbhi = (aref icbh i)
           do (fortran:write ibhi)
           do (fortran:write jbhi)
           do (fortran:write icbhi))
        (fortran:end-line)
        
        ;; write IBH, JBH, ICBH

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG BONDS_WITHOUT_HYDROGEN")
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-25-")
        (fortran:format 10 "%8d")
        (loop for i below (length ib)
           for ibi = (aref ib i)
           for jbi = (aref jb i)
           for icbi = (aref icb i)
           do (fortran:write ibi)
           do (fortran:write jbi)
           do (fortran:write icbi))
        (fortran:end-line)
        ;; write IB, JB, ICB

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG ANGLES_INC_HYDROGEN")
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-26-")
        (fortran:format 10 "%8d")
        (loop for i below (length ith)
           for ithi = (aref ith i)
           for jthi = (aref jth i)
           for kthi = (aref kth i)
           for icthi = (aref icth i)
           do (fortran:write ithi)
           do (fortran:write jthi)
           do (fortran:write kthi)
           do (fortran:write icthi))
        (fortran:end-line)
        ;; write ITH, JTH, KTH, ICTH

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG ANGLES_WITHOUT_HYDROGEN")
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-27-")
        (fortran:format 10 "%8d")
        (loop for i below (length it)
           for iti = (aref it i)
           for jti = (aref jt i)
           for kti = (aref kt1 i)
           for icti = (aref ict i)
           do (fortran:write iti)
           do (fortran:write jti)
           do (fortran:write kti)
           do (fortran:write icti))
        (fortran:end-line)
        ;; write IT, JT, KT, ICT

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG DIHEDRALS_INC_HYDROGEN")
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-28-")
        (fortran:format 10 "%8d")
        (loop for i below (length iph)
           for iphi = (aref iph i)
           for jphi = (aref jph i)
           for kphi = (aref kph i)
           for lphi = (aref lph i)
           for icphi = (aref icph i)
           do (fortran:write iphi)
           do (fortran:write jphi)
           do (fortran:write kphi)
           do (fortran:write lphi)
           do (fortran:write icphi))
         (fortran:end-line)
       ;; write IPH, JPH, KPH, LPH, ICPH

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG DIHEDRAL_WITHOUT_HYDROGEN")
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-29-")
        (fortran:format 10 "%8d")
        (loop for i below (length ip)
           for ipi = (aref ip i)
           for jpi = (aref jp i)
           for kpi = (aref kp i)
           for lpi = (aref lp i)
           for icpi = (aref icp i)
           do (fortran:write ipi)
           do (fortran:write jpi)
           do (fortran:write kpi)
           do (fortran:write lpi)
           do (fortran:write icpi))
        (fortran:end-line)
        ;; write IP, JP, KP, LP, ICP

        ;; Next
        (fortran:format 1 "%-80s")
        (fortran:write "%FLAG EXCLUDED_ATOMS_LIST")
        (fortran:write "%FORMAT(10I8)")
        (fortran:debug "-30-")
        (fortran:format 10 "%8d")
        (fortran:end-line)
        ;; write the excluded atom list

        ;; pick up at unitio.cc:4969
        ))))
