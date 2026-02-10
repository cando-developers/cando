(in-package #:rosetta.params)

(defclass rosetta-atom-attributes ()
  ((atom-name :initarg :atom-name :reader rosetta-atom-name)
   (rosetta-type :initarg :rosetta-type :reader rosetta-atom-rosetta-type)
   (mm-type :initarg :mm-type :reader rosetta-atom-mm-type)
   (charge :initarg :charge :reader rosetta-atom-charge)
   (parse-charge :initarg :parse-charge :reader rosetta-atom-parse-charge)
   (properties :initarg :properties :reader rosetta-atom-properties)))

(defmethod print-object ((obj rosetta-atom-attributes) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream ":atom-name ~s :rosetta-type ~s" (rosetta-atom-name obj) (rosetta-atom-rosetta-type obj))))

(defclass rosetta-residue ()
  ((name :initarg :name :accessor rosetta-residue-name)
   (residue :initarg :residue :accessor rosetta-residue-residue)
   (atoms-by-name :initarg :atoms-by-name :accessor rosetta-residue-atoms-by-name)
   (atom-attributes :initarg :atom-attributes :accessor rosetta-residue-atom-attributes)
   (lower-connect :initarg :lower-connect :accessor rosetta-residue-lower-connect)
   (upper-connect :initarg :upper-connect :accessor rosetta-residue-upper-connect)))

(defclass rosetta-capped-residue ()
  ((name :initarg :name :accessor rosetta-capped-residue-name)
   (molecule :initarg :molecule :accessor rosetta-capped-residue-molecule)
   (atom-attributes :initarg :atom-attributes :accessor rosetta-capped-residue-atom-attributes)))

(defmethod print-object ((obj rosetta-capped-residue) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream ":name ~a" (rosetta-capped-residue-name obj))))

(defun make-rosetta-residue (&key name residue atoms-by-name atom-attributes lower-connect upper-connect)
  (make-instance 'rosetta-residue
                 :name name
                 :residue residue
                 :atoms-by-name atoms-by-name
                 :atom-attributes atom-attributes
                 :lower-connect lower-connect
                 :upper-connect upper-connect))

(defun make-rosetta-capped-residue (&key name molecule atom-attributes)
  (make-instance 'rosetta-capped-residue
                 :name name
                 :molecule molecule
                 :atom-attributes atom-attributes))

(defun whitespacep (ch)
  (or (char= ch #\Space)
      (char= ch #\Tab)
      (char= ch #\Newline)
      (char= ch #\Return)))

(defun strip-comment (line)
  (let ((pos (position #\# line)))
    (if pos
        (subseq line 0 pos)
        line)))

(defun split-fields (line)
  (let ((len (length line))
        (pos 0)
        (fields nil))
    (loop while (< pos len) do
      (loop while (and (< pos len) (whitespacep (char line pos)))
            do (incf pos))
      (when (< pos len)
        (let ((start pos))
          (loop while (and (< pos len) (not (whitespacep (char line pos))))
                do (incf pos))
          (push (subseq line start pos) fields))))
    (nreverse fields)))

(defun parse-float (token)
  (handler-case
      (coerce (read-from-string token) 'double-float)
    (error () 0.0d0)))

(defun maybe-parse-float (token)
  (handler-case
      (let ((val (read-from-string token)))
        (when (numberp val)
          (coerce val 'double-float)))
    (error () nil)))

(defun make-keyword (name)
  (intern (string-upcase name) :keyword))

(defun make-case-sensitive-keyword (name)
  (intern name :keyword))

(defun add-rosetta-atom (residue atoms-by-name atom-attributes
                         atom-name element rosetta-type mm-type
                         charge parse-charge properties)
  (let ((atom (chem:make-atom atom-name element)))
    (chem:set-atom-type atom rosetta-type)
    (chem:set-charge atom charge)
    (chem:add-matter residue atom)
    (setf (gethash atom-name atoms-by-name) atom)
    (setf (gethash atom atom-attributes)
          (make-instance 'rosetta-atom-attributes
                         :atom-name atom-name
                         :rosetta-type rosetta-type
                         :mm-type mm-type
                         :charge charge
                         :parse-charge parse-charge
                         :properties properties))
    atom))

(defun normalize-atom-name-string (name)
  (let ((len (length name)))
    (if (and (> len 0) (digit-char-p (char name 0)))
        (concatenate 'string (subseq name 1) (string (char name 0)))
        name)))

(defun make-atom-keyword (name)
  (make-keyword (normalize-atom-name-string name)))

(defun rosetta-element-from-atom-name (name)
  (let* ((norm (normalize-atom-name-string name))
         (upper (string-upcase norm)))
    (if (and (>= (length upper) 2)
             (string= "CL" upper :end1 2 :end2 2))
        :|Cl|
        (chem:element-from-atom-name-string-case-insensitive norm))))

(defun normalize-atom-name-symbol (sym)
  (let* ((name (symbol-name sym))
         (norm (normalize-atom-name-string name)))
    (if (string= name norm)
        sym
        (intern (string-upcase norm) :keyword))))

(defun bond-key (name1 name2)
  (if (string< (symbol-name name1) (symbol-name name2))
      (list name1 name2)
      (list name2 name1)))

(defun bond-order-from-token (token)
  (let ((tok (string-upcase token)))
    (cond
      ((string= tok "1") :single-bond)
      ((string= tok "2") :double-bond)
      ((string= tok "3") :triple-bond)
      ((or (string= tok "ARO") (string= tok "AROM")) :aromatic-bond)
      (t :single-bond))))

(defparameter *rosetta-aromatic-atom-types*
  '("AROC" "NTRP" "NHIS" "NTRR" "OARO"))

(defun rosetta-aromatic-atom-type-p (rosetta-type)
  (and rosetta-type
       (member (string-upcase (symbol-name rosetta-type))
               *rosetta-aromatic-atom-types*
               :test #'string=)))

(defun rosetta-atom-type= (rosetta-type name)
  (and rosetta-type
       (string= (string-upcase (symbol-name rosetta-type)) name)))

(defun promote-aromatic-bonds (residue atom-attributes)
  (let ((to-promote nil))
    (chem:map-bonds
     nil
     (lambda (a1 a2 bond-order bond)
       (declare (ignore bond))
       (when (eq bond-order :single-bond)
         (let* ((attr1 (gethash a1 atom-attributes))
                (attr2 (gethash a2 atom-attributes))
                (t1 (and attr1 (rosetta-atom-rosetta-type attr1)))
                (t2 (and attr2 (rosetta-atom-rosetta-type attr2))))
           (when (and (rosetta-aromatic-atom-type-p t1)
                      (rosetta-aromatic-atom-type-p t2))
             (push (list a1 a2) to-promote)))))
     residue)
    (dolist (pair to-promote)
      (destructuring-bind (a1 a2) pair
        (when (eq :single-bond (chem:bond-order-to a1 a2))
          (chem:remove-bond-to a1 a2)
          (chem:bond-to a1 a2 :aromatic-bond))))))

(defun promote-carbonyl-bonds (residue atom-attributes)
  (let ((to-promote nil))
    (chem:map-bonds
     nil
     (lambda (a1 a2 bond-order bond)
       (declare (ignore bond))
       (when (eq bond-order :single-bond)
         (let* ((attr1 (gethash a1 atom-attributes))
                (attr2 (gethash a2 atom-attributes))
                (t1 (and attr1 (rosetta-atom-rosetta-type attr1)))
                (t2 (and attr2 (rosetta-atom-rosetta-type attr2))))
           (when (or (and (rosetta-atom-type= t1 "COBB")
                          (rosetta-atom-type= t2 "OCBB"))
                     (and (rosetta-atom-type= t1 "OCBB")
                          (rosetta-atom-type= t2 "COBB")))
             (push (list a1 a2) to-promote)))))
     residue)
    (dolist (pair to-promote)
      (destructuring-bind (a1 a2) pair
        (when (eq :single-bond (chem:bond-order-to a1 a2))
          (chem:remove-bond-to a1 a2)
          (chem:bond-to a1 a2 :double-bond))))))

(defun default-rosetta-home ()
  (or (uiop:getenv-absolute-directory "ROSETTA_HOME")
      (probe-file "/opt/rosetta/")))

(defun default-fa-standard-lcaa-directory ()
  (let ((home (default-rosetta-home)))
    (when home
      (merge-pathnames
       "database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/"
       home))))

(defun default-ace-params-path ()
  (let ((home (default-rosetta-home)))
    (when home
      (merge-pathnames
       "database/chemical/residue_type_sets/fa_standard/residue_types/terminal/ACE.params"
       home))))

(defun default-nme-params-path ()
  (let ((home (default-rosetta-home)))
    (when home
      (merge-pathnames
       "database/chemical/residue_type_sets/fa_standard/residue_types/terminal/NME.params"
       home))))

(defun parse-params-file (pathname)
  (let ((atoms-by-name (make-hash-table :test #'eq))
        (atom-attributes (make-hash-table :test #'eq))
        (bond-lines nil)
        (bond-type-map (make-hash-table :test #'equal))
        (virtual-shadow (make-hash-table :test #'eq))
        (pending-virtual-atoms nil)
        (pending-aliases nil)
        (res-name nil)
        (residue nil)
        (lower-connect nil)
        (upper-connect nil))
    (with-open-file (fin pathname)
      (loop for raw-line = (read-line fin nil nil)
            while raw-line do
              (let* ((line (string-trim '(#\Space #\Tab #\Newline #\Return)
                                        (strip-comment raw-line))))
                (when (> (length line) 0)
                  (let ((fields (split-fields line)))
                    (when fields
                      (let ((head (string-upcase (first fields))))
                        (cond
                          ((string= head "NAME")
                           (setf res-name (make-keyword (second fields)))
                           (setf residue (chem:make-residue res-name)))
                          ((string= head "ATOM")
                           (unless residue
                             (error "ATOM line before NAME in ~a" pathname))
                           (let* ((raw-name (second fields))
                                  (norm-name (normalize-atom-name-string raw-name))
                                  (atom-name (make-atom-keyword raw-name))
                                  (rosetta-type (make-case-sensitive-keyword (third fields)))
                                  (mm-type (make-keyword (fourth fields)))
                                  (charge-token (fifth fields))
                                  (maybe-parse (sixth fields))
                                  (parse-charge (maybe-parse-float maybe-parse))
                                  (properties-start (if parse-charge 6 5))
                                  (properties (mapcar #'make-keyword (nthcdr properties-start fields)))
                                  (element (rosetta-element-from-atom-name norm-name))
                                  (charge (parse-float charge-token))
                                  (virt-p (or (eq rosetta-type :VIRT)
                                              (eq mm-type :VIRT))))
                             (if virt-p
                                 (push (list :atom-name atom-name
                                             :element element
                                             :rosetta-type rosetta-type
                                             :mm-type mm-type
                                             :charge charge
                                             :parse-charge parse-charge
                                             :properties properties)
                                       pending-virtual-atoms)
                                 (add-rosetta-atom residue atoms-by-name atom-attributes
                                                   atom-name element rosetta-type mm-type
                                                   charge parse-charge properties))))
                          ((string= head "ATOM_ALIAS")
                           (let* ((real-name (make-atom-keyword (second fields)))
                                  (alias-name (make-atom-keyword (third fields)))
                                  (atom (gethash real-name atoms-by-name)))
                             (if atom
                                 (let ((existing (gethash alias-name atoms-by-name)))
                                   (unless (and existing (not (eq existing atom)))
                                     (setf (gethash alias-name atoms-by-name) atom)))
                                 (push (list real-name alias-name) pending-aliases))))
                          ((string= head "LOWER_CONNECT")
                           (let ((val (string-upcase (second fields))))
                             (unless (string= val "NONE")
                               (setf lower-connect (make-atom-keyword val)))))
                          ((string= head "UPPER_CONNECT")
                           (let ((val (string-upcase (second fields))))
                             (unless (string= val "NONE")
                               (setf upper-connect (make-atom-keyword val)))))
                          ((string= head "VIRTUAL_SHADOW")
                           (let ((virt (make-atom-keyword (second fields)))
                                 (shadow (make-atom-keyword (third fields))))
                             (setf (gethash virt virtual-shadow) shadow)))
                          ((string= head "BOND")
                           (push (list (make-atom-keyword (second fields))
                                       (make-atom-keyword (third fields)))
                                 bond-lines))
                          ((string= head "BOND_TYPE")
                           (let* ((a1 (make-atom-keyword (second fields)))
                                  (a2 (make-atom-keyword (third fields)))
                                  (order (bond-order-from-token (fourth fields))))
                             (setf (gethash (bond-key a1 a2) bond-type-map) order)
                             (push (list a1 a2) bond-lines)))
                          (t nil)))))))))
    (dolist (data pending-virtual-atoms)
      (let ((atom-name (getf data :atom-name)))
        (unless (gethash atom-name virtual-shadow)
          (add-rosetta-atom residue atoms-by-name atom-attributes
                            atom-name
                            (getf data :element)
                            (getf data :rosetta-type)
                            (getf data :mm-type)
                            (getf data :charge)
                            (getf data :parse-charge)
                            (getf data :properties)))))
    (dolist (alias pending-aliases)
      (destructuring-bind (real-name alias-name) alias
        (let ((atom (gethash real-name atoms-by-name)))
          (when atom
            (let ((existing (gethash alias-name atoms-by-name)))
              (unless (and existing (not (eq existing atom)))
                (setf (gethash alias-name atoms-by-name) atom)))))))
    (let ((bond-done (make-hash-table :test #'equal)))
      (dolist (bond bond-lines)
        (destructuring-bind (a1 a2) bond
          (let* ((ra1 (gethash a1 virtual-shadow a1))
                 (ra2 (gethash a2 virtual-shadow a2))
                 (key (bond-key ra1 ra2))
                 (raw-key (bond-key a1 a2)))
            (unless (gethash key bond-done)
              (setf (gethash key bond-done) t)
              (let ((atom1 (gethash ra1 atoms-by-name))
                    (atom2 (gethash ra2 atoms-by-name)))
                (unless atom1
                  (error "Unknown atom ~a in bond line of ~a" ra1 pathname))
                (unless atom2
                  (error "Unknown atom ~a in bond line of ~a" ra2 pathname))
                (chem:bond-to atom1 atom2
                                  (or (gethash key bond-type-map)
                                      (gethash raw-key bond-type-map)
                                      :single-bond))))))))
    (promote-aromatic-bonds residue atom-attributes)
    (promote-carbonyl-bonds residue atom-attributes)
    (when (gethash lower-connect virtual-shadow)
      (setf lower-connect (gethash lower-connect virtual-shadow)))
    (when (gethash upper-connect virtual-shadow)
      (setf upper-connect (gethash upper-connect virtual-shadow)))
    (make-rosetta-residue :name res-name
                          :residue residue
                          :atoms-by-name atoms-by-name
                          :atom-attributes atom-attributes
                          :lower-connect lower-connect
                          :upper-connect upper-connect)))

(defun rosetta-residue-atom (rosetta-residue atom-name)
  (gethash (cond
             ((symbolp atom-name)
              (normalize-atom-name-symbol atom-name))
             (t
              (make-atom-keyword atom-name)))
           (rosetta-residue-atoms-by-name rosetta-residue)))

(defun merge-atom-attributes (&rest tables)
  (let ((merged (make-hash-table :test #'eq)))
    (dolist (table tables)
      (when table
        (maphash (lambda (atom attrs)
                   (setf (gethash atom merged) attrs))
                 table)))
    merged))

(defun build-capped-residue (params-path
                             &key ace-params-path
                               nme-params-path)
  (let* ((ace-path (or ace-params-path (default-ace-params-path)))
         (nme-path (or nme-params-path (default-nme-params-path))))
    (unless ace-path
      (error "Could not locate ACE.params - provide :ace-params-path"))
    (unless nme-path
      (error "Could not locate NME.params - provide :nme-params-path"))
    (format t "Loading ~s~%" params-path)
    (let* ((base (parse-params-file params-path))
           (ace (parse-params-file ace-path))
           (nme (parse-params-file nme-path))
           (mol (chem:make-molecule (rosetta-residue-name base))))
      (chem:add-matter mol (rosetta-residue-residue ace))
      (chem:add-matter mol (rosetta-residue-residue base))
      (chem:add-matter mol (rosetta-residue-residue nme))
      (let ((ace-upper (rosetta-residue-upper-connect ace))
            (base-lower (rosetta-residue-lower-connect base))
            (base-upper (rosetta-residue-upper-connect base))
            (nme-lower (rosetta-residue-lower-connect nme)))
        (unless (and ace-upper base-lower base-upper nme-lower)
          (error "Missing connect atoms for ACE/BASE/NME for ~a" params-path))
        (chem:bond-to (rosetta-residue-atom ace ace-upper)
                      (rosetta-residue-atom base base-lower)
                      :single-bond)
        (chem:bond-to (rosetta-residue-atom base base-upper)
                      (rosetta-residue-atom nme nme-lower)
                      :single-bond))
      (make-rosetta-capped-residue
       :name (rosetta-residue-name base)
       :molecule mol
       :atom-attributes (merge-atom-attributes (rosetta-residue-atom-attributes ace)
                                               (rosetta-residue-atom-attributes base)
                                               (rosetta-residue-atom-attributes nme))))))

(defun build-capped-residues-in-directories (&key
                                               (directories (list
                                                             #+(or)(probe-file #P"rosetta:database;chemical;residue_type_sets;fa_standard;residue_types;l-caa;")
                                                             (list (probe-file #P"rosetta:database;chemical;residue_type_sets;fa_standard;residue_types;peptoids;")
                                                                   (list "011" "308" "309" "403" "506")))
                                                             )
                                               (ace-params-path (make-pathname :name "ACE" :type "params" :defaults (probe-file #P"rosetta:database;chemical;residue_type_sets;fa_standard;residue_types;terminal;")))
                                               (nme-params-path (make-pathname :name "NME" :type "params" :defaults (probe-file #P"rosetta:database;chemical;residue_type_sets;fa_standard;residue_types;terminal;")))
                                               (sort t))
  (loop for maybe-dir in directories
        for dir = (if (consp maybe-dir)
                      (first maybe-dir)
                      maybe-dir)
        for skip = (if (consp maybe-dir)
                       (second maybe-dir)
                       nil)
        append (let* ((pattern (merge-pathnames "*.params" (uiop:ensure-directory-pathname dir)))
                      (files (directory pattern)))
                 (when sort
                   (setf files (sort files #'string< :key #'namestring)))
                 (loop for file in files
                       unless (member (pathname-name file) skip :test 'equal)
                         collect (build-capped-residue file
                                                       :ace-params-path ace-params-path
                                                       :nme-params-path nme-params-path)))))
