(in-package :kin.fragment)

(defclass fragment ()
  ((name :initarg :name :accessor name)
   (atoms :initarg :atoms :accessor atoms)
   (stereoisomer-info :initarg :stereoisomer-info :accessor stereoisomer-info)
   (in-plug :initarg :in-plug :accessor in-plug)
   (out-plugs :initarg :out-plugs :accessor out-plugs)))

(defmethod print-object ((obj fragment) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (name obj))))

(defclass atom ()
  ((name :initarg :name :accessor name)
   (bonds :initarg :bonds :accessor bonds)))

(defmethod print-object ((obj atom) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (name obj))))

(defclass bond ()
  ((target :initarg :target :accessor target)))

(defmethod print-object ((obj bond) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (target obj))))

(defclass single-bond (bond) ())
(defclass double-bond (bond) ())
(defclass ring-single-bond (bond) ())
(defclass ring-double-bond (bond) ())

(defclass stereoisomers ()
  ((name-template :initarg :name-template :accessor name-template)
   (stereocenters :initform nil :initarg :stereocenters :accessor stereocenters)))

(defclass plug ()
  ((name :initform nil :initarg :name :accessor name)
   (bond0 :initarg :bond0 :accessor bond0)
   (bond1 :initform nil :initarg :bond1 :accessor bond1)))

(defmethod print-object ((obj plug) stream)
  (print-unreadable-object (obj stream :type t)
    (if (null (bond1 obj))
        (format stream "~a :bond0 ~a" (name obj) (bond0 obj) )
        (format stream "~a :bond0 ~a :bond1 ~a" (name obj) (bond0 obj) (bond1 obj)))))


(defclass in-plug (plug)
  ())

(defclass origin-plug (in-plug)
  ())

(defclass out-plug (plug)
  ())


(defun do-fragment (name &rest parts)
  (loop for part in parts
        with rev-atoms = nil
        with stereoisomers = nil
        with in-plug = nil
        with out-plugs = nil
        do (cond
             ((typep part 'atom)
              (push part rev-atoms))
             ((typep part 'stereoisomers)
              (setf stereoisomers part))
             ((typep part 'in-plug)
              (setf in-plug part))
             ((typep part 'out-plug)
              (push part out-plugs))
             (t (warn "do-fragment - handle ~a~%" part)))
        finally (return (make-instance 'fragment
                                       :name name
                                       :atoms (reverse rev-atoms)
                                       :stereoisomer-info (or stereoisomers (make-instance 'stereoisomers :name-template (string name)))
                                       :in-plug in-plug
                                       :out-plugs (reverse out-plugs)))))

(defparameter *fragments* nil)
(defparameter *cap-map* nil)

(defmacro fragment (name &rest args)
  (let ((frag-gs (gensym)))
    `(let ((,frag-gs (do-fragment ,name ,@args)))
       (push ,frag-gs *fragments*)
       )))


(defun atom (name &rest parts)
  (loop for part in parts
        with rev-bonds = nil
        do (cond
             ((typep part 'single-bond)
              (push part rev-bonds))
             ((typep part 'double-bond)
              (push part rev-bonds))
             ((typep part 'ring-single-bond)
              (push part rev-bonds))
             ((typep part 'ring-double-bond)
              (push part rev-bonds))
             (t (warn "Bad term ~a" part)))
        finally (return (make-instance 'atom :name name :bonds (reverse rev-bonds)))))

(defun single (name)
  (make-instance 'single-bond :target name))

(defun double (name)
  (make-instance 'double-bond :target name))

(defun ring-single (name)
  (make-instance 'ring-single-bond :target name))

(defun ring-double (name)
  (make-instance 'ring-double-bond :target name))

(defun stereoisomers (template &rest names)
  (make-instance 'stereoisomers :name-template template :stereocenters names))

(defun origin (root)
  (make-instance 'origin-plug
                 :bond0 root))

(defun in (name &key bond0 bond1)
  (make-instance 'in-plug
                 :name name
                 :bond0 bond0
                 :bond1 bond1))

(defun out (name &key bond0 bond1 cap)
  (make-instance 'out-plug
                 :name name
                 :bond0 bond0
                 :bond1 bond1))

(defun cap-map (&rest map)
  (setf *cap-map* (alexandria:plist-hash-table map)))

(defun stable-unique-hydrogen-names (valence atom-suffix-name &optional unique-atom-names)
  (let ((hnames (loop for id from 1 to valence
                      collect (intern (format nil "H~a~a" atom-suffix-name id) :keyword))))
    (when unique-atom-names
      (loop for name in hnames
            when (gethash name unique-atom-names)
              do (let ((names (loop for id below valence
                                    for name = (intern (format nil "H~a" (+ id (hash-table-count unique-atom-names))) :keyword)
                                    do (setf (gethash name unique-atom-names) t))))
                   (return names))))
    (loop for name in hnames
          do (setf (gethash name unique-atom-names) t))
    hnames))
       

(defun build-residue-from-fragment-fill-valences (fragment)
  (format t "---- Fragment: ~a~%" (name fragment))
  (format t "   (in-plug fragment) -> ~a~%" (in-plug fragment))
  (format t "   (out-plugs fragment) -> ~a~%" (out-plugs fragment))
  (let ((atom-name-to-atom (make-hash-table))
        (residue (chem:make-residue (name fragment)))
        (plug-atoms (make-hash-table))
        (back-spanning-names (make-hash-table))
        (root-atom-name (name (first (atoms fragment)))))
    (when (and (in-plug fragment) (not (typep (in-plug fragment) 'origin-plug)))
      (incf (gethash (bond0 (in-plug fragment)) plug-atoms 0))
      (when (bond1 (in-plug fragment))
        (incf (gethash (bond1 (in-plug fragment)) plug-atoms 0))))
    (loop for out-plug in (out-plugs fragment)
          do (incf (gethash (bond0 out-plug) plug-atoms 0))
          do (when (bond1 out-plug)
               (incf (gethash (bond1 out-plug) plug-atoms 0))))
    (loop for atom in (atoms fragment)
          for atom-name = (name atom)
          for element = (chem:element-from-atom-name-string (string atom-name))
          for chem-atom = (chem:make-atom atom-name element)
          do (chem:add-matter residue chem-atom)
          do (setf (gethash atom-name atom-name-to-atom) chem-atom))
    (loop for atom in (atoms fragment)
          for atom-name = (name atom)
          for chem-atom = (gethash atom-name atom-name-to-atom)
          for bonds = (bonds atom)
          do (loop for bond in bonds
                   for bonded-atom-name = (target bond)
                   for bonded-atom = (gethash bonded-atom-name atom-name-to-atom)
                   for bond-order = (typecase bond
                                      (single-bond :single-bond)
                                      (double-bond :double-bond)
                                      (ring-single-bond :single-bond)
                                      (ring-double-bond :double-bond)
                                      (otherwise (error "Illegal bond ~a" bond)))
                   do (setf (gethash atom-name back-spanning-names) bonded-atom-name)
                   do (chem:bond-to chem-atom bonded-atom bond-order)))
    (let ((unique-atom-names (make-hash-table)))
      (cando:do-atoms (atom residue)
        (setf (gethash (string (chem:get-name atom)) unique-atom-names) t))
      (cando:do-atoms (atom residue)
        (let* ((atom-valence (chem:number-of-open-valence atom))
               (in-plug-count (gethash (chem:get-name atom) plug-atoms 0))
               (valence (- atom-valence in-plug-count))
               (atom-name (chem:get-name atom))
               (atom-name-string (string atom-name))
               (atom-element-string (string (chem:get-element atom)))
               (atom-suffix-name (if (< (length atom-element-string) (length atom-name-string))
                                     (subseq atom-name-string (length atom-element-string))
                                     ""))
               (hydrogen-atom-names (stable-unique-hydrogen-names valence atom-suffix-name unique-atom-names)))
          (format t " For atom ~a atom-valence: ~a  in-plug-count: ~a valence: ~a ~%"
                  atom-name atom-valence in-plug-count valence)
          (format t "    hydrogens: ~a plug-atoms: ~a~%"
                  hydrogen-atom-names (alexandria:hash-table-alist plug-atoms))
          (format t "    bonds: ~a~%"  (chem:bonds-as-list atom))
          (loop for hydrogen-atom-name in hydrogen-atom-names
                for hydrogen = (chem:make-atom hydrogen-atom-name :h)
                do (setf (gethash hydrogen-atom-name back-spanning-names) atom-name)
                do (chem:add-matter residue hydrogen)
                do (chem:bond-to atom hydrogen :single-bond))))
      (let* ((root-atom (chem:atom-with-name residue root-atom-name))
             (all-atoms (chem:map-atoms 'list #'identity residue))
             (back-spanning-atoms (let ((bsa (make-hash-table)))
                                    (maphash (lambda (key value)
                                               (let ((atm (chem:atom-with-name residue key))
                                                     (back-span-atm (chem:atom-with-name residue value)))
                                                 (setf (gethash atm bsa) back-span-atm)))
                                             back-spanning-names)
                                    bsa)))
        (values residue root-atom all-atoms back-spanning-atoms)))))

(defun extract-prepare-topology (fragment)
  (multiple-value-bind (residue root-atom all-atoms atom-back-span)
      (build-residue-from-fragment-fill-valences fragment)
    (let* ((name (chem:get-name residue))
           (constitution-atoms (chem:make-constitution-atoms-from-residue residue))
           (stereo-information (stereoisomer-info fragment))
           (chiral-atoms (loop for stereocenter-name in (stereocenters stereo-information)
                               for chiral-atom = (chem:atom-with-name residue stereocenter-name)
                               collect chiral-atom))
           (number-of-chiral-atoms (length chiral-atoms))
           (number-of-stereoisomers (expt 2 number-of-chiral-atoms)))
      (let* ((stereoisomers (loop for stereoisomer-index below number-of-stereoisomers
                                  for stereoisomer = (kin::build-stereoisomer name chiral-atoms stereoisomer-index)
                                  collect stereoisomer))
             (stereo-information (chem:make-stereoinformation stereoisomers nil)))
        (let* ((fragment-in-plug (in-plug fragment))
               (in-plug (when fragment-in-plug
                          (let* ((bond0 (bond0 fragment-in-plug))
                                 (bond1 (bond1 fragment-in-plug))
                                 (in-plug (if (typep fragment-in-plug 'origin-plug)
                                              (if (null bond1)
                                                  (chem:make-origin-plug nil bond0)
                                                  (error "Origin plug can't have more than one atom"))
                                              (if (null bond1)
                                                  (chem:make-in-plug (name fragment-in-plug)
                                                                     nil
                                                                     bond0
                                                                     :single-bond)
                                                  (chem:make-in-plug (name fragment-in-plug)
                                                                     nil
                                                                     bond0
                                                                     :single-bond
                                                                     bond1
                                                                     :single-bond)))))
                            in-plug)))
               (out-plugs (loop for fragment-out-plug in (out-plugs fragment)
                                for name = (name fragment-out-plug)
                                for bond0 = (bond0 fragment-out-plug)
                                for bond1 = (bond1 fragment-out-plug)
                                for out-plug = (if (null bond1)
                                                   (chem:make-out-plug name nil nil
                                                                       bond0
                                                                       :single-bond)
                                                   (chem:make-out-plug name nil nil
                                                                       bond0
                                                                       :single-bond
                                                                       bond1
                                                                       :single-bond
                                                                       ))
                                collect out-plug))
               (plugs (if in-plug
                          (cons in-plug out-plugs)
                          out-plugs))
               (constitution (chem:make-constitution name
                                                     "no-comment"
                                                     constitution-atoms
                                                     stereo-information
                                                     plugs
                                                     nil)) 
               (joint-template-tree (kin::build-internal-coordinate-joint-template-tree-spanning-tree name residue plugs constitution root-atom all-atoms atom-back-span)) )
          (make-instance 'kin::prepare-topology
                         :name name
                         :constitution constitution
                         :constitution-atoms constitution-atoms
                         :plugs plugs
                         :stereo-information stereo-information
                         :joint-template joint-template-tree
                         :residue residue))))))
