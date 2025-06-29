(in-package :topology)

#|
Create topology instances using a graph described using an s-expression.

(topology:define-topology 'gly '((C C :plugs (+dkp.1)) (= O)
                                 - CA (- HA1) (- HA2)
                                 - (N N :plugs (+dkp.0))
                                 - CN (- HN1) (- HN2) - HN3))

(topology:define-topology 'pro4 '((CAC C :plugs (-dkp.0)) (= OA)
                                  - (CA :C :stereochemistry-type :both) (- HA)
                                  - CB (- HB1) (- HB2)
                                  - (CG :C :stereochemistry-type :both)
                                  (- (NG N :plugs (+side.0 +dkp.0)))
                                  (- (CGC C :plugs (+dkp.1)) = OG)
                                  - CD (- HD1) (- HD2) -
                                  (NA N :plugs (-dkp.1)) - CA ))

(topology:define-topology 'pro '((C C :plugs (-dkp.0)) (= O)
                                 - (CA C :stereochemistry-type :both) (- HA)
                                 - CB (- HB1) (- HB2)
                                 - CG (- HG1) (- HG2)
                                 - CD (- HD1) (- HD2)
                                 - (N N :plugs (-dkp.1)) - CA ))

(topology:define-topology 'bnz '((CM C :plugs (-side.0)) (- HM1) (- HM2)
                                 - C1
                                 = C2 (- H2)
                                 - C3 (- H3)
                                 = C4 (- H4)
                                 - C5 (- H5)
                                 = C6 (- H6)
                                 - C1))
|#


(defclass node ()
  ((name :initarg :name :accessor name)
   ))

(defmethod print-object ((obj node) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (name obj))))

(defclass atom-node (node)
  ((element :initarg :element :accessor element)
   (property-list :initform nil :initarg :property-list :accessor property-list)
   (children :initform nil :accessor children)
   ))

(defclass ordered-atom-node (atom-node)
  ((constitution-atom-index :initarg :constitution-atom-index :accessor constitution-atom-index)))

(defmethod print-object ((obj ordered-atom-node) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a" (name obj) (constitution-atom-index obj))))

(defun order-children (children)
  "Order the children so that hydrogen atoms follow heavy atoms.
Keep the heavy atoms and hydrogens in the order they were found."
  (when (= (length children) 1)
    (return-from order-children children))
  (let ((build-order-count (count-if (lambda (xxx)
                                       (getf (property-list xxx) :build-order))
                                     children)))
    (when (> build-order-count 1)
      (if (= build-order-count (length children))
          (error "Add a sort of the children by build order")
          (error "One child in ~s has :build-order and so all children must have it" children))))
  (multiple-value-bind (heavy-atom-nodes hydrogen-atom-nodes)
      (loop for child in children
            if (chem:element-is-hydrogen (element child))
              collect child into hydrogen-nodes
            else
              collect child into heavy-nodes
            finally (return (values heavy-nodes hydrogen-nodes)))
    (append heavy-atom-nodes hydrogen-atom-nodes)))

(defun add-children (node children)
  (check-type children list)
  (when children
    (let* ((all-children (append (children node) children))
           (ordered-children (order-children all-children)))
      (setf (children node) ordered-children))))

(defclass edge ()
  ((from-node :initarg :from-node :accessor from-node)
   (to-node :initarg :to-node :accessor to-node)
   (edge-type :initarg :edge-type :accessor edge-type)))

(defmethod print-object ((obj edge) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a ~a" (name (from-node obj)) (edge-type obj) (name (to-node obj)))))

(defclass graph ()
  ((name :initarg :name :accessor name)
   (root-node :initarg :root-node :accessor root-node)
   (in-plug :initform nil :initarg :in-plug :accessor in-plug)
   (edges :initform nil :initarg :edges :accessor edges)
   (nodes :initform (make-hash-table) :initarg :nodes :accessor nodes)
   (plugs :initform (make-hash-table) :accessor plugs)
   ))

;; A plug node
(defun parse-plug-name (keyword-atm atom-name plugs graph)
  "Parse plug names like -dkp.0, +dkp.1, +side, -side"
  (let* ((sname (symbol-name keyword-atm))
         (direction-char (schar sname 0)))
    (when (member direction-char '(#\- #\+))
      (let* ((dot-pos (position #\. sname))
             (plug-name-str (if dot-pos
                                (subseq sname 0 dot-pos)
                                sname))
             (plug-name (intern plug-name-str :keyword))
             (plug-index (if dot-pos
                             (parse-integer sname :start (1+ dot-pos) :junk-allowed t)
                             0))
             (plug (gethash plug-name plugs)))
        (unless plug
          (cond
            ((char= direction-char #\-)
             (let ((in-plug (make-instance 'in-plug :name plug-name)))
               (setf plug in-plug
                     (in-plug graph) in-plug)))
            ((char= direction-char #\+)
             (setf plug (make-instance 'out-plug :name plug-name)))
            (t (error "Illegal direction-char ~s" direction-char)))
          (setf (gethash plug-name plugs) plug))
        (let ((plug-bond (make-instance 'plug-bond
                                        :atom-name atom-name
                                        :bond-order :single-bond)))
          (setf (aref (plug-bonds plug) plug-index) plug-bond))))))

(defun node-from-symbol-entry (keyword-atm graph &optional properties element)
  (let* ((name-string (string keyword-atm))
         (dot-position (position #\. name-string))
         (name-only (if dot-position
                        (subseq name-string 0 dot-position)
                        name-string))
         (keyword-element (if element
                              element
                              (chem:element-from-atom-name-string name-only))))
    (cond
      ;; Names like:
      ;;   "CA.HA" generates CA -> HA
      ;;   "CB.2HB" generates CB -> (HB1 HB2)
      ;;   "CM.3HM" generates CM -> (HM1 HM2 HM3)
      ((and dot-position (< dot-position (length name-string)))
       (let ((hydrogen-specifier (subseq name-string (1+ dot-position) (length name-string))))
         (when (string= hydrogen-specifier "HZ")
           (break "Check other stuff"))
         (multiple-value-bind (num-hydrogens hprefix-start)
             (cond
               ((member (elt hydrogen-specifier 0) '(#\2 #\3))
                (values (- (char-code (elt hydrogen-specifier 0)) (char-code #\0)) 1))
               ((or (char= (elt hydrogen-specifier 0) #\H)
                    (char= (elt hydrogen-specifier 0) #\h))
                (values 1 0))
               (t (error "Cannot translate ~s to atom" name-string)))
           (let* ((hprefix (subseq hydrogen-specifier (if (= 1 num-hydrogens) 0 1) (length hydrogen-specifier)))
                  (node (let ((node (make-instance 'atom-node
                                                   :name (intern name-only :keyword)
                                                   :element keyword-element
                                                   :property-list properties
                                                   )))
                          (setf (gethash (name node) (nodes graph)) node)
                          node))
                  (hydrogens (loop for hindex from 1 to num-hydrogens
                                   for hname = (if (= num-hydrogens 1)
                                                   (string-upcase hprefix)
                                                   (format nil "~a~a" (string-upcase hprefix) hindex))
                                   for hydrogen = (make-instance 'atom-node
                                                                 :name (intern hname :keyword)
                                                                 :element :H
                                                                 )
                                   do (setf (gethash (name hydrogen) (nodes graph)) hydrogen)
                                   do (push (make-instance 'edge
                                                           :from-node node
                                                           :to-node hydrogen
                                                           :edge-type :-)
                                            (edges graph))
                                   collect hydrogen)))
             (add-children node hydrogens)
             node))))
      (dot-position
       (error "Can not interpret ~s" name-string))
      (t (when (bond-symbol keyword-atm)
           (error "Don't create atoms with name ~a" keyword-atm))
         (let ((node (make-instance 'atom-node
                                    :name keyword-atm
                                    :element keyword-element
                                    :property-list properties
                                    )))
           (setf (gethash (name node) (nodes graph)) node)
           node)))))

(defun ensure-keyword (name)
  (intern (symbol-name name) :keyword))

(defun parse-atom (atm graph)
  (cond
    ((symbolp atm)
     (let* ((keyword-atm (ensure-keyword atm))
            (seen-node (gethash keyword-atm (nodes graph))))
       (if seen-node
           (values seen-node t)
           (node-from-symbol-entry keyword-atm graph))))
    ((consp atm)
     (let* ((keyword-atm (ensure-keyword (car atm)))
            (element (if (cadr atm)
                         (ensure-keyword (cadr atm))
                         (chem:element-from-atom-name-string (symbol-name atm))))
            (property-list (let ((pl (cddr atm)))
                             (unless (and (evenp (length pl))
                                          (loop named keywords
                                                for (a b) on pl by #'cddr
                                                unless (keywordp a)
                                                  do (return-from keywords nil)
                                                finally (return-from keywords t))
                                          t)
                               (error "Atom sexp ~s should end with a property list" atm))
                             pl)))
       (when (bond-symbol keyword-atm)
         (error "Don't create atoms with name ~a" keyword-atm))
       (node-from-symbol-entry keyword-atm graph property-list element)))
    (t (error "Illegal atm ~s" atm))))

(defun bond-symbol (name)
  (when (symbolp name)
    (let ((kw-name (intern (symbol-name name) :keyword)))
      (when (member kw-name '(:- := :# :~))
        kw-name))))

(defun atom-sexp (sexp)
  (or (and (consp sexp) (null (bond-symbol (car sexp))))
      (and (symbolp sexp) (null (bond-symbol sexp)))))

(defun parse-atom-or-bond (sexp graph prev-node)
  (cond
    ((null sexp) nil)
    ((and (null prev-node) (atom-sexp (car sexp)))
     (let ((node (parse-atom (car sexp) graph)))
       (setf (root-node graph) node)
       (setf prev-node node))
     (parse-atom-or-bond (cdr sexp) graph prev-node))
    ((and (symbolp (car sexp)) (bond-symbol (car sexp)))
     (let* ((bond-symbol (bond-symbol (car sexp))))
       (multiple-value-bind (node seen)
           (if (cadr sexp)
               (parse-atom (cadr sexp) graph)
               (error "Missing atom at end of ~s" sexp))
         (if prev-node
             (progn
               (unless seen
                 (add-children prev-node (list node)))
               (push (make-instance 'edge
                                    :from-node prev-node
                                    :to-node node
                                    :edge-type bond-symbol)
                     (edges graph)))
             (error "Missing prev-node"))
         (parse-atom-or-bond (cddr sexp) graph node))))
    ((and (car sexp) prev-node (atom-sexp (caar sexp)))
     (let* ((bond-symbol :-))
       (multiple-value-bind (node seen)
           (if (car sexp)
               (parse-atom (car sexp) graph)
               (error "Missing atom at end of ~s" sexp))
         (if prev-node
             (progn
               (unless seen
                 (add-children prev-node (list node)))
               (push (make-instance 'edge
                                    :from-node prev-node
                                    :to-node node
                                    :edge-type bond-symbol)
                     (edges graph))))
         (parse-atom-or-bond (cdr sexp) graph node)
         node)))
    ((and (consp (car sexp)) (bond-symbol (car (car sexp))))
     (if prev-node
         (parse-atom-or-bond (car sexp) graph prev-node)
         (error "Missing prev-node"))
     (parse-atom-or-bond (cdr sexp) graph prev-node))
    (t (error "Atom node ~s allowed in first position only when no prev-node. prev-node = ~s" sexp prev-node))))

(defun interpret (name sexp)
  (let ((graph (make-instance 'graph :name name)))
    (parse-atom-or-bond sexp graph nil)
    graph))

(defun parse-atom-graph (sexp)
  "Parse an atom-graph like those used by define-topology"
  (interpret nil sexp))


(defun bond-order-from-edge-type (edge-type)
  (case edge-type
    (:- :single-bond)
    (:= :double-bond)
    (:# :triple-bond)
    (otherwise (error "Handle edge-type ~a" edge-type))))


(defun build-stereoisomer (name variable-chiral-atoms stereoisomer-index fixed-chiral-info stereoisomer-atoms)
  "Build a stereoisomer. Assemble a name for the stereiosomer using the name and
the stereoisomer-index.   The stereoisomer-index is like a bit-vector that indicates
the stereochemistry of this stereoisomer (0 is :S and 1 is :R).
The name is constructed by appending \"{chiral-atom-name/(S|R),...}\".
So if name is \"ALA\" and stereoisomer-index is 1 the name becomes ALA{CA/S}."
  (if (or variable-chiral-atoms fixed-chiral-info)
      (let* ((configurations (loop for chiral-index below (length variable-chiral-atoms)
                                   for config = (if (logbitp chiral-index stereoisomer-index) :S :R)
                                   collect config))
             (stereo-configurations (append
                                     (mapcar (lambda (atom config)
                                               (when (eq config :undefined-configuration)
                                                 (error "Don't allow :undefined-configuration for ~a" atom))
                                               (make-instance 'stereoconfiguration
                                                              :atom-name (progn
                                                                           (unless (symbolp atom)
                                                                             (error "atom ~s must be an atom-name" atom))
                                                                           atom)
                                                              :configuration config))
                                             variable-chiral-atoms configurations)
                                     (mapcar (lambda (info)
                                               (let* ((constitution-atom (constitution-atom info))
                                                      (name (atom-name constitution-atom))
                                                      (type (stype info)))
                                                 #+(or)(format t "fixed-chiral-info: ~s ~s~%" name type)
                                                 (when (eq type :undefined-configuration)
                                                   (error "Don't allow :undefined-configuration for ~a" name))
                                                 (make-instance 'stereoconfiguration
                                                                :atom-name name
                                                                :configuration type)))
                                             fixed-chiral-info)))
             (new-name-string (format nil "~A~{~A~}"
                                      name
                                      (mapcar (lambda (atom config)
                                                (declare (ignore atom))
                                                (format nil "~A"
                                                        (string config))) variable-chiral-atoms configurations)))
             (new-name (intern new-name-string (symbol-package name))))
        #+(or)
        (progn
          (format t "name -> ~a~%" name)
          (format t "variable-chiral-atoms -> ~a~%" variable-chiral-atoms)
          (format t "stereoisomer-index -> ~a~%" stereoisomer-index)
          (format t "new-name -> ~a~%" new-name))
        (make-instance 'stereoisomer
                       :name new-name
                       :pdb nil
                       :stereoconfigurations stereo-configurations
                       :stereoisomer-atoms stereoisomer-atoms))
      (progn
        ;;#+(or)(format t "Single stereoisomer name -> ~a~%" name)
        (make-instance 'stereoisomer
                       :name name
                       :stereoconfigurations nil
                       :stereoisomer-atoms stereoisomer-atoms))))

(defclass stereocenter-info ()
  ((constitution-atom :initarg :constitution-atom :accessor constitution-atom)
   (%type :initarg :type :accessor stype)))

(defmethod print-object ((obj stereocenter-info) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a" (constitution-atom obj) (stype obj))))


(defun build-stereo-information (name stereocenter-info stereoisomer-atoms)
  "Build stereoinformation from the stereocenter-info and return it."
  (multiple-value-bind (variable-chiral-atoms fixed-chiral-info)
      (loop for info in stereocenter-info
            if (eq (stype info) :both)
              collect (atom-name (constitution-atom info)) into variable-chiral-atoms
            else
              collect info into fixed-chiral-info
            finally (return (values variable-chiral-atoms fixed-chiral-info)))
    (let* ((number-of-variable-chiral-atoms (length variable-chiral-atoms))
           (number-of-stereoisomers (expt 2 number-of-variable-chiral-atoms))) ; every chiral atom multiplies # stereoisomers
      #+(or)
      (progn
        (format t "residue: ~s~%" residue)
        (format t "root atom: ~s~%" root-atom)
        (format t "chiral atoms: ~s~%" variable-chiral-atoms)
        (terpri)
        (terpri))
      (let ((stereoisomers
              (loop for stereoisomer-index below number-of-stereoisomers
                    for stereoisomer = (build-stereoisomer name variable-chiral-atoms stereoisomer-index fixed-chiral-info stereoisomer-atoms)
                    collect stereoisomer)))
        stereoisomers))))

(defun maybe-parse-charge (property-list)
  (let ((charge (getf property-list :charge)))
    (cond
      ((and charge (numberp charge))
       (round charge))
      ((null charge)
       0.0)
      (t (error "Illegal charge ~s in property list" property-list)))))

(defun parse-graph (graph plug-names residue-properties)
  (let ((constitution-atoms (make-array (hash-table-count (nodes graph)))))
    #+(or)(format t "graph = ~a~%" (name graph))
    (maphash (lambda (key node)
               (declare (ignore key))
               (when (typep node 'atom-node)
                 (let* ((property-list (property-list node))
                        (constitution-atom (make-instance 'constitution-atom
                                                          :atom-name (name node)
                                                          :element (element node)
                                                          :index (constitution-atom-index node)
                                                          :properties property-list)))
                   (setf (aref constitution-atoms (constitution-atom-index node))
                         constitution-atom))))
             (nodes graph))
    (loop for edge in (edges graph)
          for from-node = (from-node edge)
          for to-node = (to-node edge)
          for edge-type = (edge-type edge)
          for bond-order = (bond-order-from-edge-type edge-type)
          for from-index = (constitution-atom-index from-node)
          for to-index = (constitution-atom-index to-node)
          for from-constitution-atom = (aref constitution-atoms from-index)
          for to-constitution-atom = (aref constitution-atoms to-index)
          do (push (make-instance 'constitution-bond
                                  :to-atom-index from-index
                                  :order bond-order)
                   (bonds to-constitution-atom))
          do (push (make-instance 'constitution-bond
                                  :to-atom-index to-index
                                  :order bond-order)
                   (bonds from-constitution-atom)))
    (let ((plugs (make-hash-table)))
      (loop for index below (length constitution-atoms)
            for ca = (elt constitution-atoms index)
            for atom-name = (atom-name ca)
            for properties = (properties ca)
            for plug-names = (getf properties :plugs)
            when plug-names
              do (progn
                   (loop for name in plug-names
                         do (parse-plug-name name atom-name plugs graph))))
      (maphash (lambda (key plug)
                 (declare (ignore key))
                 (cond
                   ((every #'identity (plug-bonds plug)))
                   ((null (elt (plug-bonds plug) 1))
                    (setf (plug-bonds plug) (vector (elt (plug-bonds plug) 0))))
                   (t (error "Error in plug-bonds ~s" plug))))
               plugs)
      (multiple-value-bind (stereocenter-info stereoisomer-atoms)
          (loop with rev-stereocenter-info = nil
                for index below (length constitution-atoms)
                for constitution-atom = (elt constitution-atoms index)
                for prop-list = (properties constitution-atom)
                for stereochemistry-type = (getf prop-list :stereochemistry-type)
                for plug-names = (getf prop-list :plugs)
                for number-of-bonds = (+ (length (bonds constitution-atom)) (length plug-names))
                for ca-element = (element constitution-atom)
                ;;                                      do (format t "zxhc ~a~%" constitution-atom)
                collect (make-instance 'stereoisomer-atom
                                       :atom-name (atom-name constitution-atom)
                                       :constitution-atom-index index
                                       :atom-charge (maybe-parse-charge prop-list)
                                       ) into stereoisomer-atoms
                do (cond
                     ((and (eq ca-element :c) (= 4 number-of-bonds) stereochemistry-type)
                      (unless (member stereochemistry-type
                                      '(:both :s :r :left-handed :right-handed :undefined-center))
                        (error "Uknown stereochemistry-type ~a - must be :both, :S, :R, :left-handed, :right-handed" stereochemistry-type))
                      (unless (eq stereochemistry-type :undefined-center)
                        (push (make-instance 'stereocenter-info
                                             :constitution-atom constitution-atom
                                             :type stereochemistry-type) rev-stereocenter-info)))
                     ((and (eq ca-element :c) (= 4 number-of-bonds))
                      ;; prochiral atoms are always :left-handed
                      (push (make-instance 'stereocenter-info
                                           :constitution-atom constitution-atom
                                           :type :left-handed) rev-stereocenter-info))
                     (stereochemistry-type
                      (error "Do not specify stereochemistry-type ~s for ~a - it has only ~s bonds" stereochemistry-type constitution-atom number-of-bonds))
                     )
                finally (return (values (nreverse rev-stereocenter-info) stereoisomer-atoms)))
        ;; build stereochemistry here
        (let ((stereo-information (build-stereo-information (name graph) stereocenter-info stereoisomer-atoms)))
          (values (make-instance 'constitution
                                 :name (name graph)
                                 :constitution-atoms constitution-atoms
                                 :residue-properties residue-properties
                                 )
                  plugs
                  stereo-information))))))

(defparameter *topology-groups* (make-hash-table))

(defun parse-restraints (restraints)
  (loop for restraint-batch in restraints
        for kind = (car restraint-batch)
        for data = (cdr restraint-batch)
        append (cond
                 ((eq kind :dihedral-restraints)
                  (loop for restraint in data
                        collect (destructuring-bind (atom1 atom2 atom3 atom4 deg weight)
                                    restraint
                                  (make-instance 'dihedral-restraint
                                                 :atom1-name (intern (string atom1) :keyword)
                                                 :atom2-name (intern (string atom2) :keyword)
                                                 :atom3-name (intern (string atom3) :keyword)
                                                 :atom4-name (intern (string atom4) :keyword)
                                                 :dihedral-degrees deg
                                                 :weight weight))))
                 (t (error "Add support for restraint ~a" kind)))))

(defun make-abstract-topology (group-names &key residue-properties plug-names)
  (let* ((plugs (let ((ht (make-hash-table)))
                  (mapcar (lambda (plug-name)
                            (let ((plug (cond
                                          ((is-in-plug-name plug-name)
                                           (make-instance 'abstract-in-plug :name plug-name))
                                          ((is-out-plug-name plug-name)
                                           (make-instance 'abstract-out-plug :name plug-name))
                                          (t
                                           (make-instance 'abstract-out-plug :name plug-name)))))
                              (setf (gethash plug-name ht) plug)))
                          plug-names)
                  ht))
         (tops (loop for name in group-names
                     for topology = (make-instance 'topology:abstract-topology
                                                   :name name
                                                   :plugs plugs
                                                   :property-list residue-properties)
                     do (loop for group-name in group-names
                              do (pushnew name (gethash group-name *topology-groups* nil)))
                     do (if residue-properties
                            (setf (topology:property-list topology)
                                  (append (topology:property-list topology)
                                          residue-properties)))
                     do (chem:register-topology topology name)
                     collect topology)))
    tops))


(defclass rotamer-limit-data ()
  ((atom-name :initarg :atom-name :reader atom-name)
   (driver-kind :initarg :driver-kind :reader driver-kind)
   (dihedrals-deg :initarg :dihedrals-deg :reader dihedrals-deg)))


(defun ensure-dihedrals-deg (dihedral-list)
  (loop for dih-deg in dihedral-list
        collect (degrees-limit dih-deg)))

(defun parse-rotamer-limits (rotamer-limits)
  (loop for tail on rotamer-limits by #'cddr
        for kind = (first tail)
        for kind-name = (ecase kind
                          (:no-driver :no-driver)
                          (:basin-hopping :basin-hopping)
                          (:grid-search :grid-search))
        for values = (second tail)
        append (loop for entry in values
                     for atom-name = (first entry)
                     for dihedrals = (second entry)
                     collect (make-instance 'rotamer-limit-data
                                            :atom-name atom-name
                                            :driver-kind kind-name
                                            :dihedrals-deg (ensure-dihedrals-deg dihedrals)))))

(defun lookup-rotamer-limit-data (topology atom-name)
  (find atom-name (rotamer-limits topology) :key #'atom-name))


(defun validate-cluster-dihedrals (top cluster-dihedrals dihedrals constitution)
  (let ((atom-names (make-hash-table)))
    (loop for ca across (constitution-atoms constitution)
          do (setf (gethash (atom-name ca) atom-names) t))
    (loop for dih in dihedrals
          when (typep dih 'dihedral-info-atom)
            do (let ((an (atom-name dih)))
                 (unless (gethash an atom-names)
                   (error "When defining topology for ~s could not find dihedral atom name ~s~% in constitution atom names ~s"
                          (name constitution)
                          an (alexandria:hash-table-keys atom-names))))))
  (when cluster-dihedrals
    (let ((ht (make-hash-table)))
      (loop for dih in dihedrals
            for name = (cond
                         ((typep dih 'dihedral-info-atom)
                          (name dih))
                         ((typep dih 'dihedral-info-external)
                          (name dih))
                         (t (car (last dih))))
            do (setf (gethash name ht) t))
      (loop for cd in cluster-dihedrals
            for found = (gethash cd ht)
            unless found
              do (error "For topology ~s could not find cluster-dihedrals ~s in ~s " (name constitution) cd dihedrals)))))


(defun validate-constitution-vs-residue-properties (constitution residue-properties)
  (let ((total-atom-charge (round (loop with total-charge = 0
                                        for ca across (constitution-atoms constitution)
                                        for atom-charge = (getf (properties ca) :charge)
                                        when atom-charge
                                          do (incf total-charge atom-charge)
                                        finally (return total-charge)))))
    (let ((residue-net-charge (let ((net (getf residue-properties :residue-net-charge)))
                                (if (numberp net)
                                    (round net)
                                    net))))
      (cond
        ((/= total-atom-charge 0)
         (unless residue-net-charge
           (error "The total atom charges for ~s is ~d so there must be a corresponding :residue-net-charge residue-property" constitution total-atom-charge))
         (unless (= residue-net-charge total-atom-charge)
           (error "The total-atom-charge for ~s is ~d and it must match the :residue-net-charge ~d" constitution total-atom-charge residue-net-charge)))
        (t (when (and residue-net-charge (not (= residue-net-charge total-atom-charge)))
             (error ":residue-net-charge was specified as ~d and it must match total-atom-charge ~d" residue-net-charge total-atom-charge)))))))

(defun topologies-from-graph (graph group-names restraints constitution plugs stereoisomers
                              xyz-joints rotamer-limits residue-properties)
  "Create multiple topology objects and one constitution object from the graph and the arguments"
  ;; Create the topology objects
  (loop for stereoisomer in stereoisomers
        for name = (topology:name stereoisomer)
        for joint-template = (build-joint-template graph xyz-joints)
        for topology = (make-instance 'topology:topology
                                      :name name
                                      :constitution constitution
                                      :plugs plugs
                                      :joint-template joint-template
                                      :stereoisomer stereoisomer
                                      :restraints (parse-restraints restraints)
                                      :rotamer-limits (parse-rotamer-limits rotamer-limits))
        do (loop for group-name in (list* name group-names)
                 do (pushnew name (gethash group-name *topology-groups* nil)))
        do (setf (topology:property-list topology) (list* :joint-template joint-template (topology:property-list topology)))
        do (if residue-properties
               (setf (topology:property-list topology)
                     (append (topology:property-list topology)
                             residue-properties)))
        do (chem:register-topology topology name)
        collect topology))

(defun validate-properties (properties name)
  (loop for cur = properties then (cddr cur)
        for property = (car cur)
        for value = (cadr cur)
        until (null cur)
        unless (and property (member property '(:plugs :ring :stereochemistry-type :adjust :charge :build-order)))
          do (error "Invalid property ~s ~s defining ~s  properties: ~s" property value name properties)
        ))

(defun validate-residue-properties (properties name)
  (loop for cur = properties then (cddr cur)
        for property = (car cur)
        for value = (cadr cur)
        until (null cur)
        unless (and property (member property '(:residue-net-charge)))
          do (error "Invalid residue-property ~s ~s defining ~s  properties: ~s" property value name properties)
        )
  properties)


(defun recursively-order-children-atom-nodes (node counter)
  (loop for child in (children node)
        do (change-class child 'ordered-atom-node
                         :constitution-atom-index (incf counter)))
  (loop for child in (children node)
        do (setf counter (recursively-order-children-atom-nodes child counter)))
  counter)

(defun do-define-topology (name sexp &key
                                       restraints
                                       types
                                       xyz-joints
                                       residue-properties
                                       plug-names
                                       dihedrals
                                       cluster-dihedrals
                                       rings
                                       rotamer-limits)
  (when restraints
    #+(or)(format t "restraints = ~a~%" restraints))
  (let ((graph (interpret (if (consp name)
                              (first name)
                              name)
                          sexp))
        (group-names (if (consp name)
                         name
                         (list name))))
    ;;Change class of all atom-nodes to ordered-atom-nodes and set their constitution-atom-index
    ;; to their order as they appear in a depth first search
    (change-class (root-node graph) 'ordered-atom-node
                  :constitution-atom-index 0)
    (recursively-order-children-atom-nodes (root-node graph) 0)
    (unless (or (listp xyz-joints)
                (eq :all xyz-joints))
      (error "xyz-joints must be either a list of atom-names or :all"))
    (multiple-value-bind (constitution plugs stereoisomers)
        (parse-graph graph plug-names residue-properties)
      (validate-constitution-vs-residue-properties constitution residue-properties)
      (when (listp xyz-joints)
        (loop for name in xyz-joints
              unless (gethash name (nodes graph))
                do (error "The name ~s in xyz-joint is not one of the names in the topology ~s"
                          name (mapcar #'name (nodes graph)))))
      ;; Create the topology objects
      (let* ((tops (topologies-from-graph graph
                                          group-names
                                          restraints
                                          constitution
                                          plugs
                                          stereoisomers
                                          xyz-joints
                                          rotamer-limits
                                          residue-properties)))
        (when types
          (loop for name-type in types
                for name = (first name-type)
                for type = (second name-type)
                for ca = (constitution-atom-named constitution name)
                do (setf (atom-type ca) type)))
        (let ((dihedral-info (if dihedrals
                                 (parse-dihedral-info dihedrals)
                                 (create-dihedral-info-from-constitution constitution)))
              (ring-info (when rings
                           (parse-ring-info rings))))
          (setf (residue-properties constitution)
                (list* :dihedrals dihedral-info (residue-properties constitution)))
          (when ring-info
            (setf (residue-properties constitution)
                  (list* :rings ring-info (residue-properties constitution))))
          (when cluster-dihedrals
            (validate-cluster-dihedrals :name cluster-dihedrals dihedral-info constitution)
            (setf (residue-properties constitution)
                  (list* :cluster-dihedrals cluster-dihedrals (residue-properties constitution)))))
        tops))))


(defmacro define-topology (name sexp &key restraints types xyz-joints dihedrals cluster-dihedrals rings residue-properties plugs rotamer-limits)
  `(do-define-topology ',name ',sexp
     :restraints ',restraints
     :dihedrals ',dihedrals
     :residue-properties ',residue-properties
     :cluster-dihedrals ',cluster-dihedrals
     :plug-names ',plugs
     :rings ',rings
     :xyz-joints ',xyz-joints
     :rotamer-limits ',rotamer-limits
     ))

(defun do-define-abstract-topology (name &key residue-properties plug-names)
  (make-abstract-topology (list name)
                          :properties (validate-residue-properties residue-properties)
                          :plug-names plug-names
                          ))

(defmacro define-abstract-topology (name &key residue-properties plugs)
  "Abstract-topoology is for development of topology and foldamer rules.
Abstract-topology doesn't have any atoms.  Once you have an abstract-topology that
can build spiroligomers then you want to convert it into a real topology."
  `(do-define-abstract-topology ',name 
     :residue-properties ',residue-properties
     :plug-names ',plugs
     ))


(defun molecule-from-sexp (sexp)
  (let ((graph (parse-atom-graph sexp))
        (res (chem:make-residue :res))
        (mol (chem:make-molecule :mol))
        (agg (chem:make-aggregate :agg))
        (node-to-atom (make-hash-table)))
    (chem:add-matter agg mol)
    (chem:add-matter mol res)
    (maphash (lambda (name node)
               (let ((atom (chem:make-atom (name node) (element node))))
                 (chem:add-matter res atom)
                 (setf (gethash node node-to-atom) atom)))
             (nodes graph))
    (loop for edge in (edges graph)
          for node1 = (from-node edge)
          for atom1 = (gethash node1 node-to-atom)
          for node2 = (to-node edge)
          for atom2 = (gethash node2 node-to-atom)
          for bond-order = (bond-order-from-edge-type (edge-type edge))
          do (chem:bond-to atom1 atom2 bond-order))
    (values agg graph node-to-atom)))
