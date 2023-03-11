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
   (constitution-atom-index :initarg :constitution-atom-index :accessor constitution-atom-index)
   (children :initform nil :accessor children)
   ))

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
   (next-atom-index :initform 0 :accessor next-atom-index)))

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

(defun node-from-symbol-entry (keyword-atm graph)
  (let ((keyword-element (chem:element-from-atom-name-string (symbol-name keyword-atm))))
    (when (bond-symbol keyword-atm)
      (error "Don't create atoms with name ~a" keyword-atm))
    (make-instance 'atom-node
                   :name keyword-atm
                   :element keyword-element
                   :constitution-atom-index (prog1 (next-atom-index graph) (incf (next-atom-index graph))))))

(defun ensure-keyword (name)
  (intern (symbol-name name) :keyword))

(defun parse-atom (atm graph)
  (cond
    ((symbolp atm)
     (let* ((keyword-atm (ensure-keyword atm))
            (seen-node (gethash keyword-atm (nodes graph))))
       (if seen-node
           (values seen-node t)
           (let ((node (node-from-symbol-entry keyword-atm graph)))
             (setf (gethash keyword-atm (nodes graph)) node)
             node))))
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
                             pl))
            (node (progn
                    (when (bond-symbol keyword-atm)
                      (error "Don't create atoms with name ~a" keyword-atm))
                    (make-instance 'atom-node :name keyword-atm
                                              :element element
                                              :property-list property-list
                                              :constitution-atom-index (prog1 (next-atom-index graph) (incf (next-atom-index graph)))))))
       (setf (gethash keyword-atm (nodes graph)) node)
       node))
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
                 (setf (children prev-node) (append (children prev-node) (list node))))
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
                 (setf (children prev-node) (append (children prev-node) (list node))))
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
                                                      (type (type info)))
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
   (type :initarg :type :accessor type)))

(defmethod print-object ((obj stereocenter-info) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a" (constitution-atom obj) (type obj))))


(defun build-stereo-information (name stereocenter-info stereoisomer-atoms)
  "Build stereoinformation from the stereocenter-info and return it."
  (multiple-value-bind (variable-chiral-atoms fixed-chiral-info)
      (loop for info in stereocenter-info
            if (eq (type info) :both)
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

(defun parse-graph (graph)
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
          (loop with stereocenter-info = nil
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
                                       ;; Fill in other info (atom-charge atom-type) from prop-list if available
                                       ) into stereoisomer-atoms
                do (cond
                     ((and (eq ca-element :c) (= 4 number-of-bonds) stereochemistry-type)
                      (unless (member stereochemistry-type
                                      '(:both :s :r :left-handed :right-handed :undefined-center))
                        (error "Uknown stereochemistry-type ~a - must be :both, :S, :R, :left-handed, :right-handed" stereochemistry-type))
                      (unless (eq stereochemistry-type :undefined-center)
                        (push (make-instance 'stereocenter-info
                                             :constitution-atom constitution-atom
                                             :type stereochemistry-type) stereocenter-info)))
                     ((and (eq ca-element :c) (= 4 number-of-bonds))
                      ;; prochiral atoms are always :left-handed
                      (push (make-instance 'stereocenter-info
                                           :constitution-atom constitution-atom
                                           :type :left-handed) stereocenter-info))
                     (stereochemistry-type
                      (error "Do not specify stereochemistry-type ~s for ~a" stereochemistry-type constitution-atom))
                     )
                finally (return (values stereocenter-info stereoisomer-atoms)))
        ;; build stereochemistry here
        (let ((stereo-information (build-stereo-information (name graph) stereocenter-info stereoisomer-atoms)))
          (values (make-instance 'constitution
                                 :name (name graph)
                                 :constitution-atoms constitution-atoms
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
                        collect (destructuring-bind (atom1 atom2 atom3 atom4 min-deg max-deg weight)
                                    restraint
                                  (make-instance 'dihedral-restraint
                                                 :atom1-name (intern (string atom1) :keyword)
                                                 :atom2-name (intern (string atom2) :keyword)
                                                 :atom3-name (intern (string atom3) :keyword)
                                                 :atom4-name (intern (string atom4) :keyword)
                                                 :dihedral-min-degrees min-deg
                                                 :dihedral-max-degrees max-deg
                                                 :weight weight))))
                 (t (error "Add support for restraint ~a" kind)))))

(defun topologies-from-graph (graph group-names restraints &key types)
  (multiple-value-bind (constitution plugs stereoisomers)
      (parse-graph graph)
    (let* ((tops (loop for stereoisomer in stereoisomers
                       for name = (topology:name stereoisomer)
                       for joint-template = (build-joint-template graph)
                       for topology = (make-instance 'topology:topology
                                                     :name name
                                                     :constitution constitution
                                                     :plugs plugs
                                                     :joint-template joint-template
                                                     :stereoisomer stereoisomer
                                                     :restraints (parse-restraints restraints))
                       do (loop for group-name in (list* name group-names)
                                do (pushnew name (gethash group-name *topology-groups* nil)))
                       do (setf (topology:property-list topology) (list* :joint-template joint-template (topology:property-list topology)))
                       do (chem:register-topology topology name)
                       collect topology)))
      (when types
        (loop for name-type in types
              for name = (first name-type)
              for type = (second name-type)
              for ca = (constitution-atom-named constitution name)
              do (setf (atom-type ca) type)))
      tops)))

(defun do-define-topology (name sexp &key restraints types)
  (when restraints
    #+(or)(format t "restraints = ~a~%" restraints))
  (let ((graph (interpret (if (consp name)
                              (first name)
                              name)
                          sexp))
        (group-names (if (consp name)
                         name
                         (list name))))
    (topologies-from-graph graph group-names restraints :types types)))


(defmacro define-topology (name sexp &key restraints types)
  `(do-define-topology ',name ',sexp :restraints ',restraints))

