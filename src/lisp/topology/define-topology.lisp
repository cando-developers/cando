(in-package :topology)

#|
Create topology instances using a graph described using an s-expression.

(topology:define-topology 'gly '((C C :plugs (+dkp.1)) (= O)
                                 - CA (- HA1) (- HA2)
                                 - (N N :plugs (+dkp.0))
                                 - CN (- HN1) (- HN2) - HN3))

(topology:define-topology 'pro4 '((CAC C :plugs (-dkp.0)) (= OA)
                                  - (CA :C :stereochemistry-type :chiral) (- HA)
                                  - CB (- HB1) (- HB2)
                                  - (CG :C :stereochemistry-type :chiral)
                                  (- (NG N :plugs (+side.0 +dkp.0)))
                                  (- (CGC C :plugs (+dkp.1)) = OG)
                                  - CD (- HD1) (- HD2) -
                                  (NA N :plugs (-dkp.1)) - CA ))

(topology:define-topology 'pro '((C C :plugs (-dkp.0)) (= O)
                                 - (CA C :stereochemistry-type :chiral) (- HA)
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
        (setf (elt (atom-names plug) plug-index) atom-name)))))

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
              (property-list (cddr atm))
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



(defun build-stereoisomer (name chiral-atoms stereoisomer-index)
  "Build a stereoisomer. Assemble a name for the stereiosomer using the name and
the stereoisomer-index.   The stereoisomer-index is like a bit-vector that indicates
the stereochemistry of this stereoisomer (0 is :S and 1 is :R).
The name is constructed by appending \"{chiral-atom-name/(S|R),...}\".
So if name is \"ALA\" and stereoisomer-index is 1 the name becomes ALA{CA/S}."
  (if chiral-atoms
      (let* ((configurations (loop for chiral-index below (length chiral-atoms)
                                   for config = (if (logbitp chiral-index stereoisomer-index) :S :R)
                                   collect config))
             (stereo-configurations (mapcar (lambda (atom config)
                                              (make-instance 'stereoconfiguration
                                                             :atom-name (atom-name atom)
                                                             :configuration config))
                                            chiral-atoms configurations))
             (new-name-string (format nil "~A~{~A~}"
                                      name
                                      (mapcar (lambda (atom config)
                                                (declare (ignore atom))
                                                (format nil "~A"
                                                        (string config))) chiral-atoms configurations)))
             (new-name (intern new-name-string (symbol-package name))))
        #+(or)
        (progn
          (format t "name -> ~a~%" name)
          (format t "chiral-atoms -> ~a~%" chiral-atoms)
          (format t "stereoisomer-index -> ~a~%" stereoisomer-index)
          (format t "new-name -> ~a~%" new-name))
        (make-instance 'stereoisomer
                       :name new-name
                       :pdb nil
                       :stereoisomer-index stereoisomer-index
                       :configurations stereo-configurations))
      (progn
        ;;#+(or)(format t "Single stereoisomer name -> ~a~%" name)
        (make-instance 'stereoisomer
                       :name name
                       :name name
                       :stereoisomer-index 0
                       :configurations nil))))

(defun build-stereo-information (name chiral-constitution-atoms)
  "Build stereoinformation from all of this and return it."
  (let* ((chiral-atoms chiral-constitution-atoms)
         (number-of-chiral-atoms (length chiral-atoms))
         (number-of-stereoisomers (expt 2 number-of-chiral-atoms)))
    #+(or)
    (progn
      (format t "residue: ~s~%" residue)
      (format t "root atom: ~s~%" root-atom)
      (format t "chiral atoms: ~s~%" chiral-atoms)
      (terpri)
      (terpri))
    (let ((stereoisomers (loop for stereoisomer-index below number-of-stereoisomers
                               for stereoisomer = (build-stereoisomer name chiral-atoms stereoisomer-index)
                               collect stereoisomer)))
      stereoisomers)))

(defun constitution-atoms-from-graph (graph)
  (let ((constitution-atoms (make-array (hash-table-count (nodes graph)))))
    (maphash (lambda (key node)
               (declare (ignore key))
               (when (typep node 'atom-node)
                 (let* ((property-list (property-list node))
                        (stereochemistry-type (getf property-list :stereochemistry-type :undefined-center)))
                   (setf (aref constitution-atoms (constitution-atom-index node))
                         (make-instance 'constitution-atom
                                        :atom-name (name node)
                                        :element (element node)
                                        :index (constitution-atom-index node)
                                        :properties property-list
                                        :stereochemistry-type stereochemistry-type)))))
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
                 (setf (atom-names plug) (subseq (atom-names plug) 0 (position nil (atom-names plug)))))
               plugs)
      (let* ((stereocenters (loop for index below (length constitution-atoms)
                                  for ca = (elt constitution-atoms index)
                                  for prop-list = (properties ca)
                                  for stereochemistry-type = (stereochemistry-type ca)
                                  when (eq :chiral stereochemistry-type)
                                    collect ca))
             (stereo-information (build-stereo-information (name graph) stereocenters)))
        (values constitution-atoms plugs stereo-information)))))


(defun constitution-from-graph (graph)
  (multiple-value-bind (constitution-atoms plugs stereo-information)
      (constitution-atoms-from-graph graph)
    (make-instance 'constitution
                   :name (name graph)
                   :constitution-atoms constitution-atoms
                   :plugs plugs
                   :stereo-information stereo-information)))

(defun topologies-from-graph (graph)
  (let* ((constitution (topology:constitution-from-graph graph))
         (plugs (topology:plugs constitution))
         (stereo-information (topology:stereo-information constitution))
         (tops (loop for stereoisomer in stereo-information
                     for name = (topology:name stereoisomer)
                     for configurations = (topology:configurations stereoisomer)
                     for joint-template = (build-joint-template graph)
                     for topology = (make-instance 'topology:topology
                                                   :name name
                                                   :constitution constitution
                                                   :plugs plugs
                                                   :joint-template joint-template
                                                   :stereoisomer-atoms configurations)
                     do (push topology (topology:topology-list constitution))
                     do (format t "stereoisomer ~a ~s~%" name configurations)
                     do (setf (topology:property-list topology) (list* :joint-template joint-template (topology:property-list topology)))
                     do (cando:register-topology topology name)
                     collect topology)))
    (setf (topology:topology-list constitution) tops)
    tops))


(defun do-define-topology (name sexp)
  (let ((graph (interpret name sexp)))
    (topologies-from-graph graph)))


(defmacro define-topology (name sexp)
  `(do-define-topology ',name ',sexp))

