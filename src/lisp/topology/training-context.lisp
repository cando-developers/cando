(in-package #:topology)

;;; Declaration syntax only: no topology lookup, matcher construction, or
;;; changes to the existing training pipeline happen here.
(define-condition training-context-parse-error (simple-error)
  ((expression :initarg :expression :reader training-context-error-expression)))

(defclass shape-kind-properties (cando.serialize:serializable) ())

(defclass backbone-shape-kind-properties (shape-kind-properties)
  ((paired-with-sidechains :initarg :paired-with-sidechains :initform nil
                          :reader shape-kind-paired-with-sidechains)
   (backbone-shape-key-selector :initarg :backbone-shape-key-selector :initform nil
                               :reader shape-kind-backbone-shape-key-selector)))

(defclass sidechain-shape-kind-properties (shape-kind-properties)
  ((paired-with :initarg :paired-with :initform nil :reader shape-kind-paired-with)
   (selectors :initarg :selectors :initform nil :reader shape-kind-selectors)))

(defclass training-context (cando.serialize:serializable)
  ((name :initarg :name :reader training-context-name)
   (focus :initarg :focus :reader training-context-focus)
   (nodes :initarg :nodes :reader training-context-nodes)
   (connections :initarg :connections :reader training-context-connections)
   (caps :initarg :caps :reader training-context-caps)
   (shape-kind-properties :initarg :shape-kind-properties
                          :reader training-context-shape-kind-properties)
   (expression :initarg :expression :reader training-context-expression)))

;;; Convenience readers derive their answers from the properties object;
;;; the context does not duplicate its state.
(defmethod training-context-shape-kind ((context training-context))
  (etypecase (training-context-shape-kind-properties context)
    (backbone-shape-kind-properties :backbone)
    (sidechain-shape-kind-properties :sidechain)))

(defmethod training-context-paired-with-sidechains ((context training-context))
  (let ((properties (training-context-shape-kind-properties context)))
    (when (typep properties 'backbone-shape-kind-properties)
      (shape-kind-paired-with-sidechains properties))))

(defmethod training-context-paired-with ((context training-context))
  (let ((properties (training-context-shape-kind-properties context)))
    (when (typep properties 'sidechain-shape-kind-properties)
      (shape-kind-paired-with properties))))

(defmethod training-context-selectors ((context training-context))
  (let ((properties (training-context-shape-kind-properties context)))
    (when (typep properties 'sidechain-shape-kind-properties)
      (shape-kind-selectors properties))))

(defclass training-context-node (cando.serialize:serializable)
  ((monomer :initarg :monomer :reader training-context-node-monomer)
   (label :initarg :label :reader training-context-node-label)))

(defclass training-context-connection (cando.serialize:serializable)
  ((source :initarg :source :reader training-context-connection-source)
   (plug :initarg :plug :reader training-context-connection-plug)
   (target :initarg :target :reader training-context-connection-target)))

(defclass training-context-cap (cando.serialize:serializable)
  ;; NODE is the attachment parent: a context node, or another training cap.
  ((label :initarg :label :initform nil :reader training-context-cap-label)
   (node :initarg :node :reader training-context-cap-node)
   (plug :initarg :plug :reader training-context-cap-plug)
   (monomer :initarg :monomer :reader training-context-cap-monomer)))

(defclass training-context-selector (cando.serialize:serializable)
  ((node :initarg :node :reader training-context-selector-node)
   (cache-field :initarg :cache-field :reader training-context-selector-cache-field)
   (output-name :initarg :output-name :reader training-context-selector-output-name)))

(defun validate-shape-kind-options (options allowed)
  "Validate a shape-kind property list for use by foldamer-specific methods."
  (labels ((bad ()
             (error 'training-context-parse-error :expression options
                    :format-control "Invalid shape-kind options ~s; allowed keys: ~s"
                    :format-arguments (list options allowed))))
    (unless (and (listp options)
                 (handler-case (let ((length (list-length options)))
                                 (and length (evenp length)))
                   (type-error () nil)))
      (bad))
    (let ((seen nil))
      (loop for key in options by #'cddr
            do (unless (and (member key allowed) (not (member key seen))) (bad))
               (push key seen)))))

(defgeneric parse-training-shape-kind (foldamer kind options &key selectors legacy-p)
  (:documentation "Validate shape options and return a SHAPE-KIND-PROPERTIES object.
Specialize on a foldamer class and KIND to implement domain-specific options.
SELECTORS contains already resolved training-context-selector objects.
LEGACY-P allows a foldamer to preserve the semantics of old declarations."))

(defmethod parse-training-shape-kind (foldamer (kind (eql :backbone)) options &key selectors legacy-p)
  (declare (ignore foldamer kind selectors legacy-p))
  (validate-shape-kind-options options '(:paired-with-sidechains :backbone-shape-key-selector))
  (let ((paired (getf options :paired-with-sidechains)))
    (unless (member paired '(nil t))
      (error 'training-context-parse-error :expression options
             :format-control ":paired-with-sidechains must be T or NIL: ~s"
             :format-arguments (list options)))
    (make-instance 'backbone-shape-kind-properties :paired-with-sidechains paired
                   :backbone-shape-key-selector (getf options :backbone-shape-key-selector))))

(defmethod parse-training-shape-kind (foldamer (kind (eql :sidechain)) options &key selectors legacy-p)
  (declare (ignore foldamer kind legacy-p))
  (validate-shape-kind-options options '(:paired-with :shape-key-sources))
  (let ((paired (getf options :paired-with)))
    (when (and (get-properties options '(:paired-with))
               (not (and (stringp paired) (plusp (length paired)))))
      (error 'training-context-parse-error :expression options
             :format-control ":paired-with must be a nonempty trainer name: ~s"
             :format-arguments (list options)))
    (make-instance 'sidechain-shape-kind-properties :paired-with paired :selectors selectors)))

(defun parse-training-context (expression &key foldamer)
  "Parse (TRAINING-CONTEXT :NAME string :CONTEXT path ...).
A node is a non-keyword symbol or (symbol :AS label). A path starts with a
node and continues with plug/node pairs or parenthesized branches. Branches
start with a plug and resume the enclosing path at their attachment node.
:TRAINING-CAPS contains (label plug monomer) entries. :SHAPE-KIND defaults
to (BACKBONE), with optional :PAIRED-WITH-SIDECHAINS boolean. Additional
options are handled by PARSE-TRAINING-SHAPE-KIND specialized on FOLDAMER.
A sidechain uses
(SIDECHAIN :PAIRED-WITH trainer-name :SHAPE-KEY-SOURCES ((label field ...) ...)).
Sidechain options are optional; pairing semantics are validated elsewhere.
Kind symbols are recognized by name and normalized to keywords. Cache fields
are keywords and labels are local symbols compared by EQ. Pairing and source
options are accepted only inside :SHAPE-KIND, never at the top level.
This validates grammar and local references, not chemical compatibility."
  (let ((nodes nil) (connections nil) (caps nil) (selectors nil)
        (labels (make-hash-table :test #'eq))
        (ports (make-hash-table :test #'eq))
        (output-names (make-hash-table :test #'eq)))
    (labels ((bad (form control &rest args)
               (error 'training-context-parse-error :expression form
                      :format-control "Invalid training-context form ~s: ~a"
                      :format-arguments (list form (apply #'format nil control args))))
             (proper (form)
               (unless (and (listp form)
                            (handler-case (integerp (list-length form))
                              (type-error () nil)))
                 (bad form "expected a proper, non-circular list"))
               form)
             (name-p (name)
               (and (symbolp name) name (not (eq name t)) (not (keywordp name))))
             (check-options (options allowed)
               (proper options)
               (unless (evenp (length options)) (bad options "option has no value"))
               (let ((seen nil))
                 (loop for key in options by #'cddr
                       do (unless (member key allowed) (bad key "unknown or misplaced option"))
                          (when (member key seen) (bad key "duplicate option"))
                          (push key seen))))
             (parse-kind (form)
               (proper form)
               (let* ((symbol (first form))
                      (kind (and (symbolp symbol)
                                 (cond ((string= (symbol-name symbol) "BACKBONE") :backbone)
                                       ((string= (symbol-name symbol) "SIDECHAIN") :sidechain))))
                      (properties (rest form)))
                 (unless kind (bad form "expected (backbone ...) or (sidechain ...)"))
                 ;; Domain-specific keys are checked by the generic hook.
                 (unless (evenp (length properties)) (bad form "option has no value"))
                 (values kind properties)))
             (plug-p (plug)
               (and (keywordp plug) (> (length (symbol-name plug)) 1)
                    (find (char (symbol-name plug) 0) "+-")))
             (node (form)
               (let ((monomer form) (label nil))
                 (when (consp form)
                   (proper form)
                   (unless (and (= (length form) 3) (eq (second form) :as)
                                (name-p (third form)))
                     (bad form "expected (monomer :as label)"))
                   (setf monomer (first form) label (third form)))
                 (unless (name-p monomer) (bad form "expected a monomer/group symbol"))
                 (when (and label (gethash label labels))
                   (bad form "duplicate label ~s" label))
                 (let ((node (make-instance 'training-context-node
                                            :monomer monomer :label label)))
                   (push node nodes)
                   (when label (setf (gethash label labels) node))
                   node)))
             (reserve (node plug form)
               (unless (plug-p plug) (bad form "expected a plug keyword beginning with + or -"))
               (when (member plug (gethash node ports) :test #'eq)
                 (bad form "plug ~s is already used on ~s" plug
                      (or (training-context-node-label node)
                          (training-context-node-monomer node))))
               (push plug (gethash node ports)))
             (walk (tail current)
               (proper tail)
               (loop while tail
                     for form = (pop tail)
                     do (if (consp form)
                            (progn
                              (proper form)
                              (unless (plug-p (first form))
                                (bad form "a branch must start with a plug"))
                              (walk form current))
                            (progn
                              (reserve current form form)
                              (unless tail (bad form "missing destination monomer"))
                              (let ((target (node (pop tail))))
                                (push (make-instance 'training-context-connection
                                                     :source current :plug form :target target)
                                      connections)
                                (setf current target))))))
             (lookup (label form)
               (or (and (name-p label) (gethash label labels))
                   (bad form "undefined node label ~s" label))))
      (proper expression)
      (unless (and (name-p (first expression))
                   (string= (symbol-name (first expression)) "TRAINING-CONTEXT"))
        (bad expression "expected (training-context ...)"))
      (let ((options (rest expression)))
        (check-options options '(:name :context :training-caps :shape-kind))
        (multiple-value-bind (kind properties)
            (parse-kind (getf options :shape-kind '(backbone)))
        (let ((name (getf options :name)) (path (getf options :context)))
          (unless (and (stringp name) (plusp (length name)))
            (bad name ":name must be a nonempty string"))
          (proper path)
          (unless path (bad path ":context must contain a focus monomer"))
          (let ((focus (node (first path))))
            (walk (rest path) focus)
            (dolist (form (proper (getf options :training-caps)))
              (proper form)
              (unless (= (length form) 3) (bad form "expected (label plug monomer)"))
              (let ((node (lookup (first form) form)))
                (unless (name-p (third form)) (bad form "cap must name a monomer/group"))
                (reserve node (second form) form)
                (push (make-instance 'training-context-cap :node node
                                     :plug (second form) :monomer (third form)) caps)))
            (dolist (form (proper (getf properties :shape-key-sources)))
              (proper form)
              (unless (>= (length form) 2) (bad form "expected (label field ...)"))
              (let ((node (lookup (first form) form)))
                (dolist (field (rest form))
                  (unless (keywordp field) (bad form "cache fields must be keywords"))
                  (when (gethash field output-names)
                    (bad form "duplicate shape-key output name ~s" field))
                  (setf (gethash field output-names) t)
                  (push (make-instance 'training-context-selector :node node
                                       :cache-field field :output-name field) selectors))))
            (make-instance 'training-context :name name :focus focus
                           :nodes (nreverse nodes) :connections (nreverse connections)
                           :caps (nreverse caps)
                           :shape-kind-properties
                           (parse-training-shape-kind foldamer kind properties
                                                      :selectors (nreverse selectors))
                           :expression expression))))))))

(defun parse-legacy-training-context (expression &key foldamer)
  "Convert one legacy (:DESCRIPTION ... :PATTERN ... :SHAPE-KIND ...) form.
An optional PATTERN wrapper is accepted. Branches retain their attachment
point; everything downstream of a :CAP stays in the training environment.
Referenced context nodes receive FOO, BAR, QUX, WAL, then NODE-5, NODE-6, ...
in the TOPOLOGY package. Unreferenced nodes stay unlabeled. No input is mutated.
Ring-cap/wildcard declarations are rejected rather than losing their meaning."
  (let ((nodes nil) (edges nil) (caps nil) (selectors nil)
        (ports (make-hash-table :test #'eq)) (names (make-hash-table :test #'eq))
        (label-counter 0))
    (labels ((bad (form message)
               (error 'training-context-parse-error :expression form
                      :format-control "Invalid legacy training declaration ~s: ~a"
                      :format-arguments (list form message)))
             (proper (form)
               (unless (and (listp form)
                            (handler-case (integerp (list-length form))
                              (type-error () nil)))
                 (bad form "expected a proper non-circular list"))
               form)
             (monomer-p (x)
               (and (symbolp x) x (not (eq x t)) (not (keywordp x))))
             (plug-p (x)
               (and (keywordp x) (> (length (symbol-name x)) 1)
                    (find (char (symbol-name x) 0) "+-")))
             (inverse (plug)
               (let ((name (copy-seq (symbol-name plug))))
                 (setf (char name 0) (if (char= (char name 0) #\+) #\- #\+))
                 (intern name :keyword)))
             (label-node (node)
               (when (null (slot-value node 'label))
                 (setf (slot-value node 'label)
                       (or (nth label-counter '(foo bar qux wal))
                           (intern (format nil "NODE-~d" (1+ label-counter)) :topology)))
                 (incf label-counter)))
             (reserve (node plug)
               (unless (plug-p plug) (bad plug "expected a plug keyword"))
               (when (member plug (gethash node ports)) (bad plug "duplicate attachment plug"))
               (push plug (gethash node ports)))
             (new-node (name)
               (unless (monomer-p name) (bad name "expected a monomer/group symbol"))
               (let ((node (make-instance 'training-context-node :monomer name :label nil)))
                 (push node nodes) node))
             (walk (tail current)
               (proper tail)
               (loop while tail for token = (pop tail)
                     do (if (consp token)
                            (progn
                              (proper token)
                              (unless (plug-p (first token)) (bad token "branch must start with a plug"))
                              (walk token current))
                            (progn
                              (reserve current token)
                              (unless tail (bad token "missing destination"))
                              (let* ((destination (pop tail))
                                     (cap-p (and (consp destination) (eq (first destination) :cap)))
                                     (name (if cap-p
                                               (progn
                                                 (proper destination)
                                                 (unless (= (length destination) 2)
                                                   (bad destination "expected (:cap monomer)"))
                                                 (second destination))
                                               destination)))
                                (unless (monomer-p name)
                                  (bad destination "unsupported monomer form (including ring-cap)"))
                                (if (or cap-p (typep current 'training-context-cap))
                                    (let ((cap (make-instance 'training-context-cap
                                                              :node current :plug token :monomer name)))
                                      (label-node current)
                                      (push cap caps)
                                      (setf current cap))
                                    (let ((node (new-node name)))
                                      (push (make-instance 'training-context-connection
                                                           :source current :plug token :target node) edges)
                                      (setf current node)))
                                (reserve current (inverse token)))))))
             (step-node (current plug form)
               (unless (plug-p plug) (bad form "invalid source path plug"))
               (or (loop for edge in edges
                         when (and (eq current (training-context-connection-source edge))
                                   (eq plug (training-context-connection-plug edge)))
                           return (training-context-connection-target edge)
                         when (and (eq current (training-context-connection-target edge))
                                   (eq plug (inverse (training-context-connection-plug edge))))
                           return (training-context-connection-source edge))
                   (loop for cap in caps
                         when (and (eq current (training-context-cap-node cap))
                                   (eq plug (training-context-cap-plug cap)))
                           return cap
                         when (and (eq current cap)
                                   (eq plug (inverse (training-context-cap-plug cap))))
                           return (training-context-cap-node cap))
                   (bad form "shape-key path leaves the training molecule"))))
      (proper expression)
      (let* ((options (if (and (monomer-p (first expression))
                               (string= (symbol-name (first expression)) "PATTERN"))
                          (rest expression) expression)))
        (validate-shape-kind-options options '(:description :pattern :shape-kind))
        (let* ((name (getf options :description)) (path (proper (getf options :pattern)))
               (kind-form (getf options :shape-kind :backbone))
               (kind (if (consp kind-form) (first (proper kind-form)) kind-form))
               (properties (if (consp kind-form) (copy-list (rest kind-form)) nil)))
          (unless (and (stringp name) (plusp (length name))) (bad name "missing description"))
          (unless (member kind '(:backbone :sidechain)) (bad kind "unknown shape kind"))
          ;; Check structure/duplicates before GETF, leaving option ownership to the hook.
          (proper properties)
          (unless (evenp (length properties)) (bad properties "option has no value"))
          (validate-shape-kind-options properties (loop for key in properties by #'cddr collect key))
          (unless path (bad path "empty pattern"))
          (let* ((focus (new-node (first path))) (current focus)
                 (sources (getf properties :shape-key-sources)))
            (walk (rest path) focus)
            (when (get-properties properties '(:shape-key-sources))
              (proper sources)
              (unless (and (eq (first sources) :path) (rest sources))
                (bad sources "expected (:path (plug field ...) ...)"))
              (dolist (step (rest sources))
                (proper step)
                (setf current (step-node current (first step) step))
                (dolist (field (rest step))
                  (let ((cache field) (output field))
                    (when (consp field)
                      (proper field)
                      (unless (= (length field) 2) (bad field "expected (cache-field output-name)"))
                      (setf cache (first field) output (second field)))
                    (unless (and (keywordp cache) (keywordp output)) (bad field "fields must be keywords"))
                    (when (gethash output names) (bad field "duplicate shape-key output name"))
                    (setf (gethash output names) t)
                    (label-node current)
                    (push (make-instance 'training-context-selector :node current
                                         :cache-field cache :output-name output) selectors))))
              (unless selectors (bad sources "no cache fields selected")))
            (make-instance 'training-context :name name :focus focus
                           :nodes (nreverse nodes) :connections (nreverse edges) :caps (nreverse caps)
                           :shape-kind-properties
                           (parse-training-shape-kind foldamer kind properties
                                                      :selectors (nreverse selectors) :legacy-p t)
                           :expression expression)))))))
