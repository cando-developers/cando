(cl:defpackage #:language.smarts.language-server
  (:use
   #:cl))

(cl:in-package #:language.smarts.language-server)

(defclass sm***s-document (protocol.language-server::document)
  ((index  :accessor index)
   (tracer :accessor tracer)
   (ast    :accessor ast)))

(defclass smiles-document (sm***s-document)
  ())

(defclass smarts-document (sm***s-document)
  ())

(defclass source-tracking-and-indexing-builder
    (text.source-location.source-tracking-builder:source-tracking-builder)
  ((index :initarg :index
          :reader  index)))

(defmethod text.source-location.source-tracking-builder:make-location
    :around ((builder source-tracking-and-indexing-builder)
             (bounds  t))
  (let ((location (call-next-method)))
    (when location
      (text.source-location.lookup:add! location (index builder)))
    location))

(defmethod (setf protocol.language-server::%text) :after ((text string) (document smarts-document))
  (handler-case
      (let* ((index   (make-instance 'text.source-location.lookup::index))
             (tracer  (model.transform.trace:make-tracer))
             (builder (make-instance 'source-tracking-and-indexing-builder
                                     :index  index
                                     :source text
                                     :target 'list)))
        (setf (index document)  index
              (tracer document) tracer
              (ast    document) (model.transform.trace:with-tracer (tracer)
                                  (model.transform.trace:with-transform (:parse)
                                    (language.smarts.parser:parse text builder))))
        (log:warn (ast document) index tracer)
        (protocol.language-server.methods::publish-diagnostics
         document '()))
    (esrap:esrap-parse-error (condition)
      (let ((position (esrap:esrap-error-position condition)))
        (protocol.language-server.methods::publish-diagnostics
         document
         (text.source-location:make-annotation
          (text.source-location:make-location
           (text.source-location:make-source "todo" :content text)
           position (1+ position))
          (princ-to-string condition)
          :kind :error))))))

(defmethod protocol.language-server.methods:completion ((workspace t)
                                                        (document  smarts-document)
                                                        (position  t))
  (let* ((text   (protocol.language-server:text document))
         (index  (text.source-location:index position))
         (prefix (subseq text 0 index)))
    (log:info prefix)
    (multiple-value-bind (result condition) ; TODO store AST or condition in document, whatever was the result of the previous parse
        (ignore-errors (language.smarts.parser:parse prefix 'list))
      (cond
        ((typep condition 'esrap:esrap-parse-error)
         (coerce (mapcan (lambda (expected)
                           (when (typep expected '(or string character))
                             (list `((:label . ,(string expected))))))
                         (esrap::result-expected-input (esrap:esrap-parse-error-result condition)))
                 'vector))))))

(defun describe-node (node)
  (apply (lambda (kind &key symbol value class count operator label ((:kind kind*)) &allow-other-keys)
           (case kind
             (:atom
              (format nil "~(~A~)~@[ ~A~] atom" kind* symbol))
             (:charge
              (format nil "charge ~D" value))
             (:chirality
              (format nil "chirality~@[ class ~A~] count ~D"
                      class count))
             (:atom-map-class
              (format nil "atom map class ~D" label))
             ((:unary-operator :binary-operator)
              (format nil "~(~A~) operator" operator))
             (:bracketed-expression
              (format nil "atom expression"))
             (:labeled
              (format nil "ring label ~D" label))
             (:bond
              (format nil "~(~A~) bond" kind*))
             (:chain
              (format nil "a chain consisting of ~D atom~:P"
                      (length (architecture.builder-protocol:node-relation 'list :element node))))
             (t
              "a node")))
         (architecture.builder-protocol:node-kind 'list node)
         (architecture.builder-protocol:node-initargs 'list node)))

(defun nodes-at-position (document position)
  (let ((locations (text.source-location.lookup:lookup
                    position (index document))))
    (values (alexandria:mappend (lambda (location)
                                  (model.transform.trace:targets-for-source
                                   (tracer document) location))
                                locations)
            locations)))

(defun node-at-position (document position)
  (multiple-value-bind (nodes locations)
      (nodes-at-position document position)
    (log:info nodes locations)
    (values (first nodes) (first locations))))

(defmethod protocol.language-server.methods:highlight-in-document
    ((workspace t)
     (document  smarts-document)
     (version   t)
     (position  t))
  (multiple-value-bind (node location)
      (node-at-position document position)
    (log:info node location)
    `#(((:kind  . 3)
        (:range . ,(protocol.language-server.protocol:unparse-range location)))
       ,@(when (eq (architecture.builder-protocol:node-kind 'list node) :labeled)
           (let ((label (getf (architecture.builder-protocol:node-initargs 'list node) :label))
                 (result '()))
             (architecture.builder-protocol:walk-nodes
              'list
              (lambda (recurse relation relatin-args node* kind relations &rest initargs &key ((:label label*)) &allow-other-keys)
                (when (and (not (eq node node*))
                           (eq kind :labeled)
                           (eql label label*))
                  (log:info node*
                            (model.transform.trace:sources-for-target
                             (tracer document) node*))
                  (push `((:kind  . 1)
                          (:range . ,(protocol.language-server.protocol:unparse-range
                                      (first (model.transform.trace:sources-for-target
                                              (tracer document) node*)))))
                        result))
                (ignore-errors (funcall recurse)))
              (ast document))
             result)))))

(defmethod protocol.language-server.methods:hover ((workspace t)
                                                   (document  smarts-document)
                                                   (position  t))
  (multiple-value-bind (node location)
      (node-at-position document position)
    (when node
      `((:contents . ,(describe-node node))
        (:range    . ,(protocol.language-server.protocol:unparse-range location))))))

(defmethod protocol.language-server.methods:code-action ((workspace t)
                                                         (document  smarts-document)
                                                         (range     t)
                                                         (context   t))
  #())

(defclass smarts-workspace (protocol.language-server:workspace
                            protocol.language-server:document-container-mixin)
  ((root-uri  :initarg :root-uri) ; TODO client should not have to do this
   (root-path :initarg :root-path)))

(defmethod protocol.language-server:make-document ((container smarts-workspace)
                              (language  t)
                              (version   t)
                              (text      t))
  (make-instance (ecase language
                   (:smiles 'smiles-document)
                   (:smarts 'smarts-document))
                 :language language
                 :version  version
                 :text     text))

#.(progn
    #1=(ql:quickload :swank)
    '#1#)

(defun language-server (&rest args &aux (input *standard-input*) (output *standard-output*))
  (ignore-errors (delete-file "/tmp/swank-port.txt"))
  (swank:start-server "/tmp/swank-port.txt" :dont-close t)

  (let* ((connection (protocol.language-server.connection:make-connection
                      input output))
         (context    (make-instance 'protocol.language-server::context :connection      connection
                                                  :workspace-class 'smarts-workspace)))
    (trace "PROTOCOL.LANGUAGE-SERVER.PROTOCOL")
    (protocol.language-server::process-requests connection context)))
