
(in-package :smarts.graphviz)

(defvar *table*)
(defvar *rings*)

(defmethod label ((node chem:chem-info-node))
  (format nil "~s" node))

(defmethod label ((node chem:atom-test))
  (format nil "~s: ~a~%" node (chem:as-smarts node)))

(defmethod label ((node chem:logical))
  (with-output-to-string (sout)
    (format sout "~a<~d>" (chem:logical-operator node) (chem:get-id node))
    (when (not (eq (chem:get-ring-test node) :sarnone))
      (format sout " :ring-test ~a" (chem:get-ring-test node))
      (format sout " :ring-id ~a" (chem:get-ring-id node)))))

(defmethod label ((node chem:bond-to-atom-test))
  (format nil "(bond-to-atom-test ~a)~%" (chem:bond-type node)))

(defmethod style ((node chem:chem-info-node))
  (format nil "shape=box"))

(defmethod style ((node chem:bond-logical))
  (format nil "style=filled;color=\"yellow\""))

(defmethod style ((node chem:bond-test))
  (format nil "style=filled;color=\"yellowgreen\""))

(defmethod style ((node chem:atom-test))
  (format nil "shape=ellipse;style=filled;color=\"#A0FFA0\""))

(defmethod style ((node chem:chain))
  (format nil "shape=box;style=filled;color=\"#FFA0A0\""))

(defmethod style ((node chem:branch))
  (format nil "shape=box;style=filled;color=\"#A0A0FF\""))

(defmethod style ((node chem:bond-to-atom-test))
  (format nil "shape=box;style=filled;color=\"#F0A0FF\""))

(defmethod draw-node ((node null) stream)
  (let ((ident (setf (gethash (gensym) *table*) (gensym))))
    (format stream "  ~a [ shape = box, label=\"NIL\" ];~%" ident)
    ident))

(defmethod draw-node ((node chem:chem-info-node) stream)
  (if (null (gethash node *table*))
      (let ((ident (setf (gethash node *table*) (gensym))))
        (format stream "  ~a [ ~a; label=\"~a\" ];~%"
                ident (style node) (label node))
        (loop for child in (chem:chem-info-node-children node)
              for child-ident = (draw-node child stream)
              do (format stream "   ~a -> ~a~%" ident child-ident))
        ident)
      (gethash node *table*)))


(defmethod draw-node ((node chem:atom-or-bond-match-node) stream)
  (let ((node-label (call-next-method)))
    (cond
      ((eq (chem:get-ring-test node) :sarring-set)
       (setf (gethash (chem:get-ring-id node) *rings*) node-label))
      ((eq (chem:get-ring-test node) :sarring-test)
       (let ((ring-set (gethash (chem:get-ring-id node) *rings*)))
         (format stream "   ~a -> ~a [minlen=3,color=\"Red\",label=\"Ring ~a\"];" node-label ring-set (chem:get-ring-id node))))
      (t  ;;; Do nothing
       ))
    node-label))
      
(defun draw-graph (node stream)
  (format stream "digraph G {~%   ordering = out; ~%")
  (let ((*table* (make-hash-table :test #'eq))
        (*rings* (make-hash-table :test #'eql)))
    (draw-node node stream))
  (format stream "}~%"))

