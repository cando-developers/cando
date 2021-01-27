(in-package :cando-widgets)


(defparameter +graph-style/show-sketches+
  "* {
   label: data(label);
   font-size: 30;
  }
  edge {
    curve-style: bezier;
    arrow-scale: 5;
  }
  edge:selected {
   color: #1976d2;
  }
  edge.source-arrow {
    source-arrow-shape: triangle-backcurve;
  }
  edge.target-arrow {
    target-arrow-shape: triangle-backcurve;
  }
  node:unselected {
   background-color: white;
   border-width: 2px;
   border-style: solid;
   border-color: black;
  }
  node:selected {
   color: #1976d2;
   background-color: white;
   border-width: 2px;
   border-style: solid;
   border-color: #1976d2;
  }
  node {
   shape: rectangle;
   background-image: data(image);
   width: data(width);
   height: data(height);
  }")


(defparameter +graph-style/hide-sketches+
  "* {
   label: data(label);
   font-size: 10;
  }
  edge:selected {
   color: #1976d2;
  }
  node:selected {
   color: #1976d2;
  }
  edge {
    curve-style: bezier;
  }
  edge.source-arrow {
    source-arrow-shape: triangle-backcurve;
  }
  edge.target-arrow {
    target-arrow-shape: triangle-backcurve;
  }")


(defparameter *graph-layouts*
  '((:label "Circle layout"
     :args (cytoscape:circle-layout))
    (:label "Concentric layout"
     :args (cytoscape:concentric-layout))
    (:label "Breadth-First layout"
     :args (cytoscape:breadth-first-layout))
    (:label "CoSE layout"
     :args (cytoscape:cose-layout))
    (:label "Dagre layout"
     :args (cytoscape:dagre-layout))
    (:label "Grid layout"
     :args (cytoscape:grid-layout))))


(defclass molecule-map (jw:has-traits)
  ((add-edges-command
     :reader add-edges-command
     :initform (make-instance 'cytoscape:menu-command ; A command which adds edges
                              :content "<span class='fa fa-project-diagram fa-2x'></span>"))
   (delete-node-command
     :reader delete-node-command
     :initform (make-instance 'cytoscape:menu-command ; Delete node based on context.
                              :content "<span class='fa fa-trash fa-2x'></span>"))
   (delete-edge-command
     :reader delete-edge-command
     :initform
     (make-instance 'cytoscape:menu-command ; Delete edge based on context.
                    :content "<span class='fa fa-trash fa-2x'></span>"))
   (delete-elements-command
     :reader delete-elements-command
     :initform (make-instance 'cytoscape:menu-command ; Delete all selected elements.
                              :content "<span class='fa fa-trash fa-2x'></span>"))
   (graph
     :reader graph
     :initform (make-instance 'cytoscape:cytoscape-widget
                              :graph-style +graph-style/show-sketches+
                              :wheel-sensitivity 0.8 ; Detune the wheel sensitivity a bit.
                              :layout (make-instance 'jw:layout
                                                     :height "auto"
                                                     :align-self "stretch"
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :grid-area "molecule-map-content")))
   (layout-select
     :reader layout-select
     :initform (make-instance 'jw:dropdown ; A dropdown to select the layout
                              :%options-labels (mapcar (lambda (layout)
                                                         (getf layout :label)) *graph-layouts*)
                              :index 3
                              :layout (make-instance 'jw:layout :width "max-content" :margin ".5em")))
   (fit-button
     :reader fit-button
     :initform (make-instance 'jw:button ; A button to refit to the whole graph.
                              :description "Fit"
                              :layout (make-instance 'jw:layout :width "max-content" :margin ".5em")))
   (show-sketches-button
     :reader show-sketches-button
     :initform (make-instance 'jw:toggle-button ; A button to toggle the sketch visibility
                              :description "Show Sketches"
                              :value t
                              :layout (make-instance 'jw:layout :width "max-content" :margin ".5em")))
   (molecule-matcher
     :reader molecule-matcher
     :initarg :molecule-matcher
     :initform nil)
   (edge-constructor
     :reader edge-constructor
     :initarg :edge-constructor
     :initform nil)
   (molecules
     :accessor molecules
     :initarg :molecules
     :initform nil
     :trait t)
   (edges
     :accessor edges
     :initarg :edges
     :initform nil
     :trait t))
  (:metaclass jupyter-widgets:trait-metaclass))


(defun set-layout (map)
  (setf (cytoscape:graph-layouts (graph map))
        (list (apply #'make-instance (getf (elt *graph-layouts* (jw:widget-index (layout-select map))) :args)))))


(defun add-edges (map &rest ids)
  "Add edges between all nodes."
  (with-slots (edges edge-constructor)
              map
    (let ((names (mapcar (lambda (id)
                           (intern id 'keyword))
                         ids))
          new-edges)
      ;; Look over every pair of nodes, but only attempt to add an edge if the
      ;; ids are in normal order. This makes comparing edge ids easy, plus it avoids loops.
      (dolist (name1 names)
        (dolist (name2 names)
          (when (notany (lambda (edge) ; check to make sure that edge doesn't already exist.
                          (and (eq name1 (first edge))
                               (eq name2 (second edge))))
                        edges))
          (let ((edge (funcall edge-constructor name1 name2)))
            (when edge
              (push edge new-edges)))))
      ;; If new edges have been created then just add them. The layout will be updated automatically.
      (when new-edges
        (setf edges (append edges new-edges))))))


(defun is-edge-id (id)
  (and (position #\| id)
       t))


(defun delete-elements (map &rest ids)
  "Delete all nodes and edges associated with a set of ids."
  (with-slots (molecules edges)
              map
    (when (some #'is-edge-id ids)
      (setf edges
            (remove-if (lambda (edge)
                         (position (apply #'edge-id edge) ids :test #'equal))
                       edges)))
    (when (notevery #'is-edge-id ids)
      (setf molecules
            (remove-if (lambda (molecule)
                         (position (molecule-name molecule) ids :test #'equal))
                       molecules)))))


(defmethod initialize-instance :after ((instance molecule-map) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (with-slots (add-edges-command delete-elements-command delete-node-command delete-edge-command
               graph layout-select fit-button show-sketches-button edge-constructor)
              instance
    (when edge-constructor
      (setf (cytoscape:context-menus graph)
            (list (make-instance 'cytoscape:context-menu
                                 :selector "core"
                                 :commands (list add-edges-command delete-elements-command))
                  (make-instance 'cytoscape:context-menu
                                 :selector "node"
                                 :commands (list delete-node-command))
                  (make-instance 'cytoscape:context-menu
                                 :selector "edge"
                                 :commands (list delete-edge-command)))))
    ;; Delete the current node.
    (cytoscape:on-menu-command-select delete-node-command
      (lambda (inst id)
        (declare (ignore inst))
        (delete-elements instance id)))
    ;; Delete the current edge.
    (cytoscape:on-menu-command-select delete-edge-command
      (lambda (inst id)
        (declare (ignore inst))
        (delete-elements instance id)))
    ;; Add edges that connect all currently selected nodes.
    (cytoscape:on-menu-command-select add-edges-command
      (lambda (inst id)
        (declare (ignore inst id))
        (apply #'add-edges instance
                           (mapcan (lambda (element)
                                     (when (and (cytoscape:selected element)
                                                (not (cytoscape:removed element))
                                                (string= "nodes" (cytoscape:group element)))
                                       (list (cdr (assoc "id" (cytoscape:data element) :test #'string=)))))
                                   (cytoscape:elements (graph instance))))))
    ;; Delete all selected elements
    (cytoscape:on-menu-command-select delete-elements-command
      (lambda (inst id)
        (declare (ignore inst id))
        (apply #'delete-elements instance
                                 (mapcan (lambda (element)
                                           (when (and (cytoscape:selected element)
                                                      (not (cytoscape:removed element)))
                                             (list (cdr (assoc "id" (cytoscape:data element) :test #'string=)))))
                                         (cytoscape:elements (graph instance))))))
    (jw:observe layout-select :index
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value new-value source))
        (set-layout instance)))
    (jw:on-button-click fit-button
      (lambda (inst)
        (declare (ignore inst))
        (cytoscape:fit-elements graph)))
    (jw:observe show-sketches-button :value
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (setf (cytoscape:graph-style graph)
              (if new-value
                +graph-style/show-sketches+
                +graph-style/hide-sketches+))
        (cytoscape:layout graph)))
    (set-layout instance)
    (create-graph instance)
    ;; If molecules or edges change then update the graph.
    (jw:observe instance :molecules
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value new-value source))
        (update-graph instance)))
    (jw:observe instance :edges
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value new-value source))
        (update-graph instance)))))



(defun edge-id (molecule1 molecule2 &key id arrow1 arrow2 &allow-other-keys)
  "Construct an edge id with a normalized ordering for and undirected graph."
  (let* ((name1 (if (stringp molecule1) molecule1 (symbol-name molecule1)))
         (name2 (if (stringp molecule2) molecule2 (symbol-name molecule2)))
         (base-id (if (or (or (and arrow1 (not arrow2))
                              (and (not arrow1) arrow2))
                          (string< name1 name2))
                    (concatenate 'string name1 "|" name2)
                    (concatenate 'string name2 "|" name1))))
    (if id
      (concatenate 'string base-id "|" id)
      base-id)))


(defun encode-svg (svg)
  (concatenate 'string
               "data:image/svg+xml,"
               (quri:url-encode svg)))


(defun make-node (map molecule &key removed)
  "Construct a cytoscape node based on a molecule."
  (let ((name (molecule-name molecule)))
    (multiple-value-bind (svg width height)
                         (sketch-molecule molecule)
      ;; Save sketch dimension and the sketch itself into the data slot. The sketch is SVG text
      ;; so we just need do a percent encoded string for the reserved URL characters.
      (make-instance 'cytoscape:element
                     :group "nodes"
                     :removed removed
                     :on-trait-change (list (cons :selected
                                                  (lambda (instance type nm old-value new-value source)
                                                    (declare (ignore instance type nm old-value source))
                                                    (update-selected-sketches map (unless new-value (list name))))))
                     :data (list (cons "id" name)
                                 (cons "label" name)
                                 (cons "width" width)
                                 (cons "height" height)
                                 (cons "image" (encode-svg svg)))))))

(defun update-node (map element &optional molecule)
  (let ((name (molecule-name molecule)))
    (set-element-removed element (not molecule))
    (when molecule
      (multiple-value-bind (svg width height)
                           (sketch-molecule molecule)
        (setf (cytoscape:data element)
              (list (cons "id" name)
                    (cons "label" name)
                    (cons "width" width)
                    (cons "height" height)
                    (cons "image" (encode-svg svg))))))))


(defun make-edge (map molecule1 molecule2 &key arrow1 arrow2 label id &allow-other-keys)
  "Construct an cytoscape edge based on an edge definition when is a list of
  the form (source target label)."
  (let ((name-1 (symbol-name molecule1))
        (name-2 (symbol-name molecule2)))
    (make-instance 'cytoscape:element
                   :group "edges"
                   :classes (cond
                              ((and arrow1 arrow2)
                                (list "source-arrow" "target-arrow"))
                              (arrow1
                                (list "source-arrow"))
                              (arrow2
                                (list "target-arrow")))
                   :on-trait-change (list (cons :selected
                                                (lambda (instance type name old-value new-value source)
                                                  (declare (ignore instance type name old-value source))
                                                  (update-selected-sketches map (unless new-value
                                                                                  (list name-1 name-2))))))
                   :data (list (cons "id" (edge-id name-1 name-2
                                                   :id id
                                                   :arrow1 arrow1
                                                   :arrow2 arrow2))
                               (cons "label" (or label ""))
                               (cons "source" name-1)
                               (cons "target" name-2)))))


(defun create-graph (map)
  (setf (cytoscape:elements (graph map))
        (nconc (mapcar (lambda (molecule)
                         (make-node map molecule))
                       (molecules map))
               (mapcar (lambda (edge)
                         (apply #'make-edge map edge))
                       (edges map)))))


(defun set-element-removed (element removed)
  "Show or hide a cytoscape element and ensure that the element is not currently selected or
  grabbed."
  (when (or (and removed (not (cytoscape:removed element)))
            (and (not removed) (cytoscape:removed element)))
    (setf (cytoscape:selected element) nil
          (cytoscape:grabbed element) nil
          (cytoscape:removed element) removed)))


(defun element-id (element)
  (cdr (assoc "id" (cytoscape:data element) :test #'string=)))


(defun get-selected-names (map)
  (mapcan (lambda (element)
            (when (cytoscape:selected element)
              (let* ((id (element-id element))
                     (pos (position #\| id))) ; If there is vertical bar it is an edge.
                (if pos
                  (list (subseq id 0 pos) (subseq id (1+ pos)))
                  (list id)))))
          (cytoscape:elements (graph map))))


(defun generate-stylesheet (atoms)
  (when atoms
    (with-output-to-string (output-stream)
      (format output-stream "~{.~A~^,~%~} {~%  filter: url(#highlight);~%}" atoms))))


(defun update-selected-sketches (map unselected-names)
  (with-slots (molecule-matcher graph molecules)
              map
    (when molecule-matcher
      (let* ((selected-names (get-selected-names map))
             (selected (remove-if (lambda (molecule)
                                    (not (position (molecule-name molecule) selected-names :test #'equal)))
                                  molecules))
             (atom-matches (funcall molecule-matcher selected)))
        (dolist (element (cytoscape:elements graph))
          (let ((molecule (find (element-id element) molecules :key #'molecule-name :test #'equal)))
            (when (and molecule
                       (position (element-id element) unselected-names :test #'equal))
              (setf (cdr (assoc "image" (cytoscape:data element) :test #'string=))
                    (encode-svg (cando-widgets:sketch-molecule (find (element-id element) molecules :key #'molecule-name :test #'equal)))))
            (jupyter-widgets:notify-trait-change element :dict :data (cytoscape:data element) (cytoscape:data element) t)))
        (loop for molecule in selected
              for atoms in atom-matches
              for element = (find (molecule-name molecule) (cytoscape:elements graph) :test #'string= :key #'element-id)
              do (setf (cdr (assoc "image" (cytoscape:data element) :test #'string=))
                       (encode-svg (cando-widgets:sketch-molecule molecule (generate-stylesheet atoms))))
              do (jupyter-widgets:notify-trait-change element :dict :data (cytoscape:data element) (cytoscape:data element) t))))))


(defun make-molecule-map (container &rest args)
  (let ((instance (apply #'make-instance 'molecule-map args)))
    (setf (jw:widget-children container)
          (append (jw:widget-children container)
                  (list (graph instance)
                        (make-instance 'jw:box
                                       :children (list (layout-select instance)
                                                       (fit-button instance)
                                                       (show-sketches-button instance))
                                       :layout (make-instance 'jw:layout
                                                              :flex-flow "row wrap"
                                                              :justify-content "center"
                                                              :align-items "baseline"
                                                              :align-content "flex-start"
                                                              :grid-area "molecule-map-controls")))))
    instance))


(defun update-graph (map)
  (with-slots (graph molecules edges)
              map
    (let (new-elements)
      ;; Look for elements that don't have a corresponding match in molecules or edges.
      (dolist (element (cytoscape:elements graph))
        (let ((id (cdr (assoc "id" (cytoscape:data element) :test #'equal))))
                               (if (equal "nodes" (cytoscape:group element))
                                 (update-node map element
                                              (find id molecules ; the element is a node so look in molecules
                                                    :test #'string=
                                                    :key #'molecule-name))
                                 (set-element-removed element
                                                      (notany (lambda (edge) ; the element is an edge so look in edges.
                                                                (equal id (apply #'edge-id edge)))
                                                              edges)))))
      ;; Look through edges for new edges
      (dolist (edge edges)
        (unless (some (lambda (element)
                        (equal (cdr (assoc "id" (cytoscape:data element) :test #'equal))
                               (apply #'edge-id edge)))
                      (cytoscape:elements graph))
          (push (apply #'make-edge map edge) new-elements)))
      ;; Look through molecules for new nodes
      (dolist (molecule molecules)
        (unless (some (lambda (element)
                        (equal (cdr (assoc "id" (cytoscape:data element) :test #'equal))
                               (molecule-name molecule)))
                      (cytoscape:elements graph))
          (push (make-node map molecule) new-elements)))
      ;; If new elements have been created then just add them. The layout will be updated automatically.
      ;; Otherwise force the layout to update.
      (if new-elements
        (setf (cytoscape:elements graph) (append (cytoscape:elements graph) new-elements))
        (cytoscape:layout graph)))))


(defun molecule-name (molecule)
  "Get a string representation of the structure."
  (symbol-name (chem:get-name molecule)))

