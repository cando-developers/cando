(in-package :tirun-jupyter)


(defmacro def-dyn-widget (name initial-value)
  (let ((var-name (gensym)))
    `(progn
       (defparameter ,var-name nil)
       (defun ,name ()
         (or ,var-name
             (setf ,var-name ,initial-value))))))


(def-dyn-widget button-style (make-instance 'w:button-style :button-color "aquamarine"))


(defclass app (jupyter-widgets:has-traits)
  ((composer-json
    :accessor composer-json
    :initarg :composer-json
    :initform nil
    :trait t)
   (receptor-string
    :accessor receptor-string
    :initform nil
    :trait t)
   (ligands-string
    :accessor ligands-string
    :initform nil
    :trait t)
   (all-ligands
    :accessor all-ligands
    :initform nil
    :trait t)
   (selected-ligands
    :accessor selected-ligands
    :initform nil
    :trait t)
   (all-edges
    :accessor all-edges
    :initform nil
    :trait t)
   (smirks
    :accessor smirks
    :initform ""
    :trait t)
   (distributor
    :accessor distributor
    :initform "s103.thirdlaw.tech"
    :trait t)
   (job-name
    :accessor job-name
    :initform "default"
    :trait t)
   (submit-stream
    :accessor submit-stream
    :initform nil
    :trait t)
   (smirks-pattern
    :accessor smirks-pattern
    :initform nil
    :trait t))
  (:metaclass jupyter-widgets:trait-metaclass))


(define-symbol-macro all-ligands
  (all-ligands *app*))


(define-symbol-macro selected-ligands
  (selected-ligands *app*))


(define-symbol-macro all-edges
  (all-edges *app*))


(defun composer-save (composer app &key (tirun *tirun*))
  (format t "In composer-save - about to parse get-decoded-time ~s~%" (get-decoded-time))
  ;; This is where we put the logic to build the molecules?
  (let ((agg (handler-bind
                 ((error (lambda (err)
                          (format t "Hit error ~a" err)
                          (error err))))
               (structure-editor:parse-kekule-json (composer-json app)))))
    (let ((ligand-molecules (tirun::assemble-ligands *tirun* agg :verbose t)))
      (setf all-ligands ligand-molecules)
      (finish-output)
      )))

(defun structure-editor:set-composer-json (composer app json)
  (setf (composer-json app) json)
  (format t "In set-composer-json~%")
  (composer-save composer app)
  )

(defun structure-editor:get-composer-json (app)
  (composer-json app))


(defclass save-load-app ()
  ((composer-json :initarg :composer-json :accessor composer-json)))

(cando:make-class-save-load save-load-app)

   
(defvar *app* nil)


(def-dyn-widget box-layout
  (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
         
(defvar *tirun* nil)

(defun new-tirun ()
  (let ((dest-dir (let ((dd (ext:getenv "CANDO_JOBS_DIRECTORY")))
                    (if dd
                        dd
                        (pathname "~/jobs/")))))
    (format t "Jobs will be saved to ~a~%setenv CANDO_JOBS_DIRECTORY if you want it to go elsewhere~%" (namestring dest-dir))
    (ensure-directories-exist dest-dir))
  (with-output-to-string (sout)
    (let ((*standard-output* sout))
      (cando-user:setup-default-paths)
      (cando-user:load-atom-type-rules "ATOMTYPE_GFF.DEF")
      (cando-user:source "leaprc.ff14SB.redq")
      (cando-user:source "leaprc.gaff")
      (leap:add-pdb-res-map '((1 :NMA :NME)))))
  (setf *app* (make-instance 'app))
  (setf *tirun* (make-instance 'tirun:tirun-calculation)))

(defun load-app (&optional (name "default.dat"))
  (let ((pn (pathname name)))
    (if (probe-file pn)
        (let ((sl-app (cando:load-cando pn)))
          (setf *app* (make-instance 'app
                                     :composer-json (composer-json sl-app))))
        (setf *app* (make-instance 'app))))
  
  t)

(defun save-app (&optional (name "default.dat"))
  (let ((pn (pathname name)))
    (ensure-directories-exist pn)
    (let ((sl (make-instance 'save-load-app
                             :composer-json (composer-json *app*)
                             )))
      (cando:save-cando sl (pathname name)))))


(defun setup-tirun (name)
  (let ((job-name (intern name :keyword)))
    (setf (job-name *app*) job-name))
  )





#+(or)(defun new-pas ()
  (with-output-to-string (sout)
    (let ((*standard-output* sout))
      (cando-user:source "leaprc.lipid14")
      (cando-user:source "leaprc.water.tip3p")))
    (asdf:load-asd "/Users/yonkunas/Development/pas/pas.asd")
    (asdf:load-system :pas))
    

(defclass simple-input ()
  ((name :initarg :name :accessor name)
   (label-widget :initarg :label-widget :accessor label-widget)
   (input-widget :initarg :input-widget :accessor input-widget)
   (hbox :initarg :hbox :accessor hbox)))

(defun make-simple-input (label &key (default "") name)
  (let ((label-widget (make-instance 'w:label :value label
                                              :layout (make-instance 'w:layout :width "25em")))
        (input-widget (make-instance 'w:text
                                     :value default
                                     :layout (make-instance 'w:layout :width "20em"))))
    (make-instance 'simple-input
                   :name name
                   :label-widget label-widget
                   :input-widget input-widget
                   :hbox (make-instance 'w:h-box
                                        :children (list label-widget input-widget)))))


(defun parse-receptor (parameter progress-callback)
  "Called from file task page to parse receptor."
  (declare (ignore progress-callback))
  (handler-case
      (with-slots (receptor-string)
                  *app*
        (setf receptor-string (babel:octets-to-string (cdr (assoc "content" (car parameter) :test #'string=)))
              (tirun:receptors *tirun*)
              (list (with-input-from-string (input-stream receptor-string)
                      (leap.pdb:load-pdb-stream input-stream))))
        t)
    (leap.pdb:pdb-read-error (condition)
      (princ condition *error-output*)
      nil)
    (end-of-file (condition)
      (princ condition *error-output*)
      nil)))


(defun load-receptor ()
  "Present a graphical interface to load and view a receptor into the current application."
  (let* ((container (make-instance 'w:accordion :selected-index 0))
         (ngl (make-instance 'nglv:nglwidget
                              :layout (make-instance 'w:layout
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :grid-area "ngl"))))
    (cw:make-file-task-page container "Load Receptor" #'parse-receptor :accept ".pdb")
    ;; Make a simple page to display the receptor in nglview
    (setf (w:widget-%titles container) (append (w:widget-%titles container)
                                              (list "View Receptor"))
          (w:widget-children container) (append (w:widget-children container)
                                               (list ngl)))
    ;; When the view page opens give nglview a chance to resize.
    (w:observe container :selected-index
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (when (equal new-value 1)
          (nglview:handle-resize ngl))))
    ;; Add the receptor if one is already defined.
    (when (receptor-string *app*)
      (nglview:add-structure ngl (make-instance 'nglview:text-structure :text (receptor-string *app*))))
    ;; Listen for changes to the receptor and add or delete when needed.
    (w:observe *app* :receptor-string
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (nglview:remove-all-components ngl)
        (when new-value
          (nglview:add-structure ngl (make-instance 'nglview:text-structure :text new-value)))))
    container))


(defclass view-ligand-page (cw:page)
  ((ngl
     :reader view-ligand-ngl
     :initform (make-instance 'nglv:nglwidget
                              :layout (make-instance 'w:layout
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :grid-area "ngl"))
     :documentation "nglview for the selected ligand and the receptor.")
   (structure
     :reader view-ligand-structure
     :initform (make-instance 'w:html
                              :layout (make-instance 'w:layout
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :grid-area "structure"))
     :documentation "A sketch view for the selected ligand.")
   (slider
     :reader view-ligand-slider
     :initform (make-instance 'w:int-slider
                              :layout (make-instance 'w:layout
                                                     :width "100%"))
     :documentation "A slider to select a ligand.")
   (dropdown
     :reader view-ligand-dropdown
     :initform (make-instance 'w:dropdown
                              :layout (make-instance 'w:layout
                                                     :width "max-content"))
     :documentation "A dropdown of all the available ligands.")
   (receptor-toggle
     :reader view-ligand-receptor-toggle
     :initform (make-instance 'w:toggle-button
                              :description "Show Receptor"
                              :value t
                              :layout (make-instance 'w:layout
                                                     :margin "auto"
                                                     :grid-area "ngl-ctl"))
     :documentation "A toggle button to control the visibility of the receptor."))
  (:metaclass jupyter-widgets:trait-metaclass)
  (:documentation "A page to view ligands that has a sketch and an nglview widget along with controls and ligand selector.")
  (:default-initargs
    :layout (make-instance 'w:layout
                           :grid-gap "1em"
                           :grid-template-rows "1fr min-content"
                           :grid-template-columns "1fr 1fr"
                           :grid-template-areas "\"structure ngl\" \"structure-ctl ngl-ctl\"")))


(defvar *smallest-ligand-sketch*)


(defun molecule-name (molecule)
  "Get a string representation of the structure."
  (symbol-name (chem:get-name molecule)))


(defun on-ligand-select (instance index &optional previous)
  "Select a ligand on a view-ligand-page and hide a previous ligand if one was selected."
  (with-slots (ngl structure)
              instance
    (let* ((mol (elt all-ligands index))
           (id (molecule-name mol)))
      ;; Load the sketch
      (setf (w:widget-value structure)
            (format nil "<div style='display:flex;align-items:center;height:100%;'><div style='display:block;margin:auto;'>~A</div></div>"
                    (sketch2d:render-svg-to-string (cw:sketch-molecule mol))))
      ;; A resize can't hurt.
      (nglview:handle-resize ngl)
      ;; Hide previous ligand.
      (when previous
        (nglview:hide-components ngl (molecule-name (elt all-ligands previous))))
      ;; Show the ligand if it is already in nglview otherwise load it.
      (cond
        ((position id (nglview:component-ids ngl) :test #'string=)
          (nglview:show-components ngl id)
          (nglview:center ngl :component id))
        (t
          (let ((agg (chem:make-aggregate)))
            (chem:add-matter agg mol)
            (nglview:add-structure ngl
                                   (make-instance 'nglview:text-structure
                                                  :id id :ext "mol2"
                                                  :text (chem:aggregate-as-mol2-string agg t)))))))))


(defmethod initialize-instance :after ((instance view-ligand-page) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  ;; Link the slider and the dropdown states together
  (jupyter-widgets:link (view-ligand-slider instance) :value
                        (view-ligand-dropdown instance) :index)
  (setf (w:widget-children instance)
        (list (view-ligand-ngl instance)
              (view-ligand-structure instance)
              (make-instance 'w:h-box
                             :children (list (view-ligand-slider instance)
                                             (view-ligand-dropdown instance))
                             :layout (make-instance 'w:layout
                                                    :width "90%" :margin "auto"
                                                    :grid-area "structure-ctl"))
              (view-ligand-receptor-toggle instance)))
  ;; Set the receptor visibility based the toggle.
  (w:observe (view-ligand-receptor-toggle instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (nglview:set-visibility (view-ligand-ngl instance) new-value "receptor")))
  ;; When the slider index changes update the sketch and nglview.
  (w:observe (view-ligand-slider instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (on-ligand-select instance new-value old-value)))
  ;; When the view ligand page opens up make sure we give nglview an opportunity to resize.
  (w:observe (cw:container instance) :selected-index
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (nglview:handle-resize (view-ligand-ngl instance))))
  ;; If the receptor-string is already set then load the current receptor.
  (when (receptor-string *app*)
    (nglview:add-structure (view-ligand-ngl instance)
                           (make-instance 'nglview:text-structure
                                          :id "receptor"
                                          :text (receptor-string *app*))))
  ;; When the receptor-string changes reload the receptor in nglview.
  (w:observe *app* :receptor-string
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name source))
      (when old-value
        (nglview:remove-components (view-ligand-ngl instance) "receptor"))
      (nglview:add-structure (view-ligand-ngl instance)
                             (make-instance 'nglview:text-structure
                                            :id "receptor"
                                            :text new-value))))
  ;; If the all-ligands changes then reload the whole thing.
  (w:observe *app* :all-ligands
    (lambda (inst type name old-value ligands source)
      (declare (ignore inst type name old-value source))
      (with-slots (slider dropdown)
                  instance
        (setf (w:widget-max slider) (1- (length ligands))
              (w:widget-%options-labels dropdown) (mapcar #'molecule-name ligands)
              (w:widget-value slider) 0)
        (when ligands
          (cw:sketch-molecules ligands)
          (on-ligand-select instance 0)
          (ignore-errors
            (tirun:tirun-calculation-from-ligands *tirun* ligands)))))))


(defun make-view-ligand-page (container title)
  "Create an instance of a view-ligand-page and add it to the container."
  (let ((page (make-instance 'view-ligand-page :container container)))
    (cw:add-page container page title)
    (values)))


(defun parse-ligands (parameter progress-callback)
  "Called from the file task page to actually parse the ligands file."
  (declare (ignore progress-callback))
  (handler-case
      (with-slots (ligands-string)
                  *app*
        (setf ligands-string (babel:octets-to-string (cdr (assoc "content" (car parameter) :test #'string=)))
              all-ligands (with-input-from-string (input-stream ligands-string)
                               (let ((ext (pathname-type (cdr (assoc "name" (car parameter) :test #'string=)))))
                                 (cond
                                   ((string-equal ext "sdf")
                                     (sdf:parse-sdf-file input-stream))
                                   ((string-equal ext "mol2")
                                     (chem:read-mol2-list input-stream))
                                   (t
                                     (error "Unknown file type ~A." ext)))))
              selected-ligands all-ligands)
        t)
    (esrap:esrap-parse-error (condition)
      (princ condition *error-output*)
      (fresh-line *error-output*)
      nil)
    (simple-error (condition) ; This is temporary until we have mol2-parse-error
      (princ condition *error-output*)
      (fresh-line *error-output*)
      nil)
    (sdf:sdf-parse-error (condition)
      (princ condition *error-output*)
      (fresh-line *error-output*)
      nil)
    (end-of-file (condition)
      (princ condition *error-output*)
      (fresh-line *error-output*)
      nil)))


(defun load-ligands ()
  "Present a graphical interface to load a collection of ligands from a file into the current
  application and then display the ligand structure using two-dimension diagrams and three
  dimension representation."
  (let ((container (make-instance 'w:accordion :selected-index 0)))
    (cw:make-file-task-page container "Load Ligands" #'parse-ligands :accept ".mol2,.sdf")
    (make-view-ligand-page container "View Ligands")
    container))


(defun select-ligands ()
  "select-ligands in the current application and displays a molecule selection list that allows one
  to select a subset of the molecules to use in further calculations."
  (let ((sel (cw:make-molecule-select :molecules all-ligands
                                      :selected selected-ligands)))
    (w:link sel :molecules *app* :all-ligands)
    (w:link sel :selected *app* :selected-ligands)
    (make-instance 'w:accordion
                   :selected-index 0
                   :%titles (list "Select Ligands")
                   :children (list sel))))


(defun build-prototype-sketch (ligands)
  "Find the smallest ligand and build a sketch for it"
  (let ((min-atoms (chem:number-of-atoms (car ligands)))
        (min-mol (car ligands)))
    (loop for mol in (cdr ligands)
          for num-atoms = (chem:number-of-atoms mol)
          when (< num-atoms min-atoms)
            do (setf min-atoms num-atoms
                     min-mol mol))
    (sketch2d:sketch2d min-mol)))


(defun run-pas (parameter progress-callback)
  "Called by the simple-task-page to do the PAS calculation."
  (declare (ignore parameter))
  (write-line "Running position analogue scanning...")
  (finish-output)
  ;; Loop over the current selected ligands and run PAS on each one accumulating the results.
  (setf all-ligands
        (loop for template-ligand in selected-ligands
              for new-ligands = (pas:pas (smirks-pattern *app*) template-ligand)
              append new-ligands
              do (format t "Using template ligand of ~a built the ligands ~{~a~#[~;, and ~:;, ~]~}...~%"
                         (molecule-name template-ligand)
                         (mapcar #'molecule-name new-ligands))
              do (finish-output)
              do (funcall progress-callback (length selected-ligands)))
        selected-ligands all-ligands))


(defvar +smirks-patterns+
  '(("[D3&CH1:1]-[#1:2]>>[N:1]" . "Replace aromatic methine carbons with nitrogen.")))


(defclass pas-config-page (cw:page)
  ((help
     :accessor help
     :initform (make-instance 'w:html
                              :layout (make-instance 'w:layout :grid-area "help"))
     :documentation "A help page listing SMIRKS examples.")
   (snippets
     :accessor snippets
     :initform (make-instance 'w:select
                              :%options-labels (mapcar #'car +smirks-patterns+)
                              :description "SMIRKS snippets"
                              :style (make-instance 'w:description-style :description-width "12em")
                              :layout (make-instance 'w:layout :width "100%" :grid-area "snippets")))
   (smirks-pattern-editor
     :accessor smirks-pattern-editor
     :initform (make-instance 'w:text
                              :description "SMIRKS pattern"
                              :value (smirks-pattern *app*)
                              :style (make-instance 'w:description-style :description-width "12em")
                              :layout (make-instance 'w:layout :width "100%" :grid-area "smirks"))
     :documentation "The SMIRKS pattern editor"))
  (:metaclass jupyter-widgets:trait-metaclass)
  (:default-initargs
    :layout (make-instance 'w:layout
                           :width "100%"
                           :max-height "12em"
                           :grid-gap "0.1em 1em"
                           :grid-template-rows "min-content 1fr"
                           :grid-template-columns "1fr 1fr"
                           :grid-template-areas "\"smirks help\" \"snippets help\"")))


(defmethod initialize-instance :after ((instance pas-config-page) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (with-slots (snippets smirks-pattern-editor help)
              instance
    ;; Connect the SMIRKS editor to *app* smirks-pattern
    (w:link *app* :smirks-pattern smirks-pattern-editor :value)
    (setf (w:widget-children instance) (list help smirks-pattern-editor snippets))
    (w:observe snippets :index
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (destructuring-bind (snippet . text)
                            (nth new-value +smirks-patterns+)
          (setf (w:widget-value (smirks-pattern-editor instance)) snippet
                (w:widget-value (help instance)) text))))))


(defun make-pas-config-page (container title)
  "Create a pas-config-page instance and add it to the container."
  (let* ((page (make-instance 'pas-config-page :container container)))
    (cw:add-page container page title)
    (values)))


(defun pasrun ()
  "pasrun loads a ligand template and runs PAS on it to generate multiple ligands.
It will put those multiple ligands into all-ligands and selected-ligands"
  (let ((container (make-instance 'w:accordion :selected-index 0)))
    (cw:make-file-task-page container "Load Template Ligands" #'parse-ligands :accept ".mol2,.sdf")
    (make-pas-config-page container "Configure PAS")
    (cw:make-simple-task-page container "Calculate PAS" #'run-pas
                              :label "Click button to calculate PAS.")
    (make-view-ligand-page container "View Ligands")
    container))


(defun edge-id (name-1 name-2)
  "Construct an edge id with a normalized ordering for and undirected graph."
  (if (string< name-1 name-2)
    (concatenate 'string name-1 "|" name-2)
    (concatenate 'string name-2 "|" name-1)))


(defun make-node (ligand)
  "Construct a cytoscape node based on a ligand."
  (let ((name (molecule-name ligand))
        (sketch (cw:sketch-molecule ligand)))
    ;; Save sketch dimension and the sketch itself into the data slot. The sketch is SVG text
    ;; so we just need do a percent encoded string for the reserved URL characters.
    (make-instance 'cytoscape:element
                   :group "nodes"
                   :removed (not (position ligand selected-ligands))
                   :data (list (cons "id" name)
                               (cons "label" name)
                               (cons "width" (sketch2d:width sketch))
                               (cons "height" (sketch2d:height sketch))
                               (cons "image" (concatenate 'string
                                                          "data:image/svg+xml,"
                                                          (quri:url-encode (sketch2d:render-svg-to-string sketch))))))))


(defun make-edge (edge)
  "Construct an cytoscape edge based on an edge definition when is a list of
  the form (source target label)."
  (let ((name-1 (symbol-name (first edge)))
        (name-2 (symbol-name (second edge))))
    (make-instance 'cytoscape:element
                   :group "edges"
                   :data (list (cons "id" (edge-id name-1 name-2))
                               (cons "label" (or (third edge) :null))
                               (cons "source" name-1)
                               (cons "target" name-2)))))


(defun create-graph (graph)
  "Create a graph from all-ligands and all-edges."
  (setf (cytoscape:elements graph)
        (nconc (mapcar #'make-node all-ligands)
               (mapcar #'make-edge all-edges))))


(defun set-element-removed (element removed)
  "Show or hide a cytoscape element and ensure that the element is not currently selected or
  grabbed."
  (when (or (and removed (not (cytoscape:removed element)))
            (and (not removed) (cytoscape:removed element)))
    (setf (cytoscape:selected element) nil
          (cytoscape:grabbed element) nil
          (cytoscape:removed element) removed)))


(defun update-graph (graph)
  "Update a graph in reponse to a change in selected-ligands or all-edges"
  (let (new-edges)
    ;; Look for elements that don't have a corresponding match in selected-ligands or all-edges.
    (dolist (element (cytoscape:elements graph))
      (let ((id (cdr (assoc "id" (cytoscape:data element) :test #'equal))))
        (set-element-removed element
                             (if (equal "nodes" (cytoscape:group element))
                               (not (position id selected-ligands ; the element is node so look in selected-ligands.
                                              :test #'string=
                                              :key #'molecule-name))
                               (notany (lambda (edge) ; the element is an edge so look in all-edges.
                                         (equal id (edge-id (symbol-name (first edge)) ; the ids are already in normal order.
                                                            (symbol-name (second edge)))))
                                       all-edges)))))
    ;; Look through all-edges for new edges
    (dolist (edge all-edges)
      (unless (some (lambda (element)
                      (equal (cdr (assoc "id" (cytoscape:data element) :test #'equal))
                                    (edge-id (symbol-name (first edge))
                                             (symbol-name (second edge)))))
                    (cytoscape:elements graph))
        (push (make-edge edge) new-edges)))
    ;; If new edges have been created then just add them. The layout will be updated automatically.
    ;; Otherwise force the layout to update.
    (if new-edges
      (setf (cytoscape:elements graph) (append (cytoscape:elements graph) new-edges))
      (cytoscape:layout graph))))


(defun run-lomap (parameter progress-callback)
  "This function is called by the task page and updates a progress widget as lomap does the calculation."
  (declare (ignore parameter))
  (format t "Calculating similarity matrix of ~{~a~#[~;, and ~:;, ~]~}...~%"
          (mapcar #'molecule-name selected-ligands))
  (finish-output) ; Make sure the message is synced to the frontend.
  (let ((multigraph (lomap::lomap-multigraph
                      (lomap::similarity-multigraph
                        selected-ligands
                        (lomap::similarity-matrix selected-ligands
                                                  :advance-progress-callback progress-callback))))
        new-edges)
    ;; Loop over the subgraphs and the vertices and edges and create edges.
    (dolist (subgraph (lomap::subgraphs multigraph))
      (dolist (edge (lomap::edges subgraph))
        (push (nconc (sort (list (chem:get-name (lomap:molecule (lomap:vertex1 edge)))
                                 (chem:get-name (lomap:molecule (lomap:vertex2 edge))))
                           #'string<)
                     (list (format nil "S = ~3,2f" (lomap:sim-score edge))))
              new-edges)))
    (setf all-edges new-edges))
  (write-line "Similarity matrix calculation completed."))


(defparameter +similarity-graph-style/show-sketches+
  "* {
   label: data(label);
   font-size: 30;
  }
  edge:selected {
   color: #1976d2;
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


(defparameter +similarity-graph-style/hide-sketches+
  "* {
   label: data(label);
   font-size: 10;
  }
  edge:selected {
   color: #1976d2;
  }
  node:selected {
   color: #1976d2;
  }")


(defun add-edges (graph &rest ids)
  "Add edges between all nodes."
  (let ((names (mapcar (lambda (id)
                         (intern id 'keyword))
                       ids))
        new-edges)
    ;; Look over every pair of nodes, but only attempt to add an edge if the
    ;; ids are in normal order. This makes comparing edge ids easy, plus it avoids loops.
    (dolist (name-1 names)
      (dolist (name-2 names)
        (when (and (string< name-1 name-2)
                   (notany (lambda (edge) ; check to make sure that edge doesn't already exist.
                             (and (eq name-1 (first edge))
                                  (eq name-2 (second edge))))
                           all-edges))
          (push (list name-1 name-2) new-edges))))
    ;; If new edges have been created then just add them. The layout will be updated automatically.
    ;; Otherwise force the layout to update.
    (when new-edges
      (setf all-edges (append all-edges new-edges))
      (cytoscape:layout graph))))


(defun delete-elements (graph &rest ids)
  "Delete all nodes and edges associated with a set of ids."
  (let* ((names (mapcar (lambda (id) ; extract the names.
                          (let ((pos (position #\| id))) ; If there is vertical bar it is an edge.
                            (if pos
                              (list (intern (subseq id 0 pos) 'keyword)
                                    (intern (subseq id (1+ pos)) 'keyword))
                              (intern id 'keyword))))
                        ids))
         (new-selected (remove-if (lambda (ligand) ; ligands that do not appear in names
                                    (position (chem:get-name ligand) names))
                                  selected-ligands))
         (new-edges (remove-if (lambda (edge) ; edges that do not appear in or do not refer to names.
                                 (some (lambda (name)
                                         (or (and (listp name)
                                                  (eq (first edge) (first name))
                                                  (eq (second edge) (second name)))
                                             (and (symbolp name)
                                                  (or (eq (first edge) name)
                                                      (eq (second edge) name)))))
                                       names))
                               all-edges)))
    ;; Remove edges first.
    (unless (= (length all-edges)
               (length new-edges))
      (setf all-edges new-edges))
    (unless (= (length selected-ligands)
               (length new-selected))
      (setf selected-ligands new-selected))
    ;; If any update has happened then relayout the graph.
    (when (or (/= (length all-edges)
                  (length new-edges))
              (/= (length selected-ligands)
                  (length new-selected)))
      (cytoscape:layout graph))))


(defun lomap ()
  "Calculate molecule similarities using a graph theory approach based on LOMAP
  (J Comput Aided Mol Des 27, 755–770 (2013). https://doi.org/10.1007/s10822-013-9678-y)"
  (let* ((container (make-instance 'w:accordion :selected-index 0))
         (add-edges-command (make-instance 'cytoscape:menu-command ; A command which adds edges
                                          :content "<span class='fa fa-project-diagram fa-2x'></span>"))
         (delete-node-command (make-instance 'cytoscape:menu-command ; Delete node based on context.
                                             :content "<span class='fa fa-trash fa-2x'></span>"))
         (delete-edge-command (make-instance 'cytoscape:menu-command ; Delete edge based on context.
                                             :content "<span class='fa fa-trash fa-2x'></span>"))
         (delete-elements-command (make-instance 'cytoscape:menu-command ; Delete all selected elements.
                                                 :content "<span class='fa fa-trash fa-2x'></span>"))
         (graph (make-instance 'cytoscape:cytoscape-widget
                               :graph-layouts (list (make-instance 'cytoscape:cose-layout))
                               :context-menus (list (make-instance 'cytoscape:context-menu
                                                                   :selector "core"
                                                                   :commands (list add-edges-command delete-elements-command))
                                                    (make-instance 'cytoscape:context-menu
                                                                   :selector "node"
                                                                   :commands (list delete-node-command))
                                                    (make-instance 'cytoscape:context-menu
                                                                   :selector "edge"
                                                                   :commands (list delete-edge-command)))
                               :graph-style +similarity-graph-style/show-sketches+
                               :wheel-sensitivity 0.8 ; Detune the wheel sensitivity a bit.
                               :layout (make-instance 'w:layout
                                                      :height "640px"
                                                      :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                      :grid-area "graph")))
         (fit (make-instance 'w:button ; A button to refit to the whole graph.
                             :description "Fit"
                             :on-click (list (lambda (inst)
                                               (declare (ignore inst))
                                               (cytoscape:fit-elements graph)))
                             :layout (make-instance 'w:layout
                                                    :grid-area "fit")))
         (show-sketches (make-instance 'w:toggle-button ; A button to toggle the sketch visibility
                                       :description "Show Sketches"
                                       :value t
                                       :layout (make-instance 'w:layout
                                                              :grid-area "show-sketches")
                                       :on-trait-change (list (cons :value
                                                                    (lambda (inst type name old-value new-value source)
                                                                      (declare (ignore inst type name old-value source))
                                                                      (setf (cytoscape:graph-style graph)
                                                                            (if new-value
                                                                              +similarity-graph-style/show-sketches+
                                                                              +similarity-graph-style/hide-sketches+))
                                                                      (cytoscape:layout graph))))))
         (full-screen (make-instance 'w:button ; Full screen toggle button
                                     :description "Full Screen"
                                     :on-click (list (lambda (inst)
                                                       (declare (ignore inst))
                                                       (cytoscape:toggle-fullscreen graph)))
                                     :layout (make-instance 'w:layout
                                                            :grid-area "full-screen")))
         (grid (make-instance 'w:grid-box ; Grid for the graph. Control buttons are in the bottom row.
                              :children (list graph fit show-sketches full-screen)
                              :layout (make-instance 'w:layout
                                                     :grid-gap "var(--jp-widgets-container-padding)"
                                                     :grid-template-columns "1fr min-content min-content min-content min-content 1fr"
                                                     :grid-template-rows "1fr min-content"
                                                     :grid-template-areas "\"graph graph graph graph graph graph\" \". fit show-sketches add-edges full-screen .\""))))
    ;; Delete the current node.
    (cytoscape:on-menu-command-select delete-node-command
      (lambda (inst id)
        (declare (ignore inst))
        (delete-elements graph id)))
    ;; Delete the current edge.
    (cytoscape:on-menu-command-select delete-edge-command
      (lambda (inst id)
        (declare (ignore inst))
        (delete-elements graph id)))
    ;; Add edges that connect all currently selected nodes.
    (cytoscape:on-menu-command-select add-edges-command
      (lambda (inst id)
        (declare (ignore inst id))
        (apply #'add-edges graph
                           (mapcan (lambda (element)
                                     (when (and (cytoscape:selected element)
                                                (not (cytoscape:removed element))
                                                (string= "nodes" (cytoscape:group element)))
                                       (list (cdr (assoc "id" (cytoscape:data element) :test #'string=)))))
                                   (cytoscape:elements graph)))))
    ;; Delete all selected elements
    (cytoscape:on-menu-command-select delete-elements-command
      (lambda (inst id)
        (declare (ignore inst id))
        (apply #'delete-elements graph
                                 (mapcan (lambda (element)
                                           (when (and (cytoscape:selected element)
                                                      (not (cytoscape:removed element)))
                                             (list (cdr (assoc "id" (cytoscape:data element) :test #'string=)))))
                                         (cytoscape:elements graph)))))
    (cw:make-simple-task-page container "Calculate Similarities" #'run-lomap
                              :label "Click button to calculate molecular similarities.")
    (setf (w:widget-%titles container) (append (w:widget-%titles container)
                                              (list "Calculation Graph"))
          (w:widget-children container) (append (w:widget-children container)
                                               (list grid)))
    (create-graph graph)
    ;; When all-ligands changes create a completely new graph.
    (w:observe *app* :all-ligands
      (lambda (instance type name old-value new-value source)
        (declare (ignore instance type name old-value new-value source))
        (create-graph graph)))
    ;; If selected-ligands or all-edges change then just update the graph.
    (w:observe *app* :selected-ligands
      (lambda (instance type name old-value new-value source)
        (declare (ignore instance type name old-value new-value source))
        (update-graph graph)))
    (w:observe *app* :all-edges
      (lambda (instance type name old-value new-value source)
        (declare (ignore instance type name old-value new-value source))
        (update-graph graph)))
    container))


(defun configure-jobs ()
  (let ((container (make-instance 'w:accordion)))
    (dolist (schema-group tirun::*default-calculation-settings* container)
      (cw:add-page container
                   (w:make-interactive-alist (second schema-group)
                                             (tirun::settings *tirun*))
                   (first schema-group)))))


(defun ensure-write-jobs (tirun jobs-dir &optional progress-callback)
  (ensure-directories-exist jobs-dir)
  (let ((*default-pathname-defaults* (pathname jobs-dir)))
    (ext:chdir *default-pathname-defaults*)
    (tirun:build-job-nodes tirun)
    (tirun:connect-job-nodes tirun all-edges)
    (let ((worklist (tirun:generate-jobs tirun :progress-callback progress-callback)))
      (with-open-file (fout #P"conf.sh" :direction :output :if-exists :supersede)
        (format fout "pmemd_cuda=pmemd.cuda
execute_cpu_local=0
execute_gpu_local=0
execute_lisp_local=0
worker=schando
distributor=s103.thirdlaw.tech
lisp_jobs_only_on=172.234.2.1
"))
      (with-open-file (sout #P"eg5_singlestep.dot" :direction :output)
        (tirundot:draw-graph-stream (list worklist) sout)))))

(defun calculate-jobs-dir (name)
  (let* ((jobs-dir (pathname (if (ext:getenv "TIRUN_DIRECTORY")
                                 (uiop:ensure-directory-pathname (ext:getenv "TIRUN_DIRECTORY"))
                                 (merge-pathnames (make-pathname :directory (list :relative "jobs"))
                                                  (uiop:ensure-directory-pathname (ext:getenv "HOME"))))))
         (dir (merge-pathnames (make-pathname :directory (list
                                                          :relative
                                                          name))
                               jobs-dir))
         (tar-file (merge-pathnames (make-pathname :name name
                                                   :type "tgz")
                                    jobs-dir)))
    (values dir tar-file jobs-dir)))

(defparameter *tar-thread* nil)
(defun generate-tar-file (messages job-name dir tar-file jobs-dir)
  (let ((cmd (list "tar"
                   "czvf"
                   (namestring tar-file)
                   "-C"
                   (namestring jobs-dir)
                   (format nil "./~a" job-name))))
    (setf (w:widget-value messages)
          (format nil "~atar cmd: ~a~%" (w:widget-value messages) cmd))
    (setf *tar-thread*
          (bordeaux-threads:make-thread
           (lambda ()
             (multiple-value-bind (ret child-pid stream)
                 (ext:vfork-execvp cmd t)
               (if (= 0 ret)
                   (setf (w:widget-value messages)
                         (format nil "~aSuccessfully generated tar file ~a~%" (w:widget-value messages) tar-file))
                   (setf (w:widget-value messages)
                         (format nil "~aFailed to submit job errno: ~a" (w:widget-value messages) ret)))
               ))))))

(defun write-jobs (&key (tirun *tirun*))
  (let* ((job-name (make-simple-input "Job name" :default (job-name *app*)))
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (progress (make-instance 'w:int-progress :description "Progress"))
         (go-button (make-instance
                     'w:button :description "Write jobs"
                     :style (button-style)
                     :on-click (list
                                (lambda (&rest args)
                                  (let ((name (w:widget-value (input-widget job-name))))
                                    (multiple-value-bind (dir tar-file from-path)
                                        (calculate-jobs-dir name)
                                      (setf (w:widget-value messages) (format nil "Writing jobs to -> ~a ...~%" dir))
                                      (if (probe-file dir)
                                          (progn
                                            (setf (w:widget-value messages) (format nil "The directory ~a already exists - aborting writing jobs" dir)))
                                          (progn
                                            (ensure-write-jobs tirun dir (let ((index 0)
                                                                               (lock (bordeaux-threads:make-recursive-lock "progress")))
                                                                           (lambda (max-progress)
                                                                             (bordeaux-threads:with-recursive-lock-held (lock)
                                                                               (setf (w:widget-value progress) index
                                                                                     (w:widget-max progress) max-progress)
                                                                               (incf index)))))
                                            (setf (w:widget-value messages) (format nil "~a~%Generating tar file.~%"
                                                                                    (w:widget-value messages)))
                                            (generate-tar-file messages name dir tar-file from-path)))))))))
         )
    (w:observe (input-widget job-name) :value (lambda (instance type name old-value new-value source)
                                                (setf (job-name *app*) new-value)))
    (make-instance 'w:v-box
                   :children (list
                              (hbox job-name)
                              (make-instance 'w:h-box
                                             :layout (box-layout)
                                             :children (list go-button))
                              (make-instance 'w:h-box
                                             :children (list progress))
                              (make-instance 'w:h-box
                                             :layout (box-layout)
                                             :children (list messages))))))
 
(defun read-submit-stream (num)
  (let ((buffer (make-array 15 :adjustable t :fill-pointer 0)))
    (flet ((grab ()
             (loop for line = (read-line (tirun-jupyter::submit-stream tirun-jupyter::*app*) nil :eof)
                   for index below num
                   do (vector-push-extend line buffer)
                   until (eq line :eof))))
      (grab)
      buffer)))

(defvar *submit-thread* nil)
(defun submit-calc (&key (calc *tirun*))
  (let* ((desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (distributor (make-simple-input "Distributor" :default (distributor *app*)))
         (job-name (make-simple-input "Job name" :default "default"))
         (ssh-key-file (make-instance 'jupyter-widgets:file-upload :description "SSH key file"
                                                                   :style (button-style)))
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (go-button (make-instance
                     'jupyter-widgets:button
                     :description "Submit jobs"
                     :style (button-style)
                     :tooltip "Click me"
                     :on-click (list
                                (lambda (&rest args)
                                  (let* ((schando-dir (UIOP/PATHNAME:ENSURE-DIRECTORY-PATHNAME
                                                       (or (ext:getenv "CLASP_SCANDO_PATH")
                                                           "/Users/meister/Development/schando/src/")))
                                         (submit-cmd (merge-pathnames (make-pathname :name "submit-to-distributor")
                                                                      schando-dir))
                                         (name (w:widget-value (input-widget job-name)))
                                         (cmd (list (namestring submit-cmd)
                                                    (namestring (calculate-jobs-dir name))
                                                    name)))
                                    (setf (w:widget-value messages)
                                          (format nil "About to evaluate: ~a" cmd))
                                    (setf *submit-thread*
                                          (bordeaux-threads:make-thread
                                           (lambda ()
                                             (multiple-value-bind (ret child-pid stream)
                                                 (ext:vfork-execvp cmd t)
                                               #+(or)(if (= 0 ret)
                                                         (setf (w:widget-value messages)
                                                               (format nil "Successfully submitted job ~s~%" cmd)
                                                               (submit-stream *app*) stream)
                                                         (setf (w:widget-value messages)
                                                               (format nil "Failed to submit job errno: ~a" ret)))
                                               (sleep 1000000)))
                                           :name "submit"))))))))
    (jupyter-widgets:observe *app* :job-name (lambda (instance type name old-value new-value source)
                                               (setf (w:widget-value (input-widget job-name)) new-value)))
    (make-instance 'w:v-box
                   :children (append (list (hbox distributor) (hbox job-name) ssh-key-file)
                                     (list (make-instance 'w:h-box
                                                          :layout (box-layout)
                                                          :children (list go-button))
                                           (make-instance 'w:h-box
                                                          :layout (box-layout)
                                                          :children (list messages)))))))


(defun check-calc (&key (tirun *tirun*))
  (let* ((job-name (make-simple-input "Job name" :default (job-name *app*)))
         (distributor (make-simple-input "Distributor" :default (distributor *app*)))
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (go-button (make-instance 'w:button :description "Check calculation"
                                             :style (button-style)
                                             :on-click (list
                                                        (lambda (&rest args)
                                                          (let* ((schando-dir (UIOP/PATHNAME:ENSURE-DIRECTORY-PATHNAME
                                                                               (or (ext:getenv "CLASP_SCANDO_PATH")
                                                                                   "/Users/meister/Development/schando/src/")))
;;; execute-in-jobdir jobdir default ls \*.SUCCESS \| wc -l
                                                                 (submit-cmd (merge-pathnames (make-pathname :name "execute-in-jobdir")
                                                                                              schando-dir))
                                                                 (name (w:widget-value (input-widget job-name))))
                                                            (multiple-value-bind (ret child-pid stream)
                                                                (ext:vfork-execvp (list (namestring submit-cmd)
                                                                                        (namestring (calculate-jobs-dir name))
                                                                                        name
                                                                                        "ls" "-l"
                                                                                        ) t)
                                                              (if (= 0 ret)
                                                                  (let ((buffer (make-string-output-stream)))
                                                                    (sleep 2)
                                                                    (tagbody
                                                                     top
                                                                       (let ((line (read-line stream nil :eof)))
                                                                         (unless (eq line :eof)
                                                                           (format buffer "~a" line)
                                                                           (go top)))
                                                                       )
                                                                    (setf (w:widget-value messages)
                                                                          (get-output-stream-string buffer)))
                                                                  (setf (w:widget-value messages)
                                                                        (format nil "Failed to submit job errno: ~a" ret))))))))))
    (make-instance 'w:v-box
                   :children (list
                              (hbox distributor)
                              (hbox job-name)
                              (make-instance 'w:h-box
                                             :layout (box-layout)
                                             :children (list go-button))
                              (make-instance 'w:h-box
                                             :layout (box-layout)
                                             :children (list messages))))))

(defun composer (&key (app *app*))
  (let ((composer (make-instance 'structure-editor::composer :app app)))
    (w:observe *app* :composer-json
               (lambda (inst type name old-value new-value source)
                 (when new-value
                   (save-app)
                   )))
    (let ((widgets (structure-editor::composer-accordion composer)))
      (jupyter-widgets:display widgets)
      (kekule:fit (structure-editor::composer-kekule composer))
      (values))))


