(in-package :tirun-jupyter)


(defmacro def-dyn-widget (name initial-value)
  (let ((var-name (gensym)))
    `(progn
       (defparameter ,var-name nil)
       (defun ,name ()
         (or ,var-name
             (setf ,var-name ,initial-value))))))


(def-dyn-widget button-style (make-instance 'w:button-style :button-color "aquamarine"))


(defparameter +default-distributor+ "s103.thirdlaw.tech")
(defparameter +default-job-name+ "default")

(defclass workspace (jupyter-widgets:has-traits)
  ((workspace-path
    :accessor workspace-path
    :initform #P"workspace-1.tirun")
   (composer-data
    :accessor composer-data
    :initarg :composer-data
    :initform nil
    :trait t)
   (receptor
    :accessor receptor
    :initarg :receptor
    :initform nil
    :trait t)
   (template-ligand
    :accessor template-ligand
    :initarg :template-ligand
    :initform nil
    :trait t)
   (all-ligands
    :accessor all-ligands
    :initarg all-ligands
    :initform nil
    :trait t)
   (selected-ligands
    :accessor selected-ligands
    :initarg selected-ligands
    :initform nil
    :trait t)
   (similarities
    :accessor similarities
    :initarg :similarities
    :initform nil
    :trait t)
   (all-edges
    :accessor all-edges
    :initarg :all-edges
    :initform nil
    :trait t)
   (distributor
    :accessor distributor
    :initarg :distributor
    :initform +default-distributor+
    :trait t)
   (job-name
    :accessor job-name
    :initarg :job-name
    :initform +default-job-name+
    :trait t)
   (submit-stream
    :accessor submit-stream
    :initform nil
    :trait t)
   (smirks-pattern
    :accessor smirks-pattern
    :initarg :smirks-pattern
    :initform nil
    :trait t)
   (tirun-settings
    :accessor tirun-settings
    :initarg :tirun-settings
    :initform (tirun:default-calculation-settings)
    :trait t)
   (tirun-stages
    :accessor tirun-stages
    :initarg :tirun-stages
    :initform 1
    :trait t))
  (:metaclass jupyter-widgets:trait-metaclass))


(define-symbol-macro all-ligands
  (all-ligands *workspace*))


(define-symbol-macro selected-ligands
  (selected-ligands *workspace*))


(define-symbol-macro template-ligand
  (template-ligand *workspace*))


(define-symbol-macro all-edges
  (all-edges *workspace*))


(cando:make-class-save-load workspace :skip-slot-names (jupyter-widgets:on-trait-change))
   
(defvar *workspace* nil)


(defmethod make-instance ((class (eql (find-class 'workspace))) &rest initargs)
  (if *workspace*
    (apply #'reinitialize-instance *workspace* initargs)
    (setf *workspace* (call-next-method))))


(def-dyn-widget box-layout
  (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))


(defun initialize-workspace ()
  (let ((dest-dir (let ((dd (ext:getenv "CANDO_JOBS_DIRECTORY")))
                    (if dd
                      dd
                      (pathname "~/jobs/")))))
    (format t "Jobs will be saved to ~a~%setenv CANDO_JOBS_DIRECTORY if you want it to go elsewhere~%" (namestring dest-dir))
    (ensure-directories-exist dest-dir))
  (let ((*standard-output* (make-string-output-stream)))
    (cando-user:setup-default-paths)
    (cando-user:load-atom-type-rules "ATOMTYPE_GFF.DEF")
    (cando-user:source "leaprc.ff14SB.redq")
    (cando-user:source "leaprc.gaff")
    (leap:add-pdb-res-map '((1 :NMA :NME))))
  t)


(defun run-workspace-task (action parameter progress-callback)
  (declare (ignore progress-callback))
  (case action
    (:new
      (new-workspace parameter))
    (:initialize
      (initialize-workspace))
    (:open
      (open-workspace parameter))
    (:save
      (save-workspace parameter))))


(defun workspace (&rest initargs)
  (apply #'make-instance 'workspace initargs)
  (let* ((container (make-instance 'w:accordion :selected-index 0))
         (page (cw:make-workspace-task-page container "TIRUN Workspace" #'run-workspace-task :file-type "tirun")))
    (cando-widgets:run-task page :initialize nil)
    container))


(defun new-workspace (&optional (path (workspace-path *workspace*)))
  (setf (composer-data *workspace*) nil
        (receptor *workspace*) nil
        (template-ligand *workspace*) nil
        (all-ligands *workspace*) nil
        (selected-ligands *workspace*) nil
        (similarities *workspace*) nil
        (all-edges *workspace*) nil
        (distributor *workspace*) +default-distributor+
        (job-name *workspace*) +default-job-name+
        (submit-stream *workspace*) nil
        (smirks-pattern *workspace*) nil
        (tirun-settings *workspace*) (tirun:default-calculation-settings)
        (tirun-stages *workspace*) 1)
  (save-workspace path))


(defun open-workspace (&optional (path (workspace-path *workspace*)))
  (setf (workspace-path *workspace*) path)
  (cond
    ((probe-file path)
      (finish-output)
      (format t "Opening workspace ~A...~%" path)
      (finish-output)
      (cando:load-cando path)
      (write-line "Open complete.")
      (finish-output)
      t)
    (t
      (format *error-output* "Unable to find workspace file named ~A~%" path)
      (finish-output *error-output*)
      nil)))


(defun save-workspace (&optional (path (workspace-path *workspace*)))
  (setf (workspace-path *workspace*) path)
  (format t "Saving workspace ~A...~%" path)
  (finish-output)
  (ensure-directories-exist path)
  (cando:save-cando *workspace* path)
  (write-line "Save complete.")
  (finish-output)
  t)


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


(defun parse-receptor (action parameter progress-callback)
  "Called from file task page to parse receptor."
  (declare (ignore action progress-callback))
  (handler-case
      (with-slots (receptor)
                  *workspace*
        (let ((receptor-string (babel:octets-to-string (cdr (assoc "content" (car parameter) :test #'string=)))))
          (setf receptor
                (with-input-from-string (input-stream receptor-string)
                  (leap.pdb:load-pdb-stream input-stream)))
          (save-workspace)))
    (leap.pdb:pdb-read-error (condition)
      (princ condition *error-output*)
      nil)
    (end-of-file (condition)
      (princ condition *error-output*)
      nil)))


(defun parse-template-ligand (action parameter progress-callback)
  "Called from file task page to parse template-ligand."
  (declare (ignore action progress-callback))
  (handler-case
      (with-slots (template-ligand)
                  *workspace*
        (let ((template-ligand-string (babel:octets-to-string (cdr (assoc "content" (car parameter) :test #'string=)))))
          (format t "About to parse template ligand pdb file~%")
          (setf template-ligand
                (with-input-from-string (input-stream template-ligand-string)
                  (leap.pdb:load-pdb-stream input-stream :ignore-missing-topology t)))
          (format t "Done parsing ligand pdb file~%")
          (finish-output)
          (save-workspace)))
    (leap.pdb:pdb-read-error (condition)
      (princ condition *error-output*)
      nil)
    (end-of-file (condition)
      (princ condition *error-output*)
      nil)))


(defun load-receptor ()
  "Present a graphical interface to load and view a receptor into the current application."
  (let* ((container (make-instance 'w:accordion :selected-index 0))
         (ngl-page (make-instance 'resizable-box:resizable-grid-box
                                  :enable-full-screen t
                                  :layout (make-instance 'resizable-box:resizable-layout
                                                         :resize "vertical"
                                                         :grid-gap "1em"
                                                         :min-height "480px"
                                                         :overflow "hidden"
                                                         :padding "0 16px 0 0"
                                                         :grid-template-rows "1fr min-content"
                                                         :grid-template-columns "1fr min-content"
                                                         :grid-template-areas "'ngl-content ngl-vertical-toolbar'
                                                                               'ngl-horizontal-toolbar .")))
         (ngl (cw:make-ngl-structure-viewer ngl-page :auto-view "receptor")))
    (cw:make-file-task-page container "Load Receptor" #'parse-receptor :accept ".pdb")
    ;; Make a simple page to display the receptor in nglview
    (setf (w:widget-%titles container) (append (w:widget-%titles container)
                                               (list "View Receptor"))
          (w:widget-children container) (append (w:widget-children container)
                                               (list ngl-page)))
    ;; Add the receptor if one is already defined.
    (when (receptor *workspace*)
      (cw:add-receptor ngl (receptor *workspace*)))
    ;; Listen for changes to the receptor and add or delete when needed.
    (w:observe *workspace* :receptor
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (when new-value
          (cw:add-receptor ngl new-value))))
    container))


(defun load-template-ligand ()
  "Present a graphical interface to load and view a template-ligand into the current application."
  (let* ((container (make-instance 'w:accordion :selected-index 0))
         (ngl-page (make-instance 'resizable-box:resizable-grid-box
                                  :enable-full-screen t
                                  :layout (make-instance 'resizable-box:resizable-layout
                                                         :resize "vertical"
                                                         :grid-gap "1em"
                                                         :overflow "hidden"
                                                         :padding "0 16px 0 0"
                                                         :min-height "480px"
                                                         :grid-template-rows "1fr min-content"
                                                         :grid-template-columns "1fr min-content"
                                                         :grid-template-areas "'ngl-content ngl-vertical-toolbar'
                                                                               'ngl-horizontal-toolbar .")))
         (ngl (cw:make-ngl-structure-viewer ngl-page :auto-view "template")))
    (cw:make-file-task-page container "Load Template ligand" #'parse-template-ligand :accept ".pdb")
    ;; Make a simple page to display the receptor in nglview
    (setf (w:widget-%titles container) (append (w:widget-%titles container)
                                              (list "View Template Ligand"))
          (w:widget-children container) (append (w:widget-children container)
                                               (list ngl-page)))
    ;; Add the template ligand if one is already defined.
    (when (template-ligand *workspace*)
      (cw:add-template ngl (template-ligand *workspace*)))
    ;; Listen for changes to the template-ligand and add or delete when needed.
    (w:observe *workspace* :template-ligand
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (when new-value
          (cw:add-template ngl new-value))))
    container))


(defclass view-ligand-page (cw:page)
  ((ngl-structure-viewer
     :accessor view-ligand-ngl-structure-viewer)
   (structure
     :reader view-ligand-structure
     :initform (make-instance 'w:html
                              :layout (make-instance 'w:layout
                                                     :grid-area "structure"))
     :documentation "A sketch view for the selected ligand.")
   (dropdown
     :reader view-ligand-dropdown
     :initform (make-instance 'w:dropdown
                              :description "Active Ligand"
                              :style (make-instance 'w:description-style
                                                    :description-width "min-content")
                              :layout (make-instance 'w:layout
                                                     :margin ".5em"
                                                     :width "max-content"))
     :documentation "A dropdown of all the available ligands.")
   (template-dropdown
     :reader view-ligand-template-dropdown
     :initform (make-instance 'w:dropdown
                              :description "Template Ligand"
                              :style (make-instance 'w:description-style
                                                    :description-width "min-content")
                              :layout (make-instance 'w:layout
                                                     :margin ".5em"
                                                     :width "max-content"))
     :documentation "A dropdown of all the available ligands."))
  (:metaclass jupyter-widgets:trait-metaclass)
  (:documentation "A page to view ligands that has a sketch and an nglview widget along with controls and ligand selector.")
  (:default-initargs
    :enable-full-screen t
    :layout (make-instance 'resizable-box:resizable-layout
                           :min-height "480px"
                           :overflow "hidden"
                           :resize "vertical"
                           :grid-gap "1em"
                           :grid-template-rows "1fr min-content"
                           :grid-template-columns "2fr 3fr min-content"
                           :grid-template-areas "'structure     ngl-content            ngl-vertical-toolbar'
                                                 'structure-ctl ngl-horizontal-toolbar .'")))


(defun molecule-name (molecule)
  "Get a string representation of the structure."
  (symbol-name (chem:get-name molecule)))


(defun on-ligand-select (instance id &optional previous-id)
  "Select a ligand on a view-ligand-page and hide a previous ligand if one was selected."
  (with-slots (structure)
              instance
    (let ((mol (find id all-ligands :key #'molecule-name :test #'equal)))
      (when mol
        ;; Load the sketch
        (setf (w:widget-value structure) (cw:sketch-molecule mol))
        ;; Show the ligand if it is already in nglview otherwise load it.
        (let ((agg (chem:make-aggregate (chem:get-name mol))))
          (chem:add-matter agg mol)
          (cw:add-ligand (view-ligand-ngl-structure-viewer instance) (molecule-name mol) agg))))))


(defun refresh-ligands-view (instance)
  "Refresh the ligand view and the ngl view."
  (with-slots (dropdown template-dropdown)
              instance
    (cond
      (all-ligands
        (setf (w:widget-disabled dropdown) nil
              (w:widget-%options-labels template-dropdown) (append (list "None")
                                                                   (when template-ligand
                                                                     (list "Template Ligand"))
                                                                   (mapcar #'molecule-name all-ligands))
              (w:widget-index template-dropdown) (if template-ligand 1 0)
              (w:widget-%options-labels dropdown) (mapcar #'molecule-name all-ligands)
              (w:widget-index dropdown) 0)
        (cw:sketch-molecules all-ligands)
        (on-ligand-select instance (molecule-name (car all-ligands))))
      (t
        (setf (w:widget-disabled dropdown) t
              (w:widget-%options-labels dropdown) nil)))))

(defmethod initialize-instance :after ((instance view-ligand-page) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (setf (w:widget-children instance)
        (list (make-instance 'w:box
                             :children (list (view-ligand-structure instance))
                             :layout (make-instance 'w:layout
                                                    :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                    :justify-content "center"
                                                    :align-items "center"))
              (make-instance 'w:box
                             :children (list (view-ligand-dropdown instance)
                                             (view-ligand-template-dropdown instance))
                             :layout (make-instance 'w:layout
                                                    :flex-flow "row wrap"
                                                    :justify-content "center"
                                                    ;:margin "-.5em"
                                                    :align-items "baseline"
                                                    :align-content "flex-start"
                                                    :grid-area "structure-ctl"))))
  (setf (view-ligand-ngl-structure-viewer instance)
        (cw:make-ngl-structure-viewer instance :auto-view "ligand"))
  ;; When the ligand changes update the sketch and nglview.
  (w:observe (view-ligand-dropdown instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (on-ligand-select instance new-value old-value)))
  ;; When the template changes update the sketch and nglview.
  (w:observe (view-ligand-template-dropdown instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (unless (equal new-value "None")
        (cw:add-template (view-ligand-ngl-structure-viewer instance)
                         (let ((agg (chem:make-aggregate :template)))
                           (chem:add-matter agg
                                            (if (equal "Template Ligand" new-value)
                                              template-ligand
                                              (find new-value all-ligands :key #'molecule-name :test #'equal)))
                           agg)))))
  ;; When the receptor changes reload the receptor in nglview.
  (w:observe *workspace* :receptor
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name source))
      (cw:add-receptor (view-ligand-ngl-structure-viewer instance) new-value)))
  ;; When the template ligand changes reload the receptor in nglview.
  (w:observe *workspace* :template-ligand
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name source))
      (with-slots (template-dropdown)
                  instance
        (cw:add-template (view-ligand-ngl-structure-viewer instance) new-value)
        (setf (w:widget-%options-labels template-dropdown) (append (list "None")
                                                                   (when new-value
                                                                     (list "Template Ligand"))
                                                                   (mapcar #'molecule-name all-ligands))
              (w:widget-index template-dropdown) (if new-value 1 0)))))
  ;; If the all-ligands changes then reload the whole thing.
  (w:observe *workspace* :all-ligands
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name new-value source))
      (cw:clear-ligands (view-ligand-ngl-structure-viewer instance))
      (refresh-ligands-view instance))))


(defun make-view-ligand-page (container title)
  "Create an instance of a view-ligand-page and add it to the container."
  (let ((page (make-instance 'view-ligand-page :container container)))
    ;(setf (w:widget-grid-area (w:widget-layout (view-ligand-ngl-structure-viewer page))) "ngls")
    (cw:add-page container page title)
    ;; If the receptor is already set then load the current receptor.
    (when (receptor *workspace*)
      (cw:add-receptor (view-ligand-ngl-structure-viewer page) (receptor *workspace*)))
    (when (template-ligand *workspace*)
      (cw:add-template (view-ligand-ngl-structure-viewer page) (template-ligand *workspace*)))
    ;; Update the view
    (refresh-ligands-view page)
    (values)))


(defun parse-ligands (action parameter progress-callback)
  "Called from the file task page to actually parse the ligands file."
  (declare (ignore action progress-callback))
  (handler-case
      (let ((ligands-string (babel:octets-to-string (cdr (assoc "content" (car parameter) :test #'string=))))
            (ligands-format (pathname-type (cdr (assoc "name" (car parameter) :test #'string=)))))
        (setf all-ligands
              (with-input-from-string (input-stream ligands-string)
                (cond
                  ((string-equal ligands-format "sdf")
                    (sdf:parse-sdf-file input-stream))
                  ((string-equal ligands-format "mol2")
                    (chem:read-mol2-list input-stream))
                  (t
                    (error "Unknown file type ~A." ligands-format))))
              selected-ligands all-ligands
              all-edges nil)
        (save-workspace))
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
    (w:link sel :molecules *workspace* :all-ligands)
    (w:link sel :selected *workspace* :selected-ligands)
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


(defun run-pas (action parameter progress-callback)
  "Called by the simple-task-page to do the PAS calculation."
  (declare (ignore action parameter))
  (write-line "Running position analogue scanning...")
  (finish-output)
  ;; Loop over the current selected ligands and run PAS on each one accumulating the results.
  (setf all-ligands
        (loop for template-ligand in selected-ligands
              for new-ligands = (pas:pas (smirks-pattern *workspace*) template-ligand)
              append new-ligands
              do (format t "Using template ligand of ~a built the ligands ~{~a~#[~;, and ~:;, ~]~}...~%"
                         (molecule-name template-ligand)
                         (mapcar #'molecule-name new-ligands))
              do (finish-output)
              do (funcall progress-callback (length selected-ligands)))
        selected-ligands all-ligands
        all-edges nil)
  (save-workspace))


(defparameter +smirks-patterns+
  '(("[D3&CH1:1]-[#1:2]>>[N:1]" . "Replace aromatic methine carbons with nitrogen.")
    ("[D3&CH1:1]-[#1:2]>>[O:1]" . "Replace aromatic methine carbons with oxygen.")
    ("[D3&CH1:1]-[#1:2]>>[CH1:1]-[F:2]" . "Generate all mono-halo [choose F,Br,Cl,I] analouges of template.")
    ("[D3&CH1:1]-[#1:2]>>[CH1:1]-[CH3:2]" . "Perform all alkylations of template.")))


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
                              :value (smirks-pattern *workspace*)
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
    ;; Connect the SMIRKS editor to *workspace* smirks-pattern
    (w:link *workspace* :smirks-pattern smirks-pattern-editor :value)
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


(defun encode-svg (svg)
  (concatenate 'string
               "data:image/svg+xml,"
               (quri:url-encode svg)))


(defun element-id (element)
  (cdr (assoc "id" (cytoscape:data element) :test #'string=)))


(defun get-selected-names (graph)
  (mapcan (lambda (element)
            (when (cytoscape:selected element)
              (let* ((id (element-id element))
                     (pos (position #\| id))) ; If there is vertical bar it is an edge.
                (if pos
                  (list (subseq id 0 pos) (subseq id (1+ pos)))
                  (list id)))))
          (cytoscape:elements graph)))

(defun match-ligands (ligands)
  (if (= 2 (length ligands))
    (multiple-value-bind (equivs diff1 diff2)
                         (tirun::calculate-masks-for-molecules (first ligands) (second ligands))
      (list (mapcar (lambda (mol) (chem:get-name mol)) diff1)
            (mapcar (lambda (mol) (chem:get-name mol)) diff2)))
    (mapcar (lambda (ligand)) ligands)))

(defun generate-stylesheet (atoms)
  (when atoms
    (with-output-to-string (output-stream)
      (format output-stream "~{.~A~^,~%~} {~%  filter: url(#highlight);~%}" atoms))))


(defun update-selected-sketches (graph unselected-names)
  (let* ((selected-names (get-selected-names graph))
         (selected (mapcan (lambda (ligand)
                             (when (position (molecule-name ligand) selected-names :test #'string=)
                               (list ligand)))
                           all-ligands))
         (atom-matches (match-ligands selected)))
    (dolist (element (cytoscape:elements graph))
      (when (and (position (element-id element) unselected-names :test #'equal))
        ;(not (position (element-id element) selected-names :test #'equal)))
        (setf (cdr (assoc "image" (cytoscape:data element) :test #'string=))
              (encode-svg (cando-widgets:sketch-molecule (find (element-id element) all-ligands :test #'string= :key #'molecule-name))))
        (jupyter-widgets:notify-trait-change element :dict :data (cytoscape:data element) (cytoscape:data element) t)))
    (loop for ligand in selected
          for atoms in atom-matches
          for element = (find (molecule-name ligand) (cytoscape:elements graph) :test #'string= :key #'element-id)
          do (setf (cdr (assoc "image" (cytoscape:data element) :test #'string=))
                   (encode-svg (cando-widgets:sketch-molecule ligand (generate-stylesheet atoms))))
          do (jupyter-widgets:notify-trait-change element :dict :data (cytoscape:data element) (cytoscape:data element) t))))


(defun make-node (graph ligand)
  "Construct a cytoscape node based on a ligand."
  (let ((name (molecule-name ligand)))
    (multiple-value-bind (svg width height)
                         (cw:sketch-molecule ligand)
      ;; Save sketch dimension and the sketch itself into the data slot. The sketch is SVG text
      ;; so we just need do a percent encoded string for the reserved URL characters.
      (make-instance 'cytoscape:element
                     :group "nodes"
                     :removed (not (position ligand selected-ligands))
                     :on-trait-change (list (cons :selected
                                                  (lambda (instance type nm old-value new-value source)
                                                    (declare (ignore instance type nm old-value source))
                                                    (update-selected-sketches graph (unless new-value (list name))))))
                     :data (list (cons "id" name)
                                 (cons "label" name)
                                 (cons "width" width)
                                 (cons "height" height)
                                 (cons "image" (encode-svg svg)))))))


(defun make-edge (graph edge)
  "Construct an cytoscape edge based on an edge definition when is a list of
  the form (source target label)."
  (let ((name-1 (symbol-name (first edge)))
        (name-2 (symbol-name (second edge)))
        (similarity (cdr (assoc edge (similarities *workspace*) :test #'equal))))
    (make-instance 'cytoscape:element
                   :group "edges"
                   :on-trait-change (list (cons :selected
                                                (lambda (instance type name old-value new-value source)
                                                  (declare (ignore instance type name old-value source))
                                                  (update-selected-sketches graph (unless new-value
                                                                                    (list name-1 name-2))))))
                   :data (list (cons "id" (edge-id name-1 name-2))
                               (cons "label" (format nil "~@[S = ~3,2f~]" similarity))
                               (cons "source" name-1)
                               (cons "target" name-2)))))


(defun create-graph (graph)
  "Create a graph from all-ligands and all-edges."
  (setf (cytoscape:elements graph)
        (nconc (mapcar (lambda (ligand)
                         (make-node graph ligand))
                       all-ligands)
               (mapcar (lambda (edge)
                         (make-edge graph edge))
                       all-edges))))


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
        (push (make-edge graph edge) new-edges)))
    ;; If new edges have been created then just add them. The layout will be updated automatically.
    ;; Otherwise force the layout to update.
    (if new-edges
      (setf (cytoscape:elements graph) (append (cytoscape:elements graph) new-edges))
      (cytoscape:layout graph))))


(defun run-lomap (action parameter progress-callback)
  "This function is called by the task page and updates a progress widget as lomap does the calculation."
  (declare (ignore action parameter))
  (format t "Calculating similarity matrix of ~{~a~#[~;, and ~:;, ~]~}...~%"
          (mapcar #'molecule-name selected-ligands))
  (finish-output) ; Make sure the message is synced to the frontend.
  (let* ((similarity-matrix (lomap::similarity-matrix selected-ligands
                                                      :advance-progress-callback progress-callback))
         (multigraph (lomap::lomap-multigraph (lomap::similarity-multigraph selected-ligands
                                                                            similarity-matrix)))
         new-edges)
    ;; Loop over the similarity matrix and add the similarities
    (setf (similarities *workspace*)
          (loop for ligand1 in selected-ligands
                for index1 from 0
                append (loop for ligand2 in selected-ligands
                             for index2 from 0 below index1
                             collect (cons (sort (list (chem:get-name ligand1)
                                                       (chem:get-name ligand2))
                                                 #'string<)
                                           (aref similarity-matrix index1 index2)))))
    ;; Loop over the subgraphs and the vertices and edges and create edges.
    (dolist (subgraph (lomap::subgraphs multigraph))
      (dolist (edge (lomap::edges subgraph))
        (push (sort (list (chem:get-name (lomap:molecule (lomap:vertex1 edge)))
                          (chem:get-name (lomap:molecule (lomap:vertex2 edge))))
                    #'string<)
              new-edges)))
    (setf all-edges new-edges))
  (write-line "Similarity matrix calculation completed.")
  (finish-output)
  (save-workspace))


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
         (layouts (list (list "Circle layout"
                              (make-instance 'cytoscape:circle-layout))
                        (list "Concentric layout"
                              (make-instance 'cytoscape:concentric-layout))
                        (list "Breadth-First layout"
                              (make-instance 'cytoscape:breadth-first-layout))
                        (list "CoSE layout"
                              (make-instance 'cytoscape:cose-layout))
                        (list "Dagre layout"
                              (make-instance 'cytoscape:dagre-layout))
                        (list "Grid layout"
                              (make-instance 'cytoscape:grid-layout))))
         (graph (make-instance 'cytoscape:cytoscape-widget
                               :graph-layouts (cdr (fourth layouts))
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
                                                      :height "auto"
                                                      :align-self "stretch"
                                                      :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                      :grid-area "graph")))
         (layout-select (make-instance 'w:dropdown ; A dropdown to select the layout
                                       :%options-labels (mapcar #'car layouts)
                                       :index 3
                                       :on-trait-change (list (cons :index
                                                                    (lambda (inst type name old-value new-value source)
                                                                      (declare (ignore inst type name old-value source))
                                                                      (setf (cytoscape:graph-layouts graph)
                                                                            (cdr (nth new-value layouts)))
                                                                      (cytoscape:layout graph))))
                                       :layout (make-instance 'w:layout
                                                              :width "max-content"
                                                              :grid-area "layout")))
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
         (grid (make-instance 'resizable-box:resizable-grid-box ; Grid for the graph. Control buttons are in the bottom row.
                              :children (list graph fit layout-select show-sketches)
                              :enable-full-screen t
                              :layout (make-instance 'resizable-box:resizable-layout
                                                     :resize "vertical"
                                                     :overflow "hidden"
                                                     :min-height "480px"
                                                     :padding "0 24px 0 0"
                                                     :grid-gap "var(--jp-widgets-container-padding)"
                                                     :grid-template-columns "1fr min-content min-content min-content min-content 1fr"
                                                     :grid-template-rows "1fr min-content"
                                                     :grid-template-areas "\"graph graph graph graph graph graph\" \". fit layout show-sketches add-edges .\""))))
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
    (w:observe *workspace* :all-ligands
      (lambda (instance type name old-value new-value source)
        (declare (ignore instance type name old-value new-value source))
        (create-graph graph)))
    ;; If selected-ligands or all-edges change then just update the graph.
    (w:observe *workspace* :selected-ligands
      (lambda (instance type name old-value new-value source)
        (declare (ignore instance type name old-value new-value source))
        (update-graph graph)))
    (w:observe *workspace* :all-edges
      (lambda (instance type name old-value new-value source)
        (declare (ignore instance type name old-value new-value source))
        (update-graph graph)))
    container))


(defclass jobs-config-page (cw:page)
  ((job-name
     :accessor job-name
     :initform (make-instance 'w:text
                              :description "Job name"
                              :value (job-name *workspace*)
                              :style (make-instance 'w:description-style :description-width "12em")
                              :layout (make-instance 'w:layout :grid-area "name"))
     :documentation "The job name.")
   (ti-stages
     :accessor ti-stages
     :initform (make-instance 'w:radio-buttons
                              :options '(1 3)
                              :%options-labels '("one stage" "three stage")
                              :description "TI stages"
                              :value (tirun-stages *workspace*)
                              :style (make-instance 'w:description-style :description-width "12em")
                              :layout (make-instance 'w:layout :grid-area "stages"))
     :documentation "The number of stages in the TI calculation."))
  (:metaclass jupyter-widgets:trait-metaclass)
  (:documentation "A page to set job name and the number of stages.")
  (:default-initargs
    :layout (make-instance 'w:layout
                           :width "100%"
                           :max-height "12em"
                           :grid-gap "0.1em 1em"
                           :grid-template-rows "min-content min-content"
                           :grid-template-columns "1fr"
                           :grid-template-areas "\"name\" \"stages\"")))


(defmethod initialize-instance :after ((instance jobs-config-page) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (with-slots (job-name ti-stages)
              instance
    ;; Connect the job-name to *workspace* job-name
    (w:link *workspace* :job-name job-name :value)
    (setf (w:widget-children instance) (list job-name ti-stages))
    ;; ti-stages isn't a trait so we can't use link.
    (w:observe ti-stages :value
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (setf (tirun-stages *workspace*) new-value)))))


(defun make-jobs-config-page (container title)
  "Create a jobs-config-page instance and add it to the container."
  (let* ((page (make-instance 'jobs-config-page :container container)))
    (cw:add-page container page title)
    (values)))


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

(defun generate-tar-file (job-name dir tar-file jobs-dir)
  (let ((cmd (list "tar"
                   "czvf"
                   (namestring tar-file)
                   "-C"
                   (namestring jobs-dir)
                   (format nil "./~a" job-name))))
    (format t "tar cmd: ~a~%" cmd)
    (finish-output)
    (multiple-value-bind (ret child-pid stream)
                         (ext:vfork-execvp cmd nil)
      (cond
        ((or (null ret)
             (zerop ret))
          (format t "Successfully generated tar file ~a~%" tar-file)
          (finish-output)
          t)
        (t
          (format *error-output* "Failed to submit job errno: ~a~%" ret)
          (finish-output *error-output*)
          nil)))))


(defun run-write-jobs (action parameter progress-callback)
  "Called by the write task page to start the write jobs task."
  (declare (ignore action parameter))
  (multiple-value-bind (dir tar-file from-path)
      (calculate-jobs-dir (job-name *workspace*))
    (save-workspace)
    (format t "Writing jobs to -> ~a ...~%" dir)
    (finish-output)
    (cond
      ((probe-file dir)
       (format *error-output* "The directory ~a already exists - aborting writing jobs" dir)
       (finish-output *error-output*)
       nil)
      (t
       (let ((tirun-calculation (make-instance 'tirun:tirun-calculation)))
         (prepare-calculation tirun-calculation *workspace*)
         (ensure-write-jobs tirun-calculation dir progress-callback)
         (write-line "Generating tar file.")
         (finish-output)
         (generate-tar-file (job-name *workspace*) dir tar-file from-path))))))

(defun only-residue (molecule)
  (unless (= (chem:content-size molecule) 1)
    (error "The molecule ~a must have only one residue" molecule))
  (chem:content-at molecule 0))

(defun short-residue-name (index)
  (let ((label (format nil "!~36,2,'0r" index)))
    (intern label :keyword)))

(defun prepare-calculation (calc workspace)
  (let ((residue-name-to-pdb-alist nil))
    (loop for molecule in (selected-ligands workspace)
          for residue = (only-residue molecule)
          for residue-index from 0
          for residue-name = (chem:get-name residue)
          when (> (length (string residue-name)) 3)
            do (push (cons residue-name (short-residue-name residue-index))
                     residue-name-to-pdb-alist))
    (let ((ligands (mapcar
                    (lambda (molecule)
                      (let ((residue (only-residue molecule)))
                        (format t "Ligand: ~a~%" molecule)
                        (finish-output)
                        (make-instance 'tirun:tirun-structure
                                       :name (chem:get-name molecule)
                                       :drawing molecule
                                       :molecule molecule
                                       :core-residue residue
                                       :core-residue-name (chem:get-name residue)
                                       :core-atoms (chem:all-atoms-as-list molecule nil)
                                       :net-charge 0.0 #|not always!|#)))
                    (selected-ligands workspace))))
      (setf (tirun:mask-method calc) :closest
            (tirun:settings calc) (tirun-settings workspace)
            (tirun:ti-stages calc) (tirun-stages workspace)
            (tirun:residue-name-to-pdb-alist calc) residue-name-to-pdb-alist
            (tirun:receptors calc) (list (receptor workspace))
            (tirun:ligands calc) ligands))))

(defun jobs ()
  "Configure and write TI jobs to jobs directory."
  (let ((container (make-instance 'w:accordion :selected-index 0)))
    (dolist (schema-group tirun::*default-calculation-settings*)
      (cw:add-page container
                   (w:make-interactive-alist (second schema-group)
                                             (tirun-settings *workspace*)
                                             :owner *workspace*
                                             :name :tirun-settings)
                   (first schema-group)))
    (make-jobs-config-page container "Configure Jobs")
    (cw:make-simple-task-page container "Write Jobs" #'run-write-jobs
                              :label "Click button to write jobs.")
    container))


(defun read-submit-stream (num)
  (let ((buffer (make-array 15 :adjustable t :fill-pointer 0)))
    (flet ((grab ()
             (loop for line = (read-line (submit-stream *workspace*) nil :eof)
                   for index below num
                   do (vector-push-extend line buffer)
                   until (eq line :eof))))
      (grab)
      buffer)))

(defvar *submit-thread* nil)
(defun submit-calc ()
  (let* ((desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (distributor (make-simple-input "Distributor" :default (distributor *workspace*)))
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
                                                               (submit-stream *workspace*) stream)
                                                         (setf (w:widget-value messages)
                                                               (format nil "Failed to submit job errno: ~a" ret)))
                                               (sleep 1000000)))
                                           :name "submit"))))))))
    (jupyter-widgets:observe *workspace* :job-name (lambda (instance type name old-value new-value source)
                                               (setf (w:widget-value (input-widget job-name)) new-value)))
    (make-instance 'w:v-box
                   :children (append (list (hbox distributor) (hbox job-name) ssh-key-file)
                                     (list (make-instance 'w:h-box
                                                          :layout (box-layout)
                                                          :children (list go-button))
                                           (make-instance 'w:h-box
                                                          :layout (box-layout)
                                                          :children (list messages)))))))


(defun check-calc ()
  (let* ((job-name (make-simple-input "Job name" :default (job-name *workspace*)))
         (distributor (make-simple-input "Distributor" :default (distributor *workspace*)))
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



(defun composer ()
  (let ((instance (structure-editor:composer "Draw Ligands" "Parse Ligands")))
    (w:link *workspace* :composer-data instance :data t)
    (w:observe instance :aggregate
               (lambda (inst type name old-value new-value source)
                 (declare (ignore (inst type name old-value source)))
                 ;; tirun-ligands are what goes in (ligands *workspace*)
                 (multiple-value-bind (molecules tirun-ligands)
                                      (tirun:assemble-ligands new-value :verbose t)
                   (setf all-ligands (tirun:pose-molecules-using-similarity molecules template-ligand)
                         selected-ligands all-ligands
                         all-edges nil)
                   (save-workspace))))
    (make-view-ligand-page instance "View Ligands")
    instance))



(defun load-chemdraw (filename)
  "Like composer - but gets the structure sketch from a chemdraw file"
  (let ((sketch (tirun:load-chem-draw-tirun filename)))
    (let ((molecules (tirun::assemble-ligands sketch)))
      (let ((posed-molecules (tirun:pose-molecules-using-similarity molecules (template-ligand *workspace*))))
        (setf all-ligands posed-molecules
              selected-ligands all-ligands)))))
