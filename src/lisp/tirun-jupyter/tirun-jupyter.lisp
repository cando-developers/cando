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
   (loaded-ligands
    :accessor loaded-ligands
    :initform nil
    :trait t)
   (selected-ligands
    :accessor selected-ligands
    :initform nil
    :trait t)
   (similarity-matrix
    :accessor similarity-matrix
    :initform nil
    :trait t)
   (all-nodes
    :accessor all-nodes
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
   (smirk-pattern
    :accessor smirk-pattern
    :initform nil
    :trait t))
  (:metaclass jupyter-widgets:trait-metaclass))



(defun composer-save (composer app)
  (format t "In composer-save About to parse ~d~%" (random 1000))
  ;; This is where we put the logic to build the molecules?
  (let ((agg (structure-editor:parse-kekule-json (composer-json app))))
    (format t "agg2 -> ~s~%" agg)
    (let* ((assembly (tirun::build-assembly agg))
           (ligand-molecules (tirun::build-ligand-molecules assembly)))
      (format t "ligand-molecules -> ~s~%" ligand-molecules)
      (setf (loaded-ligands *app*) ligand-molecules)
      (let ((molecules-svg-str (with-output-to-string (sout)
                                 (loop for molecule in ligand-molecules
                                       for sketch2d = (sketch2d:sketch2d molecule)
                                       for __1 = (format t "composer ~a~%" (random 100))
                                       for svg = (sketch2d:svg sketch2d)
                                       for svg-str = (sketch2d:render-svg-to-string svg)
                                       do (format sout "~a" svg-str)))))
        (setf (w:widget-value (structure-editor:composer-view-grid composer)) molecules-svg-str))
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

#|
(defmethod initialize-instance ((instance app) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
#|  ; Wire up the observers
  (jupyter-widgets:observe instance :loaded-ligands (lambda ....))
  ; Or use link
  (jupyter-widgets:link instance :loaded-ligands other-instance :trait-name)
|#)
|#
   
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


(defun on-ligand-select (instance index &optional previous)
  "Select a ligand on a view-ligand-page and hide a previous ligand if one was selected."
  (with-slots (ngl structure)
              instance
    (let* ((mol (elt (loaded-ligands *app*) index))
           (id (symbol-name (chem:get-name mol))))
      ;; Load the sketch
      (setf (w:widget-value structure)
            (format nil "<div style='display:flex;align-items:center;height:100%;'><div style='display:block;margin:auto;'>~A</div></div>"
                    (sketch2d:render-svg-to-string (cw:sketch-molecule mol))))
      ;; A resize can't hurt.
      (nglview:handle-resize ngl)
      ;; Hide previous ligand.
      (when previous
        (nglview:hide-components ngl (symbol-name (chem:get-name (elt (loaded-ligands *app*) previous)))))
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
      (when (equal new-value 1)
        (nglview:handle-resize (view-ligand-ngl instance)))))
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
  ;; If the loaded-ligands changes then reload the whole thing.
  (w:observe *app* :loaded-ligands
    (lambda (inst type name old-value ligands source)
      (declare (ignore inst type name old-value source))
      (with-slots (slider dropdown)
                  instance
        (setf (w:widget-description slider) (format nil "~A ligands" (length ligands))
              (w:widget-max slider) (1- (length ligands))
              (w:widget-%options-labels dropdown) (mapcar (lambda (mol)
                                                            (symbol-name (chem:get-name mol)))
                                                          ligands)
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
      (with-slots (ligands-string loaded-ligands selected-ligands)
                  *app*
        (setf ligands-string (babel:octets-to-string (cdr (assoc "content" (car parameter) :test #'string=)))
              loaded-ligands (with-input-from-string (input-stream ligands-string)
                               (sdf:parse-sdf-file input-stream))
              selected-ligands loaded-ligands)
        t)
    (sdf:sdf-parse-error (condition)
      (princ condition *error-output*)
      nil)
    (end-of-file (condition)
      (princ condition *error-output*)
      nil)))


(defun load-ligands ()
  "Present a graphical interface to load a collection of ligands from a file into the current
  application and then display the ligand structure using two-dimension diagrams and three
  dimension representation."
  (let ((container (make-instance 'w:accordion :selected-index 0)))
    (cw:make-file-task-page container "Load Ligands" #'parse-ligands :accept ".sdf")
    (make-view-ligand-page container "View Ligands")
    container))


(defun select-ligands ()
  "select-ligands in the current application and displays a molecule selection list that allows one
  to select a subset of the molecules to use in further calculations."
  (let ((sel (make-instance 'cw:molecule-select
                            :molecules (loaded-ligands *app*)
                            :selected (selected-ligands *app*))))
    (w:link sel :molecules *app* :loaded-ligands)
    (w:link sel :selected *app* :selected-ligands)
    (make-instance 'w:accordion
                   :selected-index nil
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


(defun pasrun (&key (calc *tirun*))
  "pasrun loads a ligand template and runs PAS on it to generate multiple ligands.
It will put those multiple ligands into (loaded-ligands *app*) and (selected-ligands *app*)"
  (let* (template-sketch2d
         template-molecule
         (file-upload-msg (make-instance 'w:label :value "PAS ligand template: "))
         (file-upload (make-instance 'jupyter-widgets:file-upload :description "ligand to run PAS on"
                                                                  :style (button-style)
                                                                  :width "auto"
                                                                  :layout (make-instance 'w:layout :width "auto")))
         (file-upload-hbox (make-instance 'w:h-box
                                          :layout (box-layout)
                                          :children (list file-upload-msg file-upload)))
         (template-view (make-instance 'w:html
                                       :layout (make-instance 'w:layout :flex "1 1 0%" :width "auto" :height "auto")))
         (smirks-pattern-label (make-instance 'w:label :value "SMIRKS pattern:"))
         (smirks-box (make-instance 'w:text :value (smirks *app*)))
         (messages (make-instance 'w:text-area))
         (structure-view (make-instance 'w:html
                                        :layout (make-instance 'w:layout :flex "1 1 0%" :width "auto" :height "auto")))
         (nglview (make-instance 'nglv:nglwidget
                                 :layout (make-instance 'w:layout :width "auto" :height "auto")))
         (Ligand-selector (make-instance 'w:bounded-int-text :description "Ligand"))
         (ligand-total (make-instance 'w:label :value ""))
         (ligand-name (make-instance 'w:label :value "Molecule"))
         (ligand-h-box (make-instance 'w:h-box :children (list ligand-selector ligand-total ligand-name)))
         (go-button (make-instance 'jupyter-widgets:button
                                   :description "Run PAS"
                                   :style (button-style)
                                   :layout (make-instance 'w:layout :width "30em")
                                   :tooltip "Use the smirks pattern entered above to generate a PAS"
                                   :on-click (list
                                              (lambda (&rest args)
                                                (setf (w:widget-value messages)
                                                      (format nil "~aRunning position-analogue-scanning~%"
                                                              (w:widget-value messages)))
                                                (let* ((smirks (w:widget-value smirks-box))
                                                       (lig template-molecule)
                                                       (new-mols (pas:pas smirks lig))
                                                       (_ (setf (w:widget-value messages)
                                                                (format nil "~aGot ligand ~a with smirks: ~a new-mols: ~a~%"
                                                                        (w:widget-value messages)
                                                                        lig smirks new-mols)))
                                                       #+(or)(new-mols-svg (with-output-to-string (sout)
                                                                             (loop for mol in new-mols
                                                                                   do (format sout (sketch2d:render-svg-to-string (sketch2d:svg (sketch2d:similar-sketch2d mol template-sketch2d))))))))
                                                  (setf (w:widget-max ligand-selector) (1- (length new-mols)))
                                                  (setf (selected-ligands *app*) new-mols
                                                        (loaded-ligands *app*) new-mols
                                                        (w:widget-value messages) (format nil "~aBuilt ~d molecules~%"
                                                                                          (w:widget-value messages)
                                                                                          (length new-mols)))
                                                  (tirun:tirun-calculation-from-ligands calc (selected-ligands *app*))
                                                  
                                                  #+(or)(setf (w:widget-value structure-view)
                                                              new-mols-svg)))))))
    (w:observe smirks-box :value
               (lambda (instance type name old-value new-value source)
                 (declare (ignore instance type name old-value source))
                 (setf (smirks *app*) new-value)))
    (w:observe file-upload :data
               (lambda (instance type name old-value new-value source)
                 (declare (ignore instance type name old-value source))
                 (format t "observe fu~%")
                 (when new-value
                   (let* ((last-sdf (car (last new-value)))
                          (octets (make-array (length last-sdf) :element-type 'ext:byte8
                                                                :initial-contents last-sdf)))
                     (let ((as-text (babel:octets-to-string octets)))
                       (let ((ligands (with-input-from-string (sin as-text)
                                        (sdf:parse-sdf-file sin))))
                         (setf *smallest-ligand-sketch* (build-prototype-sketch (list (first ligands))))
                         (setf template-molecule (first ligands)
                               template-sketch2d *smallest-ligand-sketch* #+(or)(sketch2d:sketch2d (first ligands)))
                         (let ((str (sketch2d:render-svg-to-string (sketch2d:svg template-sketch2d))))
                           (setf (w:widget-value template-view) str))))))))
    (w:observe ligand-selector :value (lambda (instance type name old-value new-value source)
                                        (declare (ignore instance type name old-value source))
                                        (when new-value
                                          (observe-ligand-selector nglview structure-view messages calc ligand-name new-value))))
    #||
    

    (setf (w:widget-value structure-view)
    (sketch2d:render-svg-to-string (sketch2d:svg (sketch2d:similar-sketch2d new-mol *smallest-ligand-sketch*))))))))))

||#
    #+(or)(w:observe smirks-box :value (lambda (instance type name old-value new-value source)
                                   (declare (ignore instance type name old-value source))
                                   (when new-value
                                     (setf (smirk-pattern *app*) new-value))))    
    (make-instance 'w:v-box
                   :children (list
                              file-upload-hbox
                              template-view
                              (make-instance 'w:h-box
                                             :layout (box-layout)
                                             :children (list smirks-pattern-label smirks-box))
                              (make-instance 'w:h-box
                                             :layout (box-layout)
                                             :children (list go-button))
                              messages
                              ligand-h-box
                              structure-view
                              nglview))))


(defun create-similarity-graph (graph)
  "Create a cytoscape representation of the similarity graph in *app*. Also make an edge list for
  the job list."
  (let ((multigraph (lomap::lomap-multigraph (lomap::similarity-multigraph (selected-ligands *app*)
                                                                           (similarity-matrix *app*))
                                             :debug nil))
        nodes edges all-edges)
    ;; Loop over the subgraphs and the vertices and edges and create structure in cytoscape.
    (dolist (subgraph (lomap::subgraphs multigraph))
      (dolist (vertex (lomap::vertices subgraph))
        (let* ((molecule (lomap:molecule vertex))
               (name (symbol-name (chem:get-name molecule)))
               (sketch (cw:sketch-molecule molecule)))
          ;; Save sketch dimension and the sketch itself into the data slot. The sketch is SVG text
          ;; so we just need do a percent encoded string for the reserved URL characters.
          (push (make-instance 'cytoscape:element
                               :group "nodes"
                               :data (list (cons "id" name)
                                           (cons "label" name)
                                           (cons "width" (sketch2d:width sketch))
                                           (cons "height" (sketch2d:height sketch))
                                           (cons "image" (concatenate 'string
                                                                      "data:image/svg+xml,"
                                                                      (quri:url-encode (sketch2d:render-svg-to-string sketch))))))
                nodes)))
      (dolist (edge (lomap::edges subgraph))
        (let ((name-1 (chem:get-name (lomap:molecule (lomap:vertex1 edge))))
              (name-2 (chem:get-name (lomap:molecule (lomap:vertex2 edge))))
              (score (lomap:sim-score edge)))
          ;; Save the edge for the job configuration.
          (push (list name-1 name-2 score) all-edges)
          (push (make-instance 'cytoscape:element
                               :group "edges"
                               :data (list (cons "source" (symbol-name name-1))
                                           (cons "target" (symbol-name name-2))
                                           (cons "label" (format nil "~3,2f" score))))
                edges))))
    ;; Update the graph and the edge list.
    (setf (cytoscape:elements graph) (nconc nodes edges)
          (all-edges *app*) all-edges)))


(defun run-lomap (parameter progress-callback)
  "This function is called by the task page and updates a progress widget as lomap does the calculation."
  (declare (ignore parameter))
  (format t "Calculating similarity matrix of ~{~a~#[~;, and ~:;, ~]~}...~%"
          (mapcar (lambda (molecule)
                    (symbol-name (chem:get-name molecule)))
                  (selected-ligands *app*)))
  (finish-output) ; Make sure the message is synced to the frontend.
  (setf (similarity-matrix *app*)
        (lomap::similarity-matrix (selected-ligands *app*)
                                  :advance-progress-callback progress-callback))
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


(defun lomap ()
  "Calculate molecule similarities using a graph theory approach based on LOMAP
  (J Comput Aided Mol Des 27, 755–770 (2013). https://doi.org/10.1007/s10822-013-9678-y)"
  (let* ((container (make-instance 'w:accordion :selected-index 0))
         (graph (make-instance 'cytoscape:cytoscape-widget
                               :graph-layouts (list (make-instance 'cytoscape:cose-layout))
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
                                                     :grid-template-columns "1fr min-content min-content min-content 1fr"
                                                     :grid-template-rows "1fr min-content"
                                                     :grid-template-areas "\"graph graph graph graph graph\" \". fit show-sketches full-screen .\""))))
    (cw:make-simple-task-page container "Calculate Similarities" #'run-lomap
                              :label "Click button to calculate molecular similarities.")
    (setf (w:widget-%titles container) (append (w:widget-%titles container)
                                              (list "Similarity Graph"))
          (w:widget-children container) (append (w:widget-children container)
                                               (list grid)))
    (w:observe *app* :similarity-matrix ; Create the graph on a change to the similarity matrix.
      (lambda (instance type name old-value new-value source)
        (declare (ignore instance type name old-value new-value source))
        (create-similarity-graph graph)))
    container))


(defun configure-layout (default-settings &key testing)
  (let (all-inputs)
    (multiple-value-bind (titles vboxes)
        (loop for batch in default-settings
              for title = (first batch)
              for settings = (second batch)
              collect (make-instance 'w:v-box
                                     :children
                                     (loop for setting in settings
                                           for setting-key = (first setting)
                                           for setting-name = (second setting)
                                           for setting-value = (third setting)
                                           for input = (make-simple-input setting-name
                                                                          :default (if (consp setting-value)
                                                                                       (if testing
                                                                                           (second setting-value)
                                                                                           (first setting-value))
                                                                                       setting-value)
                                                                          :name setting-key)
                                           do (push input all-inputs)
                                           collect (hbox input)))
                into vboxes
              collect title into titles
              finally (return (values titles vboxes)))
      (values (make-instance 'jupyter-widgets:accordion
                             :%titles titles
                             :selected-index nil
                             :children vboxes)
              all-inputs))))

(defun configure-jobs (&key (calc *tirun*))
  (let* ((desc-width "12em")
         (steps (make-instance 'jupyter-widgets:dropdown :%options-labels (list "single" "3-step")
                                                         :options (list "single" "3-step")
                                                         :index 0
                                                         :description "Steps:")))
    (multiple-value-bind (accordion all-inputs)
        (configure-layout tirun::*default-calculation-settings*)
      (let* ((desc-style (make-instance 'w:description-style :description-width desc-width))
             (messages (make-instance 'w:text-area :description "Messages"
                                                   :layout (make-instance 'w:layout :width "60em")))
             (go-button (make-instance 'jupyter-widgets:button :description "Configure"
                                                               :style (button-style)
                                                               :tooltip "Click me"
                                                               :on-click (list
                                                                          (lambda (&rest args)
                                                                             (setf (w:widget-value messages)
                                                                                  (format nil "Configure~%")))))))
        (make-instance 'w:v-box
                       :children (append (list steps
                                               accordion)
                                         (list 
                                          (make-instance 'w:h-box
                                                         :layout (box-layout)
                                                         :children (list go-button))
                                          (make-instance 'w:h-box
                                                         :layout (box-layout)
                                                         :children (list messages)))))))))

(defun ensure-write-jobs (tirun jobs-dir &optional progress-callback)
  (ensure-directories-exist jobs-dir)
  (let ((*default-pathname-defaults* (pathname jobs-dir)))
    (ext:chdir *default-pathname-defaults*)
    (tirun:build-job-nodes tirun)
    (tirun:connect-job-nodes tirun (all-edges *app*))
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
    (structure-editor::composer-accordion composer)))


