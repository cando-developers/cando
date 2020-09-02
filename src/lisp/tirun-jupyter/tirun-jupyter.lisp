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
    :trait :list)
   (receptor-string
    :accessor receptor-string
    :initform nil
    :trait :string)
   (loaded-ligands
    :accessor loaded-ligands
    :initform nil
    :trait :list)
   (selected-ligands
    :accessor selected-ligands
    :initform nil
    :trait :list)
   (ligand-widget
    :accessor ligand-widget
    :initform nil)
   (all-nodes
    :accessor all-nodes
    :initform nil
    :trait :list)
   (all-edges
    :accessor all-edges
    :initform nil
    :trait :list)
   (smirks
    :accessor smirks
    :initform ""
    :trait :string)
   (distributor
    :accessor distributor
    :initform "s103.thirdlaw.tech"
    :trait :string)
   (job-name
    :accessor job-name
    :initform "default"
    :trait :string)
   (submit-stream
    :accessor submit-stream
    :initform nil
    :trait :stream)
   (cyto-observe
    :accessor cyto-observe
    :initform nil
    :trait :list)
   (cyto-widget
    :accessor cyto-widget
    :initform nil
    :trait :list)
   (smirk-pattern
    :accessor smirk-pattern
    :initform nil
    :trait :string)
   )
  (:metaclass jupyter-widgets:trait-metaclass))



(defun composer-save (composer app)
  (format t "In composer-save About to parse ~d~%" (random 1000))
  ;; This is where we put the logic to build the molecules?
  (let ((agg (structure-editor:parse-kekule-json (composer-json app))))
    (format t "agg2 -> ~s~%" agg)
    (let* ((assembly (tirun::build-assembly agg))
           (ligands (tirun::build-tiruns assembly)))
      (format t "ligands -> ~s~%" ligands)
      (let ((molecules-svg-str (with-output-to-string (sout)
                                 (loop for ligand in ligands
                                       for molecule = (tirun::drawing ligand)
                                       for sketch2d = (sketch2d:sketch2d molecule)
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
         
(defvar *fu* nil)
(defvar *nglv* nil)

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
    

(defun widget-visibility (widget visibility)
  (setf (w:widget-display (w:widget-layout widget)) visibility))
(defun widget-hide (widget)
  (widget-visibility widget "none"))
(defun widget-show (widget)
  (widget-visibility widget ""))


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


(defclass receptor-loader (loader)
  ((ngl
     :reader receptor-loader-ngl
     :initform (make-instance 'nglv:nglwidget
                              :layout (make-instance 'w:layout
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :grid-area "ngl")))
   (calc
     :reader receptor-loader-calc
     :initform *tirun*)))


(defmethod loader-parse ((instance receptor-loader) data)
  (nglview:remove-all-components (receptor-loader-ngl instance))
  (handler-case
      (let ((as-text (babel:octets-to-string data)))
        (list as-text
              (with-input-from-string (sin as-text)
                (leap.pdb:load-pdb-stream sin))))
    (leap.pdb:pdb-read-error (condition)
      (write-string (leap.pdb:messages condition))
      nil)
    (error (condition)
      (print condition *error-output*)
      nil)))


(defmethod loader-show ((instance receptor-loader) data)
  (destructuring-bind (as-text agg)
                      data
    (nglview:add-structure (receptor-loader-ngl instance)
    (make-instance 'nglview:text-structure :text as-text))
    (nglview:handle-resize (receptor-loader-ngl instance))
    (setf (receptor-string *app*) as-text)
    (setf (tirun:receptors (receptor-loader-calc instance)) (list agg))))


(defmethod initialize-instance :after ((instance receptor-loader) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (setf (w:widget-accept (loader-upload-button instance)) ".pdb")
  (setf (w:widget-%titles (loader-accordion instance)) (list "Structure Load" "Structure View"))
  (setf (w:widget-children (loader-view-grid instance))
        (list (receptor-loader-ngl instance)))
  (setf (w:widget-layout (loader-view-grid instance))
        (make-instance 'w:layout
                       :grid-template-rows "1fr"
                       :grid-template-columns "1fr"
                       :grid-template-areas "\"ngl\"")))


(defun load-receptor (&optional (calc *tirun*))
  (loader-accordion (make-instance 'receptor-loader :calc calc)))


(defclass ligand-loader (loader)
  ((ngl
     :reader ligand-loader-ngl
     :initform (make-instance 'nglv:nglwidget
                              :layout (make-instance 'w:layout
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :grid-area "ngl")))
   (structure
     :reader ligand-loader-structure
     :initform (make-instance 'w:html
                              :layout (make-instance 'w:layout
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :grid-area "structure")))
   (slider
     :reader ligand-loader-slider
     :initform (make-instance 'w:int-slider
                              :layout (make-instance 'w:layout
                                                     :width "100%")))
   (dropdown
     :reader ligand-loader-dropdown
     :initform (make-instance 'w:dropdown
                              :layout (make-instance 'w:layout
                                                     :width "max-content")))
   (receptor-toggle
     :reader ligand-loader-receptor-toggle
     :initform (make-instance 'w:toggle-button
                              :description "Show Receptor"
                              :value t
                              :layout (make-instance 'w:layout
                                                     :margin "auto"
                                                     :grid-area "ngl-ctl")))
   (calc
     :reader ligand-loader-calc
     :initform *tirun*)))


(defvar *ligands-string* nil)
(defvar *smallest-ligand-sketch*)


(defun on-ligand-select (instance index &optional previous)
  (with-slots (ngl structure calc)
              instance
    (let* ((mol (elt (selected-ligands *app*) index))
           (id (symbol-name (chem:get-name mol))))
      (setf (w:widget-value structure)
            (format nil "<div style='display:flex;align-items:center;height:100%;'><div style='display:block;margin:auto;'>~A</div></div>"
                    (sketch2d:render-svg-to-string (sketch2d:svg (sketch2d:similar-sketch2d mol *smallest-ligand-sketch*)))))
      (nglview:handle-resize ngl)
      (when previous
        (nglview:hide-components ngl (symbol-name (chem:get-name (elt (selected-ligands *app*) previous)))))
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


(defmethod loader-parse ((instance ligand-loader) data)
  (nglview:remove-all-components (ligand-loader-ngl instance))
  (handler-case
      (let ((as-text (babel:octets-to-string data)))
        (list as-text
              (with-input-from-string (sin as-text)
                (sdf:parse-sdf-file sin))))
    (error (condition)
      (print condition *error-output*)
      nil)))


(defmethod loader-show ((instance ligand-loader) data)
  (with-slots (slider dropdown ngl calc)
              instance
    (destructuring-bind (as-text ligands)
                        data
      (setf *ligands-string* as-text
            (loaded-ligands *app*) ligands
            (selected-ligands *app*) ligands
            *smallest-ligand-sketch* (build-prototype-sketch ligands)
            (w:widget-description slider) (format nil "~A ligands" (length ligands))
            (w:widget-max slider) (1- (length ligands))
            (w:widget-%options-labels dropdown) (mapcar (lambda (mol)
                                                          (symbol-name (chem:get-name mol)))
                                                        ligands)
            (w:widget-value slider) 0)
      (nglview:add-structure (ligand-loader-ngl instance)
                             (make-instance 'nglview:text-structure
                                            :id "receptor"
                                            :text (receptor-string *app*)))
      (nglview:handle-resize ngl)
      (on-ligand-select instance 0)
      (ignore-errors
        (tirun:tirun-calculation-from-ligands calc ligands)))))


(defmethod initialize-instance :after ((instance ligand-loader) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (jupyter-widgets:link (ligand-loader-slider instance) :value
                        (ligand-loader-dropdown instance) :index)
  (setf (w:widget-accept (loader-upload-button instance)) ".sdf")
  (setf (w:widget-%titles (loader-accordion instance)) (list "Ligand Load" "Ligand View"))
  (setf (w:widget-children (loader-view-grid instance))
        (list (ligand-loader-ngl instance)
              (ligand-loader-structure instance)
              (make-instance 'w:h-box
                             :children (list (ligand-loader-slider instance)
                                             (ligand-loader-dropdown instance))
                             :layout (make-instance 'w:layout
                                                    :width "90%" :margin "auto"
                                                    :grid-area "structure-ctl"))
              (ligand-loader-receptor-toggle instance)))
  (setf (w:widget-layout (loader-view-grid instance))
        (make-instance 'w:layout
                       :grid-gap "1em"
                       :grid-template-rows "1fr min-content"
                       :grid-template-columns "1fr 1fr"
                       :grid-template-areas "\"structure ngl\" \"structure-ctl ngl-ctl\""))
  (w:observe (ligand-loader-receptor-toggle instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (nglview:set-visibility (ligand-loader-ngl instance) new-value "receptor")))
  (w:observe (ligand-loader-slider instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (on-ligand-select instance new-value old-value))))


(defun load-ligands (&optional (calc *tirun*)); (show-receptor t))
  (loader-accordion (make-instance 'ligand-loader :calc calc)))


(defun select-ligands (&key (calc *tirun*))
  "select-ligands takes the ligands in (loaded-ligands *app*) and lets the user downselect a
subset that it then puts into (selected-ligands *app*)"
  (labels ((get-loaded-names ()
             (let* ((names (loop for lig in (loaded-ligands *app*)
                                 collect (string (chem:get-name lig))))
                    (sorted-names (sort names #'string<)))
               sorted-names)))
    (let* ((sorted-names (get-loaded-names))
           (messages (make-instance 'w:text-area))
           (multi-select (make-instance 'jupyter-widgets:select-multiple
                                        :%options-labels (get-loaded-names)))
           (vbox (make-instance 'w:v-box :children (list multi-select messages))))
      (w:observe *app* :loaded-ligands (lambda (&rest args)
                                         (let ((new-sorted-names (get-loaded-names)))
                                           (setf (w:widget-%options-labels multi-select) new-sorted-names
                                                 sorted-names new-sorted-names))))
      (w:observe multi-select :index (lambda (instance type name old-value new-value source)
                                       (setf (w:widget-value messages) (format nil "Starting multi-select index sorted-names: ~a" sorted-names))
                                       (let ((mols (if sorted-names
                                                       (loop for index in new-value
                                                             for name = (elt sorted-names index)
                                                             for name-key = (intern name :keyword)
                                                             collect (find name-key (loaded-ligands *app*) :key #'chem:get-name))
                                                       (loaded-ligands *app*))))
                                         (setf (selected-ligands *app*) mols)
                                         (setf (w:widget-value messages)
                                               (format nil "Molecules:~% ~a~%" (if mols (mapcar #'chem:get-name mols) "all"))))))
      vbox)))

(defun generate-nodes-and-edges (multigraph)
  (let (all-nodes all-edges)
    (loop for subgraph in (lomap::subgraphs multigraph)
          do (let ((nodes (loop for vertex in (lomap::vertices subgraph)
                                for mol = (lomap:molecule vertex)
                                for name = (chem:get-name mol)
                                collect name))
                   (edges (loop for edge in (lomap::edges subgraph)
                                for vertex1 = (lomap:vertex1 edge)
                                for vertex2 = (lomap:vertex2 edge)
                                for score = (lomap:sim-score edge)
                                for name1 = (chem:get-name (lomap:molecule vertex1))
                                for name2 = (chem:get-name (lomap:molecule vertex2))
                                collect (list name1 name2 score))))
               (setf all-nodes (append all-nodes nodes)
                     all-edges (append all-edges edges))))
    (values all-nodes all-edges)))

(defun on-selection (instance name type old-value new-value source)
  (cond
    ((string= "nodes" (cytoscape:group instance))
     (let* ((name-str (cdr (assoc "id" (cytoscape:data instance) :test #'string=)))
            (name-sym (intern name-str :keyword))
            (mol (find name-sym (selected-ligands *app*) :key #'chem:get-name)))
       (setf (jupyter-widgets:widget-value (cyto-observe *app*))
             (format nil "~a~%~a"
                     name-str
                     (sketch2d:render-svg-to-string (sketch2d:svg (sketch2d:similar-sketch2d mol *smallest-ligand-sketch*)))))))
    ((string= "edges" (cytoscape:group instance))
     (let* ((source-str (cdr (assoc "source" (cytoscape:data instance) :test #'string=)))
            (source-sym (intern source-str :keyword))
            (target-str (cdr (assoc "target" (cytoscape:data instance) :test #'string=)))
            (target-sym (intern target-str :keyword))
            (source-mol (find source-sym (selected-ligands *app*) :key #'chem:get-name))
            (target-mol (find target-sym (selected-ligands *app*) :key #'chem:get-name)))
       (setf (jupyter-widgets:widget-value (cyto-observe *app*))
             (format nil "~a~%~a~%~a~%~a"
                     source-str
                     (sketch2d:render-svg-to-string (sketch2d:svg (sketch2d:similar-sketch2d source-mol *smallest-ligand-sketch*)))
                     (sketch2d:render-svg-to-string (sketch2d:svg (sketch2d:similar-sketch2d target-mol *smallest-ligand-sketch*)))
                     target-str))))))

(defun cyto-fill-graph (&optional (all-nodes (all-nodes *app*)) (all-edges (all-edges *app*)))
  (let ((all-nodes (loop for name in all-nodes
                         collect (make-instance 'cytoscape:element
                                                :group "nodes"
                                                :data (list (cons "id" (string name)) (cons "label" (string name))))))
        (all-edges (loop for edge in all-edges
                         collect (destructuring-bind (name1 name2 score)
                                     edge
                                   (make-instance 'cytoscape:element
                                                  :group "edges"
                                                  :data (list (cons "source" (string name1))
                                                              (cons "target" (string name2))
                                                              (cons "label" (format nil "~3,2f" score))))))))
    (setf (cytoscape:elements (cyto-widget *app*)) nil)
    (setf (cytoscape:elements (cyto-widget *app*)) (append all-nodes all-edges))
    (let ((closure (lambda (&rest args)
                     (apply 'on-selection args))))
      (dolist (instance (cytoscape:elements (cyto-widget *app*)))
        (jupyter-widgets:observe instance :selected closure)))))


(defun cyto-layout-graph ()
  (let* ((observe (make-instance 'jupyter-widgets:html
                                 :layout (make-instance 'w:layout
                                                        :border "solid"
                                                        :height "auto")))
         (cyto-widget (make-instance 'cytoscape:cytoscape-widget
                                     :graph-layouts (list (make-instance 'cytoscape:cose-layout))
                                     :graph-style "* { label: data(label); font-size: 10; }"
                                     :context-menus
                                     (list
                                      (make-instance
                                       'cytoscape:context-menu 
                                       :commands (list (make-instance
                                                        'cytoscape:menu-command 
                                                        :content "add edge"
                                                        :on-select (list
                                                                    (lambda (instance id)
                                                                      (setf (jupyter-widgets:widget-value observe)
                                                                            (format nil "~A~&examine ~A" (jupyter-widgets:widget-value observe) id))))))))))
         (boxen (make-instance 'w:v-box
                               :children (list cyto-widget observe))))
    (setf (cyto-observe *app*) observe
          (cyto-widget *app*) cyto-widget)
    boxen))

(defun lomap-progress (tirun &key (molecules (selected-ligands *app*)) advance-progress-callback )
  (let* ((mat (lomap::similarity-matrix molecules :advance-progress-callback advance-progress-callback))
         (multigraph (lomap::similarity-multigraph molecules mat)))
    (setf multigraph (lomap::lomap-multigraph multigraph :debug nil))
    (generate-nodes-and-edges multigraph)))

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

(defun lomap (&key (tirun *tirun*))
  (let* ((progress (make-instance 'w:int-progress :description "Progress")) 
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (boxen (cyto-layout-graph))
         (go-button (make-instance 'jupyter-widgets:button
                                   :description "Calculate molecular similarities"
                                   :style (button-style)
                                   :layout (make-instance 'w:layout :width "30em")
                                   :tooltip "Click me"
                                   :on-click (list
                                              (lambda (&rest args)
                                                (let ((mols (if (selected-ligands *app*)
                                                                (selected-ligands *app*)
                                                                (loaded-ligands *app*))))
                                                  (setf (w:widget-value messages) (format nil "Calculating lomap on ~a" mols))
                                                  (multiple-value-bind (nodes edges)
                                                      (lomap-progress tirun
                                                                      :molecules mols
                                                                      :advance-progress-callback (let ((index 0))
                                                                                                   (lambda (max-progress)
                                                                                                     (core:check-pending-interrupts)
                                                                                                     (setf (w:widget-value progress) index
                                                                                                           (w:widget-max progress) (1- max-progress))
                                                                                                     (incf index))))
                                                    (setf (w:widget-value messages) (format nil "Done lomap on ~a" mols))
                                                    (widget-show boxen)
                                                    (loop for element in (cytoscape:elements (cyto-widget *app*))
                                                          do (setf (cytoscape:removed element) t))
                                                    (setf (all-nodes *app*) nodes
                                                          (all-edges *app*) edges)
                                                    (cyto-fill-graph)
                                                    (loop for element in (cytoscape:elements (cyto-widget *app*))
                                                          do (setf (cytoscape:removed element) nil)))))))))
    (widget-hide boxen)
    (make-instance 'w:v-box
                   :children (list
                              (make-instance 'w:h-box
                                             :layout (box-layout)
                                             :children (list go-button))
                              (make-instance 'w:h-box
                                             :layout (box-layout)
                                             :children (list progress))
                              (make-instance 'w:h-box
                                             :layout (box-layout)
                                             :children (list messages))
                              boxen))))


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


