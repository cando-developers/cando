(in-package :tirun-jupyter)

(defvar *button-style* (make-instance 'w:button-style :button-color "aquamarine"))

(defclass app (jupyter-widgets:has-traits)
  ((receptor-string
    :accessor receptor-string
    :initform nil
    :trait :string)
   (loaded-ligands
    :accessor loaded-ligands
    :initform nil
    :trait :list)
   (ligands
    :accessor ligands
    :initform nil
    :trait :list)
   (ligand-names
    :accessor ligand-names
    :initform nil
    :trait :list)
   (all-nodes
    :accessor all-nodes
    :initform nil
    :trait :list)
   (all-edges
    :accessor all-edges
    :initform nil
    :trait :list)
   (distributor
    :accessor distributor
    :initform "s103.thirdlaw.tech"
    :trait :string)
   (job-name
    :accessor job-name
    :initform "default"
    :trait :string)
   (cyto-observe
    :accessor cyto-observe
    :initform nil
    :trait :list)
   (cyto-widget
    :accessor cyto-widget
    :initform nil
    :trait :list)



   )  
  (:metaclass jupyter-widgets:trait-metaclass))

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


(defvar *box-layout*
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


(defun widget-visibility (widget visibility)
  (setf (w:widget-display (w:widget-layout widget)) visibility))
(defun widget-hide (widget)
  (widget-visibility widget "none"))
(defun widget-show (widget)
  (widget-visibility widget ""))


(defclass simple-input ()
  ((label-widget :initarg :label-widget :accessor label-widget)
   (input-widget :initarg :input-widget :accessor input-widget)
   (hbox :initarg :hbox :accessor hbox)))

(defun make-simple-input (label &key (default ""))
  (let ((label-widget (make-instance 'w:label :value label
                                              :layout (make-instance 'w:layout :width "20em")))
        (input-widget (make-instance 'w:text
                                     :value default
                                     :layout (make-instance 'w:layout :width "20em"))))
    (make-instance 'simple-input
                   :label-widget label-widget
                   :input-widget input-widget
                   :hbox (make-instance 'w:h-box
                                        :children (list label-widget input-widget)))))

(defun load-receptor (&optional (calc *tirun*))
  (let* ((desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (file-upload-msg (make-instance 'w:label :value "Receptor structure: "))
         (file-upload (make-instance 'jupyter-widgets:file-upload :description "PDB file for receptor"
                                                                  :style *button-style*
                                                                  :width "auto"
                                                                  :layout (make-instance 'w:layout :width "auto")))
         (file-upload-vbox (make-instance 'w:h-box
                                          :layout *box-layout*
                                          :children (list file-upload-msg file-upload)))
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (nglview (make-instance 'nglv:nglwidget)))
    (setf (w:widget-value messages) "Start message")
    (labels ((setup-receptor* (new-value)
               (let* ((last-pdb (car (last new-value)))
                      (octets (make-array (length last-pdb) :element-type 'ext:byte8
                                                            :initial-contents last-pdb)))
                 (let ((as-text (babel:octets-to-string octets)))
                   (setf (receptor-string *app*) nil)
                   (nglview:remove-all-components nglview)
                   (let ((agg (handler-case
                                  (with-input-from-string (sin as-text)
                                    (leap.pdb:load-pdb-stream sin))
                                (leap.pdb:pdb-read-error (condition)
                                  (setf (w:widget-value messages) (leap.pdb:messages condition))
                                  (widget-show messages)
                                  (widget-hide nglview)
                                  (return-from setup-receptor* nil)))))
                     (nglview:add-structure nglview (make-instance 'nglview:text-structure :text as-text))
                     (setf (w:widget-value messages) "Loaded successfully")
                     (widget-hide messages)
                     (widget-show nglview)
                     (setf (receptor-string *app*) as-text)
                     (setf (tirun:receptors calc) (list agg)))))))
      (widget-show messages)
      (widget-hide nglview)
      (w:observe file-upload :data
                 (lambda (instance type name old-value new-value source)
                   (declare (ignore instance type name old-value source))
                   (when new-value
                     (setup-receptor* new-value))))
      (make-instance 'w:v-box
                     :children (list
                                (make-instance 'w:h-box
                                               :layout *box-layout*
                                               :children (list file-upload-vbox))
                                (make-instance 'w:h-box
                                               :layout *box-layout*
                                               :children (list messages))
                                (make-instance 'w:h-box
                                               :layout *box-layout*
                                               :children (list nglview)))))))


(defvar *ligands-string* nil)
(defvar *current-component* nil)
(defvar *smallest-ligand-sketch*)
(defun observe-ligand-selector (nglview structure-view calc ligand-name new-value)
  (let ((mol (elt (ligands *app*) new-value)))
    (let* ((agg (chem:make-aggregate)))
      (chem:add-matter agg mol)
      (setf (w:widget-value ligand-name) (string (chem:get-name mol)))
      (let* ((mol2 (chem:aggregate-as-mol2-string agg t))
             (new-component (nglview:add-structure nglview (make-instance 'nglview:text-structure :text mol2 :ext "mol2"))))
        (when *current-component*
          (nglview:remove-components nglview *current-component*))
        (setf (w:widget-value structure-view)
              (sketch2d:render-svg-to-string (sketch2d:svg (sketch2d:similar-sketch2d mol *smallest-ligand-sketch*))))
        (setf *current-component* new-component)
        ))))

(defvar *lv* nil)
(defun load-ligands (&key (calc *tirun*) (show-receptor t))
  (let* ((fu (make-instance 'w:file-upload :description "Upload SDF file"
                            :style *button-style*))
         (Ligand-selector (make-instance 'w:bounded-int-text :description "Ligand"))
         (ligand-total (make-instance 'w:label :value ""))
         (ligand-name (make-instance 'w:label :value "Molecule"))
         (ligand-h-box (make-instance 'w:h-box :children (list ligand-selector ligand-total ligand-name)))
         (control-box (make-instance 'w:v-box :children (list fu ligand-h-box)
                                              :layout (make-instance 'w:layout :flex "1 1 0%" :width "auto" :height "auto")))
         (structure-view (make-instance 'w:html
                                        :layout (make-instance 'w:layout :flex "1 1 0%" :width "auto" :height "auto")))
         (nglview (make-instance 'nglv:nglwidget
                                 :layout (make-instance 'w:layout :width "auto" :height "auto"))))
    (setf *lv* nglview)
    (widget-hide ligand-h-box)
    (widget-hide nglview)
    (w:observe fu :data
               (lambda (instance type name old-value new-value source)
                 (declare (ignore instance type name old-value source))
                 (format t "observe fu~%")
                 (when new-value
                   (let* ((last-pdb (car (last new-value)))
                          (octets (make-array (length last-pdb) :element-type 'ext:byte8
                                                                :initial-contents last-pdb)))
                     (let ((as-text (babel:octets-to-string octets)))
                       (setf *ligands-string* as-text)
                       (let ((ligands (with-input-from-string (sin as-text)
                                        (sdf:parse-sdf-file sin))))
                         (setf (loaded-ligands *app*) ligands
                               (ligands *app*) ligands))
                       (flet ((build-prototype-sketch (ligands)
                                (let ((min-atoms (chem:number-of-atoms (car ligands)))
                                      (min-mol (car ligands)))
                                  (loop for mol in (cdr ligands)
                                        for num-atoms = (chem:number-of-atoms mol)
                                        when (< num-atoms min-atoms)
                                          do (setf min-atoms num-atoms
                                                   min-mol mol))
                                  (sketch2d:sketch2d min-mol))))
                         (setf *smallest-ligand-sketch* (build-prototype-sketch (ligands *app*))))
                       (setf (w:widget-value ligand-selector) 0
                             (w:widget-value ligand-total) (format nil "~a" (length (loaded-ligands *app*)))
                             (w:widget-max ligand-selector) (1- (length (loaded-ligands *app*)))
                             (w:widget-min ligand-selector) 0)
                       (widget-show ligand-h-box)
                       (widget-show nglview)
                       (when show-receptor
                         (nglview:add-structure nglview (make-instance 'nglview:text-structure :text (receptor-string *app*))))
                       (observe-ligand-selector nglview structure-view calc ligand-name 0)
                       (tirun:tirun-calculation-from-ligands calc (ligands *app*))
                       (format t "Done observe fu~%"))))))
    (w:observe ligand-selector :value (lambda (instance type name old-value new-value source)
                                        (declare (ignore instance type name old-value source))
                                        (when new-value
                                          (observe-ligand-selector nglview structure-view calc ligand-name new-value))))
    (make-instance 'w:v-box
                   :children (list
                              (make-instance 'w:h-box
                                             :layout *box-layout*
                                             :children (list control-box structure-view))
                              #+(or)(make-instance 'w:h-box
                                                   ;;                                             :layout *box-layout*
                                                   :children (list ligand-h-box))
                              (make-instance 'w:h-box
;;;                                             :layout *box-layout*
                                            #+(or):layout #+(or)(make-instance 'w:layout
                                                                          :width "100em"
                                                                          :flex-flow "row"
                                                                          :align-items "stretch")
                                             :children (list nglview))
                              ))))


(defun select-ligands (&key (calc *tirun*))
  (labels ((get-loaded ()
             (let* ((names (loop for lig in (loaded-ligands *app*)
                                 collect (string (chem:get-name lig))))
                    (sorted-names (sort names #'string<)))
               sorted-names)))
    (let* ((messages (make-instance 'w:text-area))
           (multi-select (make-instance 'jupyter-widgets:select-multiple
                                        :%options-labels (get-loaded)))
           (vbox (make-instance 'w:v-box :children (list multi-select messages))))
      (w:observe *app* :loaded-ligands (lambda (&rest args)
                                         (let ((sorted-names (get-loaded)))
                                           (setf (ligand-names *app*) sorted-names
                                                 (w:widget-%options-labels multi-select) sorted-names))))
      (w:observe multi-select :index (lambda (instance type name old-value new-value source)
                                       (setf (w:widget-value messages) (format nil "Starting multi-select index ligand-names: ~a" (ligand-names *app*)))
                                       (let ((mols (if (ligand-names *app*)
                                                       (loop for index in new-value
                                                             for name = (elt (ligand-names *app*) index)
                                                             for name-key = (intern name :keyword)
                                                             collect (find name-key (loaded-ligands *app*) :key #'chem:get-name))
                                                       (loaded-ligands *app*))))
                                         (setf (ligands *app*) mols)
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
            (mol (find name-sym (ligands *app*) :key #'chem:get-name)))
       (setf (jupyter-widgets:widget-value (cyto-observe *app*))
             (format nil "~a~%~a"
                     name-str
                     (sketch2d:render-svg-to-string (sketch2d:svg (sketch2d:similar-sketch2d mol *smallest-ligand-sketch*)))))))
    ((string= "edges" (cytoscape:group instance))
     (let* ((source-str (cdr (assoc "source" (cytoscape:data instance) :test #'string=)))
            (source-sym (intern source-str :keyword))
            (target-str (cdr (assoc "target" (cytoscape:data instance) :test #'string=)))
            (target-sym (intern target-str :keyword))
            (source-mol (find source-sym (ligands *app*) :key #'chem:get-name))
            (target-mol (find target-sym (ligands *app*) :key #'chem:get-name)))
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
    (setf (cytoscape:elements (cyto-widget *app*)) (append all-nodes all-edges)
          (cytoscape:graph-style (cyto-widget *app*))
          (append (cytoscape:graph-style (cyto-widget *app*))
                  (list (jupyter:json-new-obj
                         ("selector" "node[label]")
                         ("style" (jupyter:json-new-obj
                                   ("label" "data(label)")
                                   ("font-size" "10"))))
                        (jupyter:json-new-obj
                         ("selector" "edge[label]")
                         ("style" (jupyter:json-new-obj
                                   ("label" "data(label)")))))))
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
                                     :graph-layout (list (cons "name" "cose")
                                                         (cons "quality" "default"))
                                     :context-menus
                                     (list
                                      (make-instance
                                       'cytoscape:context-menu 
                                       :commands (list (make-instance
                                                        'cytoscape:menu-command 
                                                        :content "add edge"
                                                        :on-select (list
                                                                    (lambda (instance id)
                                                                      (jupyter-widgets:widget-value
                                                                       (setf (jupyter-widgets:widget-value observe)
                                                                             (format nil "~A~&fu on ~A" (jupyter-widgets:widget-value observe) id)))))
                                                        )))
                                      (make-instance 'cytoscape:context-menu 
                                                     :selector "edge"
                                                     :commands (list (make-instance
                                                                      'cytoscape:menu-command
                                                                      :content "remove"
                                                                      :on-select (list
                                                                                  (lambda (instance id)
                                                                                    (setf (jupyter-widgets:widget-value observe)
                                                                                          (format nil "~A~&remove on ~A" (jupyter-widgets:widget-value observe) id)))))
                                                                     (make-instance 'cytoscape:menu-command
                                                                                    :content "examine"
                                                                                    :on-select (list
                                                                                                (lambda (instance id)
                                                                                                  (setf (jupyter-widgets:widget-value observe)
                                                                                                        (format nil "~A~&examine ~A" (jupyter-widgets:widget-value observe) id))))))))))
         (boxen (make-instance 'w:v-box
                               :children (list cyto-widget observe))))
    (setf (cyto-observe *app*) observe
          (cyto-widget *app*) cyto-widget)
    boxen))

(defun lomap-progress (tirun &key (molecules (ligands *app*)) advance-progress-callback )
  (let* ((mat (lomap::similarity-matrix molecules :advance-progress-callback advance-progress-callback))
         (multigraph (lomap::similarity-multigraph molecules mat)))
    (setf multigraph (lomap::lomap-multigraph multigraph :debug nil))
    (generate-nodes-and-edges multigraph)))

(defun lomap (&key (tirun *tirun*))
  (let* ((progress (make-instance 'w:int-progress :description "Progress")) 
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (boxen (cyto-layout-graph))
         (go-button (make-instance 'jupyter-widgets:button
                                   :description "Calculate molecular similarities"
                                   :style *button-style*
                                   :layout (make-instance 'w:layout :width "30em")
                                   :tooltip "Click me"
                                   :on-click (list
                                              (lambda (&rest args)
                                                (let ((mols (if (ligands *app*)
                                                                (ligands *app*)
                                                                (loaded-ligands *app*))))
                                                  (setf (w:widget-value messages) (format nil "Calculating lomap on ~a" mols))
                                                  (multiple-value-bind (nodes edges)
                                                      (lomap-progress tirun
                                                                      :molecules mols
                                                                      :advance-progress-callback (let ((index 0))
                                                                                                   (lambda (max-progress)
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
                                             :layout *box-layout*
                                             :children (list go-button))
                              (make-instance 'w:h-box
                                             :layout *box-layout*
                                             :children (list progress))
                              (make-instance 'w:h-box
                                             :layout *box-layout*
                                             :children (list messages))
                              boxen))))

(defun configure-jobs (&key (calc *tirun*))
  (let* ((desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (simulation-length (make-instance 'jupyter-widgets:text :description "Simulation length"))
         (lambda-windows (make-instance 'jupyter-widgets:text :description "Lambda windows"))
         (settings (make-instance 'jupyter-widgets:v-box ;;; :layout *box-layout*
                                                         :children (list simulation-length
                                                                         lambda-windows)))
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (go-button (make-instance 'jupyter-widgets:button :description "Configure"
                                                           :style *button-style*
                                                           :tooltip "Click me"
                                                           :on-click (list
                                                                      (lambda (&rest args)
                                                                        (setf (w:widget-value messages)
                                                                              (format nil "Configure: ~a~%" args)))))))
    (setf *job-messages* messages)
    (make-instance 'w:v-box
                   :children (list
                              (make-instance 'w:h-box
                                             :layout *box-layout*
                                             :children (list settings))
                              (make-instance 'w:h-box
                                             :layout *box-layout*
                                             :children (list go-button))
                              (make-instance 'w:h-box
                                             :layout *box-layout*
                                             :children (list messages))))))

(defun ensure-write-jobs (tirun jobs-dir progress-callback)
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

(defun write-jobs (&key (tirun *tirun*))
  (let* ((jobs-dir (pathname (if (ext:getenv "TIRUN_DIRECTORY")
                                 (uiop:ensure-directory-pathname (ext:getenv "TIRUN_DIRECTORY"))
                                 (merge-pathnames (make-pathname :directory (list :relative "jobs"))
                                                  (uiop:ensure-directory-pathname (ext:getenv "HOME"))))))
         (job-name (make-simple-input "Job name" :default (job-name *app*)))
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (progress (make-instance 'w:int-progress :description "Progress"))
         (go-button (make-instance 'w:button :description "Write jobs"
                                             :style *button-style*
                                             :on-click (list
                                                        (lambda (&rest args)
                                                          (let* ((name (w:widget-value (input-widget job-name)))
                                                                 (dir (merge-pathnames (make-pathname :directory (list
                                                                                                                  :relative
                                                                                                                  name))
                                                                                       jobs-dir)))
                                                            (setf (w:widget-value messages) (format nil "Starting to write jobs dir-> ~a" dir))
                                                            (if (probe-file dir)
                                                                (progn
                                                                  (setf (w:widget-value messages) (format nil "The directory ~a already exists - aborting writing jobs" dir)))
                                                                (progn
                                                                  (ensure-directories-exist dir)
                                                                  (ensure-write-jobs tirun dir (let ((index 0)
                                                                                                     (lock (bordeaux-threads:make-recursive-lock "progress")))
                                                                                                 (lambda (max-progress)
                                                                                                   (bordeaux-threads:with-recursive-lock-held (lock)
                                                                                                     (setf (w:widget-value progress) index
                                                                                                           (w:widget-max progress) max-progress)
                                                                                                     (incf index)))))
                                                                  (setf (w:widget-value messages) "Done writing jobs"))))))))
         )
    (w:observe (input-widget job-name) :value (lambda (instance type name old-value new-value source)
                                                (setf (job-name *app*) new-value)))
    (make-instance 'w:v-box
                   :children (list
                              (hbox job-name)
                              (make-instance 'w:h-box
                                             :layout *box-layout*
                                             :children (list go-button))
                              (make-instance 'w:h-box
                                             :children (list progress))
                              (make-instance 'w:h-box
                                             :layout *box-layout*
                                             :children (list messages))))))
                              
                                             
    



(defvar *job-messages*)
(defun submit-calc (&key (calc *tirun*))
  (let* ((desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (distributor (make-simple-input "Distributor" :default (distributor *app*)))
         (job-name (make-simple-input "Job name" :default "default"))
         (ssh-key-file (make-instance 'jupyter-widgets:file-upload :description "SSH key file"
                                                                   :style *button-style*))
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (go-button (make-instance 'jupyter-widgets:button :description "Submit jobs"
                                                           :style *button-style*
                                                           :tooltip "Click me"
                                                           :on-click (list
                                                                      (lambda (&rest args)
                                                                        (setf (w:widget-value messages)
                                                                              (format nil "Submitting job ~a on ~a~%"
                                                                                      (w:widget-value (input-widget job-name))
                                                                                      (w:widget-value (input-widget distributor)))))))))
    (jupyter-widgets:observe *app* :job-name (lambda (instance type name old-value new-value source)
                                               (setf (w:widget-value (input-widget job-name)) new-value)))
    (make-instance 'w:v-box
                   :children (append (list (hbox distributor) (hbox job-name) ssh-key-file)
                                     (list (make-instance 'w:h-box
                                                          :layout *box-layout*
                                                          :children (list go-button))
                                           (make-instance 'w:h-box
                                                          :layout *box-layout*
                                                          :children (list messages)))))))


(defun check-calc (&key (tirun *tirun*))
  (let* ((job-name (make-simple-input "Job name" :default (job-name *app*)))
         (distributor (make-simple-input "Distributor" :default (distributor *app*)))
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (go-button (make-instance 'w:button :description "Check calculation"
                                             :style *button-style*
                                             :on-click (list
                                                        (lambda (&rest args)
                                                          (setf (w:widget-value messages)
                                                                (format nil "Checking job ~a on ~a~%"
                                                                        (w:widget-value (input-widget job-name))
                                                                        (w:widget-value (input-widget distributor)))))))))
    (jupyter-widgets:observe *app* :job-name (lambda (instance type name old-value new-value source)
                                               (setf (w:widget-value (input-widget job-name)) new-value)))
    (make-instance 'w:v-box
                   :children (list
                              (hbox distributor)
                              (hbox job-name)
                              (make-instance 'w:h-box
                                             :layout *box-layout*
                                             :children (list go-button))
                              (make-instance 'w:h-box
                                             :layout *box-layout*
                                             :children (list messages))))))
