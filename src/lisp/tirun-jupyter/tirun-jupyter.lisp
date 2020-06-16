(in-package :tirun-jupyter)

(defvar *fu* nil)
(defvar *nglv* nil)

(defvar *receptor-string* nil)

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
  (make-instance 'tirun:tirun-calculation))


(defun widget-visibility (widget visibility)
  (setf (w:widget-display (w:widget-layout widget)) visibility))
(defun widget-hide (widget)
  (widget-visibility widget "none"))
(defun widget-show (widget)
  (widget-visibility widget ""))

(defun setup-receptor (calc)
  (let* ((box-layout (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
         (desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (file-upload-msg (make-instance 'w:label :value "Receptor structure: "))
         (file-upload (make-instance 'w:file-upload :description "PDB file for receptor"
                                                    :width "auto"
                                                    :layout (make-instance 'w:layout :width "auto")))
         (file-upload-vbox (make-instance 'w:h-box
                                          :layout box-layout
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
                   (setf *receptor-string* nil)
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
                     (setf *receptor-string* as-text)
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
                                               :layout box-layout
                                               :children (list file-upload-vbox))
                                (make-instance 'w:h-box
                                               :layout box-layout
                                               :children (list messages))
                                (make-instance 'w:h-box
                                               :layout box-layout
                                               :children (list nglview)))))))


(defvar *ligands-string* nil)
(defvar *ligands* nil)
(defvar *current-component* nil)
(defvar *smallest-ligand-sketch*)
(defun observe-ligand-selector (nglview structure-view calc ligand-name new-value)
  (let ((mol (elt *ligands* new-value)))
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
(defun setup-ligands (calc &optional (show-receptor t))
  (let* ((box-layout (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
         (fu (make-instance 'w:file-upload :description "Upload SDF file"))
         (ligand-selector (make-instance 'w:bounded-int-text :description "Ligand"))
         (ligand-total (make-instance 'w:label :value ""))
         (ligand-name (make-instance 'w:label :value "Molecule"))
         (ligand-h-box (make-instance 'w:h-box :children (list ligand-selector ligand-total ligand-name)))
         (control-box (make-instance 'w:v-box :children (list fu ligand-h-box)
                                              :layout (make-instance 'w:layout :flex "1 1 0%" :width "auto" :height "auto")))
         (structure-view (make-instance 'w:html
                                        :layout (make-instance 'w:layout :flex "1 1 0%" :width "auto" :height "auto")))
         (nglview (make-instance 'nglv:nglwidget
                                 :layout (make-instance 'w:layout :width "auto" :height "auto")))
         )
    (setf *lv* nglview)
    (widget-hide ligand-h-box)
    (widget-hide nglview)
    (when show-receptor
      (nglview:add-structure nglview (make-instance 'nglview:text-structure :text *receptor-string*)))
    (w:observe fu :data
               (lambda (instance type name old-value new-value source)
                 (declare (ignore instance type name old-value source))
                 (when new-value
                   (let* ((last-pdb (car (last new-value)))
                          (octets (make-array (length last-pdb) :element-type 'ext:byte8
                                                                :initial-contents last-pdb)))
                     (let ((as-text (babel:octets-to-string octets)))
                       (setf *ligands-string* as-text)
                       (setf *ligands* (with-input-from-string (sin as-text)
                                         (sdf:parse-sdf-file sin)))
                       (flet ((build-prototype-sketch (ligands)
                                (let ((min-atoms (chem:number-of-atoms (car ligands)))
                                      (min-mol (car ligands)))
                                  (loop for mol in (cdr ligands)
                                        for num-atoms = (chem:number-of-atoms mol)
                                        when (< num-atoms min-atoms)
                                          do (setf min-atoms num-atoms
                                                   min-mol mol))
                                  (sketch2d:sketch2d min-mol))))
                         (setf *smallest-ligand-sketch* (build-prototype-sketch *ligands*)))
                       (setf (w:widget-value ligand-selector) 0
                             (w:widget-value ligand-total) (format nil "~a" (length *ligands*))
                             (w:widget-max ligand-selector) (1- (length *ligands*))
                             (w:widget-min ligand-selector) 0)
                       (widget-show ligand-h-box)
                       (widget-show nglview)
                       (observe-ligand-selector nglview structure-view calc ligand-name 0)
                       (tirun:tirun-calculation-from-ligands calc *ligands*))))))
    (w:observe ligand-selector :value (lambda (instance type name old-value new-value source)
                                        (declare (ignore instance type name old-value source))
                                        (when new-value
                                          (observe-ligand-selector nglview structure-view calc ligand-name new-value))))
    (make-instance 'w:v-box
                   :children (list
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list control-box structure-view))
                              #+(or)(make-instance 'w:h-box
                                                   ;;                                             :layout box-layout
                                                   :children (list ligand-h-box))
                              (make-instance 'w:h-box
;;;                                             :layout box-layout
                                            #+(or):layout #+(or)(make-instance 'w:layout
                                                                          :width "100em"
                                                                          :flex-flow "row"
                                                                          :align-items "stretch")
                                             :children (list nglview))
                              ))))


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

(defvar *all-nodes*)
(defvar *all-edges*)
(defvar *cyto-graph*)

(defun cyto-label ()
  (let ((widget *cyto-graph*))
    (setf (cytoscape:cytoscape-style widget)
          (append (cytoscape:cytoscape-style widget) (list (jupyter:json-new-obj
                                                            ("selector" "node[label]")
                                                            ("style" (jupyter:json-new-obj
                                                                      ("label" "data(label)"))))
                                                           (jupyter:json-new-obj
                                                            ("selector" "edge[label]")
                                                            ("style" (jupyter:json-new-obj
                                                                      ("label" "data(label)"))))))))
  nil)


(defun lomap-graph (&optional (all-nodes *all-nodes*) (all-edges *all-edges*))
  (let* ((all-nodes (loop for name in all-nodes
                          collect (make-instance 'cytoscape:node :data (list (cons "id" (string name)) (cons "label" (string name))))))
         (all-edges (loop for edge in all-edges
                          collect (destructuring-bind (name1 name2 score)
                                      edge
                                    (make-instance 'cytoscape:edge
                                                   :data (list (cons "source" (string name1))
                                                               (cons "target" (string name2))
                                                               (cons "label" (format nil "~3,2f" score))))))))
    (let* ((graph (make-instance 'cytoscape:graph
                                 :nodes all-nodes
                                 :edges all-edges))
           (widget (make-instance 'cytoscape:cytoscape-widget
                                  :graph graph
                                  :cytoscape-layout (list (cons "name" "cose")
                                                          (cons "quality" "default")
                                  ))))
      (setf *cyto-graph* widget)
      (cyto-label)
      widget)))

    

(defun lomap-progress (tirun &key (molecules *ligands*) advance-progress-callback)
  (let* ((mat (lomap::similarity-matrix molecules :advance-progress-callback advance-progress-callback))
         (multigraph (lomap::similarity-multigraph molecules mat)))
    (setf multigraph (lomap::lomap-multigraph multigraph :debug nil))
    (generate-nodes-and-edges multigraph)))


#|
(defun lomap-display (tirun)
  (let* ((all-nodes (loop for name in *all-nodes*
                          collect (make-instance 'cytoscape:node :data (list (cons "id" (string name)) (cons "label" (string name))))))
         (all-edges (loop for edge in *all-edges*
                          collect (destructuring-bind (name1 name2 score)
                                      edge
                                    (make-instance 'cytoscape:edge
                                                   :data (list (cons "source" (string name1))
                                                               (cons "target" (string name2))
                                                               (cons "label" (format nil "~3,2f" score)))))))
                                                         (graph (make-instance 'cytoscape:graph
                                                                               :nodes all-nodes
                                                                               :edges all-edges))
         (widget (make-instance 'cytoscape:cytoscape-widget
                                :cytoscape-layout (list (cons "name" "cose")
                                                        (cons "quality" "default")
                                                        )))
                  
    (setf (cytoscape:graph widget) graph)
                                                    (w:display widget)
                                                    (setf *all-nodes* all-nodes
                                                          *all-edges* all-edges
                                                          *cyto-graph* widget)

(make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list widget))
|#


(defun lomap (tirun &key (molecules *ligands*))
  (let* ((box-layout (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
         (progress (make-instance 'w:int-progress :description "Progress")) 
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (callback (let ((index 0))
                     (lambda (max-progress)
                       (setf (w:widget-value progress) index
                             (w:widget-max progress) max-progress)
                       (incf index))))
         (go-button (make-instance 'jupyter-widgets:button
                                   :description "Calculate molecular similarities"
                                   :layout (make-instance 'w:layout :width "30em")
                                   :tooltip "Click me"
                                   :on-click (list
                                              (lambda (&rest args)
                                                (setf (w:widget-value messages) "Calculating lomap")
                                                (multiple-value-bind (nodes edges)
                                                    (lomap-progress tirun
                                                                    :molecules molecules
                                                                    :advance-progress-callback callback)
                                                  (setf (w:widget-value messages) "Done lomap")
                                                  (setf *all-nodes* nodes
                                                        *all-edges* edges)))))))
    (make-instance 'w:v-box
                   :children (list
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list go-button))
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list progress))
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list messages))))))

(defun configure-jobs (calc)
  (let* ((box-layout (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
         (desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (simulation-length (make-instance 'jupyter-widgets:text :description "Simulation length"))
         (lambda-windows (make-instance 'jupyter-widgets:text :description "Lambda windows"))
         (settings (make-instance 'jupyter-widgets:v-box ;;; :layout box-layout
                                                         :children (list simulation-length
                                                                         lambda-windows)))
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (go-button (make-instance 'jupyter-widgets:button :description "Configure"
                                                           :tooltip "Click me"
                                                           :on-click (list
                                                                      (lambda (&rest args)
                                                                        (setf (w:widget-value messages)
                                                                              (format nil "Configure: ~a~%" args)))))))
    (setf *job-messages* messages)
    (make-instance 'w:v-box
                   :children (list
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list settings))
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list go-button))
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list messages))))))


(defun write-jobs (tirun &optional (name "full"))
  (tirun:build-job-nodes tirun)
  (tirun:connect-job-nodes tirun *all-edges*)
  (let* ((*default-pathname-defaults* (merge-pathnames (pathname (format nil "~a/" name))))
         (_ (ensure-directories-exist *default-pathname-defaults*))
         (worklist
           (progn
             (format t "generate-jobs to *default-pathname-defaults* -> ~s~%" *default-pathname-defaults*)
             (tirun:generate-jobs tirun))))
    (with-open-file (fout #P"jobs/conf.sh" :direction :output :if-exists :supersede)
      (format fout "pmemd_cuda=pmemd.cuda
execute_cpu_local=0
execute_gpu_local=0
execute_lisp_local=0
worker=schando
distributor=s103.thirdlaw.tech
lisp_jobs_only_on=172.234.2.1
"))
    (with-open-file (sout #P"eg5_singlestep.dot" :direction :output)
      (tirundot:draw-graph-stream (list worklist) sout))
    ))



(defvar *job-messages*)
(defun submit-jobs (calc)
  (let* ((box-layout (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
         (desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (login-label (make-instance 'w:label :value "ssh key password"))
         (login (make-instance 'jupyter-widgets:password :layout (make-instance 'w:layout :width "40em")))
         (login-h-box (make-instance 'w:h-box :layout box-layout
                                     :children (list login-label login)))
         (messages (make-instance 'w:text-area :description "Messages"
                                               :layout (make-instance 'w:layout :width "60em")))
         (go-button (make-instance 'jupyter-widgets:button :description "Submit jobs"
                                                           :tooltip "Click me"
                                                           :on-click (list
                                                                      (lambda (&rest args)
                                                                        (setf (w:widget-value messages)
                                                                              (format nil "Submitting job args: ~a~%" args)))))))
    (setf *job-messages* messages)
    (make-instance 'w:v-box
                   :children (list
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list login-h-box))
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list go-button))
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list messages))))))
