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

(defun setup-receptor* (nglview messages calc new-value)
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
                     (setf (jupyter-widgets:widget-value messages) (leap.pdb:messages condition))
                     (setf (jupyter-widgets:widget-display (jupyter-widgets:widget-layout messages)) "") 
                     (setf (jupyter-widgets:widget-display (jupyter-widgets:widget-layout nglview)) "none") 
                     (return-from setup-receptor* nil)))))
        (nglview:add-structure nglview (make-instance 'nglview:text-structure :text as-text))
        (setf (jupyter-widgets:widget-value messages) "Loaded successfully")
        (setf (jupyter-widgets:widget-display (jupyter-widgets:widget-layout messages)) "none") 
        (setf (jupyter-widgets:widget-display (jupyter-widgets:widget-layout nglview)) "") 
        (setf *receptor-string* as-text)
        (setf (tirun:receptors calc) (list agg))))))

(defun setup-receptor (calc)
  (let* ((box-layout (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
         (desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (fu (make-instance 'jupyter-widgets:file-upload :description "Upload PDB file"))
         (messages (make-instance 'jupyter-widgets:text-area :description "Messages"
                                  :layout (make-instance 'jupyter-widgets:layout :width "60em")))
         (nglview (make-instance 'nglv:nglwidget)))
    (setf (jupyter-widgets:widget-display (jupyter-widgets:widget-layout messages)) "") 
    (setf (jupyter-widgets:widget-display (jupyter-widgets:widget-layout nglview)) "none") 
    (setf *fu* fu
          *nglv* nglview)
    (jupyter-widgets:observe fu :data
                             (lambda (instance type name old-value new-value source)
                               (declare (ignore instance type name old-value source))
                               (when new-value
                                 (setup-receptor* nglview messages calc new-value))))
    (make-instance 'w:v-box
                   :children (list
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list fu))
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list messages))
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list nglview)))))
  )


(defvar *ligands-string* nil)
(defvar *ligands* nil)
(defvar *current-component* nil)
(defun observe-ligand-selector (nglview calc ligand-name new-value)
  (let ((mol (elt *ligands* new-value)))
    (let* ((agg (chem:make-aggregate)))
      (chem:add-matter agg mol)
      (setf (jupyter-widgets:widget-value ligand-name) (string (chem:get-name mol)))
      (let* ((mol2 (chem:aggregate-as-mol2-string agg t))
             (new-component (nglview:add-structure nglview (make-instance 'nglview:text-structure :text mol2 :ext "mol2"))))
        (when *current-component*
          (nglview:remove-components nglview *current-component*))
        (setf *current-component* new-component)
        ))))

(defvar *lv* nil)
(defun setup-ligands (calc &optional (show-receptor t))
  (let* ((box-layout (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
         (desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (fu (make-instance 'jupyter-widgets:file-upload :description "Upload SDF file"))
         (ligand-selector (make-instance 'jupyter-widgets:bounded-int-text :description "Ligand"))
         (ligand-name (make-instance 'jupyter-widgets:label :value "Molecule"))
         (nglview (make-instance 'nglv:nglwidget :gui-style "ngl")))
    (setf *lv* nglview)
    (setf (jupyter-widgets:widget-display (jupyter-widgets:widget-layout ligand-selector)) "none") 
    (setf (jupyter-widgets:widget-display (jupyter-widgets:widget-layout ligand-name)) "none") 
    (setf (jupyter-widgets:widget-display (jupyter-widgets:widget-layout nglview)) "none") 
    (when show-receptor
      (nglview:add-structure nglview (make-instance 'nglview:text-structure :text *receptor-string*)))
    (jupyter-widgets:observe fu :data
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
                                     (setf (jupyter-widgets:widget-value ligand-selector) 0
                                           (jupyter-widgets:widget-max ligand-selector) (1- (length *ligands*))
                                           (jupyter-widgets:widget-min ligand-selector) 0)
                                     (setf (jupyter-widgets:widget-display (jupyter-widgets:widget-layout ligand-selector)) "") 
                                     (setf (jupyter-widgets:widget-display (jupyter-widgets:widget-layout ligand-name)) "") 
                                     (setf (jupyter-widgets:widget-display (jupyter-widgets:widget-layout nglview)) "") 
                                     (observe-ligand-selector nglview calc ligand-name 0)
                                     (tirun:tirun-calculation-from-ligands calc *ligands*))))))
    (jupyter-widgets:observe ligand-selector :value (lambda (instance type name old-value new-value source)
                                                        (declare (ignore instance type name old-value source))
                                                      (when new-value
                                                        (observe-ligand-selector nglview calc ligand-name new-value))))
    (make-instance 'w:v-box
                   :children (list
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list fu))
                              (make-instance 'w:h-box
;;                                             :layout box-layout
                                             :children (list ligand-selector ligand-name))
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list nglview))))))


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


(defun cyto-label ()
  (let ((widget *cyto-widget*))
    (setf (cytoscape:graph-style widget)
          (append (cytoscape:graph-style widget)
                  (list (jupyter:json-new-obj
                          ("selector" "node[label]")
                          ("style" (jupyter:json-new-obj
                                     ("label" "data(label)"))))
                        (jupyter:json-new-obj
                          ("selector" "edge[label]")
                          ("style" (jupyter:json-new-obj
                                     ("label" "data(label)"))))))))
  nil)


(defvar *cyto-widget*)
(defun cyto-graph (&optional (all-nodes *all-nodes*) (all-edges *all-edges*))
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
    (setf *cyto-widget* (make-instance 'cytoscape:cytoscape-widget :elements (append all-nodes all-edges)))
    (cyto-label)
    *cyto-widget*))

    

(defvar *all-nodes*)
(defvar *all-edges*)
(defvar *cyto-graph*)

(defun lomap (tirun &optional (mols *ligands*))
  (let* ((mat (lomap::similarity-matrix mols))
         (multigraph (lomap::similarity-multigraph mols mat)))
    (setf multigraph (lomap::lomap-multigraph multigraph :debug nil))
    (multiple-value-setq (*all-nodes* *all-edges*) (generate-nodes-and-edges multigraph))
    (setq *cyto-graph* (cyto-graph *all-nodes* *all-edges*))))


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
