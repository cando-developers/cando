(in-package :tirun-jupyter)

(defvar *fu* nil)
(defvar *nglv* nil)

(defvar *receptor-string* nil)

(defun setup-receptor ()
  (let* ((box-layout (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
         (desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (fu (make-instance 'jupyter-widgets:file-upload :description "Upload PDB file"))
         (nglview (make-instance 'nglv:nglwidget)))
    (setf *fu* fu
          *nglv* nglview)
    (jupyter-widgets:observe fu :data
                             (lambda (instance type name old-value new-value source)
                               (declare (ignore instance type name old-value source))
                               (when new-value
                                 (let* ((last-pdb (car (last new-value)))
                                        (octets (make-array (length last-pdb) :element-type 'ext:byte8
                                                                              :initial-contents last-pdb)))
                                   (let ((as-text (babel:octets-to-string octets)))
                                     (setf *receptor-string* as-text)
                                     (nglview:remove-all-components nglview)
                                     (nglview:add-structure nglview (make-instance 'nglview:text-structure :text as-text)))))))
    (make-instance 'w:v-box
                   :children (list
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list fu))
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list nglview))))))


(defvar *ligands-string* nil)
(defvar *ligands* nil)
(defvar *current-component* nil)
(defun observe-ligand-selector (nglview ligand-name new-value)
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
(defun setup-ligands (&optional show-receptor)
  (let* ((box-layout (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
         (desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (fu (make-instance 'jupyter-widgets:file-upload :description "Upload SDF file"))
         (ligand-selector (make-instance 'jupyter-widgets:bounded-int-text :description "Ligand"))
         (ligand-name (make-instance 'jupyter-widgets:label :value "Molecule"))
         (nglview (make-instance 'nglv:nglwidget :gui-style "ngl")))
    (setf *lv* nglview)
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
                                     (observe-ligand-selector nglview ligand-name 0))))))
    (jupyter-widgets:observe ligand-selector :value (lambda (instance type name old-value new-value source)
                                                        (declare (ignore instance type name old-value source))
                                                      (when new-value
                                                        (observe-ligand-selector nglview ligand-name new-value))))
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


(defun cyto-graph (graph)
  (let ((nodes (loop for vertex in (lomap::vertices graph)
                     for mol = (lomap:molecule vertex)
                     for name = (string (chem:get-name mol))
                     collect (make-instance 'cytoscape:node :data (list (cons "id" name) (cons "label" name)))))
        (edges (loop for edge in (lomap::edges graph)
                     for vertex1 = (lomap:vertex1 edge)
                     for vertex2 = (lomap:vertex2 edge)
                     for score = (lomap:sim-score edge)
                     for name1 = (string (chem:get-name (lomap:molecule vertex1)))
                     for name2 = (string (chem:get-name (lomap:molecule vertex2)))
                     collect (make-instance 'cytoscape:edge
                                            :data (list (cons "source" name1)
                                                        (cons "target" name2)
                                                        (cons "label" (format nil "~f" score)))))))
    (let ((graph (make-instance 'cytoscape:graph
                                :nodes nodes
                                :edges edges)))
      (make-instance 'cytoscape:cytoscape-widget
                     :graph graph))))
  
(defun layout-calculation (mols)
  (let* ((calc (tirun::tirun-calculation-from-ligands mols))
         (mat (lomap::similarity-matrix mols))
         (graph (lomap::similarity-graph mols mat)))
    ;;;(lomap::lomap-graph graph :debug nil)
    (cyto-graph graph)))

#|
        
    (defparameter +options+ '("fu" "bar" "wibble" "quux" "baz" "gronk" "kilroy"))
    (defparameter +desc-width+ "12em")
    (defparameter +desc-style+ (make-instance 'w:description-style :description-width +desc-width+))
    (defparameter +box-layout+ (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
    (defparameter +slider-style+ (make-instance 'w:slider-style :description-width +desc-width+))
    (defparameter +widget-layout+ (make-instance 'w:layout :width "28em"))
    
    (make-instance 'w:accordion
                   :children (list (make-instance 'w:h-box
                                                  :layout +box-layout+
                                                  :children (list (make-instance 'w:color-picker :description "color-picker" :style +desc-style+)
                                                                  (make-instance 'w:date-picker :description "date-picker" :style +desc-style+)))
                                   (make-instance 'w:h-box
                                                  :layout +box-layout+
                                                  :children (list (make-instance 'w:dropdown :description "dropdown" :style +desc-style+ :layout +widget-layout+ :%options-labels +options+)
                                                                  (make-instance 'w:radio-buttons :description "radio-button" :style +desc-style+ :layout +widget-layout+ :%options-labels +options+)
                                                                  (make-instance 'w:select :description "select" :style +desc-style+ :layout +widget-layout+  :%options-labels +options+)
                                                                  (make-instance 'w:selection-slider :description "selection-slider" :style +slider-style+ :layout +widget-layout+ :%options-labels +options+)
                                                                  (make-instance 'w:selection-range-slider :description "selection-range-slider" :style +slider-style+ :layout +widget-layout+  :%options-labels +options+)
                                                                  (make-instance 'w:select-multiple :description "select-multiple" :style +desc-style+ :layout +widget-layout+  :%options-labels +options+)))
                                   (make-instance 'w:h-box
                                                  :layout +box-layout+
                                                  :children (list (make-instance 'w:button :description "button" :style +desc-style+)
                                                                  (make-instance 'w:toggle-button :description "toggle-button" :style +desc-style+))))
                   :%titles '("Picker Widgets" "Selection Widgets" "Buttons")))


|#
