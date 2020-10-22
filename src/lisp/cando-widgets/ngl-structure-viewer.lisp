(in-package :cando-widgets)


(defun get-json-obj (obj &rest keys)
  (dolist (key keys (values obj t))
    (unless (and (listp obj)
                 (eql :obj (car obj)))
      (return (values nil nil)))
    (let ((pair (assoc key (cdr obj) :test #'string=)))
      (unless pair
        (return (values nil nil)))
      (setf obj (cdr pair)))))


(defconstant +degrees-to-radians+ (/ pi 180d0))

(defun rad2deg (value)
  (/ value +degrees-to-radians+))

(defun deg2rad (value)
  (* value +degrees-to-radians+))


(defparameter *ngl-representations*
  '((:name "Backbone"
     :receptor t
     :ligand nil
     :template nil
     :definition ((("type" . "backbone"))))
    (:name "Ball & Stick"
     :receptor t
     :ligand t
     :template t
     :definition ((("type" . "ball+stick"))))
    (:name "Cartoon"
     :receptor t
     :ligand nil
     :template nil
     :definition ((("type" . "cartoon")
                   ("colorScheme" . "residueindex"))))
    (:name "Cartoon/Ball Stick"
     :receptor t
     :ligand nil
     :template nil
     :definition ((("type" . "cartoon")
                   ("colorScheme" . "residueindex"))
                  (("type" . "ball+stick")
                   ("sele" . "ligand"))))
    (:name "Cartoon/Licorice"
     :receptor t
     :ligand nil
     :template nil
     :definition ((("type" . "cartoon")
                   ("colorScheme" . "residueindex"))
                  (("type" . "licorice")
                   ("sele" . "ligand"))))
    (:name "Cartoon/Space Fill"
     :receptor t
     :ligand nil
     :template nil
     :definition ((("type" . "cartoon")
                   ("colorScheme" . "residueindex"))
                  (("type" . "spacefill")
                   ("sele" . "ligand"))))
    (:name "Licorice"
     :receptor t
     :ligand t
     :template t
     :definition ((("type" . "licorice"))))
    (:name "Line"
     :receptor t
     :ligand t
     :template t
     :definition ((("type" . "line"))))
    (:name "Surface"
     :receptor t
     :ligand t
     :template t
     :definition ((("type" . "surface")
                   ("useWorker" . :true)
                   ("colorScheme" . "residueindex"))))
    (:name "Surface/Ball & Stick"
     :receptor t
     :ligand nil
     :template nil
     :definition ((("type" . "surface")
                   ("useWorker" . :true)
                   ("colorScheme" . "residueindex"))
                  (("type" . "ball+stick")
                   ("sele" . "ligand"))))
    (:name "Surface/Licorice"
     :receptor t
     :ligand nil
     :template nil
     :definition ((("type" . "surface")
                   ("useWorker" . :true)
                   ("colorScheme" . "residueindex"))
                  (("type" . "licorice")
                   ("sele" . "ligand"))))
    (:name "Surface/Space Fill"
     :receptor t
     :ligand nil
     :template nil
     :definition ((("type" . "surface")
                   ("useWorker" . :true)
                   ("colorScheme" . "residueindex"))
                  (("type" . "spacefill")
                   ("sele" . "ligand"))))
    (:name "Ribbon"
     :receptor t
     :ligand nil
     :template nil
     :definition ((("type" . "ribbon")
                   ("colorScheme" . "residueindex"))))
    (:name "Ribbon/Ball & Stick"
     :receptor t
     :ligand nil
     :template nil
     :definition ((("type" . "ribbon")
                   ("colorScheme" . "residueindex"))
                  (("type" . "ball+stick")
                   ("sele" . "ligand"))))
    (:name "Ribbon/Licorice"
     :receptor t
     :ligand nil
     :template nil
     :definition ((("type" . "ribbon")
                   ("colorScheme" . "residueindex"))
                  (("type" . "licorice")
                   ("sele" . "ligand"))))
    (:name "Ribbon/Space Fill"
     :receptor t
     :ligand nil
     :template nil
     :definition ((("type" . "ribbon")
                   ("colorScheme" . "residueindex"))
                  (("type" . "spacefill")
                   ("sele" . "ligand"))))
    (:name "Space Fill"
     :receptor t
     :ligand t
     :template t
     :definition ((("type" . "spacefill"))))))


(defun get-representation-definition (name)
  (dolist (rep *ngl-representations*)
    (when (string= name (getf rep :name))
      (return (getf rep :definition)))))


(defclass ngl-structure-viewer (jw:grid-box)
  ((ngl
     :reader ngl
     :initform (make-instance 'nglv:nglwidget
                              :layout (make-instance 'jw:layout
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :grid-area "ngl")))
   (receptor-representation
     :reader receptor-representation
     :initform (make-instance 'jw:dropdown
                              :description "Receptor"
                              :%options-labels (cons "None"
                                                     (mapcan (lambda (def)
                                                               (when (getf def :receptor)
                                                                 (list (getf def :name))))
                                                             *ngl-representations*))
                              :value "Cartoon/Licorice"
                              :style (make-instance 'jw:description-style
                                                    :description-width "min-content")
                              :layout (make-instance 'jw:layout
                                                     :margin ".5em"
                                                     :width "max-content"
                                                     :display "none")))
   (ligand-representation
     :reader ligand-representation
     :initform (make-instance 'jw:dropdown
                              :description "Ligand"
                              :%options-labels (mapcan (lambda (def)
                                                         (when (getf def :ligand)
                                                           (list (getf def :name))))
                                                       *ngl-representations*)
                              :value "Ball & Stick"
                              :style (make-instance 'jw:description-style
                                                    :description-width "min-content")
                              :layout (make-instance 'jw:layout
                                                     :margin ".5em"
                                                     :width "max-content"
                                                     :display "none")))
   (template-representation
     :reader template-representation
     :initform (make-instance 'jw:dropdown
                              :description "Template"
                              :%options-labels (cons "None"
                                                     (mapcan (lambda (def)
                                                               (when (getf def :template)
                                                                 (list (getf def :name))))
                                                             *ngl-representations*))
                              :value "Line"
                              :style (make-instance 'jw:description-style
                                                    :description-width "min-content")
                              :layout (make-instance 'jw:layout
                                                     :margin ".5em"
                                                     :width "max-content"
                                                     :display "none")))
   (minimize-button
     :reader minimize-button
     :initform (make-instance 'jw:button
                              :description "Minimize"
                              :style (make-instance 'jw:description-style
                                                    :description-width "min-content")
                              :layout (make-instance 'jw:layout
                                                     :margin ".5em"
                                                     :width "max-content")))
   (jostle-button
     :reader jostle-button
     :initform (make-instance 'jw:button
                              :description "Jostle"
                              :style (make-instance 'jw:description-style
                                                    :description-width "min-content")
                              :layout (make-instance 'jw:layout
                                                     :margin ".5em"
                                                     :width "max-content")))
   (angle-slider
     :reader angle-slider
     :initform (make-instance 'jw:float-slider
                              :min -180d0
                              :max 180d0
                              :step 0.5d0
                              :value 0d0
                              :readout-format ".1f"
                              :orientation "vertical"
                              :disabled t
                              :layout (make-instance 'jw:layout
                                                     :height "100%"
                                                     :grid-area "angle")))
   (controls-container
     :reader controls-container
     :initform (make-instance 'jw:box
                              :layout (make-instance 'jw:layout
                                                     :flex-flow "row wrap"
                                                     :justify-content "center"
                                                     :margin "-.5em"
                                                     :align-items "baseline"
                                                     :align-content "flex-start"
                                                     :grid-area "controls")))
   (twister
     :accessor twister
     :initform nil)
   (index
     :accessor index
     :initform nil))
  (:metaclass jupyter-widgets:trait-metaclass)
  (:default-initargs
    :layout (make-instance 'jw:layout
                           :grid-gap "1em"
                           :grid-template-rows "1fr min-content"
                           :grid-template-columns "1fr min-content"
                           :grid-template-areas "'ngl angle' 'controls .'")))


(defun on-picked (instance picked)
  (setf (index instance) nil
        (twister instance) nil
        (jw:widget-value (angle-slider instance)) 0d0
        (jw:widget-disabled (angle-slider instance)) t)
  (dotimes (component (length (nglview:components (ngl instance))))
    (nglview:remove-representations-by-name (ngl instance) "dihedral" :component component))
  (multiple-value-bind (bond presentp)
                       (get-json-obj picked "bond")
    (declare (ignore bond))
    (when presentp
      (let* ((component (get-json-obj picked "component"))
             (aggregate (matter (elt (nglview:components (ngl instance)) component)))
             (atom-vec (chem:map-atoms 'vector #'identity aggregate))
             (atom1 (aref atom-vec (get-json-obj picked "atom1" "index")))
             (atom2 (aref atom-vec (get-json-obj picked "atom2" "index")))
             (quad (find-if (lambda (d)
                              (or (and (equal 1 (position atom1 d))
                                       (equal 2 (position atom2 d)))
                                  (and (equal 2 (position atom1 d))
                                       (equal 1 (position atom2 d)))))
                            (chem:map-dihedrals 'list #'list aggregate)))
             (twister (chem:make-twister)))
        (when quad
          (nglview:add-representation (ngl instance) "dihedral"
                                      :atom-quad (list (mapcar (lambda (a)
                                                                 (position a atom-vec))
                                                               quad))
                                      :sdf :false
                                      :label-size 1.0
                                      :label-color "blue"
                                      :color-value "red"
                                      :sector-opacity 0.75)
          (apply #'chem:twister-define-for-dihedral twister quad))
          (setf (index instance) component
                (twister instance) twister
                (jw:widget-value (angle-slider instance)) (rad2deg (chem:twister-current-dihedral-angle-radians twister))
                (jw:widget-disabled (angle-slider instance)) nil)))))


(defun on-twist (instance angle)
  (with-slots (twister index ngl)
              instance
    (when (and angle index twister)
      (chem:twister-rotate-absolute twister (deg2rad angle))
      (nglv:set-coordinates ngl (list (cons index (crd (matter (elt (nglview:components ngl) index)))))))))


(defun on-jostle-relax (instance jostle)
  (loop with ngl = (ngl instance)
        for component in (nglview:components (ngl instance))
        for index from 0
        when (and (nglview:shown component)
                  (string/= "receptor" (nglview:id component)))
        do (jostle-relax ngl index (matter component) jostle)))


(defun change-receptor-representation (instance name)
  (let ((representation (get-representation-definition name)))
    (when representation
      (nglview:set-representations (ngl instance) representation "receptor"))
    (nglview:set-visibility (ngl instance) (and representation t) "receptor")))


(defun change-ligand-representation (instance name)
  (let ((representation (get-representation-definition name)))
    (when representation
      (apply #'nglview:set-representations
             (ngl instance) representation
             (remove-if (lambda (id)
                          (or (string= "receptor" id)
                              (string= "template" id)))
                        (mapcar #'nglview:id (nglview:components (ngl instance))))))))


(defun change-template-representation (instance name)
  (let ((representation (get-representation-definition name)))
    (when representation
      (nglview:set-representations (ngl instance) representation "template"))
    (nglview:set-visibility (ngl instance) (and representation t) "template")))


(defmethod initialize-instance :after ((instance ngl-structure-viewer) &rest initargs &key &allow-other-keys)
  (setf (jw:widget-children (controls-container instance))
        (list (receptor-representation instance)
              (template-representation instance)
              (ligand-representation instance)
              (minimize-button instance)
              (jostle-button instance))
        (jw:widget-children instance)
        (list (ngl instance)
              (controls-container instance)
              (angle-slider instance)))
  ; hack
  (jw:notify-trait-change (controls-container instance) :widget-list :children nil (jw:widget-children (controls-container instance)) t)
  (jw:observe (ngl instance) :picked
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (on-picked instance new-value)))
  (jw:observe (angle-slider instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (on-twist instance new-value)))
  (jw:on-button-click (minimize-button instance)
    (lambda (inst)
      (declare (ignore inst))
      (on-jostle-relax instance nil)))
  (jw:on-button-click (jostle-button instance)
    (lambda (inst)
      (declare (ignore inst))
      (on-jostle-relax instance t)))
  (jw:observe (ligand-representation instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (change-ligand-representation instance new-value)))
  (jw:observe (template-representation instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (change-template-representation instance new-value)))
  ;; Set the receptor visibility based the toggle.
  (jw:observe (receptor-representation instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (change-receptor-representation instance new-value))))


(defun crd (agg)
  (let ((arr (make-array (* 3 (chem:number-of-atoms agg)) :element-type 'single-float))
        (index -1))
    (cando:do-atoms (atm agg)
      (let ((pos (chem:get-position atm)))
        (setf (aref arr (incf index)) (float (geom:vx pos) 1.0s0)
              (aref arr (incf index)) (float (geom:vy pos) 1.0s0)
              (aref arr (incf index)) (float (geom:vz pos) 1.0s0))))
    arr))


(defun jostle-relax (ngl index aggregate jostle)
  (when jostle
    (cando:jostle aggregate))
  (let* ((energy-func (chem:make-energy-function :matter aggregate :use-excluded-atoms t :assign-types t))
         (minimizer (chem:make-minimizer energy-func)))
    (flet ((show-coords (pos)
             (let ((coords (make-array (length pos) :element-type 'single-float)))
               (loop for index below (length pos)
                     do (setf (aref coords index) (float (aref pos index) 1.0s0)))
               (nglv:set-coordinates ngl (list (cons index coords))))))
      (chem:set-step-callback minimizer #'show-coords)
      (energy:minimize-minimizer minimizer))))


(defun add-receptor (instance receptor)
  (with-slots (ngl receptor-representation)
              instance
    (nglview:remove-components ngl "receptor")
    (setf (jw:widget-display (jw:widget-layout receptor-representation))
          (unless receptor
            "none"))
    (when receptor
      (let ((representation (get-representation-definition (jw:widget-value receptor-representation))))
        (nglview:add-structure ngl
                               (make-instance 'cw:cando-structure
                                              :id "receptor" :matter receptor)
                               "defaultRepresentation" :false)
        (if representation
          (nglview:set-representations ngl representation "receptor")
          (nglview:hide-components ngl "receptor"))))))


(defun add-ligand (instance id ligand)
  (with-slots (ngl ligand-representation)
              instance
    (setf (jw:widget-display (jw:widget-layout ligand-representation))
          (unless ligand
            "none"))
    (when ligand
      (let ((representation (get-representation-definition (jw:widget-value ligand-representation))))
        (nglview:add-structure ngl
                               (make-instance 'cw:cando-structure
                                              :id id :matter ligand)
                               "defaultRepresentation" :false)
        (when representation
          (nglview:set-representations ngl representation id))))))


(defun add-template (instance template)
  (with-slots (ngl template-representation)
              instance
    (setf (jw:widget-display (jw:widget-layout template-representation))
          (unless template
            "none"))
    (when template
      (let ((representation (get-representation-definition (jw:widget-value template-representation))))
        (nglview:add-structure ngl
                               (make-instance 'cw:cando-structure
                                              :id "template" :matter template)
                               "defaultRepresentation" :false)
        (when representation
          (nglview:set-representations ngl representation "template"))))))


(defun make-ngl-structure-viewer ()
  (make-instance 'ngl-structure-viewer))

