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
     :definition ((("type" . "backbone"))))
    (:name "Ball & Stick"
     :receptor t
     :ligand t
     :definition ((("type" . "ball+stick"))))
    (:name "Cartoon"
     :receptor t
     :ligand nil
     :definition ((("type" . "cartoon")
                   ("colorScheme" . "residueindex"))))
    (:name "Cartoon/Ball Stick"
     :receptor t
     :ligand nil
     :definition ((("type" . "cartoon")
                   ("colorScheme" . "residueindex"))
                  (("type" . "ball+stick")
                   ("sele" . "ligand"))))
    (:name "Cartoon/Licorice"
     :receptor t
     :ligand nil
     :definition ((("type" . "cartoon")
                   ("colorScheme" . "residueindex"))
                  (("type" . "licorice")
                   ("sele" . "ligand"))))
    (:name "Licorice"
     :receptor t
     :ligand t
     :definition ((("type" . "licorice"))))
    (:name "Line"
     :receptor t
     :ligand t
     :definition ((("type" . "line"))))
    (:name "Surface"
     :receptor t
     :ligand t
     :definition ((("type" . "surface")
                   ("useWorker" . :true)
                   ("colorScheme" . "residueindex"))))
    (:name "Surface/Ball & Stick"
     :receptor t
     :ligand nil
     :definition ((("type" . "surface")
                   ("useWorker" . :true)
                   ("colorScheme" . "residueindex"))
                  (("type" . "ball+stick")
                   ("sele" . "ligand"))))
    (:name "Surface/Licorice"
     :receptor t
     :ligand nil
     :definition ((("type" . "surface")
                   ("useWorker" . :true)
                   ("colorScheme" . "residueindex"))
                  (("type" . "licorice")
                   ("sele" . "ligand"))))
    (:name "Ribbon"
     :receptor t
     :ligand nil
     :definition ((("type" . "ribbon")
                   ("colorScheme" . "residueindex"))))
    (:name "Ribbon/Ball & Stick"
     :receptor t
     :ligand nil
     :definition ((("type" . "ribbon")
                   ("colorScheme" . "residueindex"))
                  (("type" . "ball+stick")
                   ("sele" . "ligand"))))
    (:name "Ribbon/Licorice"
     :receptor t
     :ligand nil
     :definition ((("type" . "ribbon")
                   ("colorScheme" . "residueindex"))
                  (("type" . "licorice")
                   ("sele" . "ligand"))))))


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
                                                     :grid-area "ngl1")))
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
                              :layout (make-instance 'jw:layout
                                                     :margin "auto"
                                                     :width "max-content"
                                                     :display "none"
                                                     :grid-area "receptor")))
   (ligand-representation
     :reader ligand-representation
     :initform (make-instance 'jw:dropdown
                              :description "Ligand"
                              :%options-labels (mapcan (lambda (def)
                                                         (when (getf def :ligand)
                                                           (list (getf def :name))))
                                                       *ngl-representations*)
                              :value "Ball & Stick"
                              :layout (make-instance 'jw:layout
                                                     :margin "auto"
                                                     :width "max-content"
                                                     :grid-area "ligand")))
   (minimize-button
     :reader minimize-button
     :initform (make-instance 'jw:button
                              :description "Minimize"
                              :layout (make-instance 'jw:layout
                                                     :margin "auto"
                                                     :grid-area "minimize")))
   (jostle-button
     :reader jostle-button
     :initform (make-instance 'jw:button
                              :description "Jostle"
                              :layout (make-instance 'jw:layout
                                                     :margin "auto"
                                                     :grid-area "jostle")))
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
                           :grid-template-columns "1fr min-content min-content min-content min-content 1fr min-content"
                           :grid-template-areas "'ngl1 ngl1 ngl1 ngl1 ngl1 ngl1 angle' '. receptor ligand jostle minimize . .'")))


(defun on-picked (instance picked)
  (multiple-value-bind (bond presentp)
                       (get-json-obj picked "bond")
    (declare (ignore bond))
    (if presentp
      (let* ((component (get-json-obj picked "component"))
             (aggregate (matter (elt (nglview:components (ngl instance)) component)))
             (atom-vec (chem:map-atoms 'vector #'identity aggregate))
             (twister (chem:make-twister)))
        (chem:twister-define-for-bond twister
                                      (aref atom-vec (get-json-obj picked "atom1" "index"))
                                      (aref atom-vec (get-json-obj picked "atom2" "index")))
        (setf (index instance) component
              (twister instance) twister
              (jw:widget-value (angle-slider instance)) (rad2deg (chem:twister-current-dihedral-angle-radians twister))
              (jw:widget-disabled (angle-slider instance)) nil))
      (setf (index instance) nil
            (twister instance) nil
            (jw:widget-value (angle-slider instance)) 0d0
            (jw:widget-disabled (angle-slider instance)) t))))


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
             (remove "receptor" (mapcar #'nglview:id (nglview:components (ngl instance))) :test #'string=)))))


(defmethod initialize-instance :after ((instance ngl-structure-viewer) &rest initargs &key &allow-other-keys)
  (setf (jw:widget-children instance)
        (list (ngl instance)
              (receptor-representation instance)
              (ligand-representation instance)
              (angle-slider instance)
              (minimize-button instance)
              (jostle-button instance)))
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
    (when ligand
      (let ((representation (get-representation-definition (jw:widget-value ligand-representation))))
        (nglview:add-structure ngl
                               (make-instance 'cw:cando-structure
                                              :id id :matter ligand)
                               "defaultRepresentation" :false)
        (when representation
          (nglview:set-representations ngl representation id))))))
