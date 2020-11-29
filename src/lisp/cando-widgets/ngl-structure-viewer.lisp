(in-package :cando-widgets)


(defconstant +pose-frame-rate+ 60)


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


(defun ligandp (component)
  (not (position (ngl:name component) '("receptor" "template") :test #'equal)))


(defparameter *ngl-representations*
  '((:label "Backbone"
     :receptor t
     :ligand nil
     :template nil
     :names ("backbone"))
    (:label "Ball & Stick"
     :receptor t
     :ligand t
     :template t
     :names ("ball-and-stick"))
    (:label "Cartoon"
     :receptor t
     :ligand nil
     :template nil
     :names ("cartoon"))
    (:label "Cartoon/Ball Stick"
     :receptor t
     :ligand nil
     :template nil
     :names ("cartoon" "ligand_ball-and-stick"))
    (:label "Cartoon/Licorice"
     :receptor t
     :ligand nil
     :template nil
     :names ("cartoon" "ligand_licorice"))
    (:label "Cartoon/Space Fill"
     :receptor t
     :ligand nil
     :template nil
     :names ("cartoon" "ligand_spacefill"))
    (:label "Licorice"
     :receptor t
     :ligand t
     :template t
     :names ("licorice"))
    (:label "Line"
     :receptor t
     :ligand t
     :template t
     :names ("line"))
    (:label "Surface"
     :receptor t
     :ligand t
     :template t
     :names ("surface"))
    (:label "Surface/Ball & Stick"
     :receptor t
     :ligand nil
     :template nil
     :names ("surface" "ligand_ball-and-stick"))
    (:label "Surface/Licorice"
     :receptor t
     :ligand nil
     :template nil
     :names ("surface" "ligand_licorice"))
    (:label "Surface/Space Fill"
     :receptor t
     :ligand nil
     :template nil
     :names ("surface" "ligand_spacefill"))
    (:label "Ribbon"
     :receptor t
     :ligand nil
     :template nil
     :names ("ribbon"))
    (:label "Ribbon/Ball & Stick"
     :receptor t
     :ligand nil
     :template nil
     :names ("ribbon" "ligand_ball-and-stick"))
    (:label "Ribbon/Licorice"
     :receptor t
     :ligand nil
     :template nil
     :names ("ribbon" "ligand_licorice"))
    (:label "Ribbon/Space Fill"
     :receptor t
     :ligand nil
     :template nil
     :names ("ribbon" "ligand_spacefill"))
    (:label "Space Fill"
     :receptor t
     :ligand t
     :template t
     :names ("spacefill"))))


(defun create-representations ()
  (list (make-instance 'ngl:dihedral :name "dihedral" :visible nil :sdf nil :sector-opacity 0.75
                                     :label-size 1.0 :label-color "blue" :color-value "red")
        (make-instance 'ngl:backbone :name "backbone")
        (make-instance 'ngl:ball-and-stick :name "ball-and-stick")
        (make-instance 'ngl:ball-and-stick :name "ligand_ball-and-stick" :sele "ligand")
        (make-instance 'ngl:cartoon :name "cartoon" :color-scheme "residueindex")
        (make-instance 'ngl:licorice :name "licorice")
        (make-instance 'ngl:licorice :name "ligand_licorice" :sele "ligand")
        (make-instance 'ngl:line :name "line")
        (make-instance 'ngl:ribbon :name "ribbon" :color-cheme "residueindex")
        (make-instance 'ngl:spacefill :name "ligand_spacefill" :sele "ligand")
        (make-instance 'ngl:spacefill :name "spacefill")
        (make-instance 'ngl:surface :name "surface" :use-worker t :color-scheme "residueindex")))


(defun update-representations (component value)
  (let ((definition (getf (find value *ngl-representations* :key (lambda (rep) (getf rep :label)) :test #'equal) :names)))
    (jupyter:inform :info component "~S" definition)
    (cond
      (definition
        (dolist (representation (ngl:representations component))
          (unless (equal (ngl:name representation) "dihedral")
            (setf (ngl:visible representation)
                  (and (position (ngl:name representation) definition :test #'equal)
                       t)))))
      (t
        (dolist (representation (ngl:representations component))
          (setf (ngl:visible representation) nil))))))


(defclass ngl-structure-viewer ()
  ((stage
     :reader stage
     :initform (make-instance 'ngl:stage
                              :clip-dist 0
                              :background-color "white"
                              :layout (make-instance 'jw:layout
                                                     :width "100%"
                                                     :height "auto"
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :grid-area "ngl-content")))
   (receptor-representation
     :reader receptor-representation
     :initform (make-instance 'jw:dropdown
                              :description "Receptor Representation"
                              :%options-labels (cons "None"
                                                     (mapcan (lambda (def)
                                                               (when (getf def :receptor)
                                                                 (list (getf def :label))))
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
                              :description "Ligand Representation"
                              :%options-labels (mapcan (lambda (def)
                                                         (when (getf def :ligand)
                                                           (list (getf def :label))))
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
                              :description "Template Representation"
                              :%options-labels (cons "None"
                                                     (mapcan (lambda (def)
                                                               (when (getf def :template)
                                                                 (list (getf def :label))))
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
                              :description "Align & Minimize"
                              :style (make-instance 'jw:description-style
                                                    :description-width "min-content")
                              :layout (make-instance 'jw:layout
                                                     :margin ".5em"
                                                     :width "max-content"
                                                     :display "none")))
   (jostle-button
     :reader jostle-button
     :initform (make-instance 'jw:toggle-button
                              :description "Jostle"
                              :value t
                              :style (make-instance 'jw:description-style
                                                    :description-width "min-content")
                              :layout (make-instance 'jw:layout
                                                     :margin ".5em"
                                                     :width "max-content"
                                                     :display "none")))
   (auto-view-type
     :reader auto-view-type
     :initform (make-instance 'jw:dropdown
                              :%options-labels '("All" "Ligand" "Receptor" "Template")
                              :options '("all" "ligand" "receptor" "template")
                              :value "all"
                              :description "Auto View Type"
                              :style (make-instance 'jw:description-style
                                                    :description-width "min-content")
                              :layout (make-instance 'jw:layout
                                                     :margin ".5em"
                                                     :width "max-content")))
   (auto-view-button
     :reader auto-view-button
     :initform (make-instance 'jw:button
                              :description "Auto View"
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
                                                     :grid-area "ngl-vertical-toolbar"
                                                     :display "none")))
   (controls-container
     :reader controls-container
     :initform (make-instance 'jw:box
                              :layout (make-instance 'jw:layout
                                                     :flex-flow "row wrap"
                                                     :justify-content "center"
                                                     ;:margin "-.5em"
                                                     :align-items "baseline"
                                                     :align-content "flex-start"
                                                     :grid-area "ngl-horizontal-toolbar")))
   (twister
     :accessor twister
     :initform nil)
   (index
     :accessor index
     :initform nil)
   (structures
     :accessor structures
     :initform (make-hash-table :test #'equal))))


(defun on-picked (instance picked)
  (dolist (component (ngl:components (stage instance)))
    (setf (ngl:visible (find "dihedral" (ngl:representations component) :key #'ngl:name :test #'equal)) nil))
  (setf (index instance) nil
        (twister instance) nil
        (jw:widget-value (angle-slider instance)) 0d0
        (jw:widget-disabled (angle-slider instance)) t)
  (when (and (equal (getf picked :signal) "click")
             (equal (getf picked :type) "bond"))
    (let* ((bond (getf picked :bond))
           (name (getf picked :component))
           (aggregate (gethash name (structures instance)))
           (component (find name (ngl:components (stage instance)) :test #'equal :key #'ngl:name))
           (atom-vec (chem:map-atoms 'vector #'identity aggregate))
           (atom1 (aref atom-vec (getf bond :atom-index1)))
           (atom2 (aref atom-vec (getf bond :atom-index2)))
           (quad (find-if (lambda (d)
                            (or (and (equal 1 (position atom1 d))
                                     (equal 2 (position atom2 d)))
                                (and (equal 2 (position atom1 d))
                                     (equal 1 (position atom2 d)))))
                          (chem:map-dihedrals 'list #'list aggregate)))
           (twister (chem:make-twister))
           (dihedral-representation (find "dihedral" (ngl:representations component) :test #'equal :key #'ngl:name)))
      (when quad
        (setf (ngl:atom-quad dihedral-representation)
              (list (mapcar (lambda (a)
                              (position a atom-vec))
                            quad))
              (ngl:visible dihedral-representation) t)
        (apply #'chem:twister-define-for-dihedral twister quad)
        (setf (index instance) component
              (twister instance) twister
              (jw:widget-value (angle-slider instance)) (rad2deg (chem:twister-current-dihedral-angle-radians twister))
              (jw:widget-disabled (angle-slider instance)) nil)))))


(defun on-twist (instance angle)
  (with-slots (twister index)
              instance
    (when (and angle index twister)
      (chem:twister-rotate-absolute twister (deg2rad angle))
      (ngl:update-position index (crd (gethash (ngl:name index) (structures instance)))))))


(defun on-jostle-relax (instance jostle)
  (dolist (component (ngl:components (stage instance)))
    (when (and (ngl:visible component)
               (ligandp component))
      (jostle-relax component (gethash (ngl:name component) (structures instance))
                    (gethash "template" (structures instance)) jostle))))


(defun change-receptor-representation (instance name)
  (let ((component (find "receptor" (ngl:components (stage instance)) :key #'ngl:name :test #'string=)))
    (when component
      (update-representations component name))))


(defun change-ligand-representation (instance name)
  (dolist (component (ngl:components (stage instance)))
    (when (ligandp component)
      (update-representations component name))))


(defun change-template-representation (instance name)
  (let ((component (find "template" (ngl:components (stage instance)) :key #'ngl:name :test #'string=)))
    (when component
      (update-representations component name))))


(defmethod initialize-instance :after ((instance ngl-structure-viewer) &rest initargs &key &allow-other-keys)
  (setf (jw:widget-children (controls-container instance))
        (list (receptor-representation instance)
              (ligand-representation instance)
              (template-representation instance)
              (auto-view-type instance)
              (auto-view-button instance)
              (jostle-button instance)
              (minimize-button instance)))
  (ngl:on-stage-pick (stage instance) (lambda (inst pick)
                                       (declare (ignore inst))
                                       (on-picked instance pick)))
  (jw:observe (angle-slider instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (on-twist instance new-value)))
  (jw:on-button-click (auto-view-button instance)
    (lambda (inst)
      (declare (ignore inst))
      (auto-view instance)))
  (jw:on-button-click (minimize-button instance)
    (lambda (inst)
      (declare (ignore inst))
      (on-jostle-relax instance (jw:widget-value (jostle-button instance)))))
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


(defun auto-view (instance)
  (let* ((components (ngl:components (stage instance)))
         (type (jw:widget-value (auto-view-type instance)))
         (obj (cond
                ((= 1 (length components))
                  (car components))
                ((equal "all" type)
                  (stage instance))
                ((equal "ligand" type)
                  (find-if #'ligandp components))
                (t
                  (find type components :key #'ngl:name :test #'equal)))))
    (when obj
      (ngl:auto-view obj 500))))


(defun jostle-relax (component aggregate template jostle)
  (chem:map-atoms nil #'chem:clear-restraints aggregate)
  (when jostle
    (cando:jostle aggregate))
  (when template
    (tirun::pose-one-molecule-using-similarity aggregate template))
  (bordeaux-threads:make-thread
    (lambda ()
      (let* ((energy-func (chem:make-energy-function :matter aggregate :use-excluded-atoms t :assign-types t))
             (minimizer (chem:make-minimizer energy-func))
             (frame-period (floor internal-time-units-per-second +pose-frame-rate+))
             (next-frame-time (get-internal-real-time)))
        (chem:set-step-callback minimizer
                                (lambda (pos)
                                  (sleep (max 0 (/ (- next-frame-time (get-internal-real-time)) internal-time-units-per-second)))
                                  (ngl:update-position component
                                                       (map '(vector single-float *)
                                                            (lambda (x)
                                                              (coerce x 'single-float))
                                                            pos))
                                  (setf next-frame-time (+ (get-internal-real-time) frame-period))))
        (energy:minimize-minimizer minimizer :print-intermediate-results nil)))))


(defun add-receptor (instance receptor)
  (with-slots (stage receptor-representation)
              instance
    (setf (ngl:components stage) (remove "receptor" (ngl:components stage) :key #'ngl:name :test #'string=))
    (setf (jw:widget-display (jw:widget-layout receptor-representation))
          (unless receptor
            "none"))
    (when receptor
      (let ((component (make-instance 'ngl:structure
                                      :name "receptor"
                                      :value (chem:aggregate-as-mol2-string receptor t)
                                      :ext "mol2"
                                      :representations (create-representations))))
        (update-representations component (jw:widget-value receptor-representation))
        (setf (ngl:components stage)
              (cons component
                    (ngl:components stage)))
        (auto-view instance)))))


(defun add-ligand (instance id ligand)
  (with-slots (stage ligand-representation minimize-button jostle-button angle-slider)
              instance
    (let ((display (unless ligand "none")))
      (setf (jw:widget-display (jw:widget-layout ligand-representation)) display
            (jw:widget-display (jw:widget-layout minimize-button)) display
            (jw:widget-display (jw:widget-layout jostle-button)) display
            (jw:widget-display (jw:widget-layout angle-slider)) display))
    (dolist (component (ngl:components stage))
      (when (ligandp component)
        (setf (ngl:visible component) nil)))
    (when ligand
      (let ((component (find id (ngl:components stage) :test #'equal :key #'ngl:name)))
        (cond
          (component
            (setf (ngl:visible component) t))
          (t
            (setf component (make-instance 'ngl:structure
                                           :name id
                                           :value (chem:aggregate-as-mol2-string ligand t)
                                           :ext "mol2"
                                           :representations (create-representations)))
            (update-representations component (jw:widget-value ligand-representation))
            (setf (gethash id (structures instance)) ligand
                  (ngl:components stage) (cons component
                                             (ngl:components stage)))))
        (auto-view instance)))))


(defun add-template (instance template)
  (with-slots (stage template-representation minimize-button jostle-button angle-slider)
              instance
    (setf (jw:widget-display (jw:widget-layout template-representation)) (unless template "none"))
    (let ((components (remove "template" (ngl:components stage) :key #'ngl:name :test #'equal)))
      (cond
        (template
          (let ((component (make-instance 'ngl:structure
                                          :name "template"
                                          :value (chem:aggregate-as-mol2-string template t)
                                          :ext "mol2"
                                          :representations (create-representations))))
            (update-representations component (jw:widget-value template-representation))
            (setf (gethash "template" (structures instance)) template
                  (ngl:components stage) (append components (list component)))
            (auto-view instance)))
        (t
          (remhash "template" (structures instance))
          (setf (ngl:components stage) components))))))


(defun clear-ligands (instance)
  (setf (ngl:components (stage instance))
        (remove-if #'ligandp (ngl:components (stage instance)))))


(defun make-ngl-structure-viewer (container)
  (let ((instance (make-instance 'ngl-structure-viewer)))
    (setf (jw:widget-children container)
          (append (jw:widget-children container)
                  (list (stage instance)
                        (controls-container instance)
                        (angle-slider instance))))
    instance))


