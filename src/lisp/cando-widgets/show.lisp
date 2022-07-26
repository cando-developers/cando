(in-package :cando-widgets)

(defparameter *panes* (make-hash-table :test #'equal))

(defparameter *default-pane* nil)

(defgeneric make-pane (object &rest rest &key &allow-other-keys))

(defgeneric show-on-pane (pane-instance object &rest rest &key &allow-other-keys))

(defun show (object &rest rest &key (pane *default-pane*) &allow-other-keys)
  "Display a residue, molecule, aggregate or a trajectory in a Jupyter notebook"
  (if pane
      (let ((existing (gethash pane *panes*)))
        (if existing
            (apply #'show-on-pane existing object rest)
            (let ((sidecar (jw:make-sidecar :title pane)))
              (jw:with-output sidecar
                (setf (gethash pane *panes*)
                      (apply #'show-on-pane nil object rest))))))
      (apply #'show-on-pane nil object rest))
  (values))

;;; sketch2d

(defmethod show-on-pane (pane-instance (sketch sketch2d:sketch2d) &rest rest &key &allow-other-keys)
  (apply #'show-on-pane pane-instance (sketch2d:svg sketch) rest))

(defmethod show-on-pane (pane-instance (sketch sketch2d:sketch-svg) &rest rest &key pane &allow-other-keys)
  (declare (ignore rest))
  (jupyter:svg (sketch2d:render-svg-to-string sketch) :display t :id pane :update (and pane-instance t))
  pane)

(defun make-axes (visible)
  (make-instance 'ngl:buffer-representation
                 :name nil
                 :visible visible
                 :buffer (list :type "wideline"
                               :position1 (vector 0 0 0
                                                  0 0 0
                                                  0 0 0)
                               :position2 (vector 10 0 0
                                                  0 10 0
                                                  0 0 10)
                               :color #(#xFF 0 0 0 0xFF 0 0 0 0xFF)
                               :color2 #(#xFF 0 0 0 0xFF 0 0 0 0xFF))))

(defun cylinder-buffer (position1 position2 radius color)
  (list :type "cylinder"
        :position1 (apply #'append position1)
        :position2 (apply #'append position2)
        :radius (make-list (length position1) :initial-element radius)
        :color (loop for p on position1
                     append color)
        :color2 (loop for p on position1
                      append color)))

(defun sphere-buffer (position radius color)
  (list :type "sphere"
        :position (apply #'append position)
        :radius (make-list (length position) :initial-element radius)
        :color (loop for p on position
                     append color)))

(defvar +bounding-box-color+ '(1 .65 0))
(defvar +bounding-box-radius+ 0.15)

(defun add-bounding-box (component aggregate)
  (when (chem:bounding-box-bound-p aggregate)
    (let* ((bounding-box (chem:bounding-box aggregate))
           (min (chem:min-corner bounding-box))
           (max (chem:max-corner bounding-box))
           (minx (geom:get-x min))
           (miny (geom:get-y min))
           (minz (geom:get-z min))
           (maxx (geom:get-x max))
           (maxy (geom:get-y max))
           (maxz (geom:get-z max))
           (p000 (list minx miny minz))
           (p001 (list minx miny maxz))
           (p010 (list minx maxy minz))
           (p011 (list minx maxy maxz))
           (p100 (list maxx miny minz))
           (p101 (list maxx miny maxz))
           (p110 (list maxx maxy minz))
           (p111 (list maxx maxy maxz))
           (cyl (cylinder-buffer (list p000 p000 p000 p001 p001 p010
                                       p010 p100 p100 p111 p111 p111)
                                 (list p001 p010 p100 p011 p101 p011
                                       p110 p101 p110 p110 p101 p011)
                                 +bounding-box-radius+
                                 +bounding-box-color+))
           (sph (sphere-buffer (list p000 p001 p010 p011 p100 p101 p110 p111)
                               +bounding-box-radius+
                               +bounding-box-color+)))
      (setf (ngl:representations component)
            (append (ngl:representations component)
                    (list (ngl:make-buffer-representation :name nil :visible t
                                                          :buffer cyl)
                          (ngl:make-buffer-representation :name nil :visible t
                                                          :buffer sph)))))))

(defun ngl-show-trajectory (trajectory)
  (when trajectory
    (let* ((desc-style (jw:make-description-style :description-width "min-content"))
           (play-back-2 (jw:make-toggle-button :icon "backward" :tooltip "Fast backward"
                                               :style desc-style
                                               :layout (jw:make-layout :margin ".5em .1em .5em .5em"
                                                                       :width "max-content")))
           (play-back-1 (jw:make-toggle-button :icon "caret-left" :tooltip "Backward"
                                               :style desc-style
                                               :layout (jw:make-layout :margin ".5em .1em .5em .1em"
                                                                       :width "max-content")))
           (stop-button (jw:make-button :icon "stop"
                                        :tooltip "Stop"
                                        :style desc-style
                                        :layout (jw:make-layout :margin ".5em .1em .5em .1em"
                                                                :width "max-content")))
           (pause-button (jw:make-button :icon "pause"
                                         :tooltip "Pause"
                                         :style desc-style
                                         :layout (jw:make-layout :margin ".5em .1em .5em .1em"
                                                                 :width "max-content")))
           (play-fore-1 (jw:make-toggle-button :icon "play" :tooltip "Foreward"
                                               :style desc-style
                                               :layout (jw:make-layout :margin ".5em .1em .5em .1em"
                                                                       :width "max-content")))
           (play-fore-2 (jw:make-toggle-button :icon "forward" :tooltip "Fast foreward"
                                               :style desc-style
                                               :layout (jw:make-layout :margin ".5em .1em .5em .1em"
                                                                       :width "max-content")))
           (mode-button (jw:make-toggle-button :icon "retweet"
                                               :tooltip "Loop"
                                               :value t
                                               :style desc-style
                                               :layout (jw:make-layout :margin ".5em .1em .5em .1em"
                                                                       :width "max-content")))
           (frame-slider (jw:make-int-slider :layout (jw:make-layout :align-self "center"
                                                                     :margin ".25em"))))
      (jw:on-button-click pause-button
        (lambda (inst)
          (declare (ignore inst))
          (ngl:pause trajectory)))
      (jw:on-button-click stop-button
        (lambda (inst)
          (declare (ignore inst))
          (ngl:stop trajectory)))
      (jw:link trajectory :frame frame-slider :value)
      (jw:observe trajectory :count
        (lambda (inst type name old-value new-value source)
          (declare (ignore inst type name old-value source))
          (setf (jw:widget-max frame-slider) (1- new-value))))
      (jw:observe mode-button :value
        (lambda (inst type name old-value new-value source)
          (declare (ignore inst type name old-value source))
          (setf (ngl:mode trajectory) (if new-value "loop" "once"))))
      (jw:observe trajectory :is-running
        (lambda (inst type name old-value new-value source)
          (declare (ignore inst type name old-value source))
          (unless new-value
            (setf (jw:widget-value play-back-2) nil
                  (jw:widget-value play-back-1) nil
                  (jw:widget-value play-fore-1) nil
                  (jw:widget-value play-fore-2) nil))))
      (jw:observe play-back-2 :value
        (lambda (inst type name old-value new-value source)
          (declare (ignore inst type name old-value source))
          (cond
            (new-value
              (setf (jw:widget-value play-back-1) nil
                    (jw:widget-value play-fore-1) nil
                    (jw:widget-value play-fore-2) nil
                    (ngl:%step trajectory) 10
                    (ngl:direction trajectory) "backward")
              (ngl:play trajectory))
            ((not (or (jw:widget-value play-back-1)
                      (jw:widget-value play-fore-1)
                      (jw:widget-value play-fore-2)))
              (ngl:pause trajectory)))))
      (jw:observe play-back-1 :value
        (lambda (inst type name old-value new-value source)
          (declare (ignore inst type name old-value source))
          (cond
            (new-value
              (setf (jw:widget-value play-back-2) nil
                    (jw:widget-value play-fore-1) nil
                    (jw:widget-value play-fore-2) nil
                    (ngl:%step trajectory) 1
                    (ngl:direction trajectory) "backward")
              (ngl:play trajectory))
            ((not (or (jw:widget-value play-back-1)
                      (jw:widget-value play-fore-1)
                      (jw:widget-value play-fore-2)))
              (ngl:pause trajectory)))))
      (jw:observe play-fore-1 :value
        (lambda (inst type name old-value new-value source)
          (declare (ignore inst type name old-value source))
          (cond
            (new-value
              (setf (jw:widget-value play-back-2) nil
                    (jw:widget-value play-back-1) nil
                    (jw:widget-value play-fore-2) nil
                    (ngl:%step trajectory) 1
                    (ngl:direction trajectory) "forward")
              (ngl:play trajectory))
            ((not (or (jw:widget-value play-back-2)
                      (jw:widget-value play-back-1)
                      (jw:widget-value play-fore-2)))
              (ngl:pause trajectory)))))
      (jw:observe play-fore-2 :value
        (lambda (inst type name old-value new-value source)
          (declare (ignore inst type name old-value source))
          (cond
            (new-value
              (setf (jw:widget-value play-back-2) nil
                    (jw:widget-value play-back-1) nil
                    (jw:widget-value play-fore-1) nil
                    (ngl:%step trajectory) 10
                    (ngl:direction trajectory) "forward")
              (ngl:play trajectory))
            ((not (or (jw:widget-value play-back-2)
                      (jw:widget-value play-back-1)
                      (jw:widget-value play-fore-1)))
              (ngl:pause trajectory)))))
      (list play-back-2 play-back-1 stop-button pause-button play-fore-1 play-fore-2 mode-button frame-slider))))

(defclass ngl-pane ()
  ((stage :reader ngl-pane-stage
          :initform (ngl:make-stage :clip-dist 0 ;:background-color "white"
                                    :layout (make-instance 'jw:layout :width "100%" :height "auto"
                                                           :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                           :grid-row "row1-start / last-line"
                                                           :grid-area "stage")))
   (grid :reader ngl-pane-grid
         :initform (resizable-box:make-resizable-grid-box
                     :enable-full-screen t
                     :layout (resizable-box:make-resizable-layout
                               :grid-gap "1em"
                               :overflow "hidden"
                               :padding "0 24px 0 0"
                               :grid-auto-rows "min-content"
                               :grid-auto-flow "row"
                               :grid-template-rows "1fr"
                               :grid-template-columns "min-content min-content 1fr"
                               :grid-template-areas "'stage stage stage'")))))

(defmethod initialize-instance :after ((instance ngl-pane) &rest initargs &key pane &allow-other-keys)
  (declare (ignore initargs))
  (setf (jw:widget-children (ngl-pane-grid instance))
        (list (ngl-pane-stage instance)))
  (let ((layout (jw:widget-layout (ngl-pane-grid instance))))
    (if pane
        (setf (jw:widget-height layout) "100%")
        (setf (resizable-box:resize layout) "vertical"
              (jw:widget-min-height layout) "480px"))))

(defun make-representations (representation)
  (list (make-instance 'ngl:backbone
                       :name "Backbone"
                       :visible (eq representation :backbone)
                       :lazy t)
        (make-instance 'ngl:ball-and-stick
                       :name "Ball and Stick"
                       :visible (eq representation :ball-and-stick)
                       :lazy t)
        (make-instance 'ngl:cartoon
                       :name "Cartoon"
                       :color-scheme "residueindex"
                       :visible (eq representation :cartoon)
                       :lazy t)
        (make-instance 'ngl:licorice
                       :name "Licorice"
                       :visible (eq representation :licorice)
                       :lazy t)
        (make-instance 'ngl:line
                       :name "Line"
                       :visible (eq representation :line)
                       :lazy t)
        (make-instance 'ngl:ribbon
                       :name "Ribbon"
                       :color-scheme "residueindex"
                       :visible (eq representation :ribbon)
                       :lazy t)
        (make-instance 'ngl:spacefill
                       :name "Spacefill"
                       :visible (eq representation :spacefill)
                       :lazy t)
        (make-instance 'ngl:surface
                       :name "Surface" :use-worker t
                       :color-scheme "residueindex"
                       :visible (eq representation :surface)
                       :lazy t)))

(defun ngl-show-on-pane (pane-instance object &rest rest
                         &key pane append (representation :ball-and-stick) &allow-other-keys)
  (unless pane-instance
    (setf pane-instance (make-instance 'ngl-pane :pane pane))
    (j:display (ngl-pane-grid pane-instance)))
  (unless append
    (setf (jw:widget-children (ngl-pane-grid pane-instance))
          (list (first (jw:widget-children (ngl-pane-grid pane-instance))))
          (ngl:components (ngl-pane-stage pane-instance)) nil))
  (let* ((representations (make-representations representation))
         (component (make-ngl-structure object :auto-view-duration 0 :representations representations))
         (auto-view-button (jw:make-button :description "Auto View"
                                           :style (jw:make-description-style :description-width "min-content")
                                           :layout (jw:make-layout :align-self "center"
                                                                   :width "min-content")))
         (representation-dropdown (jw:make-dropdown :description "Representation"
                                                    :%options-labels (mapcan (lambda (rep &aux (name (ngl:name rep)))
                                                                               (when name (list name)))
                                                                             representations)
                                                    :index (position-if #'ngl:visible representations)
                                                    :style (jw:make-description-style :description-width "min-content")
                                                    :layout (jw:make-layout :align-self "center"
                                                                            :width "max-content"))))
    (add-bounding-box component object)
    (jw:on-button-click auto-view-button
      (lambda (inst)
        (declare (ignore inst))
        (ngl:auto-view component 1000)))
    (jw:observe representation-dropdown :value
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (dolist (representation (ngl:representations component))
          (when (ngl:name representation)
            (setf (ngl:visible representation) (equalp new-value (ngl:name representation)))))))
    (setf (ngl:components (ngl-pane-stage pane-instance)) (append (ngl:components (ngl-pane-stage pane-instance))
                                                         (list component))
          (jw:widget-children (ngl-pane-grid pane-instance)) (append (jw:widget-children (ngl-pane-grid pane-instance))
                                                            (list (jw:make-label :value (or (chem:get-name object) "")
                                                                                 :style (jw:make-description-style :description-width "min-content")
                                                                                 :layout (jw:make-layout :align-self "center"
                                                                                                         :width "min-content"))
                                                                  auto-view-button
                                                                  representation-dropdown))))
  pane-instance)

(defun ngl-show (instance &rest rest &key &allow-other-keys)
  (let* ((axes (make-axes (getf rest :axes)))
         (component (make-ngl-structure instance
                                        :auto-view-duration 0
                                        :representations (list (make-instance 'ngl:backbone
                                                                              :name "Backbone"
                                                                              :visible nil :lazy t)
                                                               (make-instance 'ngl:ball-and-stick
                                                                              :name "Ball and Stick"
                                                                              :visible t :lazy t)
                                                               (make-instance 'ngl:cartoon
                                                                              :name "Cartoon"
                                                                              :color-scheme "residueindex"
                                                                              :visible nil :lazy t)
                                                               (make-instance 'ngl:licorice
                                                                              :name "Licorice"
                                                                              :visible nil :lazy t)
                                                               (make-instance 'ngl:line
                                                                              :name "Line"
                                                                              :visible nil :lazy t)
                                                               (make-instance 'ngl:ribbon
                                                                              :name "Ribbon"
                                                                              :color-scheme "residueindex"
                                                                              :visible nil :lazy t)
                                                               (make-instance 'ngl:spacefill
                                                                              :name "Spacefill"
                                                                              :visible nil :lazy t)
                                                               (make-instance 'ngl:surface
                                                                              :name "Surface" :use-worker t
                                                                              :color-scheme "residueindex"
                                                                              :visible nil :lazy t)
                                                               axes)))

         (stage (apply #'make-instance 'ngl:stage
                       :clip-dist 0
                       :background-color "white"
                       :components (list component
                                         (make-instance 'ngl:shape
                                                        :primitives (getf rest :shapes)
                                                        :representations (list (make-instance 'ngl:buffer-representation
                                                                                              :opacity (getf rest :shapes-opacity 1d0)))))
                       :layout (make-instance 'jw:layout
                                              :width "100%"
                                              :height "auto"
                                              :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                              :grid-area "stage")
                       rest))
         (representation-dropdown (make-instance 'jw:dropdown
                                                 :description "Representation"
                                                 :%options-labels (mapcan (lambda (rep &aux (name (ngl:name rep)))
                                                                            (when name (list name)))
                                                                          (ngl:representations component))
                                                 :index (position-if #'ngl:visible (ngl:representations component))
                                                 :style (make-instance 'jw:description-style
                                                                       :description-width "min-content")
                                                 :layout (make-instance 'jw:layout
                                                                        :margin ".5em"
                                                                        :width "max-content")))
         (auto-view-button (make-instance 'jw:button
                                          :description "Auto View"
                                          :style (make-instance 'jw:description-style
                                                                :description-width "min-content")
                                          :layout (make-instance 'jw:layout
                                                                 :margin ".5em"
                                                                 :width "max-content")))
         (axes-button (jw:make-toggle-button :icon "cube"
                                             :tooltip "Toggle Axes"
                                             :value (getf rest :axes)
                                             :style (make-instance 'jw:description-style
                                                                       :description-width "min-content")
                                             :layout (jw:make-layout :margin ".5em .1em .5em .1em"
                                                                     :width "max-content"))))
    (when (and (typep instance 'chem:aggregate)
               (chem:bounding-box-bound-p instance))
      (add-bounding-box component instance))
    (jw:observe representation-dropdown :value
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (dolist (representation (ngl:representations component))
          (when (ngl:name representation)
            (setf (ngl:visible representation) (equalp new-value (ngl:name representation)))))))
    (jw:observe axes-button :value
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (setf (ngl:visible axes) new-value)))
    (jw:on-button-click auto-view-button
      (lambda (inst)
        (declare (ignore inst))
        (ngl:auto-view stage 1000)))
    (make-instance 'resizable-box:resizable-grid-box
                   :children (list stage
                                   (make-instance 'jw:box
                                                  :children (append (list representation-dropdown
                                                                          auto-view-button
                                                                          axes-button)
                                                                    (ngl-show-trajectory (first (ngl:trajectories component))))
                                                  :layout (make-instance 'jw:layout
                                                                         :flex-flow "row wrap"
                                                                         :justify-content "center"
                                                                         ;:margin "-.5em"
                                                                         :align-items "baseline"
                                                                         :align-content "flex-start"
                                                                         :grid-area "controls")))
                   :enable-full-screen t
                   :layout (make-instance 'resizable-box:resizable-layout
                                          :resize "vertical"
                                          :grid-gap "1em"
                                          :min-height "480px"
                                          :overflow "hidden"
                                          :padding "0 24px 0 0"
                                          :grid-template-rows "1fr min-content"
                                          :grid-template-columns "1fr"
                                          :grid-template-areas "'stage' 'controls"))))

(defun isolate-residue (residue)
  (let* ((agg (chem:make-aggregate nil))
         (mol (chem:make-molecule nil))
         (new-to-old (make-hash-table))
         (res (chem:copy residue new-to-old)))
    (let (break-bonds)
      (chem:map-bonds
       'nil
       (lambda (a1 a2 order bbond)
         (declare (ignore bbond))
         (unless (and (chem:contains-atom res a1) (chem:contains-atom res a2))
           (push (list a1 a2) break-bonds)))
       res)
      (loop for bond in break-bonds
            for a1 = (first bond)
            for a2 = (second bond)
            do (chem:remove-bond-to a1 a2)))
    (chem:add-matter agg mol)
    (chem:add-matter mol res)
    agg))

(defmethod show-on-pane (pane-instance (object chem:aggregate) &rest rest &key &allow-other-keys)
  (apply #'ngl-show-on-pane pane-instance object rest))

(defmethod show-on-pane (pane-instance (object chem:molecule) &rest rest &key &allow-other-keys)
  (apply #'ngl-show-on-pane pane-instance object rest))

(defmethod show-on-pane (pane-instance (object chem:residue) &rest rest &key &allow-other-keys)
  (apply #'ngl-show-on-pane pane-instance (isolate-residue object) rest))

(defmethod show-on-pane (pane-instance (object dynamics:trajectory) &rest rest &key &allow-other-keys)
  (apply #'ngl-show-on-pane pane-instance object rest))

(defmethod show-on-pane (pane-instance (object dynamics:simulation) &rest rest &key &allow-other-keys)
  (apply #'ngl-show-on-pane pane-instance (dynamics:make-trajectory object) rest))

(defmethod show-on-pane (pane-instance (object leap.topology:amber-topology-trajectory-pair) &rest rest &key &allow-other-keys)
  (apply #'ngl-show-on-pane pane-instance object rest))
