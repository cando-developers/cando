(in-package :cando-widgets)

(defgeneric show (instance &rest kwargs &key &allow-other-keys)
  (:documentation "Display a residue, molecule, aggregate or a trajectory in a Jupyter notebook"))


(defun add-bounding-box (component aggregate)
  (unless (chem:bounding-box-bound-p aggregate)
    (error "There is no bounding-box defined for aggregate"))
  (let* ((bounding-box (chem:bounding-box aggregate))
         (center (chem:get-bounding-box-center bounding-box))
         (x-width/2 (/ (chem:get-x-width bounding-box) 2.0))
         (y-width/2 (/ (chem:get-y-width bounding-box) 2.0))
         (z-width/2 (/ (chem:get-z-width bounding-box) 2.0))
         (minx (- (geom:vx center) x-width/2))
         (miny (- (geom:vy center) y-width/2))
         (minz (- (geom:vz center) z-width/2))
         (maxx (+ (geom:vx center) x-width/2))
         (maxy (+ (geom:vy center) y-width/2))
         (maxz (+ (geom:vz center) z-width/2)))
    (setf (ngl:representations component)
          (append (ngl:representations component)
                  (list (make-instance 'ngl:buffer-representation
                                       :name "bounding-box"
                                       :buffer (list :type "wideline"
                                                     :position1 (vector minx miny minz
                                                                        minx miny minz
                                                                        minx miny minz
                                                                        minx miny maxz
                                                                        minx miny maxz
                                                                        minx maxy minz
                                                                        minx maxy minz
                                                                        maxx miny minz
                                                                        maxx miny minz
                                                                        minx maxy maxz
                                                                        maxx miny maxz
                                                                        maxx maxy minz)
                                                     :position2 (vector minx miny maxz
                                                                        minx maxy minz
                                                                        maxx miny minz
                                                                        minx maxy maxz
                                                                        maxx miny maxz
                                                                        minx maxy maxz
                                                                        maxx maxy minz
                                                                        maxx miny maxz
                                                                        maxx maxy minz
                                                                        maxx maxy maxz
                                                                        maxx maxy maxz
                                                                        maxx maxy maxz)
                                                     :color #(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
                                                     :color2 #(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0))))))))


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


(defun ngl-show (instance &rest kwargs &key &allow-other-keys)
  (let* ((component (make-ngl-structure instance
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
                                                                              :visible nil :lazy t))))
         (stage (apply #'make-instance 'ngl:stage
                       :clip-dist 0
                       :background-color "white"
                       :components (list component
                                         (make-instance 'ngl:shape :primitives (getf kwargs :shapes)))
                       :layout (make-instance 'jw:layout
                                              :width "100%"
                                              :height "auto"
                                              :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                              :grid-area "stage")
                       kwargs))
         (representation-dropdown (make-instance 'jw:dropdown
                                                 :description "Representation"
                                                 :%options-labels (mapcar #'ngl:name (ngl:representations component))
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
                                                                 :width "max-content"))))
    (when (and (typep instance 'chem:aggregate)
               (chem:bounding-box-bound-p instance))
      (add-bounding-box component instance))
    (jw:observe representation-dropdown :value
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (dolist (representation (ngl:representations component))
          (setf (ngl:visible representation) (equalp new-value (ngl:name representation))))))
    (jw:on-button-click auto-view-button
      (lambda (inst)
        (declare (ignore inst))
        (ngl:auto-view stage 1000)))
    (make-instance 'resizable-box:resizable-grid-box
                   :children (list stage
                                   (make-instance 'jw:box
                                                  :children (append (list representation-dropdown
                                                                          auto-view-button)
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


(defmethod show ((instance chem:aggregate) &rest kwargs &key &allow-other-keys)
  (apply #'ngl-show instance kwargs))

(defmethod show ((instance chem:molecule) &rest kwargs &key &allow-other-keys)
  (apply #'ngl-show instance kwargs))

(defun isolate-residue (residue)
  (let* ((agg (chem:make-aggregate nil))
         (mol (chem:make-molecule nil))
         (new-to-old (make-hash-table))
         (res (chem:copy residue new-to-old)))
    (let (break-bonds)
      (chem:map-bonds
       'nil
       (lambda (a1 a2 order)
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

(defmethod show ((residue chem:residue)  &rest kwargs &key &allow-other-keys)
  (let ((agg (isolate-residue residue)))
    (apply 'show agg kwargs)))

(defmethod show ((sketch sketch2d:sketch2d) &rest kwargs &key &allow-other-keys)
  (show (sketch2d:svg sketch)))

(defmethod show ((sketch sketch2d:sketch-svg) &rest kwargs &key &allow-other-keys)
  (jupyter:svg (sketch2d:render-svg-to-string sketch)))

(defmethod show ((instance dynamics:trajectory) &rest kwargs &key &allow-other-keys)
  (apply #'ngl-show instance kwargs))

(defmethod show ((instance dynamics:simulation) &rest kwargs &key &allow-other-keys)
  (apply #'ngl-show (dynamics:make-trajectory instance) kwargs))
