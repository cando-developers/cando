(in-package :cando-widgets)

(defparameter *panes* (make-hash-table :test #'equal))

(defparameter *default-pane* nil)

(defgeneric make-pane (object &rest rest &key &allow-other-keys))

(defgeneric show-on-pane (pane-instance object &rest rest &key &allow-other-keys))

(defun default-pane (pane)
  (setf *default-pane* pane)
  (values))

(defun show (object &rest rest &key (pane *default-pane*) &allow-other-keys)
  "Display a residue, molecule, aggregate or a trajectory in a Jupyter notebook"
  (if pane
      (let ((existing (gethash pane *panes*)))
        (if existing
            (apply #'show-on-pane existing object rest)
            (let ((sidecar (jw:make-sidecar :title pane
                                            :on-close (list (lambda (&rest rest)
                                                              (declare (ignore rest))
                                                              (remhash pane *panes*))))))
              (jw:with-output sidecar
                (setf (gethash pane *panes*)
                      (apply #'show-on-pane nil object rest)))
              (push (lambda (inst)
                      (declare (ignore inst))
                      (remhash pane *panes*))
                    (ngl:on-remove (ngl-pane-stage (gethash pane *panes*)))))))
      (apply #'show-on-pane nil object rest))
  (values))

;;; sketch2d

(defmethod show-on-pane (pane-instance (sketch sketch2d:sketch2d) &rest rest &key &allow-other-keys)
  (apply #'show-on-pane pane-instance (sketch2d:svg sketch) rest))

(defmethod show-on-pane (pane-instance (sketch sketch2d:sketch-svg) &rest rest &key pane &allow-other-keys)
  (declare (ignore rest))
  (jupyter:svg (sketch2d:render-svg-to-string sketch) :display t :id pane :update (and pane-instance t))
  pane)

(defvar +default-color+ '(1 .65 0))
(defvar +default-radius+ 0.15)
(defvar +default-text-size+ 2)

(defun cylinder-buffer (position1 position2 &key radius color)
  (list :type "cylinder"
        :position1 (apply #'append position1)
        :position2 (apply #'append position2)
        :radius (or radius (make-list (length position1) :initial-element +default-radius+))
        :color (if color
                   (apply 'append color)
                   (loop for p on position1
                         append +default-color+))
        :color2 (if color
                    (apply 'append color)
                    (loop for p on position1
                          append +default-color+))))

(defun sphere-buffer (position &key radius color)
  (list :type "sphere"
        :position (apply #'append position)
        :radius (or radius (make-list (length position) :initial-element +default-radius+))
        :color (if color
                   (apply 'append color)
                   (loop for p on position
                         append +default-color+))))

(defun text-buffer (position text &key size radius color)
  (list :type "text"
        :text text
        :position (apply #'append position)
        :size (or size (make-list (length position) :initial-element +default-text-size+))
        :color (if color
                   (apply 'append color)
                   (loop for p on position
                         append +default-color+))))

(defun add-bounding-box (component aggregate boxp axesp)
  (unless (chem:bounding-box-bound-p aggregate)
    (leap.set-box::calculate-bounding-box aggregate t))
  (let* ((bounding-box (chem:bounding-box aggregate))
         (center-of-mass (chem:center-of-mass aggregate))
         (min (chem:min-corner bounding-box))
         (max (chem:max-corner bounding-box))
         (midx (geom:get-x center-of-mass))
         (midy (geom:get-y center-of-mass))
         (midz (geom:get-z center-of-mass))
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
                                     p110 p101 p110 p110 p101 p011)))
         (sph (sphere-buffer (list p000 p001 p010 p011 p100 p101 p110 p111)))
         (axes-cyl (cylinder-buffer (list (list minx midy midz)
                                          (list midx miny midz)
                                          (list midx midy minz))
                                    (list (list maxx midy midz)
                                          (list midx maxy midz)
                                          (list midx midy maxz))
                                    :color '((.8 0 0) (0 .8 0) (0 0 .8))))
         (axes-sph (sphere-buffer (list (list minx midy midz)
                                        (list midx miny midz)
                                        (list midx midy minz)
                                        (list maxx midy midz)
                                        (list midx maxy midz)
                                        (list midx midy maxz))
                                    :color '((.8 0 0) (0 .8 0) (0 0 .8) (.8 0 0) (0 .8 0) (0 0 .8))))
         (axes-text (text-buffer (list (list (+ maxx +default-radius+) midy midz)
                                       (list midx (+ maxy +default-radius+) midz)
                                       (list midx midy (+ maxz +default-radius+)))
                                 '("x" "y" "z")
                                 :color '((.8 0 0) (0 .8 0) (0 0 .8)))))
    (setf (ngl:representations component)
          (append (ngl:representations component)
                  (list (ngl:make-buffer-representation :name "_box" :lazy t
                                                        :visible boxp :buffer cyl)
                        (ngl:make-buffer-representation :name "_box" :lazy t
                                                        :visible boxp :buffer sph)
                        (ngl:make-buffer-representation :name "_axes" :lazy t
                                                        :visible axesp :buffer axes-cyl)
                        (ngl:make-buffer-representation :name "_axes" :lazy t
                                                        :visible axesp :buffer axes-sph)
                        (ngl:make-buffer-representation :name "_axes" :lazy t
                                                        :visible axesp :buffer axes-text))))))

(defun ngl-show-trajectory (trajectory)
  (when trajectory
    (let* ((desc-style (jw:make-description-style :description-width "min-content"))
           (stop-button (jw:make-button :icon "stop"
                                        :tooltip "Stop"
                                        :style desc-style
                                        :layout (jw:make-layout :width "max-content"
                                                                :grid-area "stop")))
           (pause-button (jw:make-button :icon "pause"
                                         :tooltip "Pause"
                                         :style desc-style
                                         :layout (jw:make-layout :width "max-content"
                                                                 :grid-area "pause")))
           (play-fore-1 (jw:make-toggle-button :icon "play" :tooltip "Foreward"
                                               :style desc-style
                                               :layout (jw:make-layout :width "max-content"
                                                                       :grid-area "fore")))
           (play-fore-2 (jw:make-toggle-button :icon "forward" :tooltip "Fast foreward"
                                               :style desc-style
                                               :layout (jw:make-layout :width "max-content"
                                                                       :grid-area "fastfore")))
           (mode-button (jw:make-toggle-button :icon "retweet"
                                               :tooltip "Loop"
                                               :value t
                                               :style desc-style
                                               :layout (jw:make-layout :width "max-content"
                                                                       :grid-area "loop")))
           (frame-slider (jw:make-int-slider :style desc-style
                                             :layout (jw:make-layout :align-self "center"
                                                                     :width "18em"
                                                                     :grid-area "frame"))))
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
            (setf (jw:widget-value play-fore-1) nil
                  (jw:widget-value play-fore-2) nil))))
      (jw:observe play-fore-1 :value
        (lambda (inst type name old-value new-value source)
          (declare (ignore inst type name old-value source))
          (cond (new-value
                 (setf (jw:widget-value play-fore-2) nil
                       (ngl:%step trajectory) 1
                       (ngl:direction trajectory) "forward")
                 (ngl:play trajectory))
                ((not (jw:widget-value play-fore-2))
                 (ngl:pause trajectory)))))
      (jw:observe play-fore-2 :value
        (lambda (inst type name old-value new-value source)
          (declare (ignore inst type name old-value source))
          (cond (new-value
                 (setf (jw:widget-value play-fore-1) nil
                       (ngl:%step trajectory) 10
                       (ngl:direction trajectory) "forward")
                 (ngl:play trajectory))
                ((not (jw:widget-value play-fore-1))
                 (ngl:pause trajectory)))))
      (list stop-button pause-button play-fore-1 play-fore-2 mode-button frame-slider))))

(defclass ngl-pane ()
  ((auto-view :reader ngl-pane-auto-view
              :initform (jw:make-button :icon "compress"
                                        :tooltip "Auto View"
                                        :layout (jw:make-layout :align-self "center"
                                                                :grid-area "autoview"
                                                                :width "max-content")))
   (camera :reader ngl-pane-camera
           :initform (jw:make-button :icon "camera"
                                     :tooltip "Capture image"
                                     :layout (jw:make-layout :align-self "center"
                                                             :grid-area "pic"
                                                             :width "max-content")))
   (stage :reader ngl-pane-stage
          :initform (ngl:make-stage :fog-near 100
                                    :layout (jw:make-layout :width "100%" :height "auto"
                                                            :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                            :grid-area "stage")))
   (grid :reader ngl-pane-grid
         :initform (resizable-box:make-resizable-grid-box
                     :enable-full-screen t
                     :layout (resizable-box:make-resizable-layout
                               ;:grid-gap "1em"
                               :overflow "hidden"
                               :padding "0 24px 0 0"
                               :grid-auto-rows "min-content"
                               :grid-template-rows "min-content 1fr"
                               :grid-template-columns "repeat(12, min-content) 1fr"
                               :grid-template-areas "'label    autoview pic   .     .     .          .     .     .      .     .         .     .      .'
                                                     'stage    stage    stage stage stage stage      stage stage stage  stage stage     stage stage  stage'
                                                     'label0   visible0 box0  axes0 repr0 selelabel0 sele0 stop0 pause0 fore0 fastfore0 loop0 frame0 .'
                                                     'label1   visible1 box1  axes1 repr1 selelabel1 sele1 stop1 pause1 fore1 fastfore1 loop1 frame1 .'
                                                     'label2   visible2 box2  axes2 repr2 selelabel2 sele2 stop2 pause2 fore2 fastfore2 loop2 frame2 .'
                                                     'label3   visible3 box3  axes3 repr3 selelabel3 sele3 stop3 pause3 fore3 fastfore3 loop3 frame3 .'
                                                     'label4   visible4 box4  axes4 repr4 selelabel4 sele4 stop4 pause4 fore4 fastfore4 loop4 frame4 .'
                                                     'label5   visible5 box5  axes5 repr5 selelabel5 sele5 stop5 pause5 fore5 fastfore5 loop5 frame5 .'
                                                     'label6   visible6 box6  axes6 repr6 selelabel6 sele6 stop6 pause6 fore6 fastfore6 loop6 frame6 .'
                                                     'label7   visible7 box7  axes7 repr7 selelabel7 sele7 stop7 pause7 fore7 fastfore7 loop7 frame7 .'
                                                     'label8   visible8 box8  axes8 repr8 selelabel8 sele8 stop8 pause8 fore8 fastfore8 loop8 frame8 .'
                                                     'label9   visible9 box9  axes9 repr9 selelabel9 sele9 stop9 pause9 fore9 fastfore9 loop9 frame9 .'")))))

(defmethod initialize-instance :after ((instance ngl-pane) &rest initargs &key pane &allow-other-keys)
  (declare (ignore initargs))
  (setf (jw:widget-children (ngl-pane-grid instance))
        (list (jw:make-label :value "Stage"
                             :style (jw:make-description-style :description-width "min-content")
                             :layout (jw:make-layout :align-self "center"
                                                     :grid-area "label"
                                                     :width "min-content"))
              (ngl-pane-auto-view instance)
              (ngl-pane-camera instance)
              (ngl-pane-stage instance)))
  (jw:on-button-click (ngl-pane-auto-view instance)
    (lambda (inst)
      (declare (ignore inst))
      (ngl:auto-view (ngl-pane-stage instance) 1000)))
  (jw:on-button-click (ngl-pane-camera instance)
    (lambda (inst)
      (declare (ignore inst))
      (ngl:make-image (ngl-pane-stage instance)
                      (lambda (data type)
                        (declare (ignore type))
                        (multiple-value-bind (second minute hour date month year day daylight-p zone)
                            (decode-universal-time (get-universal-time))
                          (declare (ignore date daylight-p zone))
                          (with-open-file (stream (format nil "stage ~4,'0d-~2,'0d-~2,'0d ~2,'0d-~2,'0d-~2,'0d.png"
                                                          year month day hour minute second)
                                           :direction :output)
                            (write-sequence data stream))))
                      :transparent t :antialias t)))
  (let ((layout (jw:widget-layout (ngl-pane-grid instance))))
    (if (or pane *default-pane*)
        (setf (jw:widget-height layout) "100%")
        (setf (resizable-box:resize layout) "vertical"
              (jw:widget-min-height layout) "480px"))))

(defun visualize-named-atom-map (from to mapping)
  (ngl:make-shape :name (format nil "~a -> ~a" (chem:get-name from) (chem:get-name to))
                  :primitives (loop for (from-name . to-name) in mapping
                                    for from-atom = (chem:first-atom-with-name from from-name)
                                    for from-position = (chem:get-position from-atom)
                                    for to-atom = (chem:first-atom-with-name to to-name)
                                    for to-position = (chem:get-position to-atom)
                                    collect (list :type "arrow"
                                                  :name (format nil "~a -> ~a" from-name to-name)
                                                  :color '(.65 0 1)
                                                  :radius 0.2
                                                  :position1 (list (geom:get-x from-position)
                                                                   (geom:get-y from-position)
                                                                   (geom:get-z from-position))
                                                  :position2 (list (geom:get-x to-position)
                                                                   (geom:get-y to-position)
                                                                   (geom:get-z to-position))))))

(defun anchor-named-atom-map (from to mapping)
  (loop for (from-name . to-name) in mapping
        for from-atom = (chem:first-atom-with-name from from-name)
        for to-atom = (chem:first-atom-with-name to to-name)
        for to-position = (chem:get-position to-atom)
        do (cando:anchor-atom from-atom to-position)))

(defun make-representations (representation sele axes)
  (list (ngl:make-backbone :name "Backbone" :sele sele :lazy t
                           :visible (eq representation :backbone))
        (ngl:make-ball-and-stick :name "Ball and Stick" :lazy t :sele sele
                                 :visible (eq representation :ball-and-stick))
        (ngl:make-cartoon :name "Cartoon" :lazy t :sele sele
                          :color-scheme "residueindex"
                          :visible (eq representation :cartoon))
        (ngl:make-licorice :name "Licorice" :lazy t :sele sele
                           :visible (eq representation :licorice))
        (ngl:make-line :name "Line" :lazy t :sele sele
                       :visible (eq representation :line))
        (ngl:make-ribbon :name "Ribbon" :lazy t :sele sele
                         :color-scheme "residueindex"
                         :visible (eq representation :ribbon))
        (ngl:make-spacefill :name "Spacefill" :lazy t :sele sele
                            :visible (eq representation :spacefill))
        (ngl:make-surface :name "Surface" :lazy t :sele sele :use-worker t
                          :color-scheme "residueindex"
                          :visible (eq representation :surface))))

(defun ngl-selectable-repr-p (representation)
  (and (ngl:name representation)
       (char/= (char (ngl:name representation) 0) #\_)))

(defun show-component (pane-instance component controls &rest rest
                       &key pane append (background "white")
                       &allow-other-keys
                       &aux (display (not pane-instance)))
  (unless pane-instance
    (setf pane-instance (make-instance 'ngl-pane :pane pane)))
  (setf (ngl:background-color (ngl-pane-stage pane-instance)) background)
  (unless append
    (setf (cddddr (jw:widget-children (ngl-pane-grid pane-instance))) nil
          (ngl:components (ngl-pane-stage pane-instance)) nil))
  (let* ((view-button (jw:make-toggle-button :icon "eye"
                                             :tooltip "Toggle Visibility"
                                             :value t
                                             :layout (jw:make-layout :align-self "center"
                                                                     :grid-area "visible"
                                                                     :width "min-content"))))
    (push view-button controls)
    (loop with count = (length (ngl:components (ngl-pane-stage pane-instance)))
          for control in controls
          for layout = (jw:widget-layout control)
          do (setf (jw:widget-grid-area layout)
                   (format nil "~a~a" (jw:widget-grid-area layout) count)))
    (jw:observe view-button :value
      (lambda (inst type name old-value new-value source)
        (declare (ignore type name old-value source))
        (setf (ngl:visible component) new-value
              (jw:widget-icon inst) (if new-value "eye" "eye-slash"))))
    (setf (ngl:components (ngl-pane-stage pane-instance)) (append (ngl:components (ngl-pane-stage pane-instance))
                                                          (list component))
          (jw:widget-children (ngl-pane-grid pane-instance)) (append (jw:widget-children (ngl-pane-grid pane-instance))
                                                                     controls)))
  (when display
    (j:display (ngl-pane-grid pane-instance)))
  pane-instance)

(defun ngl-show-on-pane (pane-instance object &rest rest
                         &key append (representation :ball-and-stick)
                              box sele axes
                         &allow-other-keys
                         &aux (display (not pane-instance))
                              (representations (make-representations representation sele axes)))
  (multiple-value-bind (component agg)
      (make-ngl-structure object :auto-view-duration 0 :representations representations)
    (let* ((traj-controls (ngl-show-trajectory (first (ngl:trajectories component))))
           (box-button (jw:make-toggle-button :icon "cube"
                                              :tooltip "Toggle Bounding Box"
                                              :value box
                                              :layout (jw:make-layout :align-self "center"
                                                                      :grid-area "box"
                                                                      :width "min-content")))
           (axes-button (jw:make-toggle-button :icon "arrows"
                                               :tooltip "Toggle Axes"
                                               :value axes
                                               :layout (jw:make-layout :align-self "center"
                                                                       :grid-area "axes"
                                                                       :width "min-content")))
           (representation-dropdown (jw:make-dropdown :tooltip "Representation"
                                                      :%options-labels (mapcan (lambda (rep)
                                                                                 (when (ngl-selectable-repr-p rep)
                                                                                   (list (ngl:name rep))))
                                                                               representations)
                                                      :index (position-if #'ngl:visible representations)
                                                      :layout (jw:make-layout :align-self "center"
                                                                              :grid-area "repr"
                                                                              :width "max-content")))
           (sele-text (jw:make-text :tooltip "Selection" :value (or sele "") :continuous-update nil
                                    :layout (jw:make-layout :align-self "center"
                                                            :grid-area "sele"
                                                            :width "8em")))
           (controls (list* (jw:make-label :value (or (chem:get-name agg) "")
                                           :style (jw:make-description-style :description-width "min-content")
                                           :layout (jw:make-layout :align-self "center"
                                                                   :grid-area "label"
                                                                   :width "min-content"))
                            box-button
                            axes-button
                            representation-dropdown
                            (jw:make-html :value "<span title='Selection Filter' class='fa fa-filter'/>"
                                          :style (jw:make-description-style :description-width "min-content")
                                          :layout (jw:make-layout :align-self "center"
                                                                  :margin "0 0 0 .5em"
                                                                  :grid-area "selelabel"
                                                                  :width "min-content"))
                            sele-text
                            traj-controls)))
      (add-bounding-box component agg box axes)
      (jw:observe sele-text :value
        (lambda (inst type name old-value new-value source)
          (declare (ignore type name old-value source))
          (loop for representation in (ngl:representations component)
                when (ngl-selectable-repr-p representation)
                  do (setf (ngl:sele representation) new-value))))
      (jw:observe box-button :value
        (lambda (inst type name old-value new-value source)
          (declare (ignore inst type name old-value source))
          (loop for representation in (ngl:representations component)
                when (equal "_box" (ngl:name representation))
                  do (setf (ngl:visible representation) new-value))))
      (jw:observe axes-button :value
        (lambda (inst type name old-value new-value source)
          (declare (ignore inst type name old-value source))
          (loop for representation in (ngl:representations component)
                when (equal "_axes" (ngl:name representation))
                  do (setf (ngl:visible representation) new-value))))
      (jw:observe representation-dropdown :value
        (lambda (inst type name old-value new-value source)
          (declare (ignore inst type name old-value source))
          (dolist (representation (ngl:representations component))
            (when (ngl-selectable-repr-p representation)
              (setf (ngl:visible representation) (equalp new-value (ngl:name representation)))))))
      (apply #'show-component pane-instance component controls rest))))

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

(defmethod show-on-pane (pane-instance (object ngl:shape) &rest rest &key &allow-other-keys)
  (apply #'show-component pane-instance object
         (list (jw:make-label :value (ngl:name object)
                              :style (jw:make-description-style :description-width "min-content")
                              :layout (jw:make-layout :align-self "center"
                                                      :grid-area "label"
                                                      :width "min-content")))
         rest))

