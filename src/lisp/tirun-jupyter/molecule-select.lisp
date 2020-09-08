(in-package :tirun-jupyter)


(defparameter *molecule-sketches* (make-hash-table :test #'eql))

(defun sketch-molecule (molecule)
  (multiple-value-bind (sketch present-p)
                       (gethash molecule *molecule-sketches*)
    (if present-p
      (lparallel:force sketch)
      (setf (gethash molecule *molecule-sketches*)
            (sketch2d:svg (sketch2d:sketch2d molecule))))))

(defparameter *common-sketch* nil)

(defun sketch-molecules (molecules)
  (let (common-sketch)
    (dolist (molecule molecules (values))
      (multiple-value-bind (sketch present-p)
                           (gethash molecule *molecule-sketches*)
        (declare (ignore sketch))
        (unless present-p
          (unless common-sketch
            (setf common-sketch (build-prototype-sketch molecules)))
          (setf *common-sketch* common-sketch)
          (setf (gethash molecule *molecule-sketches*)
                (eval `(lparallel:future
                         (sketch2d:svg (sketch2d:similar-sketch2d ,molecule ,common-sketch))))))))))


(defconstant +min-sketch-width+ 50.0)


(defclass molecule-select (w:has-traits)
  ((container
     :accessor container
     :initform (make-instance 'w:grid-box
                              :layout (make-instance 'w:layout
                                                     :height "640px"
                                                     :overflow-y "scroll"
                                                     :grid-auto-rows "min-content"
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :grid-gap ".25em"
                                                     :padding ".25em"
                                                     :grid-area "list"
                                                     :grid-template-columns "repeat(auto-fill, 200px)")))
   (outer
     :accessor outer
     :initform (make-instance 'w:grid-box
                              :layout (make-instance 'w:layout
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :grid-gap "var(--jp-widgets-container-padding)"
                                                     :padding "var(--jp-widgets-container-padding)"
                                                     :grid-template-columns "1fr min-content min-content 1fr"
                                                     :grid-template-rows "min-content 1fr"
                                                     :grid-template-areas "\". all none .\"
                                                                           \"list list list list\"")))
   (select-all
     :accessor select-all
     :initform (make-instance 'w:button
                              :description "Select All"
                              :layout (make-instance 'w:layout
                                                     :grid-area "all")))
   (select-none
     :accessor select-none
     :initform (make-instance 'w:button
                              :description "Select None"
                              :layout (make-instance 'w:layout
                                                     :grid-area "none")))
   (molecules
     :accessor molecules
     :initarg :molecules
     :initform nil
     :trait :list)
   (selected
     :accessor selected
     :initarg :selected
     :initform nil
     :trait :list)
   (max-width
     :accessor max-width
     :initform +min-sketch-width+
     :trait :float)
   (selection-mode
     :reader selection-mode
     :initarg :selection-mode
     :initform :multiple)
   (task-channel
     :reader task-channel
     :initform (lparallel:make-channel)))
  (:metaclass jupyter-widgets:trait-metaclass))


(defmethod initialize-instance :after ((instance molecule-select) &rest args &key &allow-other-keys)
  (declare (ignore args))
  (with-slots (select-all select-none outer container selection-mode)
              instance
    (w:observe instance :max-width
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (setf (w:widget-grid-template-columns (w:widget-layout (container instance)))
              (format nil "repeat(auto-fill,~Fpx)" (+ 2 new-value)))))
    (setf (w:widget-children outer)
          (if (eql :multiple selection-mode)
            (list select-all select-none container)
            (list container)))
    (w:on-button-click select-all
      (lambda (inst)
        (declare (ignore inst))
        (setf (selected instance) (molecules instance))))
    (w:on-button-click select-none
      (lambda (inst)
        (declare (ignore inst))
        (setf (selected instance) nil)))))


(defmethod jupyter-widgets:%display ((instance molecule-select) &rest args &key &allow-other-keys)
  (declare (ignore args))
  (outer instance))


(defun sketch-and-update (instance html molecule)
  (ignore-errors
    (let ((sketch (sketch-molecule molecule)))
      (setf (max-width instance) (max (max-width instance) (sketch2d:width sketch))
            (w:widget-value html) (sketch2d:render-svg-to-string sketch)))))


(defun make-molecule-view (instance molecule)
  (with-slots (selection-mode)
              instance
    (let ((toggle (if (eql :none selection-mode)
                    (make-instance 'w:label
                                   :layout (make-instance 'w:layout
                                                          :margin "var(--jp-widgets-margin) auto")
                                   :value (symbol-name (chem:get-name molecule)))
                    (make-instance 'w:toggle-button
                                   :layout (make-instance 'w:layout
                                                          :margin "var(--jp-widgets-margin) auto")
                                   :description (symbol-name (chem:get-name molecule)))))
          (html (make-instance 'w:html
                               :layout (make-instance 'w:layout
                                                      :margin "auto"
                                                      :overflow "hidden"))))
      (lparallel:submit-task (task-channel instance) #'sketch-and-update instance html molecule)
      (unless (eql :none selection-mode)
        (w:observe toggle :value
          (lambda (inst type name old-value new-value source)
            (declare (ignore inst type name old-value source))
            (setf (selected instance)
                  (cond
                    ((and new-value
                          (eql :multiple selection-mode))
                      (nconc (selected instance) (list molecule)))
                    (new-value
                      (list molecule))
                    (t
                      (remove molecule (selected instance) :test #'eql))))))
        (w:observe instance :selected
          (lambda (inst type name old-value new-value source)
            (declare (ignore type name old-value source))
            (setf (w:widget-value toggle)
                  (and (position molecule (selected instance) :test #'eql)
                       t)))))
      (make-instance 'w:v-box
                     :layout (make-instance 'w:layout
                                            :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)")
                     :children (list toggle html)))))


(defun render-molecules (instance)
  (setf (max-width instance) +min-sketch-width+
        (w:widget-children (container instance)) (mapcar (lambda (molecule)
                                                           (make-molecule-view instance molecule))
                                                         (molecules instance))))


(defmethod jupyter-widgets:on-trait-change ((instance molecule-select) type (name (eql :molecules))
                                            old-value new-value source)
  (declare (ignore type name old-value new-value source))
  (render-molecules instance))

