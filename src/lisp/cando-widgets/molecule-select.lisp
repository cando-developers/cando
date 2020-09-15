(in-package :cando-widgets)


(defparameter *molecule-sketches* (make-hash-table :test #'eql))


(defun sketch-molecule (molecule)
  (multiple-value-bind (sketch present-p)
                       (gethash molecule *molecule-sketches*)
    (if present-p
      (lparallel:force sketch)
      (setf (gethash molecule *molecule-sketches*)
            (sketch2d:svg (sketch2d:sketch2d molecule))))))


(defparameter *common-sketch* nil)


(defun build-prototype-sketch (ligands)
  "Find the smallest ligand and build a sketch for it"
  (let ((min-atoms (chem:number-of-atoms (car ligands)))
        (min-mol (car ligands)))
    (loop for mol in (cdr ligands)
          for num-atoms = (chem:number-of-atoms mol)
          when (< num-atoms min-atoms)
            do (setf min-atoms num-atoms
                     min-mol mol))
    (sketch2d:sketch2d min-mol)))


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


(defclass molecule-select (jw:grid-box)
  ((container
     :accessor container
     :initform (make-instance 'jw:grid-box
                              :layout (make-instance 'jw:layout
                                                     :height "640px"
                                                     :overflow-y "scroll"
                                                     :grid-auto-rows "min-content"
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :grid-gap ".25em"
                                                     :padding ".25em"
                                                     :grid-area "list"
                                                     :grid-template-columns "repeat(auto-fill, 200px)")))
   (select-all
     :accessor select-all
     :initform (make-instance 'jw:button
                              :description "Select All"
                              :layout (make-instance 'jw:layout
                                                     :grid-area "all")))
   (select-none
     :accessor select-none
     :initform (make-instance 'jw:button
                              :description "Select None"
                              :layout (make-instance 'jw:layout
                                                     :grid-area "none")))
   (selection-mode
     :reader selection-mode
     :initarg :selection-mode
     :initform :multiple)
   (molecules
     :accessor molecules
     :initarg :molecules
     :initform nil
     :trait t)
   (selected
     :accessor selected
     :initarg :selected
     :initform nil
     :trait t)
   (max-width
     :accessor max-width
     :initform +min-sketch-width+
     :trait t)
   (task-channel
     :reader task-channel
     :initform (lparallel:make-channel)))
  (:metaclass jupyter-widgets:trait-metaclass)
  (:default-initargs
    :layout (make-instance 'jw:layout
                           ;:border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                           :grid-gap "var(--jp-widgets-container-padding)"
                           ;:padding "var(--jp-widgets-container-padding)"
                           :grid-template-columns "1fr min-content min-content 1fr"
                           :grid-template-rows "min-content 1fr"
                           :grid-template-areas "\"list list list list\" \". all none .\"")))


(defmethod initialize-instance :after ((instance molecule-select) &rest args &key &allow-other-keys)
  (declare (ignore args))
  (with-slots (select-all select-none container selection-mode)
              instance
    (jw:observe instance :max-width
      (lambda (inst type name old-value new-value source)
        (declare (ignore inst type name old-value source))
        (setf (jw:widget-grid-template-columns (jw:widget-layout (container instance)))
              (format nil "repeat(auto-fill,~Fpx)" (+ 2 new-value)))))
    (setf (jw:widget-children instance)
          (if (eql :multiple selection-mode)
            (list select-all select-none container)
            (list container)))
    (jw:on-button-click select-all
      (lambda (inst)
        (declare (ignore inst))
        (setf (selected instance) (molecules instance))))
    (jw:on-button-click select-none
      (lambda (inst)
        (declare (ignore inst))
        (setf (selected instance) nil)))))


(defun sketch-and-update (instance html molecule)
  (ignore-errors
    (let ((sketch (sketch-molecule molecule)))
      (setf (max-width instance) (max (max-width instance) (sketch2d:width sketch))
            (jw:widget-value html) (sketch2d:render-svg-to-string sketch)))))


(defun make-molecule-view (instance molecule)
  (with-slots (selection-mode)
              instance
    (let ((toggle (if (eql :none selection-mode)
                    (make-instance 'jw:label
                                   :layout (make-instance 'jw:layout
                                                          :margin "var(--jp-widgets-margin) auto")
                                   :value (symbol-name (chem:get-name molecule)))
                    (make-instance 'jw:checkbox
                                   :layout (make-instance 'jw:layout
                                                          :margin "var(--jp-widgets-margin) auto")
                                   :description (symbol-name (chem:get-name molecule)))))
          (html (make-instance 'jw:html
                               :layout (make-instance 'jw:layout
                                                      :margin "auto"
                                                      :overflow "hidden"))))
      (lparallel:submit-task (task-channel instance) #'sketch-and-update instance html molecule)
      (unless (eql :none selection-mode)
        (jw:observe toggle :value
          (lambda (inst type name old-value new-value source)
            (declare (ignore inst type name old-value source))
            (let ((current-position (position molecule (selected instance) :test #'eql)))
              (cond
                ((and (not new-value)
                      current-position)
                  (setf (selected instance)
                        (remove molecule (selected instance) :test #'eql)))
                ((not new-value))
                ((eql :single selection-mode)
                  (setf (selected instance)
                        (list molecule)))
                ((not current-position)
                  (setf (selected instance)
                        (nconc (selected instance) (list molecule))))))))
        (jw:observe instance :selected
          (lambda (inst type name old-value new-value source)
            (declare (ignore type name old-value source))
            (setf (jw:widget-value toggle)
                  (and (position molecule (selected instance) :test #'eql)
                       t)))))
      (make-instance 'jw:v-box
                     :layout (make-instance 'jw:layout
                                            :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)")
                     :children (list toggle html)))))


(defun render-molecules (instance)
  (setf (max-width instance) +min-sketch-width+
        (jw:widget-children (container instance)) (mapcar (lambda (molecule)
                                                           (make-molecule-view instance molecule))
                                                         (molecules instance))))


(defmethod jupyter-widgets:on-trait-change ((instance molecule-select) type (name (eql :molecules))
                                            old-value new-value source)
  (declare (ignore type name old-value new-value source))
  (render-molecules instance))

