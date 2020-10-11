(in-package :cando-widgets)


(defclass workspace-task-page (task-page)
  ((workspace-name
     :accessor workspace-name
     :initform (make-instance 'jw:text
                              :description "Name"
                              :value "workspace-1"
                              :style (make-instance 'jw:description-style :description-width "8em")
                              :layout (make-instance 'jw:layout :width "100%" :grid-area "name")))
   (workspace-list
     :accessor workspace-list
     :initform (make-instance 'jw:select
                              :description "Workspaces"
                              :style (make-instance 'jw:description-style :description-width "8em")
                              :layout (make-instance 'jw:layout :width "100%" :grid-area "workspaces")))
   (new-button
     :accessor new-button
     :initform (make-instance 'jw:button
                              :description "New"
                              :icon "file"
                              :layout (make-instance 'jw:layout :grid-area "new")))
   (open-button
     :accessor open-button
     :initform (make-instance 'jw:button
                              :description "Open"
                              :icon "file-import"
                              :layout (make-instance 'jw:layout :grid-area "open")))
   (duplicate-button
     :accessor duplicate-button
     :initform (make-instance 'jw:button
                              :description "Duplicate"
                              :icon "copy"
                              :layout (make-instance 'jw:layout :grid-area "duplicate")))
   (save-button
     :accessor save-button
     :initform (make-instance 'jw:button
                              :description "Save"
                              :icon "file-export"
                              :layout (make-instance 'jw:layout :grid-area "save")))
   (file-type
     :accessor file-type
     :initarg :file-type
     :initform "dat"))
  (:metaclass jw:trait-metaclass)
  (:default-initargs
    :layout (make-instance 'jw:layout
                           :width "100%"
                           :max-height "12em"
                           :grid-gap "0.1em 1em"
                           :grid-template-rows "min-content min-content min-content min-content 1fr"
                           :grid-template-columns "auto auto auto auto 1fr"
                           :grid-template-areas "'new        open       save       duplicate  messages-label'
                                                 'name       name       name       name       messages'
                                                 'workspaces workspaces workspaces workspaces messages'
                                                 'progress   progress   progress   progress   messages'")))


(defun workspace-pathname (instance &optional from-list)
  (make-pathname :name (jw:widget-value (if from-list
                                          (workspace-list instance)
                                          (workspace-name instance)))
                 :type (file-type instance)))


(defun update-workspace-list (instance)
  (setf (jw:widget-%options-labels (workspace-list instance))
        (mapcar #'pathname-name (directory (make-pathname :name :wild :type (file-type instance))))))


(defun update-workspace-ui (instance)
  (update-workspace-list instance)
  (let ((no-workspace-selected (null (jw:widget-value (workspace-list instance))))
        (no-workspace-named (zerop (length (jw:widget-value (workspace-name instance))))))
    (setf (jw:widget-disabled (open-button instance)) no-workspace-selected
          (jw:widget-disabled (new-button instance)) no-workspace-named
          (jw:widget-disabled (save-button instance)) no-workspace-named
          (jw:widget-disabled (duplicate-button instance)) (or no-workspace-selected no-workspace-named))))


(defmethod initialize-instance :after ((instance workspace-task-page) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (setf (jw:widget-children instance)
        (append (jw:widget-children instance)
                (list (workspace-name instance)
                      (workspace-list instance)
                      (open-button instance)
                      (save-button instance)
                      (duplicate-button instance)
                      (new-button instance))))
  (jw:observe (workspace-list instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value new-value source))
      (update-workspace-ui instance)))
  (jw:observe (workspace-name instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value new-value source))
      (update-workspace-ui instance)))
  (jw:on-button-click (duplicate-button instance)
    (lambda (inst)
      (declare (ignore inst))
      (run-task instance :open (workspace-pathname instance t)
      (run-task instance :save (workspace-pathname instance)))))
  (jw:on-button-click (open-button instance)
    (lambda (inst)
      (declare (ignore inst))
      (run-task instance :open (workspace-pathname instance t))
      (setf (jw:widget-value (workspace-name instance))
            (jw:widget-value (workspace-list instance)))))
  (jw:on-button-click (save-button instance)
    (lambda (inst)
      (declare (ignore inst))
      (run-task instance :save (workspace-pathname instance))))
  (jw:on-button-click (new-button instance)
    (lambda (inst)
      (declare (ignore inst))
      (run-task instance :new (workspace-pathname instance)))))


(defmethod begin-task ((instance workspace-task-page))
  (setf (jw:widget-disabled (open-button instance)) t
        (jw:widget-disabled (save-button instance)) t
        (jw:widget-disabled (new-button instance)) t))


(defmethod end-task ((instance workspace-task-page))
  (update-workspace-ui instance))


(defun make-workspace-task-page (container title task-function &key file-type)
  (let ((page (make-instance 'workspace-task-page
                             :container container
                             :task-function task-function
                             :file-type (or file-type "dat"))))
    (add-page container page title)
    (update-workspace-list page)
    page))

