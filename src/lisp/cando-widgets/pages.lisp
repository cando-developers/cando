(in-package :cando-widgets)


(defclass page (jw:grid-box)
  ((container
     :accessor container
     :initarg :container
     :initform nil))
  (:metaclass jupyter-widgets:trait-metaclass))


(defun add-page (container page title)
  (setf (jw:widget-%titles container) (append (jw:widget-%titles container) (list title))
        (jw:widget-children container) (append (jw:widget-children container) (list page)))
  (values))


(defclass task-page (page)
  ((button
     :accessor button
     :initarg :button)
   (label
     :accessor label
     :initform (make-instance 'jw:label
                              :layout (make-instance 'jw:layout
                                                     :align-self "start"
                                                     :grid-area "label")))
   (progress
     :accessor progress
     :initform (make-instance 'jw:int-progress
                              :layout (make-instance 'jw:layout
                                                     :align-self "end"
                                                     :display "none"
                                                     :grid-area "progress")))
   (messages
     :accessor messages
     :initform (make-instance 'jw:output
                              :layout (make-instance 'jw:layout
                                                     :grid-area "messages"
                                                     :padding "var(--jp-widgets-input-padding)"
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :overflow-y "scroll")))
   (task-function
     :accessor task-function
     :initarg :task-function
     :initform (lambda (parameter progress-callback)
                 (declare (ignore parameter progress-callback))
                 t)))
  (:metaclass jupyter-widgets:trait-metaclass)
  (:default-initargs
    :layout (make-instance 'jw:layout
                           :width "100%"
                           :max-height "12em"
                           :grid-gap "0.1em 1em"
                           :grid-template-rows "min-content min-content 1fr"
                           :grid-template-columns "auto 1fr"
                           :grid-template-areas "\"upload messages\" \"label messages\" \"progress messages\"")))


(defmethod initialize-instance :after ((instance task-page) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (setf (jw:widget-children instance) (list (button instance)
                                            (label instance)
                                            (progress instance)
                                            (messages instance))))


(defun run-task (instance parameter)
  (declare (ignore args))
  (with-slots (container messages progress button)
              instance
    (setf (jw:widget-outputs messages) nil
          (jw:widget-value progress) 0
          (jw:widget-disabled button) t)
    (jw:with-output messages
      (handler-bind
          ((warning
             (lambda (wrn)
               (princ wrn)
               (terpri)
               (finish-output)
               (muffle-warning))))
        (write-line "Starting task...")
        (finish-output)
        (cond
          ((funcall (task-function instance)
                    parameter
                    (lambda (maximum)
                      (setf (jw:widget-display (jw:widget-layout progress)) nil
                            (jw:widget-max progress) maximum
                            (jw:widget-value progress) (1+ (jw:widget-value progress)))))
            (write-line "Task completed successfully.")
            (setf (jw:widget-display (jw:widget-layout progress)) "none")
            (let ((index (position instance (jw:widget-children container) :test #'eql)))
              (setf (jw:widget-selected-index container)
                    (when (and index
                               (< (1+ index) (length (jw:widget-children container))))
                      (1+ index)))))
          (t
            (write-line "Task completed with failures." *error-output*)
            (finish-output *error-output*)))
        (setf (jw:widget-disabled button) nil)
        (values)))))

(defclass file-task-page (task-page)
  ()
  (:metaclass jupyter-widgets:trait-metaclass)
  (:default-initargs
    :button (make-instance 'jw:file-upload
                           :description "Select file"
                           :width "auto"
                           :layout (make-instance 'jw:layout
                                                  :align-self "start"
                                                  :grid-area "upload"))))


(defun on-data-change (instance data metadata)
  (when (= (length data) (length metadata))
    (run-task instance
              (mapcar (lambda (content metadata)
                        (acons "content" content (copy-alist metadata)))
                      data metadata))))


(defmethod initialize-instance :after ((instance file-task-page) &rest initargs &key &allow-other-keys)
  (jw:observe (button instance) :data
    (lambda (inst type name old-value new-value source)
      (declare (ignore type name old-value source))
      (on-data-change instance new-value (jw:widget-metadata inst))))
  (jw:observe (button instance) :metadata
    (lambda (inst type name old-value new-value source)
      (declare (ignore type name old-value source))
      (on-data-change instance (jw:widget-data inst) new-value))))


(defun make-file-task-page (container title task-function &key label description accept)
  (let ((page (make-instance 'file-task-page
                             :container container
                             :task-function task-function)))
    (add-page container page title)
    (when accept
      (setf (jw:widget-accept (button page)) accept))
    (setf (jw:widget-value (label page))
          (or label
              "Click to select or drag a file onto the upload button.")
          (jw:widget-description (button page))
          (or description
              "Select file"))
    (values)))


(defclass simple-task-page (task-page)
  ()
  (:metaclass jupyter-widgets:trait-metaclass)
  (:default-initargs
    :button (make-instance 'jw:button
                           :description "Start task"
                           :width "auto"
                           :layout (make-instance 'jw:layout
                                                  :align-self "start"
                                                  :grid-area "upload"))))


(defmethod initialize-instance :after ((instance simple-task-page) &rest initargs &key &allow-other-keys)
  (jw:on-button-click (button instance)
    (lambda (inst)
      (declare (ignore inst))
      (run-task instance nil))))


(defun make-simple-task-page (container title task-function &key label description)
  (let ((page (make-instance 'simple-task-page
                             :container container
                             :task-function task-function)))
    (add-page container page title)
    (when label
      (setf (jw:widget-value (label page)) label))
    (when description
      (setf (jw:widget-description (button page)) description))
    (values)))
