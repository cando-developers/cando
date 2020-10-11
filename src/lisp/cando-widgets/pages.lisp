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
  ((progress
     :accessor progress
     :initform (make-instance 'jw:int-progress
                              :layout (make-instance 'jw:layout
                                                     :align-self "end"
                                                     :display "none"
                                                     :grid-area "progress")))
   (messages-label
     :accessor messages-label
     :initform (make-instance 'jw:label
                              :value "Messages"
                              :layout (make-instance 'jw:layout
                                                     :text-align "center"
                                                     :grid-area "messages-label")))
   (messages
     :accessor messages
     :initform (make-instance 'jw:output
                              :layout (make-instance 'jw:layout
                                                     :grid-area "messages"
                                                     :padding "var(--jp-widgets-input-padding)"
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :min-height "4em"
                                                     :overflow-y "scroll")))
   (task-function
     :accessor task-function
     :initarg :task-function
     :initform (lambda (action parameter progress-callback)
                 (declare (ignore action parameter progress-callback))
                 t)))
  (:metaclass jupyter-widgets:trait-metaclass))


(defmethod initialize-instance :after ((instance task-page) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (setf (jw:widget-children instance)
        (append (jw:widget-children instance)
                (list (progress instance)
                      (messages-label instance)
                      (messages instance)))))


(defclass single-task-page (task-page)
  ((button
     :accessor button
     :initarg :button)
   (label
     :accessor label
     :initform (make-instance 'jw:label
                              :layout (make-instance 'jw:layout
                                                     :align-self "start"
                                                     :grid-area "label"))))
  (:metaclass jupyter-widgets:trait-metaclass)
  (:default-initargs
    :layout (make-instance 'jw:layout
                           :width "100%"
                           :max-height "12em"
                           :grid-gap "0.1em 1em"
                           :grid-template-rows "min-content min-content 1fr"
                           :grid-template-columns "auto 1fr"
                           :grid-template-areas "\"button messages-label\" \"label messages\" \"progress messages\"")))


(defmethod initialize-instance :after ((instance single-task-page) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (setf (jw:widget-children instance)
        (append (jw:widget-children instance)
                (list (button instance)
                      (label instance)))))


(defgeneric begin-task (instance)
  (:method (instance))
  (:method ((instance single-task-page))
    (setf (jw:widget-disabled (button instance)) t)))


(defgeneric end-task (instance)
  (:method (instance))
  (:method ((instance single-task-page))
    (setf (jw:widget-disabled (button instance)) nil)))


(defun run-task (instance action parameter)
  (declare (ignore args))
  (with-slots (container messages progress)
              instance
    (unwind-protect
        (progn
          (begin-task instance)
          (setf (jw:widget-outputs messages) nil
                (jw:widget-value progress) 0)
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
                                        action
                                        parameter
                                        (lambda (maximum)
                                          (setf (jw:widget-display (jw:widget-layout progress)) nil
                                                (jw:widget-max progress) maximum
                                                (jw:widget-value progress) (1+ (jw:widget-value progress)))))
                                (write-line "Task completed successfully.")
                                (setf (jw:widget-display (jw:widget-layout progress)) "none")
                                (let ((index (position instance (jw:widget-children container) :test #'eql)))
                                  (when (and index
                                             (< (1+ index) (length (jw:widget-children container))))
                                    (setf (jw:widget-selected-index container) (1+ index)))))
                              (t
                                (write-line "Task completed with failures." *error-output*)
                                (finish-output *error-output*)))
                            (values))))
      (end-task instance))))


(defclass file-task-page (single-task-page)
  ()
  (:metaclass jupyter-widgets:trait-metaclass)
  (:default-initargs
    :button (make-instance 'jw:file-upload
                           :description "Select file"
                           :width "auto"
                           :layout (make-instance 'jw:layout
                                                  :align-self "start"
                                                  :grid-area "button"))))


(defun on-data-change (instance data metadata)
  (when (= (length data) (length metadata))
    (run-task instance nil
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


(defclass simple-task-page (single-task-page)
  ()
  (:metaclass jupyter-widgets:trait-metaclass)
  (:default-initargs
    :button (make-instance 'jw:button
                           :description "Start task"
                           :width "auto"
                           :layout (make-instance 'jw:layout
                                                  :align-self "start"
                                                  :grid-area "button"))))


(defmethod initialize-instance :after ((instance simple-task-page) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (jw:on-button-click (button instance)
    (lambda (inst)
      (declare (ignore inst))
      (run-task instance nil nil))))


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

