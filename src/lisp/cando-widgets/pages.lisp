(in-package :cando-widgets)


(defclass page (resizable-box:resizable-grid-box)
  ((container
     :accessor container
     :initarg :container
     :initform nil))
  (:metaclass jupyter-widgets:trait-metaclass))


(defun add-page (container page title)
  (if (typep container 'jw:accordion)
    (setf (jw:widget-%titles container) (append (jw:widget-%titles container)
                                                (list title))
          (jw:widget-children container) (append (jw:widget-children container)
                                                 (list page)))
    (setf (jw:widget-children container) (append (jw:widget-children container)
                                                 (list (jw:make-html :value (format nil "<h1>~a</h1>" title)) page))))
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
                                                     :max-height "12em"
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


(defun do-run-task (instance action parameter)
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
                                        instance
                                        action
                                        parameter
                                        (lambda (&key increment value maximum)
                                          (setf (jw:widget-display (jw:widget-layout progress)) nil)
                                          (when maximum
                                            (setf (jw:widget-max progress) maximum))
                                          (when value
                                            (setf (jw:widget-value progress) value))
                                          (when increment
                                            (incf (jw:widget-value progress) increment))))
                                (write-line "Task completed successfully.")
                                (setf (jw:widget-display (jw:widget-layout progress)) "none")
                                (when (typep container 'jw:accordion)
                                  (let ((index (position instance (jw:widget-children container) :test #'eql)))
                                    (when (and index
                                               (< (1+ index) (length (jw:widget-children container))))
                                      (setf (jw:widget-selected-index container) (1+ index))))))
                              (t
                                (write-line "Task completed with failures." *error-output*)
                                (finish-output *error-output*)))
                            (values))))
      (end-task instance))))


(defgeneric run-task (instance &optional action parameter)
  (:method (instance &optional action parameter)
    (do-run-task instance action parameter)))


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


(defmethod initialize-instance :after ((instance file-task-page) &rest initargs &key &allow-other-keys)
  (jw:observe (button instance) :value
    (lambda (inst type name old-value new-value source)
      (declare (ignore type name old-value source))
      (unless (zerop (length new-value))
        (run-task instance nil new-value)))))


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
    page))


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
    page))


(defclass threaded-task-page (single-task-page)
  ((thread
     :accessor thread
     :initform nil)
   (stop-button
     :accessor stop-button
     :initform (jw:make-button :description "Stop task"
                               :width "auto"
                               :disabled t
                               :layout (make-instance 'jw:layout
                                                      :align-self "start"
                                                      :grid-area "stop-button")))
   (parameter
     :accessor parameter
     :initarg :parameter
     :initform nil))
  (:metaclass jupyter-widgets:trait-metaclass)
  (:default-initargs
    :layout (make-instance 'jw:layout
                           :width "100%"
                           :grid-gap "0.1em 1em"
                           :grid-template-rows "min-content min-content 1fr"
                           :grid-template-columns "auto auto 1fr"
                           :grid-template-areas "\"button stop-button messages-label\" \"label label messages\" \"progress progress messages\"")
    :button (make-instance 'jw:button
                           :description "Start task"
                           :width "auto"
                           :layout (make-instance 'jw:layout
                                                  :align-self "start"
                                                  :grid-area "button"))))


(defmethod run-task ((instance threaded-task-page) &optional action parameter)
  (setf (thread instance)
        (let* ((stdout (jw:make-output-widget-stream (messages instance)))
               (stderr (jw:make-output-widget-stream (messages instance) t))
               (bordeaux-threads:*default-special-bindings* `((jupyter::*kernel* . ,jupyter::*kernel*)
                                                              (jupyter::*message* . nil)
                                                              (*standard-output* . ,stdout)
                                                              (*debug-io* . ,stdout)
                                                              (*trace-output* . ,stdout)
                                                              (*error-output* . ,stderr)
                                                              ,@bordeaux-threads:*default-special-bindings*)))
          (bordeaux-threads:make-thread
            (lambda ()
              (j:handling-comm-errors
                (do-run-task instance action (or parameter (parameter instance)))))))))


(defmethod initialize-instance :after ((instance threaded-task-page) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (setf (jw:widget-children instance)
        (append (jw:widget-children instance)
                (list (stop-button instance))))
  (jw:on-button-click (stop-button instance)
    (lambda (inst)
      (declare (ignore inst))
      (bordeaux-threads:destroy-thread (thread instance))))
  (jw:on-button-click (button instance)
    (lambda (inst)
      (declare (ignore inst))
      (run-task instance nil nil))))


(defun make-threaded-task-page (container title task-function &key label description parameter)
  (let ((page (make-instance 'threaded-task-page
                             :parameter parameter
                             :container container
                             :task-function task-function)))
    (add-page container page title)
    (when label
      (setf (jw:widget-value (label page)) label))
    (when description
      (setf (jw:widget-description (button page)) description))
    page))


(defmethod begin-task :before ((instance threaded-task-page))
  (setf (jw:widget-disabled (stop-button instance)) nil
        (jw:widget-msg-id (messages instance)) nil))


(defmethod end-task :before ((instance threaded-task-page))
  (setf (jw:widget-disabled (stop-button instance)) t))

