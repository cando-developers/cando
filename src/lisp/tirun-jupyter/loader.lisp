(in-package :tirun-jupyter)


(defclass loader ()
  ((upload-button
     :reader loader-upload-button
     :initform (make-instance 'jupyter-widgets:file-upload
                              :description "Select file"
                              :width "auto"
                              :layout (make-instance 'w:layout
                                                     :align-self "start"
                                                     :grid-area "upload")))
   (upload-label
     :reader loader-upload-label
     :initform (make-instance 'w:label
                              :value "Click to select or drag a file onto the upload button."
                              :layout (make-instance 'w:layout
                                                     :align-self "start"
                                                     :grid-area "label")))

   (upload-progress
     :reader loader-upload-progress
     :initform (make-instance 'w:int-progress
                              :layout (make-instance 'w:layout
                                                     :align-self "end"
                                                     :grid-area "progress")))
   (upload-log
     :reader loader-upload-log
     :initform (make-instance 'w:output
                              :layout (make-instance 'w:layout
                                                     :grid-area "log"
                                                     :padding "var(--jp-widgets-container-padding)"
                                                     :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                     :overflow-y "scroll")))
   (upload-grid
     :reader loader-upload-grid
     :initform (make-instance 'w:grid-box
                              :layout (make-instance 'w:layout
                                                     :width "100%"
                                                     :max-height "12em"
                                                     :grid-gap "0.1em 1em"
                                                     :grid-template-rows "min-content min-content 1fr"
                                                     :grid-template-columns "auto 1fr"
                                                     :grid-template-areas "\"upload log\" \"label log\" \"progress log\"")))
   (view-grid
     :reader loader-view-grid
     :initform (make-instance 'w:grid-box))
   (accordion
     :reader loader-accordion
     :initform (make-instance 'w:accordion
                              :%titles (list "Upload" "View")
                              :selected-index 0))))


(defgeneric loader-parse (instance data))


(defgeneric loader-show (instance data))


(defmethod loader-show :before ((instance loader) data)
  (declare (ignore data))
  (setf (w:widget-selected-index (loader-accordion instance)) 1))


(defmethod initialize-instance :after ((instance loader) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))
  (setf (w:widget-children (loader-upload-grid instance)) (list (loader-upload-button instance)
                                                                (loader-upload-label instance)
                                                                (loader-upload-progress instance)
                                                                (loader-upload-log instance))
        (w:widget-children (loader-accordion instance)) (list (loader-upload-grid instance)
                                                              (loader-view-grid instance)))
  (w:observe (loader-upload-button instance) :data
    (lambda (inst type name old-value new-value source)
      (declare (ignore inst type name old-value source))
      (setf (w:widget-selected-index (loader-accordion instance)) 0)
      (setf (w:widget-outputs (loader-upload-log instance)) nil)
      (w:with-output (loader-upload-log instance)
        (when new-value
          (write-string "Parsing file...")
          (finish-output)
          (let ((data (loader-parse instance (car (last new-value)))))
            (cond
              (data
                (write-string "Parsing complete.")
                (loader-show instance data))
              (t
                (write-string "Parsing failed." *error-output*)
                (finish-output *error-output*)))))))))



