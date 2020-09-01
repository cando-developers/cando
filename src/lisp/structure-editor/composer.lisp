
(in-package :structure-editor)

(defclass composer ()
  ((app
    :accessor composer-app
    :initarg :app)
   (kekule
    :reader composer-kekule
    :initform (make-instance 'kekule:composer :format-id kekule:+kekule-json-format+
                                              :width "auto"
                                              :layout (make-instance 'w:layout
                                                                     :width "auto"
                                                                     :align-self "start"
                                                                     :grid-area "composer")))
   (save-button
    :reader composer-save-button
    :initform (make-instance 'jupyter-widgets:button
                             :description "Save"
                             :width "auto"
                             :layout (make-instance 'w:layout
                                                    :align-self "end"
                                                    :grid-area "save")
                             ))
   (composer-log
    :reader composer-log
    :initform (make-instance 'w:output
                             :layout (make-instance 'w:layout
                                                    :grid-area "log"
                                                    :padding "var(--jp-widgets-input-padding)"
                                                    :border "var(--jp-widgets-border-width) solid var(--jp-border-color1)"
                                                    :overflow-y "scroll")))
   (composer-grid
    :reader composer-grid
    :initform (make-instance 'w:grid-box
                             :layout (make-instance 'w:layout
                                                    :width "auto"
                                                    :max-height "auto"
                                                    :grid-gap "0.1em"
                                                    :grid-template-rows "1fr min-content"
                                                    :grid-template-columns "1fr 3fr"
                                                    :grid-template-areas "\"composer composer\" \"save log\"")))
   (view-grid
    :reader composer-view-grid
    :initform (make-instance 'w:html
                             :layout (make-instance 'w:layout
                                                    :width "100%"
                                                    :height "auto"
                                                    )))
   (accordion
    :reader composer-accordion
    :initform (make-instance 'w:accordion
                             :%titles (list "Composer" "View")
                             :selected-index 0))))



(defmethod initialize-instance :after ((instance composer) &rest initargs &key &allow-other-keys)
  (declare (ignore initargs))  (setf (w:widget-children (composer-accordion instance))
                                     (list (composer-grid instance) (composer-view-grid instance)))
  (setf (w:widget-children (composer-grid instance))
        (list (composer-kekule instance) (composer-save-button instance) (composer-log instance)))
  (setf (kekule:data (composer-kekule instance))
        (get-composer-json (composer-app instance)))
  (setf (w::widget-on-click (composer-save-button instance))
        (list 
         (lambda (&rest args)
           (setf (w:widget-outputs (composer-log instance)) nil)
           (w:with-output (composer-log instance)
             (let ((json (kekule:data (composer-kekule instance))))
               (set-composer-json (composer-app instance) json)
               #+(or)(format t "new observe save-button pressed~%~s~%" json)
               (finish-output)
               ))))))
  




