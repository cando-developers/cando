
(in-package :structure-editor)


(defclass composer (w:accordion)
  ((kekule
     :reader composer-kekule
     :initform (make-instance 'kekule:composer
                              :format-id kekule:+kekule-json-format+
                              :layout (make-instance 'w:layout :grid-area "kekule" :width "100%" :height "100%")))
   (data
     :accessor composer-data
     :initform nil
     :trait :json)
   (aggregate
     :accessor composer-aggregate
     :initform nil
     :trait t))
  (:metaclass jupyter-widgets:trait-metaclass)
  (:default-initargs
    :selected-index 0))


(defun run-parsing (instance progress-callback)
  (declare (ignore progress-callback))
  (handler-case
      (setf (composer-aggregate instance) (parse-kekule-json (composer-data instance)))
    (simple-error (condition) ; This is temporary until we have composer-parse-error
      (princ condition *error-output*)
      (fresh-line *error-output*)
      nil)
    (:no-error (ret)
      (declare (ignore ret))
      t)))


(defun composer (kekule-title parsing-title)
  (let* ((instance (make-instance 'composer))
         (full-screen (make-instance 'w:button ; Full screen toggle button
                                     :tooltip "Full Screen"
                                     :icon "expand"
                                     :layout (make-instance 'w:layout
                                                            :width "max-content"
                                                            :justify-self "end"
                                                            :grid-area "full-screen")))
        (kekule-page (make-instance 'resizable-box:resizable-grid-box
                                    :layout (make-instance 'resizable-box:resizable-layout
                                                           :height "400px"
                                                           :grid-template-columns "1fr"
                                                           :grid-template-rows "1fr min-content"
                                                           :grid-template-areas "'kekule' 'full-screen'"
                                                           :padding "0 0 24px 0"
                                                           :resize "vertical"))))
    (jupyter-widgets:on-button-click full-screen
      (lambda (inst)
       (declare (ignore inst))
       (resizable-box:enter-full-screen kekule-page)))
    (cw:add-page instance kekule-page kekule-title)
    (setf (jupyter-widgets:widget-children kekule-page) (list full-screen (composer-kekule instance)))
    (cw:make-simple-task-page instance parsing-title (lambda (action parameter progress-callback)
                                                       (declare (ignore action parameter))
                                                       (run-parsing instance progress-callback))
                              :label "Click button to parse Kekule output.")
    (w:link instance :data (composer-kekule instance) :data t)
    (w:observe (composer-kekule instance) :data
      (lambda (kekule-instance type name old-value new-value source)
        (declare (ignore type name source))
        (when (and new-value
                   (not old-value))
          (kekule:fit kekule-instance))))
    instance))
