
(in-package :structure-editor)


(defclass composer (jw:accordion)
  ((kekule
     :reader composer-kekule
     :initform (make-instance 'kekule:composer
                              :format-id kekule:+kekule-json-format+
                              :layout (make-instance 'jw:layout :grid-area "kekule" :width "100%" :height "100%")))
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
        (kekule-page (make-instance 'resizable-box:resizable-grid-box
                                    :enable-full-screen t
                                    :layout (make-instance 'resizable-box:resizable-layout
                                                           :height "400px"
                                                           :grid-template-columns "1fr"
                                                           :grid-template-rows "1fr min-content"
                                                           :grid-template-areas "'kekule' 'full-screen'"
                                                           :padding "0 16px 24px 0"
                                                           :resize "vertical"))))
    (cw:add-page instance kekule-page kekule-title)
    (setf (jupyter-widgets:widget-children kekule-page) (list (composer-kekule instance)))
    (cw:make-simple-task-page instance parsing-title (lambda (action parameter progress-callback)
                                                       (declare (ignore action parameter))
                                                       (run-parsing instance progress-callback))
                              :label "Click button to parse Kekule output.")
    (jw:link instance :data (composer-kekule instance) :data t)
    (jw:observe (composer-kekule instance) :data
      (lambda (kekule-instance type name old-value new-value source)
        (declare (ignore type name source))
        (when (and new-value
                   (not old-value))
          (kekule:fit kekule-instance))))
    instance))
