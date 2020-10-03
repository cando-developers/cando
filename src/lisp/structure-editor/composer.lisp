
(in-package :structure-editor)


(defclass composer (w:accordion)
  ((kekule
     :reader composer-kekule
     :initform (make-instance 'kekule:composer
                              :format-id kekule:+kekule-json-format+
                              :layout (make-instance 'w:layout :width "95%")))
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
  (let ((instance (make-instance 'composer)))
    (cw:add-page instance (composer-kekule instance) kekule-title)
    (cw:make-simple-task-page instance parsing-title (lambda (parameter progress-callback)
                                                       (declare (ignore parameter))
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
