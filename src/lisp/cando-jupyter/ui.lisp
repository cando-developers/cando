(in-package #:cando-jupyter)

(defclass jupyter-ui-client ()
  ())

(defclass progress-bar ()
  ((start-time :reader progress-bar-start-time
               :initform (get-universal-time))
   (update-time :accessor progress-bar-update-time
                :initform (get-universal-time))
   (bar :reader progress-bar-bar
        :initarg :bar)
   (label :reader progress-bar-label
          :initarg :label)
   (style :reader progress-bar-style
          :initarg :style)))

(defmethod cando:do-make-progress-bar
    ((client jupyter-ui-client) &rest initargs
     &key (message "Completed") (message-width 10) (total 100.0) (width 50) (style :time)
     &allow-other-keys)
  (finish-output)
  (finish-output *error-output*)
  (let ((bar (jw:make-float-progress :description message
                                     :style (jw:make-progress-style :description-width (format nil "~aex" message-width))
                                     :max total
                                     :layout (jw:make-layout :width (format nil "~aex" (+ width message-width)))))
        (label (jw:make-label)))
    (j:display (jw:make-h-box :children (list bar label)))
    (make-instance 'progress-bar :bar bar :label label :style style)))

(defmethod cando:progress-advance ((instance progress-bar) counter &optional message)
  (with-accessors ((start-time progress-bar-start-time)
                   (update-time progress-bar-update-time)
                   (bar progress-bar-bar)
                   (label progress-bar-label)
                   (style progress-bar-style))
      instance
    (when message
      (setf (jw:widget-description bar) message))
    (when (> (- (get-universal-time) update-time) 0.2)
      (setf (jw:widget-value bar) counter
            update-time (get-universal-time))
      (let* ((fraction (/ counter (jw:widget-max bar)))
             (elapsed-time (- (get-universal-time) start-time))
             (remaining-time (* elapsed-time (1- (/ fraction)))))
        (setf (jw:widget-value label)
              (format nil "~3d%, ~@[Elapsed: ~a, ~]Remaining: ~a"
                      (round (* 100 fraction))
                      (when (eq style :time)
                        (cando::progress-convenient-time elapsed-time))
                      (cando::progress-convenient-time remaining-time)))))))

(defmethod cando:progress-done ((instance progress-bar))
  (with-accessors ((bar progress-bar-bar)
                   (label progress-bar-label))
      instance
    (setf (jw:widget-value bar) (jw:widget-max bar)
          (jw:widget-value label) "")))
