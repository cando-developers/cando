(in-package :tirun-jupyter)

(defparameter *fu* nil)

(defun setup-receptor ()
  (format t "setup-receptor ~%")
  (let* ((box-layout (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
         (desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width))
         (fu (make-instance 'jupyter-widgets:file-upload :description "Upload PDB file"))
         (nglview (make-instance 'nglv:nglwidget)))
    (setf *fu* fu)
    (jupyter-widgets:observe fu :data
                             (lambda (instance type name old-value new-value source)
                               (declare (ignore instance type name old-value source))
                               (format t "Value changed~%")
                               (format t "old-value lengths -> ~a~%" (mapcar #'length old-value))
                               (format t "new-value lengths -> ~a~%" (mapcar #'length new-value))
                               (when new-value
                                 (let* ((last-pdb (car (last new-value)))
                                        (octets (make-array (length last-pdb) :element-type 'ext:byte8
                                                                              :initial-contents last-pdb)))
                                   (let ((as-text (babel:octets-to-string octets)))
                                     (dolist (id (nglv::ngl-component-ids nglview)) (nglv::remove-component nglview id))
                                     (nglv::add-structure nglview (make-instance 'nglv::TextStructure :text as-text)))))))
    (make-instance 'w:v-box
                   :layout box-layout
                   :children (list
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list fu))
                              (make-instance 'w:h-box
                                             :layout box-layout
                                             :children (list nglview))))))


#|
        
    (defparameter +options+ '("fu" "bar" "wibble" "quux" "baz" "gronk" "kilroy"))
    (defparameter +desc-width+ "12em")
    (defparameter +desc-style+ (make-instance 'w:description-style :description-width +desc-width+))
    (defparameter +box-layout+ (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
    (defparameter +slider-style+ (make-instance 'w:slider-style :description-width +desc-width+))
    (defparameter +widget-layout+ (make-instance 'w:layout :width "28em"))
    
    (make-instance 'w:accordion
                   :children (list (make-instance 'w:h-box
                                                  :layout +box-layout+
                                                  :children (list (make-instance 'w:color-picker :description "color-picker" :style +desc-style+)
                                                                  (make-instance 'w:date-picker :description "date-picker" :style +desc-style+)))
                                   (make-instance 'w:h-box
                                                  :layout +box-layout+
                                                  :children (list (make-instance 'w:dropdown :description "dropdown" :style +desc-style+ :layout +widget-layout+ :%options-labels +options+)
                                                                  (make-instance 'w:radio-buttons :description "radio-button" :style +desc-style+ :layout +widget-layout+ :%options-labels +options+)
                                                                  (make-instance 'w:select :description "select" :style +desc-style+ :layout +widget-layout+  :%options-labels +options+)
                                                                  (make-instance 'w:selection-slider :description "selection-slider" :style +slider-style+ :layout +widget-layout+ :%options-labels +options+)
                                                                  (make-instance 'w:selection-range-slider :description "selection-range-slider" :style +slider-style+ :layout +widget-layout+  :%options-labels +options+)
                                                                  (make-instance 'w:select-multiple :description "select-multiple" :style +desc-style+ :layout +widget-layout+  :%options-labels +options+)))
                                   (make-instance 'w:h-box
                                                  :layout +box-layout+
                                                  :children (list (make-instance 'w:button :description "button" :style +desc-style+)
                                                                  (make-instance 'w:toggle-button :description "toggle-button" :style +desc-style+))))
                   :%titles '("Picker Widgets" "Selection Widgets" "Buttons")))


|#
