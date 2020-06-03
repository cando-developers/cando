(in-package :tirun-jupyter)


(defun setup-receptor ()
  (let* ((box-layout (make-instance 'w:layout :width "auto" :flex-flow "row wrap"))
         (desc-width "12em")
         (desc-style (make-instance 'w:description-style :description-width desc-width)))
    (make-instance 'w:h-box
                   :layout box-layout
                   :children (list (make-instance 'w:file-upload :description "Upload PDB File" :style desc-style)))))


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
