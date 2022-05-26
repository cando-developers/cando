
(in-package :cando-user)

(defparameter *energy-function* nil)
(defparameter *pick-aggregate* nil)
(defparameter *ff-output* nil)

(defun identify-atom (atom-ht aggregate)
  (let ((xpos (gethash "x" atom-ht))
        (ypos (gethash "y" atom-ht))
        (zpos (gethash "z" atom-ht)))
    (cando:do-atoms (atm aggregate)
      (let ((deltax (abs (- xpos (geom:get-x (chem:get-position atm)))))
            (deltay (abs (- ypos (geom:get-y (chem:get-position atm)))))
            (deltaz (abs (- zpos (geom:get-z (chem:get-position atm))))))
        (when (and (< deltax 0.01)
                   (< deltay 0.01)
                   (< deltaz 0.01))
          (return-from identify-atom atm))))
    nil))

(defun log-interactions (picked-atoms energy-function)
  (cond
    ((and (= (length picked-atoms) 2)
          (first picked-atoms)
          (second picked-atoms))
     (destructuring-bind (a1 a2)
         picked-atoms
       (let* ((term (chem:get-stretch-component energy-function))
              (atom-table (chem:atom-table energy-function))
              (terms (chem:lookup-stretch-terms term atom-table a1 a2)))
         (format t "Stretch terms: ~%")
         (format t "~4a ~4a~%" (chem:get-name a1) (chem:get-name a2))
         (format t "~4a ~4a ~8@a ~8@a~%"
                 "T1" "T2" "kb" "r0")
         (loop for term in terms
               do (destructuring-bind (&key type1 type2 kb r0)
                      term
                    (format t "~4a ~4a ~8,3f ~8,3f~%"
                            type1 type2 kb r0))))))
    ((and (= (length picked-atoms) 3)
          (first picked-atoms)
          (second picked-atoms)
          (third picked-atoms))
     (destructuring-bind (a1 a2 a3)
         picked-atoms
       (let* ((term (chem:get-angle-component energy-function))
              (atom-table (chem:atom-table energy-function))
              (terms (chem:lookup-angle-terms term atom-table a1 a2 a3)))
         (format t "Angle terms: ~%")
         (format t "~4a ~4a ~4a~%" (chem:get-name a1) (chem:get-name a2) (chem:get-name a3))
         (format t "~4a ~4a ~4a ~8@a ~8@a~%"
                 "T1" "T2" "T3" "kt" "t0")
         (loop for term in terms
               do (destructuring-bind (&key type1 type2 type3 kt t0)
                      term
                    (format t "~4a ~4a ~4a ~8,3f ~8,3f~%"
                            type1 type2 type3 kt t0))))))
    ((and (= (length picked-atoms) 4)
          (first picked-atoms)
          (second picked-atoms)
          (third picked-atoms)
          (fourth picked-atoms))
     (destructuring-bind (a1 a2 a3 a4)
         picked-atoms
       (let* ((term (chem:get-dihedral-component energy-function))
              (atom-table (chem:atom-table energy-function))
              (dihedral-terms (chem:lookup-dihedral-terms term atom-table a1 a2 a3 a4)))
         (format t "Dihedral terms: ~%")
         (format t "~4a ~4a ~4a ~4a~%" (chem:get-name a1) (chem:get-name a2) (chem:get-name a3) (chem:get-name a4))
         (format t "~4a ~4a ~4a ~4a ~8@a ~2@a ~5@a ~8@a~%"
                 "T1" "T2" "T3" "T4" "Type" "M" "Phase" "V")
         (loop for term in dihedral-terms
               do (destructuring-bind (&key type1 type2 type3 type4 proper multiplicity phase-degrees v)
                      term
                    (format t "~4a ~4a ~4a ~4a ~[  proper~;improper~] ~2d ~5f ~8,3f~%"
                            type1 type2 type3 type4 (if proper 0 1) multiplicity phase-degrees v))))))))

(defun stage-select (out inst data)
  (declare (ignore inst))
  (let ((atoms (loop for val across data
                     for atm = (identify-atom val *pick-aggregate*)
                     collect atm)))
    (jw:with-output out
      (log-interactions atoms *energy-function*))))

(defun force-field-viewer (aggregate)
  "Interactively display bond, angle and dihedral terms of an energy function 
generated from an aggregate"
  (let* ((energy-function (chem:make-energy-function :matter aggregate :assign-types t))
         (component (cando-widgets:make-ngl-structure
                     aggregate
                     :auto-view-duration 0
                     :representations (list
                                       (make-instance 'ngl:ball-and-stick :name "Ball and Stick" :visible t))))
         (stage (apply #'make-instance 'ngl:stage
                       :clip-dist 0
                       :background-color "white"
                       :components (list component)
                       nil))
         (label (make-instance 'jw:label
                               :value "Force field terms"))
         (outs (make-instance 'jw:output
                              :layout (jw:make-layout :height "10em"
                                                      :overflow-y "scroll")))
         (all (make-instance 'jw:v-box :children (list stage label outs)))
         )
    (setf *energy-function* energy-function)
    (setf *pick-aggregate* aggregate)
    (ngl:on-stage-select stage
                         (lambda (inst data)
                           (stage-select outs inst data)))
    all))
    
