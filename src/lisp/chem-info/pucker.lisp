(in-package :chem-info)

(defun cremer-pople-ring5-positions (vector-list)
  (geom:v* (reduce #'geom:v+ vector-list)
           (/ 1.0 (length vector-list))))

(defun cremer-pople-ring5-atoms (atom-list)
  (let ((vector-list (mapcar #'chem:get-position atom-list)))
    (cremer-pople-ring5-positions vector-list)))
