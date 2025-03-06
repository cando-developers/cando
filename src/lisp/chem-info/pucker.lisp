(in-package :chem-info)



(defun cremer-pople-ring5-positions (vector-list)
  (let ((center-of-geometry (loop with acc = (geom:vec 0.0 0.0 0.0)
                                  for vec in vector-list
                                  for acc = (geom:v+ vec acc)
                                  finally (return (geom:v* acc (/ 1.0 (length vector-list)))))))
    center-of-geometry))


(defun cremer-pople-ring5-atoms (atom-list)
  (let ((vector-list (mapcar #'chem:get-position atom-list)))
    (cremer-pople-ring5-positions vector-list)))
