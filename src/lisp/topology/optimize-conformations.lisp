(in-package :foldamer)


(defun classify-dihedral-angle (angle-degrees)
  (when (< angle-degrees 0.0)
    (setf angle-degrees (+ 360.0 angle-degrees)))
  (let* ((rounded (round (/ angle-degrees 30.0)))
         (result (if (= rounded 12) 0 rounded)))
    result))


(defun classify-dihedral-angles (angles)
  "Convert the list of dihedral angles into multiples of 30.0 degrees (values 0 to 11)
and convert that into a base 12 value."
  (loop with accumulator = 0
        for angle in angles
        for classified-angle = (classify-dihedral-angle angle)
        do (setf accumulator (+ (* accumulator 12) classified-angle))
        finally (return-from classify-dihedral-angles accumulator)))


(defun optimize-fragment-conformations (conformations)
  (maphash (lambda (key value)
             (format t "~a~%" key))
           conformations))



