(in-package :geom)

(defun v3 (x y z)
  (geom:make-v3 x y z))

(defun center-of-geometry (atoms)
  (let ((sum (v3 0.0 0.0 0.0)))
    (dolist (a atoms)
      (setf sum (geom:add sum (chem:get-position a))))
    (geom:times-scalar sum (/ 1.0 (length atoms)))))
