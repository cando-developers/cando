(in-package :geom)

(defun ov3 (x y z)
  (let ((v (core:make-cxx-object 'geom:ovector3)))
    (geom:set-all3 v x y z)
    v))

(defun center-of-geometry (atoms)
  (let ((sum (ov3 0.0 0.0 0.0)))
    (dolist (a atoms)
      (setf sum (geom:add sum (chem:get-position a))))
    (geom:times-scalar sum (/ 1.0 (length atoms)))))
