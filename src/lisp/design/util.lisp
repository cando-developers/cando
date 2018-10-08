
(in-package :design)

(defun vec-to-array (vec)
  (vector (float (geom:vx vec) 1s0)
          (float (geom:vy vec) 1s0)
          (float (geom:vz vec) 1s0)))

(defun stub-as-coordinate-system (stub)
  (let* ((len 5.0)
         (o (geom:m*v stub (geom:vec 0.0 0.0 0.0)))
         (x (geom:m*v stub (geom:vec len 0.0 0.0)))
         (y (geom:m*v stub (geom:vec 0.0 len 0.0)))
         (z (geom:m*v stub (geom:vec 0.0 0.0 len))))
    (values o x y z)))

(defun generate-coordinate-system-arrows (stub)
  (multiple-value-bind (o x y z)
      (stub-as-coordinate-system stub)
    (let* ((thick 0.3)
           (shape (vector (vector "arrow" (vec-to-array o) (vec-to-array x) #(1 0 0) thick)
                          (vector "arrow" (vec-to-array o) (vec-to-array y) #(0 1 0) thick)
                          (vector "arrow" (vec-to-array o) (vec-to-array z) #(0 0 1) thick))))
      shape)))


(defun measure-distance (a b)
  (let ((va (apply #'geom:vec a))
        (vb (apply #'geom:vec b)))
    (geom:vlength (geom:v- va vb))))

(defun measure-distance (a b c)
  (let ((va (apply #'geom:vec a))
        (vb (apply #'geom:vec b))
        (vc (apply #'geom:vec c)))
    (geom:vlength (geom:v- va vb))))


(defun find-joint-using-atom (atom-tree input-atom atom-id-map)
  "Return the joint in the atom-tree that corresponds to input-atom using the atom-id-map"
  (kin:walk atom-tree
            (lambda (joint) 
              (let* ((atom (chem:lookup-atom atom-id-map (kin:atom-id joint)))
                     (pos (chem:get-position atom))
                     (delta (geom:v- pos (chem:get-position input-atom)))
                     (delta-len (geom:vlength delta)))
                (when (< delta-len 0.01)
                  (return-from find-joint-using-atom joint))))))

