
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



(defun topology-constitution-order (topology)
  (let* ((constitution (chem:get-constitution topology))
         (constitution-atoms (chem:get-constitution-atoms constitution))
         (atoms-list (chem:constitution-atoms-as-list constitution-atoms))
         (names (make-array (length atoms-list) :initial-contents (mapcar (lambda (atom) (chem:atom-name atom)) atoms-list))))
    names))

(defun walk-joint-template (node func)
  (funcall func node)
  (loop for child in (kin:children node)
        do (walk-joint-template child func)))

(defun topology-joint-template-order (topology)
  (let* ((joint-template (chem:get-property topology :joint-template))
         (max-index 0))
    (walk-joint-template joint-template
                         (lambda (node)
                           (setf max-index (max max-index (kin:id node)))))
    (let ((vec (make-array (1+ max-index))))
      (walk-joint-template joint-template
                           (lambda (node)
                             (setf (elt vec (kin:id node)) (kin:name node))))
      vec)))

(defun coordinate-order (coordinates)
  (let ((order-size (loop for index = 0 then (+ index 5)
                          while (< index (length coordinates))
                          maximize (aref coordinates index))))
    (let ((order (make-array (1+ order-size))))
      (loop for index = 0 then (+ index 5)
            while (< index (length coordinates))
            do (let ((atom-index (aref coordinates index))
                     (atom-name (aref coordinates (+ index 1))))
                 (setf (aref order atom-index) atom-name)))
      order)))


(defun topology-for-context (context)
  "Find a topology that matches the context"
  (let* (result
         (focus-name (third context))
         (context-coupling-name (second context))
         (context-out-plug-name (chem:out-plug-name context-coupling-name)))
    (cando:walk-topologys
     (lambda (topology)
       (chem:walk-stereoisomer-atoms topology
                                     (lambda (stereoisomer-atoms)
                                       (when (eq (chem:get-name stereoisomer-atoms)
                                                 focus-name)
                                         (if (chem:has-in-plug topology)
                                             (let* ((in-plug (chem:get-in-plug topology))
                                                    (in-plug-name (chem:get-name in-plug)))
                                               (when (or
                                                      (chem:out-plug-name-matches-in-plug-name context-out-plug-name in-plug-name)
                                                      (and (eq in-plug-name :-origin)
                                                           (null context-coupling-name)))
                                                 (push topology result)))))))))
    (case (length result)
      (0 (error "Could not find topology for ~a" context))
      (1 (first result))
      (otherwise (error "There were multiple topologys ~s that match the context ~s" result context)))))

(defun equal-order (order1 order2)
  (loop for index from 0 below (length order1)
        for o1 = (aref order1 index)
        for o2 = (aref order2 index)
        for match = (or (null o1)
                        (null o2)
                        (eq o1 o2))
        unless match
          do (return-from equal-order nil))
  t)

(defun dump-coordinate-order (context-coordinates-ht)
  (flet ((dump-context-order (stream context names)
           (format stream "~50s " context)
           (loop for x across names
                 do (format stream "~7s" x))
           (terpri)))
    (let ((context-coords nil))
      (maphash (lambda (key coords)
                 (push (cons key coords) context-coords))
               context-coordinates-ht)
      (let ((sorted-context-coords (sort context-coords
                                         (lambda (x y)
                                           (string< (concatenate 'string (string (third (car x))) (string (second (car x))))
                                                    (concatenate 'string (string (third (car y))) (string (second (car y)))))))))
        (loop for entry in sorted-context-coords
              for context = (car entry)
              for context-focus = (third context)
              for topology = (topology-for-context context)
              for coords = (cdr entry)
              for topology-order = (topology-constitution-order topology)
              for joint-order = (topology-joint-template-order topology)
              for coordinate-order = (coordinate-order (car coords))
              unless (and (equal-order topology-order coordinate-order)
                          (equal-order topology-order joint-order))
                do (error "For context ~s there is a mismatch between the topology atom order, the joint-tree atom order or the coordinate atom order" context)
              do (let ((*print-pretty* nil))
                   (dump-context-order t (format nil "Top: ~s" (chem:get-name topology)) (topology-constitution-order topology)))
              do (dump-context-order t context (coordinate-order (car coords))))
        sorted-context-coords
        (format t "~%~%!!!!GOOD!!!   The orders all match!!!!!")))))
