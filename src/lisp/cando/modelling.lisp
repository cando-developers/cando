(in-package :cando)

(defun unbuilt-vec3-p (v)
  (= (geom:vx v) (geom:vy v) (geom:vz v) 0.0))

(defun unbuilt-atoms (matter)
  (let (unbuilt)
    (chem:map-atoms 'nil (lambda (atom)
                      (when (unbuilt-vec3-p (chem:get-position atom))
                        (push atom unbuilt)))
                    matter)
    unbuilt))

(defun count-unbuilt-atoms (matter)
  (let ((unbuilt 0))
    (chem:map-atoms
     'nil
     (lambda (atom)
       (when (unbuilt-vec3-p (chem:get-position atom))
         (incf unbuilt)))
     matter)
    unbuilt))

(defun random-unit-vector ()
  "Use answer from here to get a random vector on a sphere.
https://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d
"
  (let* ((theta (random (* 2.0 PI)))
         (z (- (random 2.0) 1.0))
         (sqrt-1-z2 (sqrt (- 1.0 (* z z))))
         (x (* sqrt-1-z2 (cos theta)))
         (y (* sqrt-1-z2 (sin theta)))
         (point (geom:vec x y z)))
    point))

(defun random-3d-offset (radius center)
  "Use answer from here to get a random vector on a sphere.
https://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d
"
  (geom:v+ center (geom:v* (random-unit-vector) radius)))


(defun simple-build-unbuilt-atoms (matter)
  "Build coordinates for the unbuilt atoms of matter by moving them one angstrom
away from their bonded neighbor that has coordinates in a random direction."
  (let ((remaining (count-unbuilt-atoms matter)))
    (loop (let ((cur-remaining remaining)
                (built 0))
            (chem:map-atoms
             'nil
             (lambda (atom)
               (when (unbuilt-vec3-p (chem:get-position atom))
                 (let ((bonds (chem:bonds-as-list atom)))
                   (if (> (length bonds) 0) ; we have at least one bond
                       (let ((built-neighbor (loop for b in bonds
                                         for o = (chem:get-other-atom b atom)
                                         when (not (unbuilt-vec3-p (chem:get-position o)))
                                           return o)))
                         (when built-neighbor
                           (unless (unbuilt-vec3-p (chem:get-position built-neighbor))
                             (let ((pos (random-3d-offset 1.0 (chem:get-position built-neighbor))))
                               (chem:set-position atom pos)
                               (incf built)
                               (decf cur-remaining)))))
                       (progn
                         (warn "There is a single atom with no neighbors ~a" atom)
                         (decf cur-remaining) ; remove this from the remaining atoms
                         )))))
             matter)
            (when (= cur-remaining 0)
              (return-from simple-build-unbuilt-atoms))
            (when (= cur-remaining remaining)
              (error "No atoms were built - there are molecules with no built atoms"))
            (setf remaining cur-remaining)) ; repeat the loop
          )))
  
