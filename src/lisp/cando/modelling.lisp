(in-package :cando)

(defun move-geometric-center-to-origin (matter)
  "Move the geometric center of the matter to the origin"
  (let* ((geometric-center (chem:geometric-center matter))
         (delta (geom:v* geometric-center -1.0)))
    (chem:apply-transform-to-atoms matter (geom:make-m4-translate delta))))

(defun move-geometric-center-to-point (matter point)
  "Move the geometric center of the matter to the point"
  (let* ((geometric-center (chem:geometric-center matter))
         (delta (geom:v* geometric-center -1.0)))
    (chem:apply-transform-to-atoms matter (geom:make-m4-translate delta))
    (chem:apply-transform-to-atoms matter (geom:make-m4-translate point))))

(defun unbuilt-vec3-p (v)
  (not (geom:is-defined v)))

(defun heavy-atom-p (a)
  (not (member (chem:get-element a) '(:h :d :t))))

(defun unbuilt-atoms (matter)
  (let (unbuilt)
    (chem:map-atoms 'nil (lambda (atom)
                      (when (unbuilt-vec3-p (chem:get-position atom))
                        (push atom unbuilt)))
                    matter)
    unbuilt))

(defun count-unbuilt-atoms (matter &optional heavy)
  (let ((unbuilt 0))
    (chem:map-atoms
     'nil
     (lambda (atom)
       (when (and (or (null heavy) (and heavy (heavy-atom-p atom)))
                  (unbuilt-vec3-p (chem:get-position atom))
         (incf unbuilt))))
     matter)
    unbuilt))

(defun random-unit-vector ()
  "Use answer from here to get a random vector on a sphere.
https://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d
"
  (let* ((theta (random (* 2.0 (coerce PI 'double-float))))
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


(defun simple-build-unbuilt-heavy-atoms (matter)
  "Build coordinates for the unbuilt atoms of matter by moving them one angstrom
away from their bonded neighbor that has coordinates in a random direction."
  (let ((remaining-heavy (count-unbuilt-atoms matter t)))
    (when (> remaining-heavy 0)
      (format t "There are ~a unbuilt heavy atoms of a total of ~a~%" remaining-heavy (chem:number-of-atoms matter))
      (loop (let ((cur-remaining remaining-heavy)
                  (built 0))
              (chem:map-atoms
               'nil
               (lambda (atom)
                 (when (and (heavy-atom-p atom) (unbuilt-vec3-p (chem:get-position atom)))
                   (format t "Unbuilt atom ~a~%" atom)
                   (let ((bonds (chem:bonds-as-list atom)))
                     (if (> (length bonds) 0) ; we have at least one bond
                         (let ((built-neighbor (loop for b in bonds
                                                     for o = (chem:bond/get-other-atom b atom)
                                                     when (and (heavy-atom-p o)
                                                               (not (unbuilt-vec3-p (chem:get-position o))))
                                                       return o)))
                           (when built-neighbor
                             (let ((pos (random-3d-offset 1.0 (chem:get-position built-neighbor))))
                               (chem:set-position atom pos)
                               (incf built)
                               (decf cur-remaining))))
                         (progn
                           (warn "There is a single atom with no neighbors ~a" atom)
                           (decf cur-remaining) ; remove this from the remaining-heavy atoms
                           )))))
               matter)
              (format t "Just built ~a heavy atoms and there are ~a remaining~%" built cur-remaining)
              (when (= cur-remaining 0)
                (return-from simple-build-unbuilt-heavy-atoms))
              (when (= cur-remaining remaining-heavy)
                (error "No atoms were built - there are molecules with no built atoms"))
              (setf remaining-heavy cur-remaining)) ; repeat the loop
            ))))

