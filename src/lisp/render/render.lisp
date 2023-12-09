(in-package :render)


(defvar *atom-colors* (make-hash-table))


(defun vector-subtract (vec1 vec2)
  "Subtract vector vec2 from vec1."
  (vector (- (aref vec1 0) (aref vec2 0))
          (- (aref vec1 1) (aref vec2 1))
          (- (aref vec1 2) (aref vec2 2))))

(defun vector-magnitude (vec)
  "Calculate the magnitude of the vector."
  (sqrt (+ (expt (aref vec 0) 2)
           (expt (aref vec 1) 2)
           (expt (aref vec 2) 2))))

(defun vector-normalize (vec)
  "Normalize the vector."
  (let ((mag (vector-magnitude vec)))
    (vector (/ (aref vec 0) mag)
            (/ (aref vec 1) mag)
            (/ (aref vec 2) mag))))

(defun vector-dot (vec1 vec2)
  "Dot product of vec1 and vec2."
  (+ (* (aref vec1 0) (aref vec2 0))
     (* (aref vec1 1) (aref vec2 1))
     (* (aref vec1 2) (aref vec2 2))))

(defun vector-cross (vec1 vec2)
  "Cross product of vec1 and vec2."
  (vector (- (* (aref vec1 1) (aref vec2 2)) (* (aref vec1 2) (aref vec2 1)))
          (- (* (aref vec1 2) (aref vec2 0)) (* (aref vec1 0) (aref vec2 2)))
          (- (* (aref vec1 0) (aref vec2 1)) (* (aref vec1 1) (aref vec2 0)))))

(defun bond-rotation-axis (atom1 atom2)
  "Calculate the rotation axis for the bond."
  (let* ((bond-vector (vector-subtract (pos atom2) (pos atom1)))
         (y-axis (vector 0 1 0))
         (rotation-axis (vector-cross y-axis bond-vector)))
    (vector-normalize rotation-axis)))

(defun bond-rotation-angle (atom1 atom2)
  "Calculate the rotation angle for the bond."
  (let* ((bond-vector (vector-subtract (pos atom2) (pos atom1)))
         (y-axis (vector 0 1 0)))
    (acos (vector-dot (vector-normalize bond-vector) y-axis))))
(defun rotation-axis-x (atom1 atom2)
  "Calculate the x component of the rotation axis for the bond."
  (let ((axis (bond-rotation-axis atom1 atom2)))
    (aref axis 0)))

(defun rotation-axis-y (atom1 atom2)
  "Calculate the y component of the rotation axis for the bond."
  (let ((axis (bond-rotation-axis atom1 atom2)))
    (aref axis 1)))

(defun rotation-axis-z (atom1 atom2)
  "Calculate the z component of the rotation axis for the bond."
  (let ((axis (bond-rotation-axis atom1 atom2)))
    (aref axis 2)))

(defun rotation-angle (atom1 atom2)
  "Calculate the rotation angle for the bond."
  (bond-rotation-angle atom1 atom2))

;; Definitions for bond-rotation-axis and bond-rotation-angle
;; along with other necessary vector functions stay the same

(defun distance (atom1 atom2)
  "Calculate the distance between two atoms."
  (let ((dx (- (x-pos atom2) (x-pos atom1)))
        (dy (- (y-pos atom2) (y-pos atom1)))
        (dz (- (z-pos atom2) (z-pos atom1))))
    (sqrt (+ (* dx dx) (* dy dy) (* dz dz)))))

(defun midpoint (atom1 atom2)
  "Calculate the midpoint between two atoms."
  (vector (/ (+ (x-pos atom1) (x-pos atom2)) 2.0)
          (/ (+ (y-pos atom1) (y-pos atom2)) 2.0)
          (/ (+ (z-pos atom1) (z-pos atom2)) 2.0)))

(defun midpoint-x (atom1 atom2)
  "Calculate the x-coordinate of the midpoint between two atoms."
  (aref (midpoint atom1 atom2) 0))

(defun midpoint-y (atom1 atom2)
  "Calculate the y-coordinate of the midpoint between two atoms."
  (aref (midpoint atom1 atom2) 1))

(defun midpoint-z (atom1 atom2)
  "Calculate the z-coordinate of the midpoint between two atoms."
  (aref (midpoint atom1 atom2) 2))


;; The rest of the molecule rendering code remains the same.

(defclass bond ()
  ((atom1 :initarg :atom1 :accessor atom1)
   (atom2 :initarg :atom2 :accessor atom2)
   (bond-order :initarg :bond-order :accessor bond-order)))

(defclass atom ()
  ((element :initarg :element :accessor element)
   (pos :initarg :pos :accessor pos)))

(defclass molecule ()
  ((atoms :initarg :atoms :initform (make-array 16 :adjustable t :fill-pointer 0)
          :accessor atoms)
   (bonds :initarg :bonds :initform (make-array 16 :adjustable t :fill-pointer 0)
          :accessor bonds)))


(defvar *atom-colors* (make-hash-table))

(defun initialize-atom-colors ()
  (setf (gethash :C *atom-colors*) '(0.3 0.3 0.3)   ; Black for Carbon
        (gethash :N *atom-colors*) '(0.0 0.0 1.0)   ; Blue for Nitrogen
        (gethash :O *atom-colors*) '(1.0 0.0 0.0)   ; Red for Oxygen
        (gethash :H *atom-colors*) '(1.0 1.0 1.0)   ; White for Hydrogen
        (gethash :S *atom-colors*) '(1.0 1.0 0.0)   ; Yellow for Sulfur
        (gethash :|Br| *atom-colors*) '(0.6 0.3 0.0) ; Brown for Bromine
        (gethash :F *atom-colors*) '(0.5 1.0 0.5)   ; Light Green for Fluorine
        ;; Add more elements and their colors here
        ))

(defun get-color (element)
  (gethash element *atom-colors* '(1.0 1.0 1.0))) ; Default color (white) if not found

(defun x-pos (atom)
  (aref (pos atom) 0))

(defun y-pos (atom)
  (aref (pos atom) 1))

(defun z-pos (atom)
  (aref (pos atom) 2))

(defun render-atom (atom sphere-radius)
  (let ((color (get-color (element atom))))
    (format nil 
      "Transform {
        translation ~A ~A ~A
        children [
          Shape {
            geometry Sphere { radius ~A }
            appearance Appearance {
              material Material { diffuseColor ~A ~A ~A }
            }
          }
        ]
      }" 
      (x-pos atom) (y-pos atom) (z-pos atom)
      sphere-radius
      (first color) (second color) (third color))))

(defun render-bond (bond cylinder-radius)
  (let ((atom1 (atom1 bond))
        (atom2 (atom2 bond)))
    (format nil 
      "Transform {
        children [
          Shape {
            geometry Cylinder { 
              radius ~A 
              height ~A 
            }
            appearance Appearance {
              material Material { diffuseColor 0.5 0.5 0.5 }
            }
          }
        ]
        translation ~A ~A ~A
        rotation ~A ~A ~A ~A
      }"
      cylinder-radius
      (distance atom1 atom2)
      (midpoint-x atom1 atom2) (midpoint-y atom1 atom2) (midpoint-z atom1 atom2)
      (rotation-axis-x atom1 atom2) (rotation-axis-y atom1 atom2) (rotation-axis-z atom1 atom2) (rotation-angle atom1 atom2))))

(defun render-molecule (molecule sphere-radius cylinder-radius file)
  (let ((atoms (atoms molecule))
        (bonds (bonds molecule)))
    (with-open-file (stream file
                            :direction :output 
                            :if-exists :supersede
                            :if-does-not-exist :create)
      (format stream "#VRML V2.0 utf8~%")
      (format stream "Group { children [~%")
      (dolist (atom (coerce atoms 'list))
        (format stream "~A~%" (render-atom atom sphere-radius)))
      (dolist (bond (coerce bonds 'list))
        (format stream "~A~%" (render-bond bond cylinder-radius)))
      (format stream "] }"))))


(defun translate-aggregate (aggregate)
  (let ((atoms (make-array 16 :adjustable t :fill-pointer 0))
        (bonds (make-array 16 :adjustable t :fill-pointer 0))
        (atom-map (make-hash-table)))
    (format t "In translate-aggregate~%")
    (chem:do-atoms (atm aggregate)
      (let* ((pos (chem:get-position atm))
             (natm (make-instance 'atom
                                 :pos (vector (geom:vx pos) (geom:vy pos) (geom:vz pos))
                                 :element (chem:get-element atm))))
        (setf (gethash atm atom-map) natm)
        (vector-push-extend natm atoms)))
    (chem:map-bonds nil
                    (lambda (a1 a2 bo aggregate)
                      (let* ((na1 (gethash a1 atom-map))
                             (na2 (gethash a2 atom-map))
                             (nbond (make-instance 'bond
                                                   :atom1 na1
                                                   :atom2 na2
                                                   :bond-order bo)))
                        (vector-push-extend nbond bonds)))
                      aggregate)
    (make-instance 'molecule :atoms atoms
                             :bonds bonds)))

(defun vrml-render (aggregate file &key (sphere 0.2) (cylinder 0.2))
  (let ((mmol (translate-aggregate aggregate)))
    (render-molecule mmol sphere cylinder file)))

(initialize-atom-colors)
;; (render-molecule your-molecule-data sphere-radius cylinder-radius)
 
