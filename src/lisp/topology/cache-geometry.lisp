(in-package #:topology)

(defclass resolved-cache-component ()
  ((declaration :initarg :declaration :reader resolved-cache-declaration)
   (monomer :initarg :monomer :reader resolved-cache-monomer)
   (kind :initarg :kind :reader resolved-cache-kind))
  (:documentation "A cache declaration bound to one assembled conformation's geometry."))

(defclass resolved-ring-cache-component (resolved-cache-component)
  ((driver :initarg :driver :reader resolved-cache-driver)))

(defclass resolved-dihedral-cache-component (resolved-cache-component)
  ((joint :initarg :joint :reader resolved-cache-joint)))

(defun cache-ring-drivers (manipulator)
  "Collect individual rings without changing their atom-tree ordering."
  (remove-duplicates
   (loop for driver in (append (other-drivers manipulator)
                               (when (typep manipulator 'focused-manipulator)
                                 (focus-drivers manipulator)))
         append (typecase driver
                  (ring-driver (list driver))
                  (fused-ring-driver
                   (list (ring1 (rings driver)) (ring2 (rings driver))))))
   :test #'eq))

(defun resolve-cache-component (declaration monomer joint ring-drivers)
  "Bind DECLARATION to its selected JOINT, validating geometry before measurement.
RING-DRIVERS are individual ring drivers. The selector never rotates RING-JOINTS."
  (unless joint
    (error "Cache field ~s: atom ~s is missing in monomer ~s"
           (cache-component-name declaration) (cache-component-atom-name declaration) monomer))
  (ecase (cache-component-kind declaration)
    (:ring
     (let ((matches (remove-duplicates
                     (remove-if-not (lambda (driver)
                                      (member joint (ring-joints driver) :test #'eq))
                                    ring-drivers)
                     :test #'eq)))
       (unless (= 1 (length matches))
         (error "Cache field ~s: atom ~s in monomer ~s selects ~d rings; expected exactly one"
                (cache-component-name declaration) (cache-component-atom-name declaration)
                monomer (length matches)))
       (let* ((driver (first matches))
              (kind (etypecase driver (ring5-driver :ring5) (ring6-driver :ring6))))
         (make-instance 'resolved-ring-cache-component :declaration declaration
                        :monomer monomer :kind kind :driver driver))))
    (:dihedral
     ;; Measurement uses Cartesian coordinates, including Amber's XYZ joints.
     ;; Driver eligibility is stricter because drivers modify internal coordinates.
     (unless (and (typep joint '(or kin:bonded-joint kin:xyz-joint))
                  (loop for ancestor = joint then (kin:joint/parent ancestor)
                        repeat 3
                        always (kin:joint/parent-bound-p ancestor)))
       (error "Cache field ~s: atom ~s in monomer ~s needs a dihedral with three ancestors"
              (cache-component-name declaration) (cache-component-atom-name declaration) monomer))
     (make-instance 'resolved-dihedral-cache-component :declaration declaration
                    :monomer monomer :kind :dihedral :joint joint))))

(defun resolve-topology-cache-components (manipulator monomer topology)
  "Resolve the topology's declarations against this manipulator's assembled monomer."
  (let ((residue (assembler-atresidue (assembler manipulator) monomer))
        (rings (cache-ring-drivers manipulator)))
    (loop for declaration in (topology-cache-declarations topology)
          collect (resolve-cache-component
                   declaration monomer
                   (joint-with-name residue (cache-component-atom-name declaration) nil)
                   rings))))

(defgeneric measure-cache-component (component externals)
  (:documentation "Return physical parameter values, in the binning API's units and order."))

(defmethod measure-cache-component ((component resolved-ring-cache-component) externals)
  (let ((indexes (get-indexes3 (ring-joints (resolved-cache-driver component)))))
    (ecase (resolved-cache-kind component)
      (:ring5 (multiple-value-list (chem:cpring-5-coordinates-to-cremer-pople externals indexes)))
      (:ring6 (multiple-value-list (chem:cpring-6-coordinates-to-cremer-pople externals indexes))))))

(defmethod measure-cache-component ((component resolved-dihedral-cache-component) externals)
  (let* ((jd (resolved-cache-joint component))
         (jc (kin:joint/parent jd))
         (jb (kin:joint/parent jc))
         (ja (kin:joint/parent jb)))
    (list (geom:calculate-dihedral-array
           (kin:joint/position-index-x3 ja) (kin:joint/position-index-x3 jb)
           (kin:joint/position-index-x3 jc) (kin:joint/position-index-x3 jd) externals))))

(defun measure-binned-shape-key-cache (foldamer components externals)
  "Measure resolved components and return a named cache using FOLDAMER's binning policy."
  (make-binned-shape-key-cache
   (loop for component in components
         collect (cons (cache-component-name (resolved-cache-declaration component))
                       (bin-shape-component foldamer (resolved-cache-kind component)
                                            (measure-cache-component component externals))))))

(defun decode-cache-component (foldamer component bins)
  "Decode one component with its declared kind and foldamer policy."
  (unbin-shape-component foldamer (resolved-cache-kind component) bins))

(defun write-ring-cache-component (foldamer component bins cpinternals)
  "Write decoded CP parameters, preserving supporting bond lengths and angles.
The existing ring reconstruction helpers can then use CPINTERNALS unchanged."
  (check-type component resolved-ring-cache-component)
  (let* ((parameters (decode-cache-component foldamer component bins))
         (driver (resolved-cache-driver component))
         (offset (begin (cremer-pople-range driver))))
    (check-shape-component-values (resolved-cache-kind component) parameters)
    (loop for value in parameters for index from offset
          do (setf (aref cpinternals index) (coerce value 'double-float))))
  cpinternals)
