(in-package :cando-widgets)

(defgeneric show (instance &rest kwargs &key &allow-other-keys)
  (:documentation "Display a residue, molecule, aggregate or a trajectory in a Jupyter notebook"))


(defun add-bounding-box (component aggregate)
  (unless (chem:bounding-box-bound-p aggregate)
    (error "There is no bounding-box defined for aggregate"))
  (let* ((bounding-box (chem:bounding-box aggregate))
         (center (chem:get-bounding-box-center bounding-box))
         (x-width/2 (/ (chem:get-x-width bounding-box) 2.0))
         (y-width/2 (/ (chem:get-y-width bounding-box) 2.0))
         (z-width/2 (/ (chem:get-z-width bounding-box) 2.0))
         (minx (- (geom:vx center) x-width/2))
         (miny (- (geom:vy center) y-width/2))
         (minz (- (geom:vz center) z-width/2))
         (maxx (+ (geom:vx center) x-width/2))
         (maxy (+ (geom:vy center) y-width/2))
         (maxz (+ (geom:vz center) z-width/2)))
    (setf (ngl:representations component)
          (append (ngl:representations component)
                  (list (make-instance 'ngl:buffer-representation
                                       :name "bounding-box"
                                       :buffer (list :type "wideline"
                                                     :position1 (vector minx miny minz
                                                                        minx miny minz
                                                                        minx miny minz
                                                                        minx miny maxz
                                                                        minx miny maxz
                                                                        minx maxy minz
                                                                        minx maxy minz
                                                                        maxx miny minz
                                                                        maxx miny minz
                                                                        minx maxy maxz
                                                                        maxx miny maxz
                                                                        maxx maxy minz)
                                                     :position2 (vector minx miny maxz
                                                                        minx maxy minz
                                                                        maxx miny minz
                                                                        minx maxy maxz
                                                                        maxx miny maxz
                                                                        minx maxy maxz
                                                                        maxx maxy minz
                                                                        maxx miny maxz
                                                                        maxx maxy minz
                                                                        maxx maxy maxz
                                                                        maxx maxy maxz
                                                                        maxx maxy maxz)
                                                     :color #(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
                                                     :color2 #(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0))))))))


(defmethod show ((instance chem:aggregate) &rest kwargs &key &allow-other-keys)
  (let ((component (make-instance 'ngl:structure
                                  :name (symbol-name (chem:get-name instance))
                                  :value (chem:aggregate-as-mol2-string instance t)
                                  :ext "mol2")))
    (when (chem:bounding-box-bound-p instance)
      (add-bounding-box component instance))
    (jw:display (apply #'make-instance 'ngl:stage
                       :clip-dist 0
                       :background-color "white"
                       :components (list component
                                         (make-instance 'ngl:shape :primitives (getf kwargs :shapes)))
                       kwargs))
    (ngl:auto-view component)
    (values)))

(defmethod show ((molecule chem:molecule) &rest kwargs &key &allow-other-keys)
  (let ((agg (chem:make-aggregate nil)))
    (chem:add-matter agg molecule)
    (apply 'show agg kwargs)))

(defun isolate-residue (residue)
  (let* ((agg (chem:make-aggregate nil))
         (mol (chem:make-molecule nil))
         (new-to-old (make-hash-table))
         (res (chem:copy residue new-to-old)))
    (let (break-bonds)
      (chem:map-bonds
       'nil
       (lambda (a1 a2 order)
         (unless (and (chem:contains-atom res a1) (chem:contains-atom res a2))
           (push (list a1 a2) break-bonds)))
       res)
      (loop for bond in break-bonds
            for a1 = (first bond)
            for a2 = (second bond)
            do (chem:remove-bond-to a1 a2)))
    (chem:add-matter agg mol)
    (chem:add-matter mol res)
    agg))

(defmethod show ((residue chem:residue)  &rest kwargs &key &allow-other-keys)
  (let ((agg (isolate-residue residue)))
    (apply 'show agg kwargs)))

(defmethod show ((sketch sketch2d:sketch2d) &rest kwargs &key &allow-other-keys)
  (show (sketch2d:svg sketch)))

(defmethod show ((sketch sketch2d:sketch-svg) &rest kwargs &key &allow-other-keys)
  (jupyter:svg (sketch2d:render-svg-to-string sketch)))

(defmethod show ((trajectory dynamics:trajectory) &rest kwargs &key &allow-other-keys)
  (change-class trajectory 'cando-trajectory)
  (apply #'nglv:make-nglwidget :structure trajectory kwargs))

(defmethod show ((dynamics dynamics:simulation) &rest kwargs &key &allow-other-keys)
  (let ((trajectory (dynamics:make-trajectory dynamics)))
    (change-class trajectory 'cando-trajectory)
    (apply 'show trajectory kwargs)))

