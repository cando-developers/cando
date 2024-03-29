(in-package :cando-widgets)


(defgeneric make-ngl-trajectories (instance &rest initargs &key &allow-other-keys)
  (:documentation "Create a ngl-clj trajectories from an object")
  (:method (instance &rest initargs &key &allow-other-keys)
    (cond
      ((getf initargs :as-trajectory) ; Create an trajectory for the structure data
        (list (apply #'make-instance 'ngl:trajectory initargs)))
      ((getf initargs :trajectory)
        (apply #'make-ngl-trajectories (getf initargs :trajectory) initargs)))))


(defmethod make-ngl-trajectories ((instance netcdf:netcdf) &rest initargs &key &allow-other-keys)
  (let* ((atom-count (netcdf::dimension instance "atom"))
         (coords (static-vectors:make-static-vector (* 3 atom-count)
                                                    :element-type 'single-float))
         (frame-count (netcdf:get-dimension instance "frame")))
    (list (apply #'make-instance 'ngl:trajectory
                                 :value (concatenate 'string "jupyter:" (getf initargs :name ""))
                                 :end frame-count
                                 :trajectory-count (lambda (inst) ; This function returns the frame count
                                                     (declare (ignore inst))
                                                     frame-count)
                                 :trajectory-frame (lambda (inst i atom-indices) ; This function returns a specific frame. Atom subsets not implemented yet.
                                                     (declare (ignore inst atom-indices))
                                                     (netcdf:get-vara-float instance "coordinates"
                                                                            (vector i 0 0)
                                                                            (vector 1 atom-count 3)
                                                                            coords)
                                                     (values i nil coords frame-count))
                                 initargs))))


(defmethod make-ngl-trajectories ((instance dynamics:trajectory) &rest initargs &key &allow-other-keys)
  (list (apply #'make-instance 'ngl:trajectory
                               :value (concatenate 'string "jupyter:" (getf initargs :name ""))
                               :end (length (dynamics:coordinates instance))
                               :trajectory-count (lambda (inst) ; This function returns the frame count
                                                   (declare (ignore inst))
                                                   (1- (length (dynamics:coordinates instance))))
                               :trajectory-frame (lambda (inst i atom-indices) ; This function returns a specific frame. Atom subsets not implemented yet.
                                                   (declare (ignore inst atom-indices))
                                                   (values i
                                                           nil
                                                           (elt (dynamics:coordinates instance) i)
                                                           (1- (length (dynamics:coordinates instance)))))
                               initargs)))

(defmethod make-ngl-trajectories ((instance chem:trajectory) &rest initargs &key &allow-other-keys)
  (list (apply #'make-instance 'ngl:trajectory
                               :value (concatenate 'string "jupyter:" (getf initargs :name ""))
                               :end (chem:trajectory/number-of-trajectory-frames instance)
                               :trajectory-count (lambda (inst) ; This function returns the frame count
                                                   (declare (ignore inst))
                                                   (1- (chem:trajectory/number-of-trajectory-frames instance)))
                               :trajectory-frame (lambda (inst i atom-indices) ; This function returns a specific frame. Atom subsets not implemented yet.
                                                   (declare (ignore inst atom-indices))
                                                   (format *error-output* "Sending frame ~a~%" i)
                                                   (values i
                                                           nil
                                                           (chem:trajectory-frame/get-coordinates
                                                            (chem:trajectory/get-trajectory-frame instance i))
                                                           (1- (chem:trajectory/number-of-trajectory-frames instance))))
                               initargs)))


(defgeneric make-ngl-structure (instance &rest initargs &key &allow-other-keys)
  (:documentation "Create a ngl-clj structure from an object"))


(defmethod make-ngl-structure ((instance chem:aggregate) &rest initargs &key &allow-other-keys)
  (values (apply #'ngl:make-structure :ext "mol2"
                 :value (chem:aggregate-as-mol2-string instance t)
                 :trajectories (apply #'make-ngl-trajectories instance initargs)
                 initargs)
          instance))

(defmethod make-ngl-structure ((instance chem:molecule) &rest initargs &key &allow-other-keys)
  (let ((agg (chem:make-aggregate (chem:get-name instance))))
    (chem:add-matter agg instance)
    (values (apply #'ngl:make-structure :ext "mol2"
                   :value (chem:aggregate-as-mol2-string agg t)
                   :trajectories (apply #'make-ngl-trajectories instance initargs)
                   initargs)
            agg)))

(defmethod make-ngl-structure ((instance dynamics:trajectory) &rest initargs &key &allow-other-keys)
  (let ((agg (dynamics:matter instance)))
    (values (apply #'ngl:make-structure :ext "mol2"
                   :value (chem:aggregate-as-mol2-string agg t)
                   :trajectories (apply #'make-ngl-trajectories instance initargs)
                   initargs)
            agg)))

(defmethod make-ngl-structure ((instance chem:trajectory) &rest initargs &key &allow-other-keys)
  (let ((agg (chem:trajectory/get-matter instance)))
    (values (apply #'ngl:make-structure :ext "mol2"
                   :value (chem:aggregate-as-mol2-string agg t)
                   :trajectories (apply #'make-ngl-trajectories instance initargs)
                   initargs)
            agg)))

(defmethod make-ngl-structure ((instance leap.topology:amber-topology-trajectory-pair) &rest initargs &key &allow-other-keys)
  (let ((agg (leap.topology:read-frame-into-atoms instance 0)))
    (values (apply #'ngl:make-structure :ext "mol2"
                   :value (chem:aggregate-as-mol2-string agg t)
                   :trajectories (apply #'make-ngl-trajectories (leap.topology:netcdf instance) initargs)
                   initargs)
            agg)))
