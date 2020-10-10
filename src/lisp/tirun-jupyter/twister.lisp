(in-package :tirun-jupyter)

(defun json-keys (json)
  (mapcar #'car json))


(defun json-lookup (json key &optional context)
  (unless (eq (car json) :obj)
    (if context
        (error "In json context: ~s~%Looking for key: ~a - the json  must start with :obj but it doesn't json keys: ~s" context key (json-keys json))
        (error "Looking for key: ~a - the json  must start with :obj but it doesn't json keys: ~s" key (json-keys json))))
  (cdr (assoc key (cdr json) :test #'string=)))


(defun picked-bond-atoms (picked agg)
  (let ((atom-vec (chem:map-atoms 'vector #'identity agg)))
    (let* ((bond (json-lookup picked "bond"))
           (atom1-json (json-lookup picked "atom1"))
           (atom1-index (json-lookup atom1-json "index"))
           (atom2-json (json-lookup picked "atom2"))
           (atom2-index (json-lookup atom2-json "index"))
           (atom1 (aref atom-vec atom1-index))
           (atom2 (aref atom-vec atom2-index)))
      (values atom1 atom2))))


(defun crd (agg)
  (let ((arr (make-array (* 3 (chem:number-of-atoms tl)) :element-type 'single-float))
        (index -1))
    (cando:do-atoms (atm agg)
      (let ((pos (chem:get-position atm)))
        (setf (aref arr (incf index)) (float (geom:vx pos) 1.0s0)
              (aref arr (incf index)) (float (geom:vy pos) 1.0s0)
              (aref arr (incf index)) (float (geom:vz pos) 1.0s0))))
    arr))

(defun drive (widget agg twister angle)
  (chem:twister-rotate-absolute twister angle)
  (nglv:set-coordinates widget (list (cons 0 (crd agg)))))

;; This should probably be moved into the application?  Or make it a closed over variable in connect-nglview-slider?

(defparameter *twister* nil)

;; connect-nglview-slider is a demonstration - to be part of the tirun-jupyter app it
;; needs to listen to an app slot that represents the ligand that we are currently manipulating
;;l
(defun connect-nglview-slider (agg nglview slider)
  (jupyter-widgets:observe nglview :picked
                           (lambda (instance type name old-value new-value source)
                             (when new-value
                               (let* ((picked new-value)
                                      (bond (json-lookup picked "bond")))
                                 (if bond
                                     (multiple-value-bind (bond-a bond-b)
                                         (picked-bond-atoms picked agg)
                                         (setf *twister* (chem:make-twister))
                                       (chem:twister-define-for-bond *twister* bond-a bond-b))
                                     (setf *twister* nil))))))
  (jupyter-widgets:observe slider :value
                           (lambda (instance type name old-value new-value source)
                             (jupyter-widgets:with-output
                               *log*
                               #+(or)(format t "slider value: ~a  *twister* -> ~a~%" new-value *twister*)
                               (when (and new-value *twister*)
                                 (drive nglview agg *twister* (* 0.1 new-value))
                                 #+(or)(format t "Called drive~%"))))))
