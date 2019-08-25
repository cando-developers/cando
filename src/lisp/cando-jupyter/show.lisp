(in-package :cando-jupyter)

(defgeneric show (matter-or-trajectory &rest kwargs &key &allow-other-keys)
  (:documentation "Display a residue, molecule, aggregate or a trajectory in a Jupyter notebook"))

(defmethod show ((aggregate chem:aggregate) &rest kwargs &key &allow-other-keys)
  (let ((structure (make-instance 'cando-structure :matter aggregate)))
    (apply #'nglv:make-nglwidget :structure structure kwargs)))

(defmethod show ((molecule chem:molecule) &rest kwargs &key &allow-other-keys)
  (let ((agg (chem:make-aggregate nil)))
    (chem:add-matter agg molecule)
    (apply 'show agg kwargs)))

(defmethod show ((residue chem:residue)  &rest kwargs &key &allow-other-keys)
  (let ((agg (chem:make-aggregate nil))
        (mol (chem:make-molecule nil)))
    (chem:add-matter agg mol)
    (chem:add-matter mol residue)
    (apply 'show agg kwargs)))
  
(defun repr (widget representation &optional (selection "all"))
  (funcall (find-symbol "ADD-REPRESENTATION" :nglv) widget representation :selection selection))

(defun cl-jupyter-kernel-start (&optional connection-file-name)
  (let ((amber-home
          (namestring (or (if (ext:getenv "AMBERHOME")
                              (probe-file (ext:getenv "AMBERHOME"))
                              "/usr/local/amber/")
                          (probe-file "/home/app/amber16-data/")
                          "/dev/null"))))
    (setf (logical-pathname-translations "amber")
          (list (list "**;*.*" (concatenate 'string amber-home "/**/*.*"))))
    (format t "Setting amber host pathname translation -> ~a~%" amber-home))
  (let ((*readtable* (copy-readtable)))
    (let ((cl-jup (find-symbol "KERNEL-START" "CL-JUPYTER")))
      (if cl-jup
          (if connection-file-name
              (funcall cl-jup connection-file-name)
              (funcall cl-jup))
          (error "cl-jupyter is not installed")))))


