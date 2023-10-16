(in-package :topology-jupyter)

(defun sketch-svg (topology &rest args)
  "Generate an svg sketch of the topology - send args to sketch2d"
  (let ((mol (topology:build-one-molecule-for-topology topology)))
    (sketch2d:svg (apply 'sketch2d:sketch2d mol args))))
