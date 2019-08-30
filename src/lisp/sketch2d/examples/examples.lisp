(ql:quickload :sketch2d)

(defparameter *agg* (load-mol2 "017.mol2"))

(defparameter *mol* (chem:content-at *agg* 0))

(defparameter *d* (sketch2d::sketch2d-dynamics *mol*))o
