;;; Common code
(setf *default-pathname-defaults* #P"~/Development/cando/extensions/cando/src/tests/aromaticity/")
(setup-default-paths)
(easy-gaff)
(load-smirnoff-params #P"~/Development/openforcefield/openforcefield/data/forcefield/smirnoff99Frosst.offxml")


;;; Try a compilicated example
(defparameter *agg* (load-chem-draw-aggregate "ligand8a.cdxml"))
(chem:fill-in-implicit-hydrogens *agg*)
(set-stereoisomer-func (cando:gather-stereocenters *agg*) (lambda (c) :R))
(defparameter *mol* (chem:content-at *agg* 0))


(set-force-field *mol* :smirnoff)
(save-amber-parm *agg* "/tmp/test.parm" "/tmp/test.crd")



(set-force-field *mol* :smirnoff)
(:= *e-smirnoff* (chem:make-energy-function *agg* :use-excluded-atoms t :assign-types t))




(load-chem-draw-aggregate "ligand8a.cdxml")


;;; Try benzene
;; Evaluate common code
(defparameter *agg* (load-chem-draw-aggregate "benzene.cdxml"))
(chem:fill-in-implicit-hydrogens *agg*)
(set-stereoisomer-func (cando:gather-stereocenters *agg*) (lambda (c) :R))
(defparameter *mol* (chem:content-at *agg* 0))
(set-force-field *mol* :smirnoff)
(save-amber-parm *agg* "/tmp/test.parm" "/tmp/test.crd")

(defparameter *g* (chem:make-chem-info-graph (chem:compile-smarts "[C:1]1=[C:2]-[C:3]=[C:4]-[C:5]=[C:6]1")))
