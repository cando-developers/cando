(setf *default-pathname-defaults* #P"~/Development/cando/extensions/cando/src/tests/alanine-demo/")
(setup-default-paths)
(easy-gaff)
(load-smirnoff-params #P"~/Development/openforcefield/openforcefield/data/forcefield/smirnoff99Frosst.offxml")
(:= *agg* (load-chem-draw-aggregate "alanine2.cdxml"))
(chem:fill-in-implicit-hydrogens *agg*)
(set-stereoisomer-func (cando:gather-stereocenters *agg*) (lambda (c) :R))
(:= *mol* (chem:content-at *agg* 0))

(set-force-field *mol* :smirnoff)
(save-amber-parm *agg* "/tmp/agg.parm" "/tmp/agg.crd")



(set-force-field *mol* :smirnoff)
(:= *e-smirnoff* (chem:make-energy-function *agg* :use-excluded-atoms t :assign-types t))

(set-force-field *mol* :default)
(:= *e-gaff* (chem:make-energy-function *agg* :use-excluded-atoms t :assign-types t))


