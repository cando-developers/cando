(setf *default-pathname-defaults* #P"~/Development/cando/extensions/cando/src/tests/alanine-demo/")
(setup-default-paths)
(easy-gaff)
(clear-force-field :smirnoff)
(load-smirnoff-params #P"~/Development/openforcefield/openforcefield/data/forcefield/smirnoff99Frosst.offxml")
;;(:= *agg* (load-chem-draw-aggregate "alanine2.cdxml"))
(:= *agg* (load-chem-draw-aggregate "/Users/shiho/Downloads/test-sdf.cdxml"))
(chem:fill-in-implicit-hydrogens *agg*)
(set-stereoisomer-func (cando:gather-stereocenters *agg*) (lambda (c) :R))
(:= *mol* (chem:content-at *agg* 0))
(:= *ee* (energy:minimize *agg*))


(set-force-field *mol* :smirnoff)
(save-amber-parm *agg* "/tmp/test.parm" "/tmp/test.crd")



(set-force-field *mol* :smirnoff)
(:= *e-smirnoff* (chem:make-energy-function *agg* :use-excluded-atoms t :assign-types t))

(set-force-field *mol* :default)
(save-amber-parm *agg* "/tmp/testdefault.parm" "/tmp/testdefault.crd")
(:= *e-gaff* (chem:make-energy-function *agg* :use-excluded-atoms t :assign-types t))


leap.core:*force-fields*

(:= *a* (make-array (chem:get-nvector-size *e-smirnoff*) :element-type 'double-float))
(chem:load-coordinates-into-vector *e-smirnoff* *a*)
(:= *e* (chem:evaluate-energy *e-smirnoff* *a*))

*a*

(chem:save-pdb *agg* "/Users/shiho/Downloads/test-sdf.pdb")
(:= *agg* (load-chem-draw-aggregate "/Users/shiho/Downloads/ligand8a.cdxml"))
