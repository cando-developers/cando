

(ql:quickload :smirnoff)

(:= *sf* (smirnoff:load-smirnoff #P"~/Development/openforcefield/openforcefield/data/forcefield/smirnoff99Frosst.offxml"))

(setf *default-pathname-defaults* #P"~/Development/cando/extensions/cando/src/lisp/cando/smirnoff/")

(:= *agg* (load-chem-draw-aggregate "rings.cdxml"))

(smirnoff:build-energy-function *agg* *sf*)
