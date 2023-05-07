startSwank

(load-atom-type-rules "ATOMTYPE_GFF.DEF")

(load-amber-params "gaff.dat")

sm = smiles "O=C1[C@H](C)/C=C(CN(C([C@H](C)/C=C(CN(C([C@H](C)/C=C(CN(C([C@H](C)/C=C(CN1CC2=CC=CC=C2)\\C)=O)CC3=CC=CC=C3)\\C)=O)CC4=CC=CC(NC(NC)=O)=C4)\\C)=O)CC5=CC=NC=C5)\\C"

show (sketch2d:sketch2d sm)

stereoisomers = (loop for config from 0 to #b0001 for newmol = (chem:matter-copy sm) for chirals = (gather-stereocenters newmol) do (chem:set-name newmol (intern (format nil "MOL~a" config) :keyword)) do (cando:set-stereoisomer chirals config :show t) collect newmol)

(lparallel:pmap nil (lambda (agg) (build-good-geometry-from-random agg) (chem:calculate-stereochemistry-from-structure agg)) stereoisomers)

(map nil (lambda (si) (format t "Stereoisomer: ~a~%" si) (let ((chirals (gather-stereocenters si))) (dump-stereocenters chirals))) stereoisomers)

show (elt stereoisomers 0)

source "leaprc.water.tip3p"

(map nil (lambda (agg) (leap:resolvate-box agg CHCL3BOX 8.0)) stereoisomers)

show (elt stereoisomers 1)

sim = (amber:simulation)

(amber:add-aggregates sim stereoisomers)

(amber:minimize sim)

(amber:heat sim :temperature 300.0)

(amber:pressurize sim)

(amber:dynamics sim :time-ps 100)

(ql:quickload :amber-jupyter)

stream = (amber:start-simulation sim :name "data/ambertools/cyclo-many/")

(amber.jupyter::monitor sim)


