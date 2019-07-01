(setup-default-paths)
(load-atom-type-rules "ATOMTYPE_GFF.DEF")
(source "leaprc.ff14SB.redq")
(source "leaprc.gaff")
(defparameter feps (fep:make-fep))
(defparameter sk (load-sketch "ligands.cdxml"))
(fep:setup-ligands feps sk)
(defparameter tests (list (cons :c1 (lambda (a) (eq (chem:get-name a) :c1))) (cons :c3 (lambda (a) (eq (chem:get-name a) :c3))) (cons :c5 (lambda (a) (eq (chem:get-name a) :c5)))))
(defparameter pick (chem:compile-smarts "[C:6]1~[C<c1>:1]~[C:2]~[C<c3>:3]~[C:4]~[C<c5>:5]1" :tests tests))
lysozyme = loadPdb "181L_mod.pdb"
(cando:build-unbuilt-hydrogens |lysozyme|)
(simple-build-unbuilt-atoms |lysozyme|)
(fep:add-receptor |feps| |lysozyme|)
loadOff phen.lib
loadOff benz.lib
ligs = loadPdb bnz_phn.pdb
(simple-build-unbuilt-atoms |ligs|)
(build-unbuilt-hydrogens |ligs|)
(fep:pose-ligands-using-pattern |feps| |pick| |ligs|)
(fep:build-initial-jobs |feps| :connections 2 :stages 3 :windows 11)
(fep:jobs |feps|)
worklist = (fep:generate-jobs |feps|)
(with-open-file (sout "/tmp/graph.dot" :direction :output) (fepdot:draw-graph-stream (list |worklist|) sout))
(ext:system "dot -Tpdf -O /tmp/graph.dot")
(ext:system "open -n /tmp/graph.dot.pdf")
quit
