(setup-default-paths)
(load-atom-type-rules "ATOMTYPE_GFF.DEF")
(source "leaprc.ff14SB.redq" :verbose t)
(source "leaprc.gaff")
(defparameter tiruns (tirun:make-tirun))
(defparameter sk (load-chem-draw "ligands.cdxml"))
(multiple-value-bind (molecules ligands)
    (tirun:assemble-ligands sk)
  (setf (tirun:ligands tiruns) ligands))
(defparameter tests (list (cons :c1 (lambda (a) (eq (chem:get-name a) :c1))) (cons :c3 (lambda (a) (eq (chem:get-name a) :c3))) (cons :c5 (lambda (a) (eq (chem:get-name a) :c5)))))
(defparameter pick (chem:compile-smarts "[C:6]1~[C<c1>:1]~[C:2]~[C<c3>:3]~[C:4]~[C<c5>:5]1" :tests tests))
(defparameter lysozyme (loadPdb "181L_mod.pdb"))
(cando:build-unbuilt-hydrogens lysozyme)
(simple-build-unbuilt-atoms lysozyme)
(format t "(tirun:add-receptor tiruns lysozyme)~%")
(tirun:add-receptor tiruns lysozyme)
(format t "(load-off \"phen.lib\")~%")
(load-off "phen.lib")
(format t "(load-off \"benz.lib\")~%")
(load-off "benz.lib")
(format t "(defparameter ligs (loadPdb \"bnz_phn.pdb\"))~%")
(defparameter ligs (loadPdb "bnz_phn.pdb"))
(format t "(simple-build-unbuilt-atoms ligs)~%")
(simple-build-unbuilt-atoms ligs)
(format t "(build-unbuilt-hydrogens ligs)~%")
(build-unbuilt-hydrogens ligs)
(format t "(tirun:pose-ligands-using-pattern tiruns pick ligs)~%")
(tirun:pose-ligands-using-pattern tiruns pick ligs)
(format t "(tirun:build-job-nodes tiruns)~%")
(tirun:build-job-nodes tiruns)
(format t "(tirun:connect-job-nodes tiruns :simple :connections 2 :stages 3 :windows 11)~%")
(tirun:connect-job-nodes tiruns :simple :connections 2 :stages 3 :windows 11)
(format t "(defparameter worklist (tirun:generate-jobs tiruns))~%")
(trace tirun:generate-jobs tirun:setup-am1-calculations tirun:save-tiruns)
(defparameter worklist (tirun:generate-jobs tiruns))
#+(or)
(progn
  (with-open-file (sout "/tmp/graph.dot" :direction :output) (tirundot:draw-graph-stream (list worklist) sout))
  (ext:system "dot -Tpdf -O /tmp/graph.dot")
  (ext:system "open -n /tmp/graph.dot.pdf"))
(ext:quit)
