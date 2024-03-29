(load "sys:extensions;cando;src;lisp;start-cando.lisp")
(in-package :cando-user)
*default-pathname-defaults*

(require :asdf)
(asdf:load-asd (merge-pathnames "~/Development/tirun/src/tirun.asd"))
(ql:quickload :tirun)
(leap.core::clear-path)
(add-path "amber:dat;leap;prep;")
(add-path "amber:dat;leap;lib;")
(add-path "amber:dat;leap;parm;")
(add-path "amber:dat;leap;cmd;")
(add-path "amber:dat;antechamber;")
(source "leaprc.protein.ff14SB")
(:= *tiruns* (make-instance 'tirun:calculation))
*tiruns*
(:= *tests* (list (cons :c1 (lambda (a) (eq (chem:get-name a) :c1)))
                  (cons :c3 (lambda (a) (eq (chem:get-name a) :c3)))
                  (cons :c5 (lambda (a) (eq (chem:get-name a) :c5)))))
(:= *pick* (chem:compile-smarts 
             "[C:6]1~[C<c1>:1]~C~[C<c3>:3]~C~[C<c5>:5]~C1" :tests *tests*))
(:= *sk* (handler-bind ((warning #'muffle-warning))
             (with-open-file (fin (open "ligands.cdxml" :direction :input))
             (chem:make-chem-draw fin :add-hydrogens nil))))
(tirun:setup-ligands *tiruns* *sk*)
(:= *lysozyme* (loadPdb "181L_mod.pdb"))
(pushnew *lysozyme* (tirun:receptors *tiruns*))
(directory "*.*")
(load-off (merge-pathnames "phen.lib"))
(load-off (merge-pathnames "benz.lib"))
(:= *ligs* (loadPdb "bnz_phn.pdb"))
*pick*
(alexandria:hash-table-alist (tirun::pattern-atoms *pick* *ligs*))
(clear-force-field)
(load-atom-type-rules "ATOMTYPE_GFF.DEF")
(load-amber-params "gaff.dat")
(tirun:build-ligands *tiruns*)
(show (tirun:layout-ligands *tiruns* :accessor 'tirun::molecule))
(tirun::pose-ligands-using-pattern *tiruns* *pick* *ligs*)
