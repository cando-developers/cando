(in-package :cando-user)
(asdf:load-system "amber")
(use-package :amber :cando-user)
(load "~/Downloads/chromatinfiber-11million/setup.lisp")


(defparameter *1mil-layout* (amber.pdb:scanpdb "~/Downloads/chromatin-1mil/chromsmall.pdb" :progress t))
(defparameter *1mil* (amber.pdb:loadpdb "~/Downloads/chromatin-1mil/chromsmall.pdb" :layout *1mil-layout* :progress t))

(defparameter *1mil-layout* (amber.pdb:scanpdb "~/Downloads/chromatin-1mil/chrom40.pdb" :progress t))
(defparameter *1mil* (amber.pdb:loadpdb "~/Downloads/chromatin-1mil/chrom40.pdb" :layout *1mil-layout* :progress t))




(defparameter *all-layout* (amber.pdb:scanpdb "~/Downloads/chromatinfiber-11million/chromnn_400.pdb" :progress t))
(defparameter *all* (amber.pdb:loadpdb "~/Downloads/chromatinfiber-11million/chromnn_400.pdb" :layout *all-layout* :progress t))



(defparameter *small-layout* (amber.pdb:scanpdb "~/Downloads/4y0y.pdb" :progress t))
*small-layout*

(cando:chimera *small*)
(cando:save-mol2 *s* "~/Downloads/chromatinfiber-11million/chromnn.mol2")
(format t "Done~%")


amber::*objects*
(defparameter *s* (amber.pdb:scanpdb "~/Downloads/4y0y.pdb" :progress t))
(defparameter *a* (amber.pdb::loadpdb "~/Downloads/4y0y.pdb"))
(defparameter *s* 1234)
*s*
*package*

cando-user::*s*

amber.pdb::*bar*

(apropos "make-progress-bar")

;; wipe out hydrogens
(let ((hydrogens (chem:all-atoms-of-element-as-list *a* :h)))
  (mapc (lambda (h) (chem:setf-needs-build h t)
                (chem:set-position h (geom:vec 0.0 0.0 0.0)))
        hydrogens)
  (format t "There are ~a hydrogens~%" (length hydrogens)))
(cando::build-unbuilt-hydrogens *a*)
(cando:chimera *a*)

(/ 1 0.0174533)57.295757



amber.pdb::*serial-to-atoms*
(print "Done")
*a*
(defparameter *ab* nil)
(dolist (a (chem:all-atoms-of-element-as-list *a* :h))
  (when (chem:needs-build a)
    (setq *ab* a)
    (format t "Atom that needs build: ~a~%" a)))

(defparameter *ht* (cando::group-unbuilt-hydrogens *a*))
(maphash (lambda (k v) (print (list k v))) *ht*)
(defparameter *sl* (chem:make-spanning-loop *ab*))
(chem:advance-loop-and-process *sl*)
(chem:get-atom *sl*)
(apropos "spanning")


*seqs*
(cando:chimera *a*)
amber:*topologies*



(defparameter *a* '(1 2 3 4 5))
(setf (car *a*) 10)
(setf (car (last *a*)) 99)
*a*

(with-input-from-string (sin "A B C D        ")
  (loop for res = (read sin nil :eof)
     unless (eq res :eof)
     collect res))




amber.off::*res*
amber.off::*in-plug-idx*
*t*
(progn
  (defparameter *a* (chem:make-aggregate))
  (defparameter *m* (chem:make-molecule))
  (defparameter *r* (chem:build-residue *t*))
  (chem:add-matter *a* *m*)
  (chem:add-matter *m* *r*))

(chem:build-residue *t*)
*r*

*t*
(cando:as-string (chem:get-constitution *t*))"


(let ((*print-readably* t))
  (core:print-cxx-object *m* t))

(cando:chimera *a*)
(chem:atom-info *t*)
  (defparameter *a* (gethash :ALA *off*))
  (defparameter *r* (chem:content-at (chem:content-at *a* 0) 0)))

(core:encode *r*)


(cando::make-topology-from-residue *r*)



(defparameter *pr* (make-instance 'pdb-reader :stream (open "~/Downloads/chromnn_400.cpptraj.pdb" :direction :input)))

(defparameter *l* (readline (stream *pr*)))
(parse-line "ATOM  80059  N   ALA    14     -53.432 -44.1623941.690  1.00  0.00           N  " *pr*)(:ATOM 80059 :N :ALA 14 -53.432 -44.162 3941.69)
(read-from-string (subseq "1 2 3 4" 2 3))
(read-from-string "1 2 3 4" :start 2 :end 3)




(chem:make-topology-from-residue *r*)
(loop for v being the hash-value of *off*
   for x from 0
     for trans = (geom:make-m4-translate (geom:vec (* x 5) 0 0))
     do (chem:apply-transform-to-atoms v trans))
(maphash (lambda (k v) (cando:chimera v) (sleep 0.1)) *off*)
(defparameter *agg* (gethash :TRP *off*))
(defparameter *mol* (chem:content-at *agg* 0))
(defparameter *res* (chem:content-at *mol* 0))
(chem:atom-with-id (cando:constitution-from-residue *res*) 1)

(chem:c

(apropos "chimera")
(in-package :amber.off)
(defparameter *agg* (assemble-hierarchy *read-hierarchy* *residues* *atoms*))

(chem:contents-as-cons *agg*)
(chem:molecules *agg*)


(loop for e across #(10 20 30 40)
     for i from 0
     do (print (list i e)))


(chem:atomic



(defvar *data*)
(with-open-file (fin "~/Development/amber/dat/leap/lib/solvents.lib" :direction :input)
  (setf *data* (amber.off:read-off-data-block fin nil nil)))


(setq *data*
      (with-open-file (fin "~/Development/amber/dat/leap/lib/all_amino03.lib" :direction :input)
        (let ((names (amber.off::read-off-lib-header fin)))
          (amber.off::read-off-unit-parts fin))))

(gethash "atoms" (gethash "PHE" *data*))


(make-array 
