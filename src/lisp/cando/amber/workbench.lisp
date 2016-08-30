
(setf cl:*features* (cons :debug-require *features*))

(asdf:load-system :cando)

(list-all-packages)
(asdf:load-system "amber")

(with-open-file (fin "~/Development/amber/dat/leap/lib/all_amino94.lib" :direction :input)
  (defparameter *off* (amber.off:read-off-unit-lib fin)))
(loop for v being the hash-value of *off*
   for x from 0
     for trans = (geom:make-m4-translate (geom:vec (* x 5) 0 0))
     do (chem:apply-transform-to-atoms v trans))
(maphash (lambda (k v) (cando:chimera v) (sleep 0.1)) *off*)

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
