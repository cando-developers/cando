;;;Build C6H12N2SO2

(defparameter *mol1* (chem:make-molecule))
(defparameter *mol2* (chem:make-molecule))
(defparameter *res1* (chem:make-residue :foo))
(defparameter *res2* (chem:make-residue :bar))

(defparameter *a1* (chem:make-atom :CL :|Cl|))
(defparameter *a2* (chem:make-atom :CL :|Cl|))

(chem:set-position *a1* '(0.0 0.0 0.0))
(chem:set-position *a2* '(3.0 0.0 0.0))
(chem:add-matter *res1* *a1*)
(chem:add-matter *res1* *a2*)
(chem:add-matter *mol1* *res1*)
(chem:add-matter *mol2* *res2*)

(defparameter *agg* (chem:make-aggregate))
(chem:add-matter *agg* *mol1*)
(chem:add-matter *agg* *mol2*)                                                                         
(chem:set-property *agg* :solvent-box '(20.0 30.0 40.0))
(leap:setup-amber)

energy::*ff*

;(cando:jostle *mol1* 80)
;(chem:all-atoms *agg*)

;(defparameter *a1-pos* (chem:get-position *a1*))
;(defparameter *a2-pos* (chem:get-position *a2*))
(defparameter *me* (chem:make-energy-function *agg* energy::*ff* :use-excluded-atoms t))

energy::*ff*

(leap.topology::save-amber-parm-format *agg* "clcl_cando.top" "clcl_cando.crd" energy::*ff*)


(unwind-protect
     (progn
       (chem:enable-debug (chem:get-nonbond-component *me*))
       (chem:calculate-energy *me*))
  (chem:disable-debug (chem:get-nonbond-component *me*)))



(setf *default-pathname-defaults* #P"/Users/tuj24515/Development/clasp/extensions/cando/src/lisp/tests/tiny/")
(progn
  (defparameter *fif* (fortran:make-fortran-input-file :stream (open "clcl_cando.top" :direction :input)))
  (defparameter *ee* (leap.topology::read-amber-parm-format *fif*)))

(progn
  (defparameter *topologyfif* (fortran:make-fortran-input-file :stream (open "clcl_cando.crd" :direction :input)))
  (defparameter *coords* (leap.topology::read-amber-coordinate-file *topologyfif*)))


(chem:evaluate-energy *ee* *coords*)
