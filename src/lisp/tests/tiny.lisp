;;;Build C6H12N2SO2

(defparameter *mol1* (chem:make-molecule))
(defparameter *mol2* (chem:make-molecule))
(defparameter *res1* (chem:make-residue :foo))
(defparameter *res2* (chem:make-residue :bar))

(defparameter *a1* (chem:make-atom :NA :|Na|))
(defparameter *a2* (chem:make-atom :CL :|Cl|))

(chem:add-matter *res1* *a1*)
(chem:add-matter *res1* *a2*)
(chem:add-matter *mol1* *res1*)
(chem:add-matter *mol2* *res2*)

(defparameter *agg* (chem:make-aggregate))
(chem:add-matter *agg* *mol1*)
(chem:add-matter *agg* *mol2*)

(energy:setup-amber)
;
;(cando:jostle *mol1* 80)
;(chem:all-atoms *agg*)

;(defparameter *a1-pos* (chem:get-position *a1*))
;(defparameter *a2-pos* (chem:get-position *a2*))
;(defparameter *me* (chem:make-energy-function *agg* energy::*ff* :use-excluded-atoms t))
energy::*ff*

;(leap.topology::save-amber-parm-format *agg* "/tmp/top.dat" "/tmp/crd.dat" energy::*ff*)
;(chem:calculate-energy *me*)

