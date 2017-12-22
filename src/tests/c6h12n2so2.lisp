;;;Build C6H12N2SO2

(defparameter *met* (chem:make-molecule))
(defparameter *res* (chem:make-residue :foo))
(defparameter *res1* (chem:make-residue :foo1))

(chem:add-matter *met* *res*)
(chem:add-matter *met* *res1*)

(defparameter *c1* (chem:make-atom :C1 :C))
(defparameter *c2* (chem:make-atom :C2 :C))
(defparameter *c3* (chem:make-atom :C3 :C))
(defparameter *c4* (chem:make-atom :C4 :C))
(defparameter *c5* (chem:make-atom :C5 :C))
(defparameter *c6* (chem:make-atom :C6 :C))
(defparameter *o1* (chem:make-atom :O1 :O))
(defparameter *o2* (chem:make-atom :O2 :O))
(defparameter *o3* (chem:make-atom :O3 :O))
(defparameter *n1* (chem:make-atom :N1 :N))
(defparameter *n2* (chem:make-atom :N2 :N))
(defparameter *h1* (chem:make-atom :H1 :h))
(defparameter *h2* (chem:make-atom :H2 :h))
(defparameter *h3* (chem:make-atom :H3 :h))
(defparameter *h4* (chem:make-atom :H4 :h))
(defparameter *h5* (chem:make-atom :H5 :h))
(defparameter *h6* (chem:make-atom :H6 :h))
(defparameter *h7* (chem:make-atom :H7 :h))
(defparameter *h8* (chem:make-atom :H8 :h))
(defparameter *h9* (chem:make-atom :H9 :h))
(defparameter *h10* (chem:make-atom :H10 :h))
(defparameter *h11* (chem:make-atom :H11 :h))
(defparameter *h12* (chem:make-atom :H12 :h))
(defparameter *s1* (chem:make-atom :S1 :s))

(chem:add-matter *res* *c1*)
(chem:add-matter *res* *c2*)
(chem:add-matter *res* *c3*)
(chem:add-matter *res1* *c4*)
(chem:add-matter *res1* *c5*)
(chem:add-matter *res1* *c6*)
(chem:add-matter *res* *o1*)
(chem:add-matter *res* *o2*)
(chem:add-matter *res1* *o3*)
(chem:add-matter *res* *n1*)
(chem:add-matter *res1* *n2*)
(chem:add-matter *res* *h1*)
(chem:add-matter *res* *h2*)
(chem:add-matter *res* *h3*)
(chem:add-matter *res* *h4*)
(chem:add-matter *res* *h5*)
(chem:add-matter *res* *h6*)
(chem:add-matter *res1* *h7*)
(chem:add-matter *res1* *h8*)
(chem:add-matter *res1* *h9*)
(chem:add-matter *res1* *h10*)
(chem:add-matter *res1* *h11*)
(chem:add-matter *res1* *h12*)
(chem:add-matter *res* *s1*)

(chem:bond-to *h1* *o1* :single-bond)
(chem:bond-to *o1* *c1* :single-bond)
(chem:bond-to *c1* *o2* :double-bond)
(chem:bond-to *c1* *c2* :single-bond)
(chem:bond-to *c2* *h3* :single-bond)
(chem:bond-to *c2* *c3* :single-bond)
(chem:bond-to *c3* *h4* :single-bond)
(chem:bond-to *c3* *h5* :single-bond)
(chem:bond-to *c3* *s1* :single-bond)
(chem:bond-to *s1* *h2* :single-bond)
(chem:bond-to *c2* *n1* :single-bond)
(chem:bond-to *n1* *h6* :single-bond)
(chem:bond-to *n1* *c4* :single-bond)
(chem:bond-to *c4* *o3* :double-bond)
(chem:bond-to *c4* *c5* :single-bond)
(chem:bond-to *c5* *h7* :single-bond)
(chem:bond-to *c5* *c6* :single-bond)
(chem:bond-to *c6* *h8* :single-bond)
(chem:bond-to *c6* *h9* :single-bond)
(chem:bond-to *c6* *h10* :single-bond)
(chem:bond-to *c5* *n2* :single-bond)
(chem:bond-to *n2* *h11* :single-bond)
(chem:bond-to *n2* *h12* :single-bond)
(defparameter *agg* (chem:make-aggregate))
(chem:add-matter *agg* *met*)
(chem:set-property *agg* :solvent-box '(20.0 30.0 40.0))
*agg*

(leap:setup-amber)
(time (defparameter *glycam* (leap:load-amber-params (probe-file #P"~/Development/amber16/dat/leap/parm/GLYCAM_06EPb.dat"))))
*glycam*

(probe-file #P"~/Development/amber16/dat/leap/parm/GLYCAM_06EPb.dat")
;
(cando:jostle *met* 80)
(defparameter *me* (energy:minimize *agg*))

(leap.topology::save-amber-parm-format *agg* "/tmp/top.dat" "/tmp/crd.dat" *ff*)
(chem:calculate-energy *me*)


;(defun foo ()
;  (Chem:make-energy-function *agg* energy::*ff* :use-excluded-atoms t))

;(foo)
					
;(defparameter *e* (Chem:make-energy-function *agg* energy::*ff* :use-excluded-atoms t))
;(print (chem:calculate-energy *e*))




;(leap.topology::prepare-amber-energy-stretch *e*)

;(defun foo (x) (bar x))
;(defun bar (x) (baz x))
;(defun baz (x) (error "In baz found error"))

;(foo 1)
;(defparameter *e* (Chem:make-energy-function *agg* energy::*ff* :use-excluded-atoms nil))

;(defparameter *c1-pos* (chem:get-position *c1*))

;(leap.topology::prepare-amber-energy-stretch *e*)

;(leap.topology::prepare-amber-energy-angle *e*)

;(chem:is-bonded-to *h5* *c2*)

;;(leap.topology::prepare-residue *e*)
;(leap.topology::save-amber-parm-format *agg* "/Users/tuj24515/Development/tests/c6h12n2so3.dat" "Users/tuj24515/c6h12v2so3-coordinate.dat" energy::*ff*)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;q
;;; Test the reader
;;;

#|
(setf *default-pathname-defaults* #P"/tmp/")
(progn
  (defparameter *fif* (fortran:make-fortran-input-file :stream (open "top.dat" :direction :input)))
  (defparameter *ee* (leap.topology::read-amber-parm-format *fif*)))

(progn
  (defparameter *topologyfif* (fortran:make-fortran-input-file :stream (open "crd.dat" :direction :input)))
  (defparameter *coords* (leap.topology::read-amber-coordinate-file *topologyfif*)))


(chem:evaluate-energy *ee* *coords*)2415.141413817592d0

(chem:evaluate-energy *me* *coords*)2416.849376233198d0



(chem:calculate-energy *me*) 

(chem:get-energy (chem:get-nonbond-component *me*))
(chem:get-energy (chem:get-nonbond-component *ee*))

(chem:dump-terms (chem:get-stretch-component *ee*))

(unwind-protect
     (progn
       (chem:enable-debug (chem:get-nonbond-component *ee*))
       (chem:calculate-energy *ee*))
  (chem:disable-debug (chem:get-nonbond-component *ee*)))




(setf *default-pathname-defaults* #P"/Users/tuj24515/Development/tutorial/B0/")
(progn
  (defparameter *fif* (fortran:make-fortran-input-file :stream (open "prmtop" :direction :input)))
  (defparameter *ae* (leap.topology::read-amber-parm-format *fif*)))

*ae*

(progn
  (defparameter *crdfif* (fortran:make-fortran-input-file :stream (open "01_Min.rst" :direction :input)))
  (defparameter *coords* (leap.topology::read-amber-coordinate-file *crdfif*)))

(chem:evaluate-energy *ae* *coords*)

(chem:get-energy (chem:get-stretch-component *ae*))
(chem:get-energy (chem:get-angle-component *ae*))
(chem:get-energy (chem:get-dihedral-component *ae*))
(chem:get-energy (chem:get-nonbond-component *ae*))



|#
