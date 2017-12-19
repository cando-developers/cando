;;;Build ACE-ALA-NME

;;Load alanine from chemdraw
(ext:chdir "~/Development/clasp/extensions/cando/src/lisp/tests/alanine-demo/" t)
(defparameter *chemagg* (cando:load-chemdraw-aggregate "alanine2.cdxml"))
(defparameter *stereocenters* (cando:gather-stereocenters *chemagg*))
(cando:set-stereoisomer-func *stereocenters* (lambda (c) :S))
(chem:set-property *chemagg* :solvent-box '(50.0 50.0 50.0))
(cando:save-mol2 *chemagg* "alanine2.mol2")
(energy::setup-amber)
(cando:jostle *chemagg*)
(energy:minimize *chemagg*)
(leap.topology::save-amber-parm-format *chemagg* "alanine_cando.top" "alanine_cando.crd" energy::*ff*)
(ext:system "~/Development/amber16/bin/sander -O -i 01_Min.in -o 01_Min_cando.out -p alanine_cando.top -c alanine_cando.crd")

(defun read-sander-output (filename)
  (with-open-file (fin filename)
    (let ((all (make-string (file-length fin))))
      (read-sequence all fin)
      (with-input-from-string (sin all)
	(loop for line = (read-line  sin nil :eof)
	   until (eq line :eof)
	   for line-len = (if (eq line :eof) 0 (length line))
	   do (progn
;;;		(format t "line: ~a~%" line)
		(cond
		  ((and (> line-len #.(length " BOND"))
			(string= line " BOND" :start1 0 :end1 #.(length " BOND")))
		   (format t "~a~%" line))
		  ((and (> line-len #.(length " VDWAALS"))
			(string= line " VDWAALS" :start1 0 :end1 #.(length " VDWAALS")))
		   (format t "~a~%" line))
		  ((and (> line-len #.(length " 1-4 VDW"))
			(string= line " 1-4 VDW" :start1 0 :end1 #.(length " 1-4 VDW")))
		   (format t "~a~%" line)))))))))

(read-sander-output "01_Min_cando.out")
			   

;;; Use (ext::start-swank)    to start
;;; Then in emacs:     slime-connect
(progn
  (defparameter *started-swank* nil)
  (defun ext::start-swank ()
    ;; Bad!  This is hard-coded to work with docker
    (if *started-swank*
        (format t "Swank is already running~%")
        (progn
          (let ((swank-loader (probe-file "~/Development/slime/swank-loader.lisp")))
            (format t "swank-loader -> ~a~%" swank-loader)
            (load swank-loader))
          (let ((swank-loader-init (find-symbol "INIT" "SWANK-LOADER")))
            (funcall swank-loader-init :delete nil :reload nil :load-contribs nil))
          (let ((swank-create-server (find-symbol "CREATE-SERVER" "SWANK")))
            (mp:process-run-function 'swank-main
                                     (lambda () (funcall swank-create-server
                                                         :port 4005
                                                         :interface "0.0.0.0")))
            (setf *started-swank* t))))))


   

;;; Build it like leap
#|
(add-path "/Users/tuj24515/Development/amber16/dat/antechamber/")
(add-path "/Users/tuj24515/Development/amber16/dat/leap/prep/")
(add-path "/Users/tuj24515/Development/amber16/dat/leap/lib/")
(add-path "/Users/tuj24515/Development/amber16/dat/leap/parm/")
(add-path "/Users/tuj24515/Development/amber16/dat/leap/cmd/")

(leap.commands::load-amber-type-rules "ATOMTYPE_AMBER.DEF")






(ext:chdir "~/Development/clasp/extensions/cando/src/lisp/tests/alanine/" t)
(defparameter *agg* (cando:load-mol2 "alanine.mol2"))

*agg*

|#
(defparameter *met* (chem:make-molecule))
(defparameter *ace* (chem:make-residue :ACE))
(defparameter *ala* (chem:make-residue :ALA))
(defparameter *nme* (chem:make-residue :NME))

(chem:add-matter *met* *ace*)
(chem:add-matter *met* *ala*)
(chem:add-matter *met* *nme*)

(defparameter *h1* (chem:make-atom :HH31 :H))
(defparameter *c1* (chem:make-atom :CH3 :C))
(defparameter *h2* (chem:make-atom :HH32 :H))
(defparameter *h3* (chem:make-atom :HH33 :H))
(defparameter *c2* (chem:make-atom :C :C))
(defparameter *o1* (chem:make-atom :O :O))
(defparameter *n1* (chem:make-atom :N :N))
(defparameter *h4* (chem:make-atom :H :H))
(defparameter *c3* (chem:make-atom :CA :C))
(defparameter *h5* (chem:make-atom :HA :H))
(defparameter *c4* (chem:make-atom :CB :C))
(defparameter *h6* (chem:make-atom :HB1 :H))
(defparameter *h7* (chem:make-atom :HB2 :H))
(defparameter *h8* (chem:make-atom :HB3 :H))
(defparameter *c5* (chem:make-atom :C :C))
(defparameter *o2* (chem:make-atom :O :O))
(defparameter *n2* (chem:make-atom :N :N))
(defparameter *h9* (chem:make-atom :H :H))
(defparameter *c6* (chem:make-atom :CH3 :C))
(defparameter *h10* (chem:make-atom :HH31 :H))
(defparameter *h11* (chem:make-atom :HH32 :H))
(defparameter *h12* (chem:make-atom :HH33 :H))

(chem:add-matter *ace* *h1*)
(chem:add-matter *ace* *c1*)
(chem:add-matter *ace* *h2*)
(chem:add-matter *ace* *h3*)
(chem:add-matter *ace* *c2*)
(chem:add-matter *ace* *o1*)
(chem:add-matter *ala* *n1*)
(chem:add-matter *ala* *h4*)
(chem:add-matter *ala* *c3*)
(chem:add-matter *ala* *h5*)
(chem:add-matter *ala* *c4*)
(chem:add-matter *ala* *h6*)
(chem:add-matter *ala* *h7*)
(chem:add-matter *ala* *h8*)
(chem:add-matter *ala* *c5*)
(chem:add-matter *ala* *o2*)
(chem:add-matter *nme* *n2*)
(chem:add-matter *nme* *h9*)
(chem:add-matter *nme* *c6*)
(chem:add-matter *nme* *h10*)
(chem:add-matter *nme* *h11*)
(chem:add-matter *nme* *h12*)

(chem:bond-to *h1* *c1* :single-bond)
(chem:bond-to *h2* *c1* :single-bond)
(chem:bond-to *h3* *c1* :single-bond)
(chem:bond-to *c1* *c2* :single-bond)
(chem:bond-to *c2* *o1* :double-bond)
(chem:bond-to *c2* *n1* :single-bond)
(chem:bond-to *n1* *h4* :single-bond)
(chem:bond-to *n1* *c3* :single-bond)
(chem:bond-to *c3* *h5* :single-bond)
(chem:bond-to *c3* *c4* :single-bond)
(chem:bond-to *c4* *h6* :single-bond)
(chem:bond-to *c4* *h7* :single-bond)
(chem:bond-to *c4* *h8* :single-bond)
(chem:bond-to *c3* *c5* :single-bond)
(chem:bond-to *c5* *o2* :double-bond)
(chem:bond-to *c5* *n2* :single-bond)
(chem:bond-to *n2* *h9* :single-bond)
(chem:bond-to *n2* *c6* :single-bond)
(chem:bond-to *c6* *h10* :single-bond)
(chem:bond-to *c6* *h11* :single-bond)
(chem:bond-to *c6* *h12* :single-bond)
(defparameter *agg* (chem:make-aggregate))
(chem:add-matter *agg* *met*)
(chem:set-property *agg* :solvent-box '(20.0 30.0 40.0))

(energy::setup-amber)

(apropos "*ff*")

(cando:jostle *met* 20)
(defparameter *me* (energy:minimize *agg*))
(setf *default-pathname-defaults* #P"/Users/tuj24515/Development/clasp/extensions/cando/src/lisp/tests/alanine/")
(leap.topology::save-amber-parm-format *agg* "alanine_cando.top" "alanine_cando.crd" energy::*ff*)


;(defparameter *leapagg* (cando:load-mol2 "/Users/tuj24515/Development/clasp/extensions/cando/src/lisp/tests/alanine/alanine.mol2"))
;*leapagg*

(unwind-protect
     (progn
       (chem:enable-debug (chem:get-nonbond-component *me*))
       (chem:calculate-energy *me*))
  (chem:disable-debug (chem:get-nonbond-component *me*)));


(setf *default-pathname-defaults* #P"/Users/tuj24515/Development/clasp/extensions/cando/src/lisp/tests/alanine/")
(progn
  (defparameter *fif* (fortran:make-fortran-input-file :stream (open "alanine_cando.top" :direction :input)))
  (defparameter *ee* (leap.topology::read-amber-parm-format *fif*)))

(progn
  (defparameter *topologyfif* (fortran:make-fortran-input-file :stream (open "alanine_cando.crd" :direction :input)))
  (defparameter *coords* (leap.topology::read-amber-coordinate-file *topologyfif*)))


(chem:evaluate-energy *ee* *coords*)

(chem:evaluate-energy *me* *coords*)

(chem:get-energy (chem:get-nonbond-component *ee*))

(chem:get-energy (chem:get-stretch-component *ee*))
;(chem:get-energy (chem:get-angle-component *me*))
;(chem:get-energy (chem:get-dihedral-component *me*))
;(chem:get-energy (chem:get-nonbond-component *me*))


;(defparameter *chemagg* (cando:load-chemdraw-aggregate "untitled.cdxml"))
;(defparameter *stereocenters* (cando:gather-stereocenters *chemagg*))
;(cando:set-stereoisomer-func *stereocenters* (lambda (c) :S))

;(leap.topology::save-amber-parm-format *chemagg* "alanine_chem_cando.top" "alanine_chem_cando.crd" energy::*ff*)
|#
