(progn
  (require :asdf)
  (let ((central-registry (find-symbol "*CENTRAL-REGISTRY*" :asdf))
        (load-system (find-symbol "LOAD-SYSTEM" :asdf)))
    (setf (symbol-value central-registry)
          (cons (translate-logical-pathname "cando:cando;")
                (symbol-value central-registry)))
    (funcall load-system "cando"))
  (format t "Done initialization pid = ~a~%" (getpid))
  (setq *print-array* t)
  (setq *default-pathname-defaults*
        #P"~/Development/Clasp/projects/cando/work/molecules/")
  (defmethod cleavir-environment::symbol-macro-expansion (sym (env core:value-environment))
    (cleavir-environment::symbol-macro-expansion sym nil))
  (defmethod cleavir-environment::macro-function (sym (env core:value-environment))
    (cleavir-environment::macro-function sym nil))
  (format t "Done with initialization~%"))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Work with a molecule
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(in-package :cando)
;;;
;;; Load Antechamber type rules and AMBER parameters
;;;
(let ((*default-pathname-defaults*
       #P"~/Development/clasp/projects/cando/work/forceField/"))
  (defparameter *parms*
    (let ((parms (chem:make-read-amber-parameters)))
      (with-open-file (fin "ATOMTYPE_GFF.DEF" :direction :input)
        (chem:read-types parms fin))
      (with-open-file (fin "gaff.dat" :direction :input)
        (chem:read-parameters parms fin)
        parms)))
  (defparameter *ff* (chem:get-force-field *parms*)))


;;;
;;; Load the chemdraw file containing a molecule
;;;
(let ((*default-pathname-defaults* #P"~/Development/Clasp/projects/cando/work/molecules/"))
  (defparameter *cd*
    (with-open-file
        (fin "chanmon0.cdxml" :direction :input)
      (chem:make-chem-draw fin)))
  (defparameter *agg* (chem:as-aggregate *cd*)))

(chem:map-residues nil
                   (lambda (a)
                     (when (chem:has-property a :constitution)
                       (format t "Residue: ~a  constitution: ~a~%" (chem:get-name a) (chem:matter-get-property a :constitution))))
                *agg*)
(chem:get-properties (chem:content-at *agg* 1))

(chem:map-atoms nil (lambda (a) (print a)) (chem:content-at *agg* 1))
(print (chem:content-at (chem:content-at *agg* 1) 0))

(progn
  (format t "Residue: ~a~%" (chem:get-name (chem:content-at (chem:content-at *agg* 1) 0)))
  (chem:map-atoms nil (lambda (a) (print a)) (chem:content-at (chem:content-at *agg* 1) 0)))

(apropos "atoms")

;;;
;;; Gather the atoms that are stereocenters into a vector
(progn
  (defparameter *stereocenters*
    (gather-stereocenters *agg*))
  (sort *stereocenters* #'string< :key #'chem:get-name)
  (format t "~a~%" *stereocenters*))

;;; Set the stereochemistry using a binary pattern
(set-stereoisomer *stereocenters* #b11111111111111 :show t)


;;; Build starting geometry for one stereoisomer
(build-good-geometry-from-random *agg* *ff*)
(save-mol2 *agg* #P"test.mol2")
(chimera *agg*)




;;;
;;; Build all 32 stereoisomers
;;;
(dotimes (i 32)
  (set-stereoisomer *stereocenters* i :show t)
  (build-good-geometry-from-random *agg* *ff*)
  (save-mol2 *agg*
             (indexed-pathname #P"fuzi/run.mol2" i)))

;;;
;;; Load all of the 32 stereoisomers
;;;


(defparameter *all-aggs*
  (loop for i below 32
     collect (let ((x (chem:load-mol2 (indexed-pathname #P"fuzi/run.mol2" i))))
               (chem:set-property x :stereoisomer-index i)
               x)))


(defparameter *all-aggs*
  (loop for i to 31
     for x = (chem:load-mol2 (indexed-pathname #P"fuzi/run.mol2" i))
     do (chem:set-property x :stereoisomer-index i)
     collect x))

(chimera (find-aggregate *all-aggs* (lambda (a) (= 28 (chem:matter-get-property a :stereoisomer-index)))))
*all-aggs*
(loop for x in *all-aggs* do (chimera x))

(dolist (x *all-aggs*)
  (chimera x))

(defparameter *3o3x* (chem:load-mol2 "fuzi/3o3x.mol2"))
(defparameter *3o3x-atomspec* "#0:43.A@CA #0:43.A@CB #0:46.A@CA #0:46.A@CB #0:47.A@CA #0:47.A@CB #0:50.A@CA #0:50.A@CB")

(chimera *3o3x*)
(chimera (car *all-aggs*))


(coordinate-array-for-atomspec *3o3x* *3o3x-atomspec*)

(defparameter *agg-atomspec*
  "#0:1.A@N.z #0:1.A@C.z2
#0:1.A@C.a #0:1.A@C.a2
#0:1.A@N.b #0:1.A@C.b2
#0:1.A@N.d #0:1.A@C.d2")

(chem:atoms-with-chimera-specifications (car *all-aggs*) *agg-atomspec*)
(coordinate-array-for-atomspec (car *all-aggs*) *agg-atomspec*)
(dump-atoms (car *all-aggs*))

(defparameter *res*
  (superpose-against-template :fixed-aggregate *3o3x*
                              :fixed-atomspec "#0:43.A@CA #0:43.A@CB #0:46.A@CA #0:46.A@CB #0:47.A@CA #0:47.A@CB #0:50.A@CA #0:50.A@CB"
                              :moveable-aggregates *all-aggs*
                              :moveable-atomspec "#0:1.A@N.z #0:1.A@C.z2    #0:1.A@C.a #0:1.A@C.a2    #0:1.A@N.b #0:1.A@C.b2    #0:1.A@N.d #0:1.A@C.d2"
                              :key (lambda (x) (chem:matter-get-property x :stereoisomer-index))))


(defparameter *data* (sort (let (res) (maphash (lambda (k v) (push (list k v) res)) *res*) res) (lambda (x y) (< (car x) (car y)))))
(gnuplot-data *data*)

(defparameter *a* (find-aggregate *all-aggs* (lambda (a) (= 6 (chem:matter-get-property a :stereoisomer-index)))))
(chimera *a*)


(format t "~b~%" 6)

(ext:system "gnuplot -e 'set terminal png; set output \"/tmp/out.png\"; plot sin(x); quit'")
(ext:system "open /tmp/out.png")
(print "Hello")
(gnuplot-data '(1 2 3 4 3 2 1))
(gnuplot-xy '((1.0 1.0) (2.0 2.0)))

(chimera *agg*)


(chimera (car *all-aggs*))
(chimera (cadr *all-aggs*))
(chimera (caddr *all-aggs*))
(chimera (cadddr *all-aggs*))

(defparameter *ca* (chem:make-coordinate-array-from-atom-list *a*))

(setq *se* (core:make-cxx-object 'chem:superpose-engine))
(chem:set-fixed-all-points *se* *ca*)
(chem:set-moveable-all-points *se* *ca*)
(chem:superpose *se*)


(cando::array-stereoisomers #P"fuzi/run.mol2" 32)

(cando::chimera #P"fuzi/run_grid.mol2")

(defparameter *energy-function* (chem:make-energy-function *agg* *ff*))
(dump-atoms *agg*)
(chem:calculate-energy *energy-function*)
(defparameter *min* (chem:make-minimizer :energy-function *energy-function*))
;;(defparameter *min* (chem:make-minimizer :energy-function *energy-function*))
(chem:minimize *min*)
(chem:save-mol2 *agg* #P"/tmp/after.mol2")



(dotimes (stereoisomer (number-of-stereoisomers *stereocenters*))
  (set-stereoisomer *stereocenters* stereoisomer 31)
  ;; Do something with the *agg*
  )


(print "Hello")
;; Do minimization using energy function





(setq *print-array* t)
(make-array 5)



(let ((p (chem:pdb-writer "/tmp/mol.pdb")))
  (chem:write p *a*)
  (chem:close
(chem:map-atoms (lambda (a) (print a)) *a*)

(defparameter *pattern* (make-cxx-object 'chem:chem-info))
(chem:compile-smarts *pattern* "C")
(defparameter *atoms* nil)
(chem:map-atoms (lambda (a) (push a *atoms*))

(defparameter *parms*
  (let ((parms (chem:make-read-amber-parameters)))
    (with-open-file (fin "~/Development/clasp/projects/cando/work/forceField/ATOMTYPE_GFF.DEF" :direction :input)
      (chem:read-types parms fin))
    (with-open-file (fin "~/Development/clasp/projects/cando/work/forceField/gaff.dat" :direction :input)
      (chem:read-parameters parms fin)
      parms)
    ))
(defparameter *ff* (chem:get-force-field *parms*))
(defparameter *ff-in-str* (as-string (chem:get-force-field *parms*)))
(defparameter *ff* (cando-sys:from-string *ff-in-str*))


(defparameter *types-in-str* (as-string (chem:get-types *parms*)))

(defparameter *types-from-str* (cando-sys:from-string *types-in-str*))


(progn
  (defparameter *h* (make-cxx-object 'chem:atom :name 'H1 :element 'chem-keyword::H))
  (defparameter *cl* (make-cxx-object 'chem:atom :name 'CL1 :element 'chem-keyword::|Cl|))
  (defparameter *r* (make-cxx-object 'chem:residue))
  (defparameter *m* (make-cxx-object 'chem:molecule))
  (chem:add-matter *m* *r*)
  (chem:add-atom *r* *h*)
  (chem:add-atom *r* *cl*)
  (chem:bond-to-single *h* *cl*)
  )

(defparameter *s* (as-string *m*))
(defparameter *m1* (cando-sys:from-string *s*))





(untrace)
(trace clos:ensure-generic-function-using-class)
(clos:ensure-generic-function-using-class nil 'foo)

(trace (clos:ensure-generic-function-using-class :break t))
(apropos "ensure-generic-function-u

(print "Hello")
