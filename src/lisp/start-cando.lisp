(progn
  (require :asdf)
  (let ((central-registry (find-symbol "*CENTRAL-REGISTRY*" :asdf))
        (load-system (find-symbol "LOAD-SYSTEM" :asdf)))
    (setf (symbol-value central-registry)
          (cons (translate-logical-pathname "cando:cando;") (symbol-value central-registry)))
    (funcall load-system "cando"))
  (format t "Done initialization pid = ~a~%" (getpid)))


(setq *default-pathname-defaults* #P"~/Development/Clasp/projects/cando/work/molecules/")
(defparameter *cd* (with-open-file (fin "onemol2.cdxml" :direction :input) (chem:make-chem-draw fin)))
(defparameter *a* (chem:as-aggregate *cd*))
(defparameter *res* (chem:map-atoms 'vector (lambda (a) (chem:randomize-atom-position a) a) *a*))
(chem:save-mol2 *a* #P"/tmp/thing.mol2")



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




