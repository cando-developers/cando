;;;
;;;    File: start-cando.lisp
;;;

;; Open Source License
;; Copyright (c) 2016, Christian E. Schafmeister
;; Permission is hereby granted, free of charge, to any person obtaining a copy
;; of this software and associated documentation files (the "Software"), to deal
;; in the Software without restriction, including without limitation the rights
;; to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
;; copies of the Software, and to permit persons to whom the Software is
;; furnished to do so, subject to the following conditions:
;; The above copyright notice and this permission notice shall be included in
;; all copies or substantial portions of the Software.
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
;; IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
;; OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
;; THE SOFTWARE.
;;  
;; This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.

;; -^-
(format t "Running start-cando.lisp script~%")

(progn
  (setf (logical-pathname-translations "cando")
        '(("**;*.*" "source-dir:extensions;cando;src;lisp;**;*.*"))))
  (require :asdf)
  (let ((central-registry (find-symbol "*CENTRAL-REGISTRY*" :asdf))
        (load-system (find-symbol "LOAD-SYSTEM" :asdf)))
    (setf (symbol-value central-registry)
          (cons (translate-logical-pathname "cando:cando;") (symbol-value central-registry)))
    #+(or)(funcall load-system "cando"))
  (format t "Done initialization pid = ~a~%" (getpid)))
(push (translate-logical-pathname "cando:cando;amber;") asdf:*central-registry*)

#+(or)
(progn
  (setq *default-pathname-defaults* #P"~/Development/Clasp/projects/cando/work/molecules/")
  (defparameter *cd* (with-open-file (fin "onemol2.cdxml" :direction :input) (chem:make-chem-draw fin)))
  (defparameter *a* (chem:as-aggregate *cd*))
  (defparameter *res* (chem:map-atoms 'vector (lambda (a) (chem:randomize-atom-position a) a) *a*))
  (chem:save-mol2 *a* #P"/tmp/thing.mol2")



  (chem:map-atoms (lambda (a) (print a)) *a*)

  (defparameter *pattern* (make-cxx-object 'chem:chem-info))
  (chem:compile-smarts *pattern* "C")
  (defparameter *atoms* nil)
  (chem:map-atoms (lambda (a) (push a *atoms*)))




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
  )
