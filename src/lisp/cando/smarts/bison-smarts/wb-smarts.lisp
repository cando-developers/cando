;;;
;;;    File: wb-smarts.lisp
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

;;;
;;; Initialize everything
;;;
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
(progn
  (setq *default-pathname-defaults*
        #P"~/Development/Clasp/projects/cando/src/lisp/cando/smarts/")
  (defparameter *cd*
    (with-open-file
        (fin (probe-file "small-mol.cdxml") :direction :input)
      (chem:make-chem-draw fin)))
  (defparameter *agg* (chem:as-aggregate *cd*))
  (defparameter *mol* (chem:content-at *agg* 0)))

(progn
  (defparameter *stereocenters*
    (gather-stereocenters *agg*))
  (sort *stereocenters* #'string< :key #'chem:get-name)
  (set-stereoisomer *stereocenters* #b1 :show t)
  (format t "~a~%" *stereocenters*))


(cando:chimera *agg*)
(cando:scramble-positions *agg*)
(cando:chimera *agg*)
(cando:build-good-geometry-from-random *agg* *ff*)
(cando:chimera *agg*)

(defparameter *n* (let (n)
                    (chem:map-atoms
                     nil
                     (lambda (a)
                       (if (eq (chem:get-name a) ckw:N)
                           (setq n a)))
                     *agg*)
                    n))


(defparameter *p* (core:make-cxx-object 'chem:chem-info))
(chem:compile-smarts *p* "N")
(format t "Matches: ~a~%" (chem:matches *p* *n*))
(chem:get-code *p*)


(getpid)49930
