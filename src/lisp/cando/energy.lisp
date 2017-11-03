;;;
;;;    File: energy.lisp
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
(in-package :energy)




(defparameter *ff* nil
  "Store the force-field")

(defun setup-amber ()
  (let ((*default-pathname-defaults*
         (translate-logical-pathname #P"cando:data;force-field;")))
    (defparameter *parms*
      (let ((parms (chem:make-read-amber-parameters)))
        (with-open-file (fin "ATOMTYPE_GFF.DEF" :direction :input)
          (chem:read-types parms fin))
        (with-open-file (fin "gaff.dat" :direction :input)
          (chem:read-parameters parms fin)
          parms)))
    (defparameter *ff* (chem:get-force-field *parms*)))
  *ff*)



(defun minimize-energy-function (energy-function &key (restraints-on t)
                                                   (force-field *ff*)
                                                   (max-sd-steps 1000)
                                                   (max-cg-steps 50000)
                                                   (max-tn-steps 0)
                                                   (sd-tolerance 5000.0)
                                                   (cg-tolerance 0.5)
                                                   (tn-tolerance 0.00001))
  "Minimize the conformational energy for an energy-function"
  (let ((minimizer (chem:make-minimizer energy-function)))
    (unless restraints-on
      (let ((restraint-term (chem:get-anchor-restraint-component energy-function)))
        (chem:disable restraint-term)))
    (cando:configure-minimizer minimizer
                               :max-sd-steps max-sd-steps
                               :max-cg-steps max-cg-steps
                               :max-tn-steps max-tn-steps
                               :sd-tolerance sd-tolerance
                               :cg-tolerance cg-tolerance
                               :tn-tolerance tn-tolerance)
    (chem:enable-print-intermediate-results minimizer)
    (chem:set-option energy-function 'chem:nonbond-term nil)
    (cando:minimize-no-fail minimizer)
    (chem:set-option energy-function 'chem:nonbond-term t)
    (cando:minimize-no-fail minimizer)))

(defun minimize (agg &rest args
                 &key (restraints-on t)
                   (force-field *ff*)
                   (max-sd-steps 1000)
                   (max-cg-steps 50000)
                   (max-tn-steps 0)
                   (sd-tolerance 5000.0)
                   (cg-tolerance 0.5)
                   (tn-tolerance 0.00001)
                   (use-excluded-atoms t))
  "Minimize the conformational energy for an aggregate"
  (let ((energy-func (chem:make-energy-function agg force-field :use-excluded-atoms use-excluded-atoms)))
    (apply #'minimize-energy-function energy-func args)
    energy-func))


;;;------------------------------------------------------------
;;;
;;; Helper routines for energy-functions
;;;

(defun excluded-atom-name (excluded-atom-index energy-function)
  (let ((atom-table (chem:atom-table energy-function)))
    (chem:elt-atom-name atom-table excluded-atom-index)))

;;; Exported
(defun dump-excluded-atom-list (energy-function)
  "Dump the excluded atom list for the nonbond term of the energy-function."
  (chem:dump-terms (chem:atom-table energy-function))
  (let* ((nonbond (chem:get-nonbond-component energy-function))
	 (number-excluded-atoms (chem:number-excluded-atoms nonbond))
	 (excluded-atom-list (chem:excluded-atom-list nonbond))
	 (excluded-atom-index-index 0))
    (loop for numi below (core:native-vector<int>-size number-excluded-atoms)
       for num = (core:native-vector<int>-elt number-excluded-atoms numi)
       do (progn
	    (format t "atom#~a ~a number-excluded-atoms ~a -----------~%" numi (excluded-atom-name numi energy-function) num)
	    (loop for ei below num
	       for eaindex = (core:native-vector<int>-elt excluded-atom-list excluded-atom-index-index)
	       do (progn
                    (format t "  excluded-atom-list[~a] -> ~a name: ~a~%" excluded-atom-index-index eaindex (excluded-atom-name eaindex energy-function))
                    (incf excluded-atom-index-index)))))))
