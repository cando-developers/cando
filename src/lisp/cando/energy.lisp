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


(defun minimize (agg &key (restraints-on t)
                       (force-field *ff*)
                       (max-sd-steps 1000)
                       (max-cg-steps 50000)
                       (max-tn-steps 0)
                       (sd-tolerance 5000.0)
                       (cg-tolerance 0.5)
                       (tn-tolerance 0.00001))
  "Minimize the conformational energy"
  (format t "Got minimizer~%")
  (let* ((energy-func (chem:make-energy-function agg force-field))
         (minimizer (chem:make-minimizer :energy-function energy-func)))
    (unless restraints-on
      (let ((restraint-term (chem:get-anchor-restraint-component energy-func)))
        (chem:disable restraint-term)))
    (cando:configure-minimizer minimizer
                               :max-sd-steps max-sd-steps
                               :max-cg-steps max-cg-steps
                               :max-tn-steps max-tn-steps
                               :sd-tolerance sd-tolerance
                               :cg-tolerance cg-tolerance
                               :tn-tolerance tn-tolerance)
    (chem:enable-print-intermediate-results minimizer)
    (chem:set-option energy-func 'chem:nonbond-term nil)
    (cando:minimize-no-fail minimizer)
    (chem:set-option energy-func 'chem:nonbond-term t)
    (cando:minimize-no-fail minimizer)))
