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

(defparameter *max-sd-steps* 100000)
(defparameter *sd-tolerance* 100000.0)
(defparameter *max-cg-steps* 100000)
(defparameter *cg-tolerance* 100.0)
(defparameter *max-tn-steps* 10000)
(defparameter *tn-tolerance* 0.00000001)

(defparameter *ff* nil)
(export '*ff*)

#|
(defun setup-gaff ()
  (let ((*default-pathname-defaults*
          (translate-logical-pathname #P"source-dir:extensions;cando;src;data;force-field;"))
        (parms (chem:make-read-amber-parameters)))
    (let ((*
    (with-open-file (fin "ATOMTYPE_GFF.DEF" :direction :input)
      (chem:read-types parms fin))
    (with-open-file (fin "gaff.dat" :direction :input)
      (chem:read-parameters parms fin (symbol-value (find-symbol "*AMBER-SYSTEM*" :leap))))
    (setf *ff* (chem:get-force-field parms)))
  *ff*)
      |#


(defun minimize-minimizer (minimizer &key (restraints-on t)
                                       (max-sd-steps *max-sd-steps*)
                                       (sd-tolerance *sd-tolerance*)
                                       (max-cg-steps *max-cg-steps*)
                                       (cg-tolerance *cg-tolerance*)
                                       (max-tn-steps *max-tn-steps*)
                                       (tn-tolerance *tn-tolerance*)
                                       (print-intermediate-results t)
                                       (report-every-n-steps 10)
                                       resignal-error
                                       &allow-other-keys)
  "Minimize the conformational energy for an energy-function"
  (let ((energy-function (chem:get-energy-function minimizer)))
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
    (if print-intermediate-results
        (chem:enable-print-intermediate-results minimizer report-every-n-steps)
        (chem:disable-print-intermediate-results minimizer))
    (cando:minimize-no-fail minimizer :resignal-error resignal-error)))

(defun save-minimizer-coordinates (coordinates minimizer-trajectory)
  (let ((single-float-coordinates (make-array (length coordinates) :element-type 'single-float :adjustable nil)))
    (loop for index from 0 below (length coordinates)
          for sf-val = (float (elt coordinates index) 1.0s0)
          do (setf (elt single-float-coordinates index) sf-val))
    (vector-push-extend single-float-coordinates minimizer-trajectory)))


(defun minimize-energy-function (energy-function &rest args
                                 &key (restraints-on t)
                                   (max-sd-steps *max-sd-steps*)
                                   (sd-tolerance *sd-tolerance*)
                                   (max-cg-steps *max-cg-steps*)
                                   (cg-tolerance *cg-tolerance*)
                                   (max-tn-steps *max-tn-steps*)
                                   (tn-tolerance *tn-tolerance*)
                                   (frozen nil)
                                   resignal-error
                                   (save-trajectory nil)
                                   (print-intermediate-results t)
                                   (report-every-n-steps 10))
  "Minimize the conformational energy for an energy-function"
  (let ((minimizer (chem:make-minimizer energy-function)))
    (chem:minimizer-set-frozen minimizer frozen)
    (if save-trajectory
        (let ((minimizer-trajectory (make-array 16 :adjustable t :fill-pointer 0)))
          (chem:set-step-callback minimizer (lambda (coords) (save-minimizer-coordinates coords minimizer-trajectory)))
          (apply #'minimize-minimizer minimizer args)
          (let ((matter (chem:get-matter energy-function)))
            (format t "matter ~a~%" matter)
            (unless matter (error "There is no matter defined for energy-function ~a" energy-function))
            (values energy-function
                    (make-instance 'dynamics:trajectory
                                   :matter matter
                                   :number-of-atoms (chem:number-of-atoms matter)
                                   :coordinates minimizer-trajectory))))
        (progn
          (apply #'minimize-minimizer minimizer args)
          energy-function))))

(defun minimize (agg &rest args
                 &key (restraints-on t)
                   (max-sd-steps *max-sd-steps*)
                   (sd-tolerance *sd-tolerance*)
                   (max-cg-steps *max-cg-steps*)
                   (cg-tolerance *cg-tolerance*)
                   (max-tn-steps *max-tn-steps*)
                   (tn-tolerance *tn-tolerance*)
                   (use-excluded-atoms t)
                   (assign-types t)
                   (resignal-error nil)
                   (print-intermediate-results t)
                   (report-every-n-steps 10)
                   (save-trajectory nil))
  "Minimize the conformational energy for an aggregate"
  (let ((energy-func (chem:make-energy-function :matter agg :use-excluded-atoms use-excluded-atoms :assign-types assign-types)))
    (apply #'minimize-energy-function energy-func args)))

(defun minimize-energy-function-from-bad-geometry (energy-function &key (restraints-on t)
                                                                     (max-sd-steps *max-sd-steps*)
                                                                     (sd-tolerance *sd-tolerance*)
                                                                     (max-cg-steps *max-cg-steps*)
                                                                     (cg-tolerance *cg-tolerance*)
                                                                     (max-tn-steps *max-tn-steps*)
                                                                     (tn-tolerance *tn-tolerance*)
                                                                     (report-every-n-steps 10)
                                                                     )
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
    (chem:enable-print-intermediate-results minimizer report-every-n-steps)
    (format t "Minimizing with nonbond terms disabled~%")
    (chem:set-option energy-function 'chem:nonbond-term nil)
    (cando:minimize-no-fail minimizer)
    (format t "Minimizing with nonbond terms enabled~%")
    (chem:set-option energy-function 'chem:nonbond-term t)
    (cando:minimize-no-fail minimizer)))

(defun minimize-from-bad-geometry (agg &rest args
                                   &key (restraints-on t)
                                     (max-sd-steps *max-sd-steps*)
                                     (sd-tolerance *sd-tolerance*)
                                     (max-cg-steps *max-cg-steps*)
                                     (cg-tolerance *cg-tolerance*)
                                     (max-tn-steps *max-tn-steps*)
                                     (tn-tolerance *tn-tolerance*)
                                     (use-excluded-atoms t)
                                     (assign-types t))
  "Minimize the conformational energy for an aggregate"
  (let ((energy-func (chem:make-energy-function :matter agg :use-excluded-atoms use-excluded-atoms :assign-types assign-types)))
    (apply #'minimize-energy-function-from-bad-geometry energy-func args)
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



