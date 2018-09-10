;;;
;;;    File: conditions.lisp
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
(in-package :chem)

(define-condition chem::linear-atoms-error ()
  ((atoms :initarg :atoms :reader chem:atoms)
   (coordinates :initarg :coordinates :reader coordinates)
   (indices :initarg :indices :reader chem::indices)))

(define-condition chem:linear-angle-error (linear-atoms-error) ())
(define-condition chem:linear-dihedral-error (linear-atoms-error) ())
(define-condition chem:linear-improper-restraint-error (linear-atoms-error) ())


(define-condition minimizer-error (error)
  ((minimizer :initarg :minimizer
              :reader minimizer)
   (coordinates :initarg :coordinates :reader coordinates)))

(define-condition minimizer-exceeded-max-steps (minimizer-error)
  ((number-of-steps :initarg :number-of-steps :accessor number-of-steps)))

(defun make-minimizer-exceeded-max-steps (minimizer coordinates number-of-steps)
  (make-condition 'minimizer-exceeded-max-steps
                  :minimizer minimizer
                  :coordinates coordinates
                  :number-of-steps number-of-steps))

(define-condition minimizer-stuck (minimizer-error) ())

(defun make-minimizer-stuck (minimizer coordinates number-of-steps)
  (make-condition 'minimizer-stuck
                  :minimizer minimizer
                  :coordinates coordinates))
