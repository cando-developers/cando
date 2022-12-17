;;;
;;;    File: geom.lisp
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
(in-package :geom)

(defun center-of-geometry (atoms)
  (let ((sum (geom:vec 0.0 0.0 0.0)))
    (dolist (a atoms)
      (setf sum (geom:v+ sum (chem:get-position a))))
    (geom:v* sum (/ 1.0 (length atoms)))))

#+(or)
(define-compiler-macro geom::+ (&rest numbers)
  (core:expand-associative 'geom::+ 'geom:v+ numbers 0))

(defmethod print-object ((vec geom:v3) stream)
  (if *print-readably*
      (call-next-method)
      (format stream "#<geom:v3 ~a ~a ~a>" (geom:vx vec) (geom:vy vec) (geom:vz vec))))
