;;;
;;;    File: print-read.lisp
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
(in-package :cando)

;; use
(defmethod cl:print-object ((object core:cxx-object) stream)
  (if (and *print-readably* (core:fieldsp object))
      (core:print-cxx-object object stream)
      (call-next-method)))


(defun read-cxx-object (stream char n)
  (declare (ignore char))
  (let ((description (read stream t nil t)))
    (apply #'core:load-cxx-object (car description) (cdr description))))

(defun as-string (obj)
  (let ((*print-readably* t)
        (*print-pretty* nil)
        (*print-circle* t))
    (with-output-to-string (sout)
      (print obj sout))))


(defun from-string (str)
  (with-input-from-string (sin str)
    (read sin)))


(set-dispatch-macro-character #\# #\I #'read-cxx-object)


(defun save-cando (obj pathname)
  "* Arguments
- obj :: An object.
- pathname :: A pathname.
* Description
Save the object to the file PATHNAME as an s-expression."
  (with-open-file (fout pathname :direction :output)
    (let ((*print-readably* t)
          (*print-pretty* nil)
          (*print-circle* t))
      (print obj fout))))

(defun load-cando (pathname)
  (with-open-file (fin pathname :direction :input)
    (read fin)))
