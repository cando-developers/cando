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

(defun as-string (obj)
  (let ((*print-readably* t)
        (*print-pretty* nil)
        (*print-circle* t))
    (with-output-to-string (sout)
      (print obj sout))))

(defun from-string (str)
  (with-input-from-string (sin str)
    (read sin)))

(defun save-cando (obj pathname)
  "* Arguments
- obj :: An object.
- pathname :: A pathname.
* Description
Save the object to the file PATHNAME as an s-expression."
  (with-open-file (fout pathname :direction :output)
    (let ((*print-readably* t)
          (*print-pretty* nil)
          (*print-circle* t)
          (*package* (find-package :keyword)))
      (cl:print obj fout))))

(defun load-cando (pathname)
  (with-open-file (fin pathname :direction :input)
    (read fin)))

(defun sharp-$-reader (stream subchar arg)
  (declare (ignore subchar))
  (when (and arg (null *read-suppress*))
        (error "~S is an extra argument for the #s readmacro." arg))
  (let ((l (read stream t nil t)))
    (when *read-suppress*
      (return-from sharp-$-reader nil))
    ;; Do the construction.
    (let ((circle-subst-args (mapcar (lambda (arg) (core::circle-subst (make-hash-table :test 'eq) arg)) l)))
      (apply 'make-instance circle-subst-args))))

(set-dispatch-macro-character #\# #\$ 'sharp-$-reader)


(defun print-object-readably-with-slots (obj stream)
  (format stream "#$(~s " (class-name (class-of obj)))
  (loop for slot in (clos:class-slots (class-of obj))
        for slot-name = (clos:slot-definition-name slot)
        for initargs = (clos:slot-definition-initargs slot)
        if (and (car initargs) (slot-boundp obj slot-name))
          do (format stream "~s ~s " (car initargs) (slot-value obj slot-name)))
  (format stream ") "))

(defmacro make-class-save-load (class-name)
  `(defmethod print-object ((obj ,class-name) stream)
     (if *print-readably*
         (progn
           (print-object-readably-with-slots obj stream))
         (print-unreadable-object (obj stream)
           (format stream "~a" (class-name (class-of obj)))))))
