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
(in-package :cando.serialize)

(defgeneric skip-slot-names (obj)
  (:documentation "Specialize this method on classes that have slots that we don't want to serialize.
Return a list of slot names that shouldn't be serialized."))
(defmethod skip-slot-names ((obj t))
  "All slots that can be serialized will be serialized"
  nil)

(defun print-object-readably-with-slots (obj stream)
  (format stream "#$(~s " (class-name (class-of obj)))
  (let ((skip-slot-names (skip-slot-names obj)))
    (loop for slot in (clos:class-slots (class-of obj))
          for slot-name = (clos:slot-definition-name slot)
          for initargs = (clos:slot-definition-initargs slot)
          if (and (car initargs)
                  (not (position slot-name skip-slot-names))
                  (slot-boundp obj slot-name))
            do (format stream "~s ~s " (car initargs) (slot-value obj slot-name))))
  (format stream ") "))

(export '(print-object-readably-with-slots))

(defun as-string (obj)
  (let ((*print-readably* t)
        (*print-pretty* nil)
        (sys:*print-dense* t)
        (*print-circle* t))
    (with-output-to-string (sout)
      (prin1 obj sout))))

(defun from-string (str)
  (with-input-from-string (sin str)
    (read sin)))


(defun save-cando (obj pathname)
  "* Arguments
- obj :: An object.
- pathname :: A pathname.
* Description
Save the object to the file PATHNAME as an s-expression."
  (let ((temp-pathname (pathname (sys:mkstemp (format nil "~a-" (namestring pathname)))))
        (new-pathname (merge-pathnames pathname *default-pathname-defaults*)))
    (with-open-file (fout temp-pathname :direction :output)
      (let ((*print-readably* t)
            (*print-pretty* nil)
            (sys:*print-dense* t)
            (*print-circle* t)
            (*package* (find-package :keyword)))
        (cl:prin1 obj fout)))
    (rename-file temp-pathname new-pathname :if-exists :supersede)))


(defun sharp-$-reader (stream subchar arg)
  (declare (ignore subchar))
  (when (and arg (null *read-suppress*))
        (error "~S is an extra argument for the #s readmacro." arg))
  (let ((l (read stream t nil t)))
    (when *read-suppress*
      (return-from sharp-$-reader nil))
    ;; Do the construction.
    (apply 'make-instance l)))

(defparameter *cando-reader* (copy-readtable nil))
(set-dispatch-macro-character #\# #\$ 'sharp-$-reader *cando-reader*)

(defun load-cando (pathname)
  (with-open-file (fin pathname :direction :input)
    (let ((*readtable* *cando-reader*))
      (read fin))))

(defmacro make-class-save-load (class-name &key skip-slot-names
                                             (print-unreadably
                                              `(lambda (obj stream)
                                                 (print-unreadable-object (obj stream)
                                                   (format stream "~a" (class-name (class-of obj)))))))
  "Create a serializer for class-name. Slots that have :initarg defined and do not appear in
  skip-slot-names will be serialized."
  (unless (and (listp skip-slot-names)
               (every #'symbolp skip-slot-names))
    (error "The skip-slot-names must be a list of symbols - instead it is ~s" skip-slot-names))
  `(defmethod print-object ((obj ,class-name) stream)
     (if *print-readably*
         (progn
           (print-object-readably-with-slots obj stream))
         (funcall ,print-unreadably obj stream))))


(defclass serializable ()
  ())

(defmethod print-object ((obj serializable) stream)
  (if *print-readably*
      (progn
        (print-object-readably-with-slots obj stream))
      (print-unreadable-object (obj stream :type t))))

#+use-mpi
(eval-when (:load-toplevel :execute)
  (setf mpi:*decode-object-hook*
        (lambda (str)
          (with-input-from-string (sin str)
            (let ((*readtable* *cando-reader*))
              (read sin)))))

  (setf mpi:*encode-object-hook*
        (lambda (obj)
          (let ((str (make-array 256 :element-type 'base-char :adjustable t :fill-pointer 0)))
          (with-output-to-string (fout str)
            (let ((*print-readably* t)
                  (*print-pretty* nil)
                  (sys:*print-dense* t)
                  (*print-circle* t)
                  (*package* (find-package :keyword)))
              (cl:prin1 obj fout))
            (copy-seq str))))))


