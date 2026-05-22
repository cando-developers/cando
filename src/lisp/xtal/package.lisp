;;;; package.lisp

(defpackage #:xtal
  (:use #:cl)
  (:export
   ))

(defpackage #:xtal.symop-parser
  (:use :cl :esrap)
  (:export
   #:symbolic-to-row-major-matrix
   #:symop))
