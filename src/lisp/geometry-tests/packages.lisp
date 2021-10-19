;; -^-
(cl:in-package #:common-lisp-user)

(defpackage #:geometry-tests
  (:use #:cl )
  (:nicknames :gt)
  (:export
   #:plane-from-atoms
   )
  (:documentation
   "Run tests on geometry to determine if there are problems"))
