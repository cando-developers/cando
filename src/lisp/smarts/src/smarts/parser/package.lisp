;;;; package.lisp --- Package definition for the smarts.parser module.
;;;;
;;;; Copyright (C) 2018 Jan Moringen
;;;;
;;;; Author: Jan Moringen <jmoringe@techfak.uni-bielefeld.de>

(cl:in-package #:common-lisp-user)

;;(cl:defpackage #:language.smarts.parser
(defpackage #:language.smarts.parser
  (:use
   #:common-lisp
   #:cl
   #:alexandria
   #:esrap)

  (:shadow
   #:parse)

  ;; Exported rules
  (:export
   #:atom-pattern
   #:recursive)

  ;; Parser protocol
  (:export
   #:parse))
