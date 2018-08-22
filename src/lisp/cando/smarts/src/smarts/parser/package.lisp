;;;; package.lisp --- Package definition for the smarts.parser module.
;;;;
;;;; Copyright (C) 2018 Jan Moringen
;;;;
;;;; Author: Jan Moringen <jmoringe@techfak.uni-bielefeld.de>

(cl:defpackage #:language.smarts.parser
  (:use
   #:cl
   #:alexandria
   #:esrap)

  (:shadow
   #:parse)

  ;; Parser protocol
  (:export
   #:parse))
