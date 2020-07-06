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

  (:local-nicknames
   (#:bp     #:architecture.builder-protocol)
   (#:smiles #:language.smiles.parser))

  ;; Exported rules
  (:export
   #:atom-pattern
   #:recursive)

  ;; Parser protocol
  (:export
   #:parse)
  ;; Smirks
  (:export
   #:parse-smirks
))
