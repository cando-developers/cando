(in-package :asdf-user)

(asdf:defsystem #:cando-user
  :description "The CANDO user system"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (#:cando-primop
               #:cando
               #:charges
               #+(or)#:design
               #:smiles
               #:leap
               #:leap-commands
               #:sdf
               #:tirun
               #:sketch2d
               #:cando-build
               #:lparallel)
  :serial t
  :components ((:file "packages")
               (:file "start-swank")
               (:file "graphviz")
               (:file "cl-jupyter")
               (:file "extension")))
