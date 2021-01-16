(in-package :asdf-user)

(defsystem "cando-build"
    :description "Build atoms, residues, molecules, aggregates in convenient ways"
    :version "0.0.1"
    :author "Michael Yonkunas"
    :licence "LGPL-3.0"
    :depends-on ()
    :serial t
    :components
    ((:file "packages")
     (:file "build")
     ))
