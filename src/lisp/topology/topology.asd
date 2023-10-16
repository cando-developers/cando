(in-package :asdf-user)

(defsystem "topology"
  :description "Build topologys using s-expressions"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:chem-extras :cando-serialize :cl-conspack :netcdf :named-readtables)
  :serial t
  :components ((:file "packages")
               (:file "readtable")
               (:file "serial")
               (:file "topology-classes")
               (:file "constitution")
               (:file "joint-templates")
               (:file "define-topology")
               (:file "oligomer")
               (:file "joint-tree")
               (:file "fold")
               (:file "conformation")
               (:file "internals")
               (:file "design")
               (:file "oligomer-space-dag")
               (:file "context")
               (:file "badgeom")
               (:file "shape")
               ))
