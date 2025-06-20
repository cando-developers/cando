(in-package :asdf-user)

(defsystem "topology"
  :description "Build topologys using s-expressions"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:kinematics-extras :chem-extras :cando-serialize :cl-conspack :netcdf :named-readtables)
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
               (:file "assembler")
               (:file "internals")
               (:file "design")
               (:file "oligomer-space-dag")
               (:file "context")
               (:file "badgeom")
               (:file "shape")
               (:file "steppers")
               (:file "manipulator")
               ))



(defsystem "topology/documentation"
  :description "Documentation for topology system"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:topology :mgl-pax)
  :serial t
  :components ((:file "documentation")))
