(in-package :asdf-user)

(defsystem "topology"
  :description "Build topologys using s-expressions"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on ()
  :serial t
  :components ((:file "packages")
               (:file "topology-classes")
               (:file "joint-templates")
               (:file "define-topology")
               (:file "oligomer")
               (:file "joint-tree")
               (:file "fold")
               (:file "conformation")
               (:file "optimize-conformations")
               (:file "internals")
               (:file "design")
               (:file "context")
               (:file "badgeom")
               (:file "foldamer")
               (:file "graphviz")))
