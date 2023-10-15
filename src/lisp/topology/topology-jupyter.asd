(in-package :asdf-user)

(defsystem "topology-jupyter"
  :description "Build topologys using s-expressions jupyter tools"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:topology :sketch2d)
  :serial t
  :components ((:file "packages-jupyter")
               (:file "jupyter")
               ))
