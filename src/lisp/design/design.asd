(in-package :asdf-user)

(defsystem "design"
  :description "Design related code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (cando )
  :serial t
  :components
  ((:file "packages")
   (:file "build-internal-coordinate-joint-template-tree")
   ))
