(in-package :asdf-user)

(defsystem "design-jupyter"
  :description "Design related code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:design)
  :serial t
  :components
  ((:file "packages-debug")
   (:file "debug-graphs")
   ))
