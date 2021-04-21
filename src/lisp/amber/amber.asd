(in-package :asdf-user)

(defsystem "amber"
  :description "Amber related code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on ( :cl-ppcre)
  :serial t
  :components
  ((:file "packages")
   (:file "amber-jobs")
   (:file "graphviz-amber")))
