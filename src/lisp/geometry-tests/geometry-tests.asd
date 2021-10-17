(in-package :asdf-user)

(defsystem "geometry-tests"
  :description "Code to determine if geometry is good"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:cando :amber)
  :serial t
  :components
  ((:file "packages")
   (:file "utils")
   ))
