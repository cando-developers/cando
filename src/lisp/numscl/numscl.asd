(in-package :asdf-user)

(defsystem "numscl"
  :description "Numerical code for cando"
  :version "0.0.1"
  :author "Christian Schafmeister <meister@thirdlaw.tech>"
  :licence "LGPL-3.0"
  :depends-on (cl-mathstats)
  :serial t
  :components
  ((:file "packages")
   (:file "polyfit")
   (:file "integrate")
   ))
