(in-package :asdf-user)

(defsystem "smirnoff"
  :description "Assign force-field parameters using smirnoff"
  :version "0.0.1"
  :author "Christian Schafmeister"
  :licence "LGPL-2.0"
  :depends-on (:aromaticity
               :smarts
               :cando
               :plump)
  :serial t
  :components
  ((:file "packages")
   (:file "load")
   (:file "smirnoff")
   ))
