(in-package :asdf-user)

(defsystem "smirnoff"
  :description "Assign force-field parameters using smirnoff"
  :version "0.0.1"
  :author "Christian Schafmeister"
  :licence "LGPL-2.0"
  :depends-on (:aromaticity
               :smarts
               :cando
               :plump
               :esrap
               :parser.common-rules)
  :serial t
  :components
  ((:file "packages")
   (:file "units")
   (:file "load")
   (:file "smirnoff")
   ))
