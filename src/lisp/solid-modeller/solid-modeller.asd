
(in-package :asdf-user)

(defsystem "solid-modeller"
  :description "Programmatically build models"
  :version "0.0.1"
  :author "Christian Schafmeister"
  :licence "LGPL-2.0"
  :depends-on ()
  :serial t
  :components
  ((:file "packages")
   (:file "modeller")
   ))

