
(in-package :asdf-user)

(defsystem "spanning"
  :description "Generate spanning-trees from graphs"
  :version "0.0.1"
  :author "Michael Yonkunas and Christian Schafmeister"
  :licence "LGPL-2.0"
  :depends-on (:cando)
  :serial t
  :components (
               (:file "packages")
               (:file "spanning")
               ))

