
(in-package :asdf-user)

(defsystem "lomap"
  :description "Generate a map for FEP calculations"
  :version "0.0.1"
  :author "Michael Yonkunas and Christian Schafmeister"
  :licence "LGPL-2.0"
  :depends-on (:cando
               :molecule-graph
               :smarts)
  :serial t
  :components ((:module src
                :pathname "src"
                :serial t
                :components (
                             (:file "packages")
                             (:file "lomap")
                             (:file "graphviz")
                             (:file "spanning")
                             ))))

