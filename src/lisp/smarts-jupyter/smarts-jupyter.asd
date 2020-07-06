
(in-package :asdf-user)

(defsystem "smarts-jupyter"
    :description "Generate cytoscape graphs for SMARTS and SMIRKS"
    :version "0.0.1"
    :author "Christian Schafmeister <>"
    :licence "LGPL-2.0"
  :depends-on (:smarts
    :cytoscape-clj)
    :serial t
    :components
    ((:file "packages")
     (:file "smarts")))

