(in-package :asdf-user)

(defsystem "cando-user"
    :description "The CANDO compiler front end"
    :version "0.0.1"
    :author "Christian Schafmeister <chris.schaf@verizon.net>"
    :licence "LGPL-3.0"
  :depends-on (:cando-primop
               :cando
               :charges
               #+(or):design
               :leap
               :leap-commands
               #+(or):cando-jupyter
               :sdf
               :tirun
               :sketch2d
               :cando-widgets
               :cando-build)
    :serial t
  :components ((:file "packages")
                 (:file "start-swank")
                 (:file "graphviz")
                 (:file "cl-jupyter")
                 (:file "interactive-tools")))
