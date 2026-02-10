(in-package :asdf-user)

(defsystem "rosetta-params-load"
  :description "Rosetta .params parsing utilities"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:cando :uiop :rosetta-params)
  :serial t
  :components ((:module "src"
                :components ((:module "load"
                                      :components ((:file "packages")
                                                   (:file "load")))))))
