(in-package :asdf-user)

(defsystem "rosetta-params"
  :description "Rosetta .params parsing utilities"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:cando :uiop)
  :serial t
  :components ((:module "src"
                :components ((:module "params"
                                      :components ((:file "packages")
                                                   (:file "params")
                                                   (:file "atom-properties")
                                                   (:file "lk-type-assign")))))))

