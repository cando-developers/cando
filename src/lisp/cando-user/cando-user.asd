(in-package :asdf-user)

(defsystem "cando-user"
    :description "The CANDO compiler front end"
    :version "0.0.1"
    :author "Christian Schafmeister <chris.schaf@verizon.net>"
    :licence "LGPL-3.0"
    :depends-on (:cando :leap)
    :serial t
    :components ((:file "packages")))
