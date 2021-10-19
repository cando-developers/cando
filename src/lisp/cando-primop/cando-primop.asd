(in-package :asdf-user)

(defsystem "cando-primop"
    :description "The primitive operations for Cando"
    :version "0.0.1"
    :author "Christian Schafmeister <chris.schaf@verizon.net>"
    :licence "LGPL-3.0"
  :depends-on ()
    :serial t
  :components ((:file "primops")))
