(in-package :asdf-user)

(defsystem "inet"
    :description "Code for interacting with the internet"
    :version "0.0.1"
    :author "Christian Schafmeister <chris.schaf@verizon.net>"
    :licence "LGPL-3.0"
    :depends-on (:trivial-http)
    :serial t
    :components
    ((:file "packages")
     (:file "inet")
     ))
