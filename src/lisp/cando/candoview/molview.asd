(in-package :asdf-user)

(defsystem "molview"
    :description "Aggregate viewer"
    :version "0.0.1"
    :author "Christian Schafmeister <chris.schaf@verizon.net> - stolen from Doug Hoyte molview"
    :licence "LGPL-3.0"
    :depends-on ()
    :serial t
    :components
    ((:file "packages")
     (:file "candoview")
     ))
