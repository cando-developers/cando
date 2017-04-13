(in-package :asdf-user)

(defsystem "leap"
    :description "Amber related code"
    :version "0.0.1"
    :author "Christian Schafmeister <chris.schaf@verizon.net>"
    :licence "LGPL-3.0"
    :depends-on ()
    :serial t
    :components
    ((:file "packages")
     (:file "core")
     (:file "off")
     (:file "cando-database")
;     (:file "commands")     ;;; moved to leap
     (:file "leap")
     (:file "pdb")
     ))
