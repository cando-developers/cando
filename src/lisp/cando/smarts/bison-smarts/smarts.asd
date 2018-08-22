(in-package :asdf-user)

(defsystem "smarts"
    :description "Parse and process SMARTS and SMIRKS"
    :version "0.0.1"
    :author "Shiho Nagai <>"
    :licence "LGPL-2.0"
    :depends-on ()
    :serial t
    :components
    ((:file "packages")
     (:file "smarts")
     ))
