(in-package :asdf-user)

(defsystem "charges"
    :description "Calculate atomic charges by a variety of approaches"
    :version "0.0.1"
    :author "Shiho Nagai <>"
    :licence "LGPL-2.0"
    :depends-on (:aromaticity)
    :serial t
    :components
    ((:file "packages")
     (:file "am1bcc")
     (:file "sqm-interface")
     ;;     (:file "bcc")
     ))
