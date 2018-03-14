(in-package :asdf-user)

(defsystem "aromaticity"
    :description "Calculate aromaticity using a variety of approaches"
    :version "0.0.1"
    :author "Shiho Nagai <>"
    :licence "LGPL-2.0"
    :depends-on (:smarts)
    :serial t
    :components
    ((:file "packages")
     (:file "aromaticity")
     ))
