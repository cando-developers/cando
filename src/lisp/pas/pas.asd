
(in-package :asdf-user)

(defsystem "pas"
    :description "Positional Analogue Scanning"
    :version "0.0.1"
    :author "Shiho Nagai <>"
    :licence "LGPL-2.0"
  :depends-on ((:version :esrap "0.15")
               :language.smarts
               :parser.common-rules
               :parser.common-rules.operators
               :architecture.builder-protocol)
    :serial t
    :components
    ((:file "packages")
     (:file "pas")
     ))
