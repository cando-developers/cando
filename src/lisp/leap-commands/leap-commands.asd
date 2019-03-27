(in-package :asdf-user)

(defsystem "leap-commands"
  :description "Amber related code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on ( :leap
                (:version :esrap "0.15")
                :parser.common-rules
                :PARSER.COMMON-RULES.OPERATORS
                :architecture.builder-protocol)
  :serial t
  :components
  ((:file "packages")
   (:file "grammar")
   (:file "commands")
   ))
