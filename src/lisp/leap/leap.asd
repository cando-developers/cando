(in-package :asdf-user)

(defsystem "leap"
    :description "Amber related code"
    :version "0.0.1"
    :author "Christian Schafmeister <chris.schaf@verizon.net>"
    :licence "LGPL-3.0"
    :depends-on ( :alexandria
                  (:version :esrap "0.15")
                  :parser.common-rules
                  :PARSER.COMMON-RULES.OPERATORS
                  :architecture.builder-protocol)
    :serial t
    :components
    ((:file "packages")
     (:file "core")
     (:file "chem")
     (:file "antechamber-type-definition-parser")
     (:file "off")
     (:file "fortran")
     (:file "topology")
     (:file "cando-database")
;     (:file "commands")     ;;; moved to leap
     (:file "pdb")
     (:file "grammar")
     (:file "leap")
;     (:file "interpreter")
     ))
