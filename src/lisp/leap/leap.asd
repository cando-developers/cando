(in-package :asdf-user)

(defsystem "leap"
  :description "Amber related code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on ( :cando
                :inet
                :alexandria
                #+(or)(:version :esrap "0.15")
                #+(or):parser.common-rules
                #+(or):PARSER.COMMON-RULES.OPERATORS
                #+(or):architecture.builder-protocol)
  :serial t
  :components
  ((:file "packages")
   (:file "core")
   (:file "chem")
   (:file "antechamber-type-definition-parser")
   (:file "off")
   (:file "topology")
   (:file "solvate")
   (:file "add-ions")
   (:file "set-box")
   (:file "cando-database")
   (:file "pdb")
   (:file "easy")
   (:file "leap")
                                        ;     (:file "interpreter")
   ))
