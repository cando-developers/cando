(in-package :asdf-user)

(defsystem "leap"
  :description "Amber related code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on ( :cando
                :smirnoff
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
   (:file "antechamber-type-definition-parser")
   (:file "chem")
   (:file "grammar")
   (:file "off")
   (:file "topology")
   (:file "solvate")
   (:file "add-ions")
   (:file "set-box")
   (:file "cando-database")
   (:file "pdb")
   (:file "leap")
   (:file "commands")
   (:file "easy")
                                        ;     (:file "interpreter")
   ))
