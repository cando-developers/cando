(in-package :asdf-user)

(defsystem "cando-sdf"
    :description "Parse sdf files"
    :version "0.0.1"
    :author "Christian Schafmeisterd>"
    :licence "LGPL-2.0"
  :depends-on ((:version :esrap "0.15")
               :language.smarts
               :parser.common-rules
               :parser.common-rules.operators
               :architecture.builder-protocol)
    :components
    ((:file "packages")
     (:file "sdf" :depends-on ("packages"))
     (:file "sdf-writer" :depends-on ("packages"))
     ))

