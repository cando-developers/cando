(asdf:defsystem #:cando-serialize
  :description "The CANDO serialization framework"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on ()
  :serial t
  :components ((:file "packages")
               (:file "print-read")))
