(asdf:defsystem #:cando-jupyter
  :description "The CANDO compiler front end jupyter interface"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (#:cando
               #:common-lisp-jupyter
               #:cando-user
               #:cando-widgets
               #:smarts-jupyter
               #:tirun-jupyter)
  :serial t
  :components ((:file "packages")
               (:file "kernel")))
