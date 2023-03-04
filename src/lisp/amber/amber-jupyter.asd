(in-package :asdf-user)

(defsystem "amber-jupyter"
  :description "Amber related code that requires jupyterlab"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (#:cando
               #:amber
               #:common-lisp-jupyter
               #:cando-widgets)
  :serial t
  :components
  ((:file "packages-jupyter")
   (:file "amber-jupyter")))

