(asdf:defsystem #:cando-user-install
  :description "User install for CANDO"
  :version "1.0.0"
  :author "Tarn W. Burton <twburton@gmail.com>"
  :licence "LGPL-3.0"
  :serial t
  :depends-on (#:cando-jupyter)
  :components ((:file "packages")
               (:file "install")))
