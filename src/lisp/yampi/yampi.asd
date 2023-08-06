(in-package :asdf-user)

(defsystem "yampi"
  :description "Yet Another (Yitzi's Awesome) Message Passing Interface"
  :version "0.1.0"
  :author "Tarn W. Burton"
  :licence "LGPL-3.0"
  :depends-on ("bordeaux-threads"
               "pzmq")
  :serial t
  :components ((:file "packages")
               (:file "message")))
