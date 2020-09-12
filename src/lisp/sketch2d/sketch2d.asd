(in-package :asdf-user)

(defsystem "sketch2d"
  :description "Sketch 2d figures of molecules"
  :version "0.0.1"
  :author "Christian Schafmeister"
  :licence "LGPL-2.0"
  :depends-on (:cando :leap :cl-svg :spanning)
  :serial t
  :components
  ((:file "packages")
   (:file "infix")
   (:file "sketch2d")
   (:file "render-svg")
   ))
