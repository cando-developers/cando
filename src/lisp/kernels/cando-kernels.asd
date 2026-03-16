

(asdf:defsystem #:cando-kernels
  :description "Optimize mathematical expressions using a compiler implemented in common lisp "
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:mathkernel)
  :serial t
  :components ((:file "energy-kernels")))
