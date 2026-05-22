;;;; crystallography.asd

(asdf:defsystem #:xtal
  :description "Describe crystallography here"
  :author "Christian Schafmeister"
  :license  "Specify license here"
  :version "0.0.1"
  :serial t
  :depends-on (:esrap)
  :components ((:file "package")
               (:file "symmetry-parser")
               (:file "xtal")))
