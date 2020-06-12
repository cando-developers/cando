(in-package :asdf-user)

(asdf:defsystem load-cando-jupyter
  :version "1.2.3"
  :license "Artistic"
  :author "Christian Schafmeister <meister@temple.edu>"
  :maintainer "Christian Schafmeister <meister@temple.edu>"
  :description "A library providing stuff notebooks"
  :homepage ""
  :serial T
  :components ()
  :depends-on (:cffi :cando-user :common-lisp-jupyter :nglview-cl))
