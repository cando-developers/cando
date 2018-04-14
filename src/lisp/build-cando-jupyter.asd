(in-package :asdf-user)

(asdf:defsystem build-cando-jupyter
  :version "1.2.3"
  :license "Artistic"
  :author "Christian Schafmeister <meister@temple.edu>"
  :maintainer "Christian Schafmeister <meister@temple.edu>"
  :description "A library providing stuff notebooks"
  :homepage ""
  :serial T
  :components ()
  :build-operation asdf:monolithic-compile-bundle-op
  :build-pathname #.(merge-pathnames #P"cando-jupyter.fasb" (translate-logical-pathname "app-fasl:"))
  :depends-on (:cffi :cando-user :cl-jupyter :nglview))
