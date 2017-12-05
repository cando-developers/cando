(in-package :asdf-user)

(asdf:defsystem build-cando
  :version "1.2.3"
  :license "Artistic"
  :author "Christian Schafmeister <meister@temple.edu>"
  :maintainer "Christian Schafmeister <meister@temple.edu>"
  :description "A library providing all of the Cando Common Lisp source code"
  :homepage ""
  :serial T
  :components ()
  :build-operation asdf:monolithic-compile-bundle-op
  :build-pathname #.(merge-pathnames "cando.fasb" (translate-logical-pathname "app-fasl:"))
  :depends-on ("cando-user"))
