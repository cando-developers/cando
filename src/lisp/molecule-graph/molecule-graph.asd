;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :asdf-user)

(defsystem "molecule-graph"
  :description "Work with molecule-graphs"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>, Nagai Shiho"
  :licence "Private"
  :depends-on ()
  :serial t
;;;  :build-operation asdf:monolithic-compile-bundle-op
;;;  :build-pathname #P"/tmp/fep.fasb"
  :components
  ((:file "packages")
   (:file "molecule-graph")
   (:file "graphviz")
   ))

