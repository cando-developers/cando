;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :asdf-user)

(defsystem "tirun-jupyter"
  :description "TIRUN setup code for jupyter widgets"
  :version "0.0.2"
  :author "Christian Schafmeister <chris.schaf@verizon.net>, Nagai Shiho"
  :licence "Private"
  :depends-on (:tirun :lomap :pas :cytoscape-clj :structure-editor :cando-widgets :quri)
  :serial t
;;;  :build-operation asdf:monolithic-compile-bundle-op
;;;  :build-pathname #P"/tmp/tirun.fasb"
  :components
  ((:file "packages")
   (:file "tirun-jupyter")
   ))

