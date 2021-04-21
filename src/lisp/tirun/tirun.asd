;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :asdf-user)

(defsystem "tirun"
  :description "TIRUN setup code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>, Nagai Shiho"
  :licence "Private"
  :depends-on ( :cl-markup
                :cando
                :charges
                :leap
                :amber
                :alexandria
                :split-sequence
                #+(or):common-lisp-jupyter
                (:version :esrap "0.15")
                :parser.common-rules
                :PARSER.COMMON-RULES.OPERATORS
                :architecture.builder-protocol
                :cl-ppcre
                :molecule-graph)
  :serial t
;;;  :build-operation asdf:monolithic-compile-bundle-op
;;;  :build-pathname #P"/tmp/tirun.fasb"
  :components
  ((:file "packages")
   (:file "assemble")
   (:file "graph2")
   (:file "tirun-version")
   (:file "tirun")
   (:file "tirun-from-structures")
   (:file "script")
   (:file "ti")
   (:file "graphviz")
   ))

#+(or)
(defsystem "tirun/build"
  :description "TIRUN build code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>, Nagai Shiho"
  :licence "Private"
  :depends-on (:tirun)
  :serial t
  :build-operation asdf:monolithic-compile-bundle-op
  :build-pathname #P"/tmp/tirun.fasb")
