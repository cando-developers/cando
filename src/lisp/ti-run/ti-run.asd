;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :asdf-user)

(defsystem "ti-run"
  :description "TI-RUN setup code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>, Nagai Shiho"
  :licence "Private"
  :depends-on ( :cl-markup
                :cando
                :charges
                :leap
                :amber
                :alexandria
                #+(or):cl-jupyter
                (:version :esrap "0.15")
                :parser.common-rules
                :PARSER.COMMON-RULES.OPERATORS
                :architecture.builder-protocol
                :cl-ppcre
                :molecule-graph)
  :serial t
;;;  :build-operation asdf:monolithic-compile-bundle-op
;;;  :build-pathname #P"/tmp/ti-run.fasb"
  :components
  ((:file "packages")
   (:file "jsme")
   (:file "graph2")
   (:file "ti-run-version")
   (:file "ti-run")
   (:file "ti-run-from-structures")
   (:file "script")
   (:file "ti")
   (:file "graphviz")
   ))

#+(or)
(defsystem "ti-run/build"
  :description "TI-RUN build code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>, Nagai Shiho"
  :licence "Private"
  :depends-on (:ti-run)
  :serial t
  :build-operation asdf:monolithic-compile-bundle-op
  :build-pathname #P"/tmp/ti-run.fasb")
