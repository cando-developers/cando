;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :asdf-user)

(defsystem "membrane-packer"
  :description "Create membranes around proteins"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>, Nagai Shiho"
  :licence "Private"
  :depends-on ( :cl-markup
                :cando
                :charges
                :leap
                :amber
)
  :serial t
;;;  :build-operation asdf:monolithic-compile-bundle-op
;;;  :build-pathname #P"/tmp/fep.fasb"
  :components
  ((:file "packages")
   (:file "membrane-packer")
   ))

#+(or)
(defsystem "fep/build"
  :description "FEP build code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>, Nagai Shiho"
  :licence "Private"
  :depends-on (:fep)
  :serial t
  :build-operation asdf:monolithic-compile-bundle-op
  :build-pathname #P"/tmp/fep.fasb")
