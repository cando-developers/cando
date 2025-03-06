;;; Copyright (c) 2020, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(asdf:defsystem #:chem-info
  :description "Chem-informatics"
  :version "0.0.1"
  :author ("Christian Schafmeister <chris.schaf@verizon.net>")
  :licence "GPL2"
  :depends-on nil
  :serial t
  :components ((:file "packages")
               (:file "pucker")))

