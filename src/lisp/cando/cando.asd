(in-package :asdf-user)

(defsystem "cando"
    :description "The CANDO compiler front end"
    :version "0.0.1"
    :author "Christian Schafmeister <chris.schaf@verizon.net>"
    :licence "LGPL-3.0"
    :depends-on ()
    :serial t
    :components
    ((:file "packages")
     (:file "print-read")
     (:file "conditions")
     (:file "select")
     (:file "geom")
     (:file "anchor")
     (:file "struct")
     (:file "molecules")
     (:file "atom-tree")
     (:file "ring-builder")
     (:file "chemdraw")
     (:file "energy")
     (:file "psf-files")
     (:file "remove-overlaps")
     (:file "matter")
     ))
