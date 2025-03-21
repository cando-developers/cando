(asdf:defsystem #:cando
  :description "The CANDO compiler front end"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (#:smarts #:kinematics-extras #:chem-extras #:cando-serialize #:topology)
  :serial t
  :components ((:file "packages")
               (:file "basic")
               (:file "infix")
               (:file "mechanics")
               (:file "environment")
               (:file "conditions")
               (:file "chemistry")
               (:file "select")
               (:file "convenience")
               (:file "geom")
               (:file "build")
               (:file "progress")
               (:file "anchor")
               (:file "struct")
               (:file "fortran")
               (:file "dynamics")
               (:file "modelling")
               (:file "constitution")
               (:file "molecules")
               (:file "stereochemistry")
               (:file "cif")
               (:file "atom-tree")
               (:file "ring-builder")
               (:file "chemdraw")
               (:file "energy")
               (:file "formula-weight")
               (:file "matter-hierarchy")
               (:file "units")
               (:file "align")
               (:file "psf-files")))

(asdf:defsystem #:cando/documentation
  :description "The CANDO compiler front end documentation"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (#:cando #:mgl-pax)
  :serial t
  :components (
               (:file "documentation")
               )
  )
