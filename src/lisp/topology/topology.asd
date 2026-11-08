(in-package :asdf-user)

(defsystem "topology"
  :description "Build topologys using s-expressions"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:kinematics-extras :chem-extras :cando-serialize :cl-conspack :netcdf :named-readtables)
  :serial t
  :in-order-to ((test-op (test-op "topology/tests")))
  :components ((:file "packages")
               (:file "readtable")
               (:file "serial")
               (:file "topology-classes")
               (:file "training-context")
               (:file "shape-key-api")
               (:file "shape-key")
               (:file "shape-key-map")
               (:file "shape-key-binning")
               (:file "constitution")
               (:file "cache-declarations")
               (:file "joint-templates")
               (:file "define-topology")
               (:file "oligomer")
               (:file "joint-tree")
               (:file "fold")
               (:file "assembler")
               (:file "internals")
               (:file "derive-topology-group")
               (:file "design")
               (:file "oligomer-space-dag")
               (:file "oligomer-to-sexp")
               (:file "context")
               (:file "badgeom")
               (:file "shape")
               (:file "piece-assembler")
               (:file "steppers")
               (:file "rotamers-api")
               (:file "manipulator")
               (:file "cache-geometry")
               (:file "blueprint")
               ))



(defsystem "topology/tests"
  :description "Regression tests for topology"
  :depends-on ("topology")
  :in-order-to ((test-op (test-op "topology/shape-key-tests")
                        (test-op "topology/training-context-tests")))
  :serial t
  :components ((:file "new-adjustment-tests")
               (:file "tests/adjustments")
               (:file "tests/derive-topology-group")
               (:file "tests/manipulator"))
  :perform (test-op (operation component)
             (declare (ignore operation component))
             (uiop:symbol-call :topology-adjustment-tests :run-tests)
             (uiop:symbol-call :topology-derivation-tests :run-tests)
             (uiop:symbol-call :topology-manipulator-tests :run-tests)))

(defsystem "topology/shape-key-tests"
  :description "Regression checks for the new shape-key foundation"
  :depends-on ("topology")
  :serial t
  :components ((:file "tests/shape-key")
               (:file "tests/cache-declarations")
               (:file "tests/shape-key-binning")
               (:file "tests/shape-key-map")
               (:file "tests/cache-geometry")
               (:file "tests/cache-storage")
               (:file "tests/shape-key-consumers"))
  :perform (test-op (operation component)
             (declare (ignore operation component))
             (uiop:symbol-call :topology-shape-key-tests :run-tests)
             (uiop:symbol-call :topology-shape-key-tests :test-cache-declarations)
             (uiop:symbol-call :topology-shape-key-tests :test-shape-key-binning)
             (uiop:symbol-call :topology-shape-key-tests :test-shape-key-map)
             (uiop:symbol-call :topology-shape-key-tests :test-cache-geometry)
             (uiop:symbol-call :topology-shape-key-tests :test-cache-storage)
             (uiop:symbol-call :topology-shape-key-tests :test-shape-key-consumers)
             (uiop:symbol-call :topology-shape-key-tests :test-cp-amplitude-boundaries)))

(defsystem "topology/training-context-tests"
  :description "Grammar and reference checks for training-context declarations"
  :depends-on ("topology")
  :components ((:file "tests/training-context"))
  :perform (test-op (operation component)
             (declare (ignore operation component))
             (uiop:symbol-call :cl-user :test-training-context-parser)))

(defsystem "topology/documentation"
  :description "Documentation for topology system"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:topology :mgl-pax)
  :serial t
  :components ((:file "documentation")))
