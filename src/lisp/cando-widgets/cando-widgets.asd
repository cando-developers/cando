;;; Copyright (c) 2020, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(asdf:defsystem "cando-widgets"
  :description "Common widgets for CANDO Jupyter"
  :version "0.0.1"
  :author ("Christian Schafmeister <chris.schaf@verizon.net>"
           "Tarn W. Burton <twburton@gmail.com>")
  :licence "GPL2"
  :depends-on (:common-lisp-jupyter
               :lparallel
               :netcdf
               :nglview-cl
               :sketch2d
               :static-vectors)
  :serial t
  :components ((:file "packages")
               (:file "show")
               (:file "nglview")
               (:file "molecule-select")
               (:file "pages")
               (:file "workspace-task-page")
               (:file "ngl-structure-viewer")))

