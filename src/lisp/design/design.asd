(in-package :asdf-user)

(defsystem "design"
  :description "Design related code"
  :version "0.0.1"
  :author "Christian Schafmeister <chris.schaf@verizon.net>"
  :licence "LGPL-3.0"
  :depends-on (:cando :leap :charges :smarts :amber)
  :serial t
  :components
  ((:file "packages")
   (:file "kin");
   (:file "fold")
   (:file "joint-templates")
   (:file "joint-tree")
   (:file "debug-graphs")
   #+(or)(:file "fold-tree")
   (:file "conformation")
   (:file "load")
   (:file "oligomer")
   (:file "design")
   (:file "trainers")
   #+(or)(:file "graphviz-draw-joint-template")
   #+(or)(:file "graphviz-draw-joint-tree")
   ))
