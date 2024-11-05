(in-package :topology)

(mgl-pax:defsection @topology-manual (:title "Topology manual")
  "Topology is a library for defining topology objects, assembling them into oligomer-space, oligomer, oligomer-shape and assemblers"
  (define-topology macro)
  )



(mgl-pax:defsection @topology/assembler (:title "Assemblers")
  "An ASSEMBLER is a complex object that calculates external coordinates from
an atom-tree of internal coordinates and sets the internal coordinates from
OLIGOMER-SHAPE objects."
  (orientation class)
  (assembler class)
  (make-coordinates-for-assembler function)
  (update-externals function)
  (update-internals function)
  )
