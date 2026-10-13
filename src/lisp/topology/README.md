# topology

Resolve a monomer's selected rotamer through the oligomer-shape:

```lisp
(topology:lookup-rotamer oligomer-shape monomer)
(topology:lookup-rotamer oligomer-shape monomer rotamer-index)
```

Ordinary indexes address the monomer-context's catalogue vector directly.
`topology:+optimized-rotamer-index+` (`#xffffffff`, the final unsigned 32-bit
index) selects an ordinary `topology:rotamer` stored in the shape's
`optimized-rotamers` hash table. A missing embedded rotamer or an invalid
catalogue index signals an error.

`topology:oligomer-shape-from-coordinates` captures the final coordinates in a
private copy of a shape. Its rotamer-shapes initially select the embedded
geometry. Fixed residue-shapes retain their stored Cartesian positions.
`catalogue-shape-keys` preserves the inherited sidechain classification, so
`topology:catalogue-shape-key` can resolve allowed catalogue indexes even when
backbone selections use optimized geometry. A saved `nil` key is valid.

Captured shapes permit catalogue sidechains and the embedded optimized
sidechain, with the optimized choice last in each allowed-index vector.
Their permissible backbones contain only the optimized choice, supporting
fixed-backbone sampling. Copying a shape copies its mutable selections and
optional maps; the optimized rotamers and shared catalogue remain immutable.
Serialization stores only the embedded geometry and inherited keys. Ordinary
shapes omit these optional slots and older serialized shapes load with `nil`.
