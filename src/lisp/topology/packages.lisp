(cl:in-package #:common-lisp-user)

(defpackage #:serial
  (:use #:cl)
  (:export #:serializable))


(defpackage #:topology
  (:use #:common-lisp)
  (:nicknames #:ts)
  (:shadow #:type #:debug)
  (:export
   #:monomer-shape-locus
   #:assembler-check-joints
   #:make-topology-from-residue
   #:topologyp
   #:constitution-atom-index
   #:atom-name
   #:element
   #:root-map
   #:bonds
   #:trainer-index
   #:to-atom-index
   #:order
   #:children
   #:name
   #:plug-bonds
   #:atom-name
   #:bond-order
   #:plug-bond
   #:root-node
   #:constitution-from-graph
   #:stereo-information
   #:stereoisomer
   #:make-in-plug
   #:make-out-plug
   #:add-plug
   #:plugs
   #:name
   #:configurations
   #:make-constitution-atoms-from-residue
   #:topology-list
   #:property-list
   #:topology
   #:in-plug
   #:has-plug-named
   #:plug-named
   #:has-in-coupling-p
   #:define-topology
   #:constitution
   #:constitution-atoms
   #:stereoisomer-name
   #:stereoisomer-atoms
   #:build-residue-single-name
   #:connect-residues
   #:is-in-plug-name
   #:other-plug-name
   #:find-in-plug
   #:atom-type
   #:joint-template
   #:bonded-joint-template
   #:in-plug-bonded-joint-template
   #:complex-bonded-joint-template
   #:jump-joint-template
   #:in-plug
   #:name
   #:id
   #:source-monomer
   #:target-monomer
   #:source-plug-name
   #:target-plug-name
   #:children
   #:constitution
   #:atom-name
   #:parent
   #:input-stub-joints
   #:plugs
   #:add-child
   #:sibling
   #:aggregate
   #:ordered-monomers

   #:seen-fragment-internals
   #:seen-index
   #:bad-fragment-internals
   #:context-rotamers
   #:merge-context-rotamers
   #:fragment-internals
   #:internal
   #:jump-internal
   #:bonded-internal
   #:complex-bonded-internal
   #:internals
   #:fragments
   #:foldamer-monomer-context
   #:foldamer-rotamers-database
   #:save-clusterable-context-rotamers
   #:load-clusterable-context-rotamers
   #:context-rotamers
   #:context-rotamers-map
   #:monomer-context-to-context-rotamers
   #:dump-fragment-internals
   #:index
   #:next-index
   #:bad-geometry-p
   
   #:make-in-plug-bonded-joint-template
   #:make-bonded-joint-template
   #:make-jump-joint-template
   #:make-complex-bonded-joint-template
   #:stereo-information
   #:build-internal-coordinate-joint-template-tree
   #:build-joint-template
   #:extract-prepare-topologys
   #:target-monomer
   #:define-cap
   #:build-molecule
   #:oligomer-space
   #:oligomers
   #:monomer
   #:directional-coupling
   #:out-plug-name
   #:in-plug-name
   #:monomers
   #:couplings
   #:number-of-sequences
   #:current-stereoisomer-name
   #:other-monomer
   #:goto-sequence
   #:monomer-positions
   #:molecule-index
   #:residue-index
   #:residue
   #:copy-externals
   #:ataggregate
   #:atmolecules
   #:atresidues
   #:joints
   #:zero-all-atom-tree-external-coordinates
   #:build-all-atom-tree-external-coordinates
   #:build-atresidue-atom-tree-external-coordinates
   #:oligomer-space-directional-coupling-iterator-factory

   #:make-joint-tree
   #:joint-tree
   #:children
   #:root
   #:atom-id-map
   #:make-assembler
   #:make-oligomer
   #:topologys-in-oligomer-space

   #:in-plug-name-p
   #:out-plug-name-p
   #:coupling-name

   ;;;
   #:update-ataggregate-joint-tree-internal-coordinates
   #:update-atmolecule-joint-tree-internal-coordinates
   #:update-atresidue-joint-tree-internal-coordinates
   #:build-all-atom-tree-external-coordinates
   #:define-part
   #:*topology-groups*
   #:restraints

   #:missing-fragment-match
   #:connected-rotamers-map
   #:rotamer-context-connections
   #:make-rotamer-context-connections
   #:set-rotamer-context-connections
   #:get-rotamer-context-connections
   #:map-rotamer-context-connections
   #:rotamer-context-connections-count

   #:oligomer-monomer-name-at-index
   #:oligomer-monomer-name-for-monomer
   #:missing-rotamer-context-connections-count
   #:matched-context-rotamers-summary

   #:fragment-match-key-before-monomer-context-index
   #:fragment-match-key-after-monomer-context-index
   #:missing-fragment-match-key-before-monomer-context-index
   #:missing-fragment-match-key-after-monomer-context-index
   #:monomer-contexts-vector
   #:focus-monomer-name

   #:build-shapes
   #:monomer-shapes
   #:monomer-shape-vector
   #:root-monomer
   #:random-rotamer-index
   #:make-oligomer-shape
   #:build-one-molecule-for-topology
   #:sketch-svg
   #:make-oligomer-space
   #:calculate-number-of-sequences
   #:energy-function
   #:copy-fragment-internals
   #:copy-internal
   #:deg-to-rad
   #:rad-to-deg
   #:radians-add
   #:radians-incf
   #:radians-sub
   #:angle-sub
   #:walk-joint-template
   #:properties
   #:degrees-add
   #:probability
   #:oligomer-force-field-name
   #:oligomer
   #:build-residue-for-isomer
   #:stereoisomer-atom
   #:atom-charge
   #:copy-oligomer
   #:monomer-indexes
   #:internals-values
   #:internals-types
   #:constitution-atom-named
   #:constitution-atom-names-with-dihedrals
   #:monomer-topology
   #:calculate-cluster-dihedral-names
   #:cluster-dihedral-names
   #:cluster-dihedral-vector
   #:cluster-dihedral-rad-vector
   #:convert-dihedral-rad-vector-to-xy-vector
   #:build-all-molecules
   #:energy
   #:monomer-context
   #:clusterable-context-rotamers
   #:fill-internals-from-oligomer-shape
   #:monomer-shape-map
   #:in-monomers
   #:out-monomers
   #:rotamer-index
   #:rotor3
   #:translation
   #:monomer-shape
   #:backbone-monomer-shape
   #:rotamer-shape
   #:coordinates
   #:joint-with-name
   #:residue-properties
   #:merge-clusterable-context-rotamers
   #:foldamer-name
   #:monomer-on-other-side
   #:in-coupling-plug-name
   #:make-shape-info
   #:foldamer
   #:shape-info
   #:shape-key
   #:shape-kind
   #:keys
   #:kind
   #:kind-keys
   #:monomer-shape-build-order
   #:assign-monomer-shapes
   #:make-rotamer-shape-connections
   #:lookup-rotamer-context-connections
   #:lookup-rotamer-shape-connections
   #:out-coupling-plug-names
   #:append-rotamer-shape-connections
   #:bin-dihedral-deg
   #:shape-key-dihedral-name-p
   #:*phi*
   #:*phi-1*
   #:*psi*
   #:*psi-1*
   #:+phi+
   #:+phi-1+
   #:+psi+
   #:+psi-1+
   #:rotamer-indexes
   #:make-backbone-rotamer-shape-connections
   #:make-sidechain-rotamer-shape-connections
   #:assign-backbone-allowed-rotamers
   #:assign-backbone-random-rotamer-indexes
   #:allowed-rotamer-indexes
   #:assign-backbone-shape-keys
   #:focus-monomer
   #:focused-assembler-atresidue
   #:assembler-atresidue
   #:dihedral-info
   #:dihedral-info-atom
   #:dihedral-info-external
   #:plug-name
   #:external-dihedral-name
   #:oligomer-containing-monomer
   #:delta-energy
   #:lookup-rotamers-for-context
   #:monomer-shape-kind
   #:only-other-out-coupling-plug-name
   #:sorted-build-order
   #:adjustment
   #:external-adjustment
   #:internal-adjustment
   #:apply-fragment-internals-to-atresidue
   #:extract-bond-angle-rad-dihedral-rad
   #:fill-joint-internals
   #:make-coordinates-for-assembler
   #:initialize-adjustment
   #:do-adjustment
   #:find-joint-for-atom
   #:adjust-all-atom-tree-external-coordinates
   #:external-adjust
   #:internal-adjust
   #:adjust-internals
   #:radians-limit
   #:degrees-limit
   #:degrees-sub
   #:build-all-atom-tree-external-coordinates-and-adjust
   #:rotamers-database
   #:to-origin
   #:from-origin
   #:oligomer-monomer-name
   #:lookup-dihedral-cache
   #:shape-key-values
   #:rotamer
   #:sidechain-rotamer
   #:rotamers
   #:sidechain-rotamers
   #:shape-key-to-index
   #:backbone-rotamers-base
   #:backbone-without-sidechain-rotamers
   #:backbone-with-sidechain-rotamers
   #:rotamers-database
   #:context-to-rotamers
   #:backbone-dihedral-cache-deg
   #:sidechain-monomer-stepper
   #:backbone-monomer-stepper
   #:is-oligomer-space-supported
   #:abstract-topology
   #:define-abstract-topology
   #:is-out-plug-name
   #:joint
   #:internal-planar-adjustment
   #:read-rotamers
   #:write-rotamers
   #:read-oligomer-shape-rotamers
   #:write-oligomer-shape-rotamers
   #:labeled-monomers
   #:at-position
   #:no-rotamers
   #:permissible-rotamers
   #:oligomer-shape
   #:permissible-sidechain-rotamers
   #:permissible-backbone-rotamers
   #:make-permissible-backbone-rotamers
   #:make-permissible-sidechain-rotamers
   #:max-rotamers
   #:zero-rotamers
   #:random-rotamers
   #:build-all-sidechains
   #:make-phi-psi
   #:first-rotamers
   #:number-of-rotamers
   #:goto-rotamer
   #:make-orientation
   #:make-orientations
   #:orientation
   #:assembler-aggregate
   #:make-training-assembler
   #:make-focused-training-assembler
   #:training-assembler
   #:dihedral-rad
   #:build-atom-tree-external-coordinates
   #:adjust-atom-tree-external-coordinates
   #:build-atom-tree-external-coordinates-and-adjust
   #:copy-all-joint-positions-into-atoms
   #:assembler
   #:assembler-base
   #:oligomer-shapes
   #:walk-atmolecule-joints
   #:build-all-oligomer-shapes-in-coordinates
   #:build-external-coordinates
   #:aggregate*
   #:make-monomer-subset

   #:update-internals-for-atresidue
   #:update-internals
   #:update-externals
   #:build-oligomer-shape-externals
   #:the-root-monomer
   #:build-all-externals
   #:lookup-phi-psi-monomer-shape
   #:make-permissible-backbone-monomer-indexes
   #:goto-permissible-monomer-sequence
   #:number-of-permissible-monomer-sequences
   #:make-permissible-sidechain-monomer-indexes
   #:write-oligomer-monomer-indexes
   #:build-initial-oligomer-shape-externals
   #:allowed-index-ends
   #:monomer-index-loci
   #:copy-orientation
   #:make-coordinates-for-number-of-atoms
   #:analyze-oligomer-shape
   #:analyze-assembler
   #:monomer-subset-with-assembler
   #:monomer-set-from-permissible-rotamers
   #:analyze-atresidue
   #:monomers-of-shape-kind
   #:make-monomer-subset-with-assembler
   #:copy-oligomer-shape
   #:monomer-index
   #:copy-monomer-shape
   #:monomer-shape-info-vector
   #:apply-monomer-shape-to-atresidue-internals
   #:internal-adjustments
   #:adjustments
   #:enumerated-rotamer-state
   #:do-atresidue-residue
   #:do-joint-atom
   #:copy-xyz-joint-pos-to-atom-positions
   #:copy-atom-positions-to-xyz-joint-pos
   #:merge-internals-relative-to-heavy-atoms
   #:extract-internals
   #:write-internals
   #:ensure-oligomer-shape-is-consistent-with-permissible-rotamers
   #:find-named-fragment-line-segments
   #:cluster-dihedral-line-segments
   #:find-named-fragment-canonical-transform
   #:backbone-without-sidechain-rotamer
   #:backbone-with-sidechain-rotamer
   #:make-backbone-without-sidechain-rotamer
   #:make-backbone-with-sidechain-rotamer
   #:make-sidechain-rotamer
   #:make-backbone-with-sidechain-rotamer-from-fragment-internals
   #:make-sidechain-rotamer-from-fragment-internals
   #:make-backbone-without-sidechain-rotamer-from-fragment-internals
   #:make-fragment-internals-with-shape-key-from-fragment-internals
   #:lookup-backbone-shape-key
   #:rotamer-vector
   #:verify-oligomer-space
   #:ring-coupling
   #:plug1
   #:plug2
   #:monomer1
   #:monomer2
   #:allowed-rotamer-indexes
   #:copy-joint-positions-into-atoms
   #:permissible-rotamer-vector
   #:make-permissible-sidechain-rotamer-for-monomer
   #:parse-atom-graph
   #:molecule-from-sexp
   #:set-rotamers-state
   #:set-rotamers
   #:apply-rotamers
   #:make-rotamer-indexes
   #:make-indexes-into-rotamer-indexes
   #:rotamer-indexes
   #:indexes-into-rotamer-indexes
   #:ensure-permissible-rotamers-equal
   #:with-orientation
   #:ensure-complete-orientations
   #:make-monomer-indexes
   #:first-monomers
   #:adjust
   #:make-rotamer-shape-convenience
   #:make-oligomer-shape-of-rotamer-shapes
   #:residue-shape
   #:assembler-dump-internals
   #:randomly-mutate-sidechain-rotamer-index
   #:randomly-mutate-backbone-rotamer-index))

(defpackage #:topology.dag
  (:use #:common-lisp)
  (:nicknames #:tsd)
  (:export
   #:make-topology-from-residue
   #:parse-oligomer-space-labeled-dag
   #:parse-oligomer-space-dag
   #:validate-dag
   #:ring-cap-node
   #:label))

(defpackage #:monomer-context
  (:use #:common-lisp)
  (:export
   #:parse-sexp
   #:matches
   #:match-as-symbol
   #:match
   #:match-iterator
   #:copy-specialized
   #:node-count))

(defpackage #:topology.graphviz
  (:use #:common-lisp)
  (:nicknames #:topg)
  (:export
   #:nodes
   #:node-label
   #:directed-edges
   #:directed-edge-from
   #:directed-edge-to
   #:directed-edge-label
   #:render-dag
   #:undirected-edges
   #:undirected-edge-from
   #:undirected-edge-to
   #:undirected-edge-label
   #:graph-label
   #:graph-name
   #:render-foldamer-joint-trees
   #:make-graph
   #:node-id
   #:dot-svg-foldamer-joint-trees))

