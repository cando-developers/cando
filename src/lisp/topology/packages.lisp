(cl:in-package #:common-lisp-user)

(defpackage #:topology
  (:use #:common-lisp)
  (:nicknames #:ts)
  (:export
   #:constitution-atom-index
   #:atom-name
   #:element
   #:bonds
   #:index
   #:to-atom-index
   #:order
   #:children
   #:name
   #:atom-names
   #:root-node
   #:constitution-from-graph
   #:stereo-information
   #:plugs
   #:name
   #:configurations
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
   #:joint-template
   #:bonded-joint-template
   #:in-plug-bonded-joint-template
   #:complex-bonded-joint-template
   #:jump-joint-template
   #:in-plug
   #:name
   #:id
   #:children
   #:constitution
   #:atom-name
   #:parent
   #:input-stub-joints
   #:plugs
   #:add-child
   #:sibling
   #:aggregate

   #:seen-fragment-internals
   #:good-fragment-internals
   #:load-fragment-pool
   #:fragment-pool
   #:fragment-internals
   #:internal
   #:jump-internal
   #:bonded-internal
   #:complex-bonded-internal
   #:fragments
   #:internals
   #:bond
   #:angle
   #:dihedral
   #:save-fragment-pool
   #:load-fragment-pool
   #:fragment-pool-map
   #:monomer-context-to-fragment-pool
   #:assemble-fragment-pool-map
   #:dump-fragment-internals
   #:index
   #:total-count
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
   #:oligomer
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
   #:ataggregate
   #:atmolecules
   #:atresidues
   #:joints
   #:zero-all-atom-tree-external-coordinates
   #:build-all-atom-tree-external-coordinates
   #:copy-joint-positions-into-atoms
   #:copy-atom-positions-into-joints
   

   #:make-joint-tree
   #:joint-tree
   #:children
   #:root
   #:atom-id-map
   #:make-conformation
   #:make-oligomer
   #:topologys-in-oligomer-space

   #:in-plug-name-p
   #:out-plug-name-p
   #:coupling-name
   
   ;;;
   #:update-joint-tree-internal-coordinates
   #:build-all-atom-tree-external-coordinates
   #:define-part
   #:*topology-groups*
   #:restraints))


(defpackage #:foldamer
  (:use #:common-lisp)
  (:nicknames #:fd)
  (:export
   #:define-foldamer
   #:parse
   #:oligomer-space
   #:generate-training-oligomers
   #:find-oligomer-for-monomer-context
   #:build-trainer
   #:prepare-to-build-trainer
   #:assemble-fragment-pool-map
   #:foldamer-monomer-context
   ))

(defpackage #:monomer-context
  (:use #:common-lisp)
  (:export
   #:parse
   #:matches
   #:match-as-string
   #:match
   ))

