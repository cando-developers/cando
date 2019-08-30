;;;
;;;    File: packages.lisp
;;;

;; Open Source License
;; Copyright (c) 2016, Christian E. Schafmeister
;; Permission is hereby granted, free of charge, to any person obtaining a copy
;; of this software and associated documentation files (the "Software"), to deal
;; in the Software without restriction, including without limitation the rights
;; to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
;; copies of the Software, and to permit persons to whom the Software is
;; furnished to do so, subject to the following conditions:
;; The above copyright notice and this permission notice shall be included in
;; all copies or substantial portions of the Software.
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
;; IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
;; OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
;; THE SOFTWARE.
;;  
;; This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.

;; -^-
(cl:in-package #:common-lisp-user)


(defpackage #:geom
  (:use #:common-lisp)
  (:export
   #:center-of-geometry
   ))

(defpackage #:select
  (:use #:common-lisp)
  (:export
   #:atoms-inside-sphere
   #:atoms-outside-sphere
   #:atoms-with-property
   #:atoms-with-property-sorted
   #:atom-lookup
   #:first-atom-with-name
   #:atoms-with-unique-names
   #:atom-lookup
   ))

(defpackage #:antique
  (:use #:common-lisp)
  (:export
   #:rename-atoms-with-unique-short-names
   #:restore-original-name)
  )

(defpackage #:atom-tree
  (:use #:common-lisp)
  (:export
   #:atom-tree
   #:build-atom-tree-for-molecule
   #:build-atom-tree
   #:nodes
   #:build-atom-pos
   #:build-geometry-for-atom-tree
   #:write-geometry-to-atoms
   #:find-coupled-dihedrals
   #:assign-coupled-dihedral-angles
   #:extract-internal-coordinates-for-active-set
   )
  )

(defpackage #:linked-rings
  (:use #:common-lisp)
  (:export
   #:identify-all-ring-chains
   #:atom-set-from-linked-rings
   #:expanded-atom-set-from-linked-rings
   #:scramble-linked-rings
   #:coerce-to-atom-set
   )
  )

(defpackage #:fortran
  (:use #:common-lisp)
  (:shadow #:format #:write)
  (:export
   #:with-fortran-output-file
   #:with-fortran-input-file
   #:fformat
   #:fwrite
   #:end-line
   #:debug-on
   #:debug
   #:fortran-input-file-look-ahead
   #:fread-line
   #:fread-line-or-error
   #:parse-fortran-format-line
   #:parse-double-float
   #:fread-double-float-vector
   #:fread-vector
   ))


(defpackage #:energy
  (:use #:common-lisp)
  (:export
   #:minimize
   #:minimize-minimizer
   #:minimize-energy-function
   #:minimize-from-bad-geometry
   #:minimize-energy-function-from-bad-geometry
   #:dump-excluded-atom-list
   ))

(defpackage #:anchor
  (:use #:common-lisp)
  (:export
   #:circle-points
   #:zig-zag-cylinder-points
   #:top-bottom-cylinder-points
   #:anchor-atoms
   #:on-circle
   #:on-points
   #:clear-all-atom-restraints
   ))

(defpackage #:struct
  (:use #:common-lisp)
  (:export
   #:find-residue-that-contains-atom
   #:remove-atom
   ))

(defpackage #:antechamber
  (:use #:common-lisp)
  (:export
   #:load-antechamber))


(defpackage #:cando
  (:use #:common-lisp
        #:cando-utility
        )
  (:shadowing-import-from :chem "ATOM")
  (:import-from :ext #:chdir #:getcwd)
  (:export
   #:|:=|
   #:<-
   #:nfx
   #:rename-atoms-using-smirks
   #:register-topology
   #:lookup-topology
   #:register-entity
   #:lookup-entity
   #:protein
   #:ligand
   #:ion
   #:solvent
   #:molecule-type
   #:load-chem-draw
   #:load-chem-draw-aggregate
   #:maybe-join-molecules-in-aggregate
   #:maybe-split-molecules-in-aggregate
   #:from-string
   #:as-string
   #:scramble-positions
   #:jostle
   #:walk-topologys
   #:list-topologys
   #:dump-stereocenters
   #:calculate-all-stereochemistry
   #:stereocenters-sorted-by-name
   #:gather-stereocenters
   #:set-stereoisomer
   #:set-stereoisomer-func
   #:set-stereoisomer-mapping
   #:set-stereoisomer-using-number
   #:set-all-stereocenters-to
   #:number-of-stereoisomers
   #:optimize-structure
   #:configure-minimizer
   #:minimize-no-fail
   #:build-good-geometry-from-random
   #:build-molecules-from-atom-connectivity
   #:chimera
   #:bad-geometry-p
   #:indexed-pathname
   #:save-mol2
   #:load-mol2
   #:make-class-save-load
   #:save-cando
   #:load-cando
   #:dump-atoms
   #:load-cdxml-aggregate
   #:load-psf-pdb
   #:remove-overlaps
   #:merge-into-one-aggregate
   #:constitution-from-residue
   #:clear-all-atom-restraints
   #:anchor-atom
   #:calculate-hybridization
   #:build-unbuilt-hydrogens
   #:progress-advance
   #:progress-bar
   #:make-progress-bar
   #:progress-done
   #:home/
   #:make-simple-topology-from-residue
   ;;; Stuff for nglview demo
   #:unbuilt-atoms
   #:count-unbuilt-atoms
   #:simple-build-unbuilt-atoms
   ;;; directory utilities
   #:set-current-directory
   #:current-directory
   #:do-molecules
   #:do-residues
   #:do-atoms
   #:remove-molecules
   #:combine
   #:assign-solvent-molecules-using-residue-name
   #:identify-solvent-residue-name
   ))


(defpackage #:dynamics
  (:use #:common-lisp)
  (:export
   #:trajectory
   #:make-trajectory
   #:matter
   #:number-of-atoms
   #:simulation
   #:make-atomic-simulation
   #:scoring-function
   #:coordinates
   #:velocity
   #:forces
   #:temp-forces
   #:delta-t
   #:velocity-scale
   #:delta-t-over-mass
   #:current-time
   #:write-coordinates-back-to-matter
   #:accumulate-coordinates
   #:accumulated-coordinates
   #:velocity-verlet-step))

(defpackage #:infix
  (:use #:common-lisp)
  (:export #:infix
           #:definfix))
