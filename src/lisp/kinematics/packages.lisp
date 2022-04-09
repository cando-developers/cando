(cl:in-package #:common-lisp-user)

(defpackage #:kin.work
  (:use #:cl )
  (:shadow #:atom)
  (:export
   #:*build-conformation*
   #:build-conformation
   #:!bad-geometry-p
   ))


(defpackage #:kin.fragment
  (:use #:cl )
  (:shadow #:atom)
  (:export
   #:*fragments*
   #:*cap-map*
   #:atoms
   #:name
   #:bonds
   #:target
   #:single-bond
   #:double-bond
   #:ring-single-bond
   #:ring-double-bond
   #:extract-prepare-topology
   )
  (:documentation
   "Load fragment sexp files and build topologies"))


(defpackage #:kin.tasks
  (:use #:cl )
  (:shadow #:atom #:directory)
  (:export
   #:sanitize-filename
   #:relative-node-pathname
   #:project-node-pathname
   #:with-project
   #:make-grammar-node
   )
  (:documentation
   "Load fragment sexp files and build topologies"))


(eval-when (:compile-toplevel :load-toplevel :execute)
  (let ((package (find-package :kin)))
    (export
     (mapcar
      (lambda (name) (intern (string name) package))
      '(
        #:trainer-node-pathnames
        #:apply-rotamer-to-atresidue
        #:out-plug
        #:plug-name
        #:smarts
        #:cap
        #:in-plug
        #:make-out-plug
        #:make-in-plug
        #:atmatter
        #:ataggregate
        #:atmolecules
        #:atmolecule
        #:root-atresidue
        #:atresidues
        #:atresidue
        #:residue
        #:parent
        #:parent-plug-name
        #:id
        #:children
        #:stereoisomer-name
        #:topology
        #:conformation-index
        #:focus-residue-sequence-number
        #:joints
        #:put-joint
        #:walk-ataggregate-joints
        #:update-internal-coordinates
        #:zero-external-coordinates
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
        #:stereo-information
        #:build-internal-coordinate-joint-template-tree
        #:extract-prepare-topologys
        #:build-internal-coordinate-joint-template-tree-spanning-tree
        #:build-stereoisomer
        #:joint-tree children
        #:root
        #:atom-id-map
        #:residue-graph
        #:matter-graph
        #:joint-template-graph
        #:aggregate
        #:ataggregate
        #:atmolecule-aref
        #:atresidue-aref
        #:atjoint-aref
        #:joint-tree
        #:make-conformation
        #:walk-atoms-joints
        #:copy-atom-positions-into-joints
        #:copy-joint-positions-into-atoms
        #:build-all-atom-tree-external-coordinates
        #:update-joint-tree-internal-coordinates
        #:zero-all-atom-tree-external-coordinates
        #:conformation
        #:make-conformation
        #:build-molecule-using-oligomer
        #:load-cdxml
        #:configure
        #:load-fragments
        #:render-fragment
        #:configure
        #:make-design
        #:grammar
        #:cap-name-map
        #:build-molecule
        #:*complex-plugs*
        #:complex-plug
        #:name
        #:smarts
        #:compiled-smarts
        #:complex-plug-or-nil
        #:dkp
        #:define-part
        #:make-oligomer
        #:context
        #:design
        #:chemdraw
        #:complex-plugs
        #:topologys
        #:cap-name-map
        #:trainers
        #:oligomer
        #:canonical-sequence
        #:build-trainers
        #:save-design
        #:load-design
        #:augment-trainer-with-superposable-conformation-collection
        #:jostle-trainer
        #:build-charge-solvate
        #:build-all-trainers
        #:canonical-sequence
        #:joint-tree
        #:root
        #:build-training-oligomers
        #:trainer-filename
        #:maybe-save-trainer
        #:load-trainer
        #:save-rotamer-list
        #:load-rotamer-list
        #:rotamer-list
        #:rotamer-library
        #:extract-rotamer
        #:entries
        #:internal-match
        #:internals-match
        ))
     package)))
        

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
