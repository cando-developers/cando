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

(defpackage #:chemdraw
  (:use :cl)
  (:export #:load-only-aggregate
           )
  )

(defpackage #:cando
  (:use #:common-lisp #:core)
  (:export
   #:from-string
   #:as-string
   #:scramble-positions
   #:jostle
   #:gather-stereocenters
   #:set-stereoisomer
   #:set-stereoisomer-func
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
   #:save-cando
   #:load-cando
   #:dump-atoms
   #:load-cdxml-aggregate
   #:load-psf-pdb
   #:remove-overlaps
   #:merge-into-one-aggregate
   #:constitution-from-residue))

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


(defpackage #:energy
  (:use #:common-lisp)
  (:export
   #:minimize
   #:setup-amber
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


(defpackage #:cando-user
  (:shadowing-import-from :chem "ATOM")
  (:shadowing-import-from :geom "BOUNDING-BOX")
  (:shadowing-import-from :common-lisp "+" "-" "/" "*" ">" "<" ">=" "<=" "SQRT")
  (:shadowing-import-from :energy "MINIMIZE")
  (:shadowing-import-from :chem "SET-ELEMENT" "GET-ELEMENT" "SIZE")
  (:shadowing-import-from :cando "AS-STRING" "LOAD-MOL2" "SAVE-MOL2")
  (:import-from :core "QUIT")
  (:use :common-lisp :chem :geom :chemdraw :cando)
  )


