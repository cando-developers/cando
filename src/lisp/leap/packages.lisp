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

(defpackage #:leap.antechamber-type-definition-parser
  (:use #:common-lisp #:esrap #:parser.common-rules #:PARSER.COMMON-RULES.OPERATORS)
  (:export
   #:type-name
   #:type-name/s
   #:type-name/?s
   #:antechamber-line
   #:read-antechamber-type-rules))

(defpackage #:leap.parser
  (:use #:cl #:alexandria #:esrap #:parser.common-rules #:PARSER.COMMON-RULES.OPERATORS)
  (:export
   #:*function-names/alist*
   #:parse-sub-matter
   #:leap)
  (:documentation
   "Embryonic parser for the leap file format.  Written by Jan Moringen"))

(defpackage #:leap.interpreter
  (:use #:cl #:alexandria #:esrap #:parser.common-rules)
  (:export #:evaluate
           #:*leap-env*)
  (:documentation
   "Embryonic parser for the leap file format.  Written by Jan Moringen"))

(defpackage #:leap.atom-types
  (:use #:common-lisp)
  (:export
   #:add-one-atom-type
   ))


(defpackage #:leap.sequence
  (:use #:common-lisp)
  (:export
   #:sequence
   ))

(defpackage #:leap.cmap
  (:use #:common-lisp)
  (:shadow #:count)
  (:export
   #:cmap
   ))

(defpackage #:leap.commands
  (:use #:common-lisp)
  (:export
   #:create-atom
   #:leap
   #:leap-set
   #:leap-repl
   #:leap-repl-then-exit
   #:load-sketch
   ))

(defpackage #:leap.topology
  (:use #:common-lisp)
  (:export
   #:write-topology
   #:save-amber-parm-format
   #:generate-aggregate-for-energy-function
   #:amber-topology-coord-pair
   #:load-amber-topology-restart-pair
   #:load-amber-topology-trajectory-pair
   #:change-coordinate-file
   #:cell-lengths
   #:number-of-frames
   #:topology-filename
   #:coordinate-filename
   #:coordinate-netcdf
   #:energy-function
   #:aggregate
   #:number-of-atoms
   #:read-frame-into-atoms
   ))

(defpackage #:leap.cando-database
  (:use #:common-lisp)
  (:export 
   ))

(defpackage #:leap.pdb
  (:use #:common-lisp)
  (:export
   #:scan-pdb
   #:scan-pdb-stream
   #:load-pdb
   #:load-pdb-stream
   #:sequences
   #:build-sequence
   #:add-pdb-atom-map
   #:add-pdb-res-map
   #:pdb-read-error
   #:messages
  ))


(defpackage #:leap.off
  (:use #:common-lisp)
  (:export
   #:bug
   #:read-off-data-block
   #:name
   #:type
   #:fields
   #:read-off-unit-lib
   #:read-off-lib
   #:load-off)
  )

(defpackage #:leap.core
  (:use #:common-lisp)
  (:export
   #:with-path
   #:ensure-path
   #:nonbond-force-field-component
   #:evaluate
   #:*leap-env*
   #:lookup-variable
   #:all-variable-names
   #:function-lookup
   #:clear-path
   #:add-path
   #:search-path
   #:residue-name
   #:chain-position
   #:merged-force-field
   #:force-fields
   #:add-combined-force-field
   #:register-variable #:lookup-variable
   #:convert-leap-case
   #:clear-force-field
   #:add-force-field-or-modification
   #:find-force-field
   #:*force-fields*)
  )


(defpackage #:leap.solvate
  (:use #:common-lisp)
  (:export
   #:tool-solvate-and-shell
   #:tool-solvate-in-sphare)
  )

(defpackage #:leap.add-ions
  (:use #:common-lisp)
  (:export
   #:add-ions
   #:add-ions-2
   #:add-ions-rand)
  )

(defpackage #:leap.align-axes
  (:use #:common-lisp)
  (:export
   #:tool-orient-principle-axis-along-coordinate-axis)
  )

(defpackage #:leap.set-box
  (:use #:common-lisp)
  (:export
   #:set-box)
  )


(defpackage #:leap
  (:shadowing-import-from :chem "ATOM")
  (:shadowing-import-from :common-lisp "+" "-" "/" "*" ">" "<" ">=" "<=" "SQRT")
  (:shadowing-import-from :energy "MINIMIZE")
  (:shadowing-import-from :chem "SET-ELEMENT" "GET-ELEMENT" "SIZE")
  (:shadowing-import-from :cando "AS-STRING" "LOAD-MOL2" "SAVE-MOL2")
  (:import-from :core #:quit)
  (:import-from :leap.core
                #:add-path
                #:clear-force-field
                #:force-fields
                #:*force-fields*)
  (:import-from :leap.off
                #:load-off)
;;;  (:import-from :cando-utility #:mkdir #:set-current-directory #:current-directory #:directory-files)
  (:import-from :inet #:download-pdb)
  (:import-from :leap.pdb
                #:load-pdb
                #:add-pdb-res-map
                #:add-pdb-atom-map
                )
  (:import-from :leap.commands
                #:load-sketch
                #:leap-set)
  (:export
   #:load-sketch
   #:load-off
   #:process-command-line-options
   #:setup-amber-paths
   #:setup-default-paths
   #:list-force-fields
   #:load-smirnoff-params
   #:show-paths
   #:source
   #:set-force-field
   #:*amber-system*
   #:setup-amber
   #:save-amber-parm
   #:assign-atom-types
   #:clear-force-field
   #:desc
   #:add-path
   #:force-fields
   #:add-pdb-res-map
   #:add-pdb-atom-map
   #:atom
   #:bounding-box
   #:minimize
   #:as-string
   #:solvate-box
   #:solvate-oct
   #:solvate-shell
   #:solvate-cap
   #:quit
;;;   #:mkdir #:set-current-directory #:current-directory #:directory-files
   #:download-pdb
;   #:load-pdb
;   #:scan-pdb
   #:list-objects
   #:object
   #:load-amber-params
   #:load-atom-type-rules
   #:add-gaff
   #:easy-gaff)
  (:use :common-lisp #+(or):chem :geom))

(eval-when (:execute :load-toplevel)
  (loop for lib-name in '(:leap.set-box :leap.add-ions :leap.solvate
                          :leap.topology)
        for lib = (find-package lib-name)
        do (do-external-symbols (sym lib)
             (import sym :leap)
             (export sym :leap))))
