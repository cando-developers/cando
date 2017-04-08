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

(defpackage #:fortran
  (:use #:common-lisp)
  (:shadow #:format #:write)
  (:export
   #:with-fortran-output-file
   #:format
   #:write
   #:end-line
   #:debug-on
   #:debug))

(defpackage #:leap.topology
  (:use #:common-lisp)
  (:export
   #:write-topology))

(defpackage #:leap.commands
  (:use #:common-lisp)
  (:export
   #:load-script
   #:load-amber-params
   ))

(defpackage #:leap.pdb
  (:use #:common-lisp)
  (:export
   #:scanpdb
   #:loadpdb
   #:add-pdb-atom-map
   #:add-pdb-res-map)
  )


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
  (:import-from :leap.pdb
                #:add-pdb-atom-map
                #:add-pdb-res-map)
  (:export
   #:add-path #:search-path
   #:residue-name
   #:chain-position
   #:register-variable #:lookup-variable
   #:add-force-field-or-modification
   #:*force-fields*
   #:*objects*
   #:add-pdb-atom-map
   #:add-pdb-res-map)
  )

(defpackage #:leap
  (:shadowing-import-from :chem "ATOM")
  (:shadowing-import-from :geom "BOUNDING-BOX")
  (:shadowing-import-from :common-lisp "+" "-" "/" "*" ">" "<" ">=" "<=" "SQRT")
  (:shadowing-import-from :energy "MINIMIZE")
  (:shadowing-import-from :chem "SET-ELEMENT" "GET-ELEMENT" "SIZE")
  (:shadowing-import-from :cando "AS-STRING" "LOAD-MOL2" "SAVE-MOL2")
  (:import-from :core "QUIT")
  (:import-from :cando-utility #:mkdir #:set-current-directory #:current-directory #:directory-files)
  (:import-from :inet #:download-pdb)
;;  (:import-from :leap.off #:load-off)
  (:import-from :leap.pdb
                #:add-pdb-res-map
                #:add-pdb-atom-map
                )
  (:export
   #:atom
   #:bounding-box
   #:+ #:- #:/ #:* #:< #:> #:>= #:<= #:sqrt
   #:minimize
   #:set-element #:get-element #:size
   #:as-string #:load-mol2 #:save-mol2
   #:quit
   #:mkdir #:set-current-directory #:current-directory #:directory-files
   #:download-pdb
;;;   #:load-off
   #:add-pdb-res-map
   #:add-pdb-atom-map
   #:load-amber-params
   #:load-off
   #:list-force-fields
   #:list-objects
   #:desc
   #:object
   )
  (:use :common-lisp :chem :geom :chemdraw :cando))
