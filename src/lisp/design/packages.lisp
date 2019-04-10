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

(defpackage #:design.load
  (:use #:cl )
  (:export
   load-cdxml
   )
  (:documentation
   "Load cdxml files and build topologies for doing design"))

(defpackage #:design.joint-tree
  (:use #:cl )
  (:export
   name
   constitution
   plugs
   stereo-information
   build-internal-coordinate-joint-template-tree
   extract-prepare-topologys)
  (:documentation
   "Embryonic joint-tree code"))

(defpackage #:design
  (:use #:cl )
  (:export
   define-part
   make-oligomer
   design
   topologys
   cap-name-map
   trainers
   oligomer
   build-trainers
   save-design
   load-design
   augment-trainer-with-superposable-conformation-collection
   jostle-trainer
   )
  (:documentation
   "Embryonic joint-tree code"))


(defpackage #:design.graphviz-draw-joint-template
  (:use #:cl)
  (:export
   #:draw-joint-template))

(defpackage #:design.graphviz-draw-joint-tree
  (:use #:cl)
  (:export
   #:draw-joint-tree))

(defpackage #:design.graphviz-draw-fold-tree
  (:use #:cl)
  (:export
   #:draw-fold-tree))
