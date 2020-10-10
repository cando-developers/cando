;;;
;;;    File: packages.lisp
;;;

;;; Copyright (c) 2016-2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

;;  
;; This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.

;; -^-
(cl:in-package #:common-lisp-user)

(defpackage #:molecule-graph
  (:use #:common-lisp)
  (:export
   #:atom-cluster
   #:atom-cluster-p
   #:atom-cluster-atoms
   #:compare-molecules
   #:element-match-with-frontier
   #:element-match
   #:element-and-name-match
   ))

(defpackage #:molecule-graph.max-clique
  (:use #:common-lisp)
  (:export
   #:compare-molecules
   #:largest-connected-cluster-of-equivalent-atoms
   #:compare-molecules-return-largest-connected-cluster-of-equivalent-atoms
   ))

(defpackage #:graphviz-molecule-graph
  (:use #:common-lisp)
  (:local-nicknames (#:mg #:molecule-graph))
  (:export
   #:draw-graph-to-file
   #:draw-graph-to-stream
   ))


