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

(defpackage #:amber.params
  (:use #:common-lisp)
  (:export
   ))

(defpackage #:amber
  (:use #:common-lisp)
  (:shadow #:stream)
  (:export
   #:simulation
   #:add-aggregates
   #:add-aggregate
   #:save-simulation
   #:start-simulation
   #:convert-to-script
   #:job
   #:job-file
   #:node-file
   #:script-mixin
   #:substitutions
   #:node-pathname
   #:make-node-file
   #:argument
   #:node
   #:option
   #:definers
   #:users
   #:inputs
   #:outputs
   #:script
   #:name
   #:filename
   #:jupyter-job
   #:simple-jupyter-job
   #:minimize
   #:heat
   #:pressurize
   #:dynamics
   #:mdcrd
   #:generate-all-code
   #:connect-graph
   #:arguments
   #:updated
   #:make-dashboard
   #:current-info-file
   #:panels
   #:update
   #:job-name
   #:panel-summary-percent-done
   #:summary
   #:target
   #:action
   #:panel-summary-estimated-time-remaining
   #:replica-exchange-equilibrate
   #:replica-exchange
   #:single-point-energy
   #:single-point-energy-debug))

(defpackage #:graphviz-amber
  (:use #:common-lisp)
  (:export
   #:draw-graph-stream
   ))


