;;;; package.lisp --- Package definition for the smiles.parser module.
;;;;
;;;; Copyright (C) 2018 Jan Moringen
;;;;
;;;; Author: Jan Moringen <jmoringe@techfak.uni-bielefeld.de>

(cl:defpackage #:language.smiles.parser
  (:use
   #:cl
   #:alexandria
   #:esrap)

  (:local-nicknames
   (#:bp #:architecture.builder-protocol))

  (:shadow
   #:parse)

  ;; Variables
  (:export
   #:*atom-maps?*)

  ;; Entry point rule
  (:export
   #:smiles)

  ;; Chain rules
  (:export
   #:chain
   #:branch
   #:atom)

  ;; Atom rules
  (:export
   #:atom-symbol
   #:organic-atom-symbol
   #:aromatic-atom-symbol
   #:inorganic-atom-symbol)

  ;; Modified atom rules
  (:export
   #:modified-atom
   #:modified-atom-body)

  ;; Atom modifier rules
  (:export
   #:atom-weight
   #:atom-modifier

   #:hydrogen-count
   #:charge
   #:chirality

   #:atom-map-class)

  ;; Bond rules
  (:export
   #:bond

   #:bond-none
   #:bond-single
   #:bond-double
   #:bond-quadruple
   #:bond-aromatic
   #:bond-up
   #:bond-down)

  ;; Parser protocol
  (:export
   #:parse))
