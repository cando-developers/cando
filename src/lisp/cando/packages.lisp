(cl:in-package #:common-lisp-user)

(defpackage #:chem
  (:use #:common-lisp)
  (:export
   )
  )


(defpackage #:cando
  (:use #:common-lisp #:core)
  (:export
   #:from-string
   #:as-string
   #:scramble-positions
   #:gather-stereocenters
   #:set-stereoisomer
   #:build-good-geometry-from-random
   #:chimera
   #:bad-geometry-p
   #:indexed-pathname
   #:save-mol2
   ))

