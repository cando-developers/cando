(cl:in-package #:common-lisp-user)


(defpackage #:foldamer
  (:use #:common-lisp)
  (:nicknames #:fd)
  (:export
   #:define-foldamer
   #:parse
   #:oligomer-space
   #:unused-trainer-contexts
   #:maybe-remove-unused-trainers
   #:generate-training-oligomers
   #:find-oligomer-for-monomer-context
   #:build-trainer
   #:prepare-to-build-trainer
   #:foldamer-monomer-context
   #:calculate-files
   #:verify-foldamer-describes-oligomer-space
   #:valid-trainer-contexts
   #:assemble-fragment-conformations-map
   ))
