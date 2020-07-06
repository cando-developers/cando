;;;; protocol.lisp --- Protocol provided by the smarts.parser module.
;;;;
;;;; Copyright (C) 2018 Jan Moringen
;;;;
;;;; Author: Jan Moringen <jmoringe@techfak.uni-bielefeld.de>

(cl:in-package #:language.smarts.parser)

(defgeneric parse (input builder)
  (:documentation
   "Parse INPUT as a SMARTS expression, build result using BUILDER."))

(defmethod parse ((input string) (builder t))
  (let ((smiles:*atom-maps?* t))
    (bp:with-builder (builder)
      (esrap:parse 'smarts input))))


(defun parse-smirks (input builder)
  (let ((smiles:*atom-maps?* t))
    (bp:with-builder (builder)
      (let ((cons-result (esrap:parse 'smirks input)))
        (values (car cons-result) (cdr cons-result))))))
