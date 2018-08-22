;;;(cl:in-package #:language.smiles.parser)
(in-package :language.smiles.parser)

(defgeneric parse (input builder &key atom-maps?)
  (:documentation
   "Parse INPUT as a SMILES expression, build result using BUILDER.

    ATOM-MAPS? controls whether the […:N] construct may be used in
    INPUT to assign atom map classes to individual atoms."))

(defmethod parse ((input string) (builder t) &key atom-maps?)
  (let ((*atom-maps?* atom-maps?))
    (architecture.builder-protocol:with-builder (builder)
      (esrap:parse 'smiles input))))
