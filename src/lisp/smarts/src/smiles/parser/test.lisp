;;;(cl:in-package #:language.smiles.parser)
(in-package :language.smiles.parser)

(do-external-symbols (s '#:language.smiles.parser)
  (format t "~36A -> ~A~%" s (not (null (esrap:find-rule s)))))

(parse "B1=C(-P)-P1" 'list)

(parse "B1=C(-[30Al+5@OH3])-[P--]1" 'list)

(parse "C1=CC=CC=C1" 'list)

#+not-sure-how-this-works (parse "C:1:C:C:C:C:C1" 'list)

(parse "n1c[nH]cc1" 'list)

(parse "[CH3:1][C:2](=[O:3])[O-:4].[N+:5]" 'list :atom-maps? t)
