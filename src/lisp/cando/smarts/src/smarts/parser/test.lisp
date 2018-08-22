(cl:in-package #:language.smarts.parser)

;;; Symbols

(do-external-symbols (s '#:language.smarts.parser)
  (format t "~36A -> ~A~%" s (not (null (esrap:find-rule s)))))

;;; Test

(parse "[CH1D3:3]" 'list)

(parse "B1=C(-P)-P1" 'list)

(parse "B1=C(-[30Al+5@OH3])-[P--]1" 'list)

(parse "C1=CC=CC=C1" 'list)

#+not-sure-how-this-works (architecture.builder-protocol:with-builder ('list)
                            (esrap:parse 'smarts "C:1:C:C:C:C:C1"))

(parse "n1c[nH]cc1" 'list)

(parse "[c;C&na]" 'list)

(progn
  (terpri)
  (architecture.builder-protocol.print-tree:serialize
   'list (parse "[C&!5C@@;BrH3]-C+3/?[3Cl]~C[*]C[R7]" 'list) *standard-output*))

(progn
  (terpri)
  (architecture.builder-protocol.print-tree:serialize
   'list (parse "[$([OH2][C,S,P]=O),$([10O]1nnnc1)]" 'list) *standard-output*))
