
(in-package :leap.antechamber-type-definition-parser)


(esrap:defrule skippable
    parser.common-rules:whitespace+)

(esrap:defrule skippable?
    parser.common-rules:whitespace*)

(parser.common-rules:defrule/s (atd)
  "ATD")

(parser.common-rules:defrule/s (ap-tag)
  (and (esrap:character-ranges (#\a #\z) (#\A #\Z))
       (* (or (esrap:character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9)))))
  (:text t))

(parser.common-rules:defrule/s (type-name)
    (and (esrap:character-ranges (#\a #\z) (#\A #\Z))
         (* (or (esrap:character-ranges (#\! #\~))))) ; entire ascii printable range
  (:text t))
(parser.common-rules:defrule/s (wild-card)
  "*")

(parser.common-rules:defrule/s (number)
    (+ (esrap:character-ranges (#\0 #\9)))
  (:lambda (s) (parse-integer (esrap:text s))))

(esrap:defrule element
    (or "He" "Li" "Be" "Ne" "Na" "Mg" "Al" "Si" "Cl" "Ar" "Ca" "Sc" "Ti" "Cr" "Mn" "Fe" "Co"
        "Ni" "Cu" "Zn" "Ga" "Ge" "As" "Se" "Br" "Kr" "Rb" "Sr" "Zr" "Nb" "Mo" "Tc" "Ru" "Rh"
        "Pd" "Ag" "Cd" "In" "Sn" "Sb" "Te" "Xe" "Cs" "Ba" "La" "Ce" "Pr" "Nd" "Pm" "Sm" "Eu"
        "Gd" "Tb" "Dy" "Ho" "Er" "Tm" "Yb" "Lu" "Hf" "Ta" "Re" "Os" "Ir" "Pt" "Au" "Hg" "Tl"
        "Pb" "Bi" "Po" "At" "Rn" "Fr" "Ra" "Ac" "Th" "Pa" "Np" "Pu" "Am" "Cm" "Bk" "Cf" "Es"
        "Fm" "Md" "No" "Lr" "H" "W" "U" "B" "C" "N" "O" "F" "P" "S" "K" "V" "I" "Y"))


(parser.common-rules:defrule/s (terminator)
  "&")


#|
(esrap:defrule antechamber-line.-.bond-definitions
    (and type-name residue-list number number number-or-* number-or-* bracketed-atomic-property chemical-environment bond-definitions ap-end-of-line)
  (:destructure (type-name residue-list atomic-number attached-atoms attached-hydrogens electron-withdrawing-atoms property environment bond-definitions)
                (let ((focus (chem:make-antechamber-focus-atom-match residue-list atomic-number
|#

(defmacro defrule-wild-card-or-number (name)
  (let ((option1-name (gensym))
        (option2-name (gensym))
        (num (gensym)))
  `(progn
     (esrap:defrule ,name
       (or ,option1-name ,option2-name))
     (esrap:defrule ,option1-name
       wild-card (:constant -1))
     (esrap:defrule ,option2-name
       number (:destructure (,num) ,num)))))

(defrule-wild-card-or-number focus-atomic-number)
(defrule-wild-card-or-number focus-number-of-attached-atoms)
(defrule-wild-card-or-number focus-number-of-attached-hydrogen-atoms)

(esrap:defrule bracketed-atomic-property
  (or bracketed-atomic-property.wild-card
      bracketed-atomic-property.log-op))

(esrap:defrule bracketed-atomic-property.wild-card
    (and wild-card)
  (:lambda (x)
    (declare (ignore x))
    (:constant t)))

(esrap:defrule bracketed-atomic-property.log-op
    (and "[" log-op "]")
  (:lambda (x)
    (second x)))

(esrap:defrule chemical-environment
  (or chemical-environment.wild-card
      chemical-environment.branches))

(esrap:defrule chemical-environment.wild-card
  wild-card
  (:lambda (s)
    (core:make-cxx-object 'chem:bond-list-match-node)))

(esrap:defrule chemical-environment.branches
  (and branches)
  (:destructure (b)))

(esrap:defrule bond-definitions
    (and ap-tag ":" ap-tag ":" ap-bond)
  (:destructure (tag1 tag2 bond)
                (chem:make-after-match-bond-test (intern tag1 :keyword) (intern tag2 :keyword) bond)))

(esrap:defrule log-op
  (or log-op.and
      log-op.or))

(esrap:defrule log-op.and
    (and log-op "," log-op)
  (:destructure (left-op and right-op)
                (declare (ignore and))
                (chem:make-logical :log-high-precedence-and left-op right-op)))

(esrap:defrule log-op.or
    (and log-op "." log-op)
  (:destructure (left-op or right-op)
                (declare (ignore or))
                (chem:make-logical :log-or left-op right-op)))

(esrap:defrule atomic-test
  (or atomic-test.bond
      atomic-test.bond
      atomic-test.number.bond
      atomic-test.bond.not-bonded-to-previous
      atomic-test.bond.bonded-to-previous
      atomic-test.antechamber-arlevel
      atomic-test.antechamber-ring-membership
      atomic-test.antechamber-ring-membership.number
      atomic-test.number.antechamber-ring-membership.number
      atomic-test.antechamber-no-ring-membership))

(esrap:defrule atomic-test.bond
    (and ap-bond)
  (:destructure (bond)
                (chem:make-atom-test-2arg :sapin-bond bond 1)))

(esrap:defrule atomic-test.number.bond
    (and number ap-bond)
  (:destructure (num bond)
                (chem:make-atom-test-2arg :sapin-bond bond num)))

(esrap:defrule atomic-test.bond.not-bonded-to-previous
    (and ap-bond ap-not-bonded-to-previous)
  (:destructure (bond not-bonded-to-previous)
                (declare (ignore not-bonded-to-previous))
                (chem:make-atom-test-1arg :sapnot-bonded-to-previous bond)))

(esrap:defrule atomic-test.bond.bonded-to-previous
    (and ap-bond ap-bonded-to-previous)
  (:destructure (bond bonded-to-previous)
                (declare (ignore bonded-to-previous))
                (chem:make-atom-test-1arg :sapbonded-to-previous bond)))

(esrap:defrule atomic-test.antechamber-arlevel
    (and ap-antechamber-arlevel)
  (:lambda (x)
    (chem:make-atom-test-1arg :sapar-level (first x))))

(esrap:defrule atomic-test.antechamber-ring-membership
    (and ap-antechamber-ring-membership)
  (:destructure (dummy num2)
                (declare (ignore dummy))
                (chem:make-atom-test-1arg :sapring-size 1)))

(esrap:defrule atomic-test.antechamber-ring-membership.number
    (and ap-antechamber-ring-membership number)
  (:destructure (dummy num2)
                (declare (ignore dummy))
                (chem:make-atom-test-1arg :sapring-size num2)))

(esrap:defrule atomic-test.number.antechamber-ring-membership.number
    (and number ap-antechamber-ring-membership number)
  (:destructure (num1 dummy num3)
                (declare (ignore dummy))
                (chem:make-atom-test-2args :sapring-size num3 num1)))

(esrap:defrule atomic-test.antechamber-no-ring-membership
    (and ap-antechamber-no-ring-membership)
  (:lambda (x)
    (chem:make-atom-test :sapno-ring)))

(esrap:defrule branches
    (and "(" branches-list ")")
  (:lambda (x)
    (second x)))

(esrap:defrule branches-list
  (or chain
      branches-list.and.chain))

(esrap:defrule branches-list.and.chain
    (and branches-list ap-operator-and chain)
  (:lambda (x)
    (chem:make-branch.left.right (first x) (third x))))

(esrap:defrule chain
  (or chain.chemical-environment-atom
      chemical-environment-atom.branches))

(esrap:defrule chain.chemical-environment-atom
  chemical-environment-atom
  (:lambda (x)
    (chem:make-chain.head (first x))))

(esrap:defrule chain.chemical-environment-atom.branches
    (and chemical-environment-atom branches)
  (:lambda (x)
    (chem:make-chain.head.tail (first x) (second x))))

(esrap:defrule chemical-environment
  (or chemical-environment-atom.element.number.bracketed-atom-property-or-null.tag-name-or-null
      chemical-environment-atom.element.bracketed-atom-property-or-null.tag-name-or-null))

(esrap:defrule chemical-environment-atom.element.number.bracketed-atom-property-or-null.tag-name-or-null
    (and element number bracketed-atom-property-or-null tag-name-or-null)
  (:destructure (element number property tag)
                (chem:make-antechamber-bond-test (intern element :keyword) number property (intern tag :keyword))))

(esrap:defrule chemical-environment-atom.element.bracketed-atom-property-or-null.tag-name-or-null
    (and element bracketed-atom-property-or-null tag-name-or-null)
  (:destructure (element property tag)
                (chem:make-antechamber-bond-test (intern element :keyword) -1 property (intern tag :keyword))))

(esrap:defrule element
  (or ap-element
      ap-element-wild-card)
  (:lambda (x)
    (first x)))

(esrap:defrule bracketed-atom-property-or-null
  (or bracketed-atom-property-or-null.null
      bracketed-atom-property-or-null.log-op))

(esrap:defrule bracketed-atom-property-or-null.null
  skippable?
  (:constant t))

(esrap:defrule bracketed-atom-property-or-null.log-op
    (and "[" log-op "]")
  (:lambda (x)
    (second x)))

(esrap:defrule tag-name-or-null
  (or tag-name-or-null.null
      tag-name-or-null.ap-tag))

(esrap:defrule tag-name-or-null.null
  (and)
  (:constant ""))

(esrap:defrule tag-name-or-null.ap-tag
    (and "<" ap-tag ">")
  (:lambda (x)
    (second x)))

;;; add more options for residue-list if we ever encounter a non * residue-list
(parser.common-rules:defrule/s (residue-list)
  wild-card)

(esrap:parse 'tag-name-or-null "")
