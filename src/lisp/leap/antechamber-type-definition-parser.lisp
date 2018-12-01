
(in-package :leap.antechamber-type-definition-parser)

;;; Utilities

;;; Lexical stuff and reserved words

(defvar *wild-dict* (core:make-cxx-object 'chem:wild-element-dict))
(defvar *next-chem-info-id*)


(esrap:defrule skippable
    (+ (or parser.common-rules:shell-style-comment
           parser.common-rules:whitespace)))

(esrap:defrule skippable?
    (* (or parser.common-rules:shell-style-comment
           parser.common-rules:whitespace)))

(parser.common-rules:defrule/s wildatom
  "WILDATOM")

(esrap:defrule wild-atom-line
    (and wildatom/s wild-element/s (+ wild-element/?s))
  (:destructure (wa wild-element elements)
                (declare (ignore wa))
                (list* wild-element elements)))
  
(esrap:defrule antechamber-line
    (and atd/s antechamber-line-no-amp "&")
  (:destructure (atd result ampersand)
                (declare (ignore atd ampersand))
                result))

(esrap:defrule antechamber-line-no-amp
    (or antechamber-line...bond-definitions
        antechamber-line...chemical-environment
        antechamber-line...bracketed-atomic-property
        antechamber-line...focus-number-of-electron-withdrawing-atoms
        antechamber-line...focus-number-of-attached-hydrogen-atoms
        antechamber-line...focus-number-of-attached-atoms
        antechamber-line...focus-atomic-number
        antechamber-line...residue-list
        antechamber-line...type-name))

(esrap:defrule antechamber-line...bond-definitions
    (and type-name/s residue-list/s focus-atomic-number focus-number-of-attached-atoms focus-number-of-attached-hydrogen-atoms focus-number-of-electron-withdrawing-atoms bracketed-atomic-property/s chemical-environment/s bond-definitions/s)
  (:destructure (type-name residue-list focus-atomic-number focus-number-of-attached-atoms focus-number-of-attached-hydrogen-atoms focus-number-of-electron-withdrawing-atoms bracketed-atomic-property chemical-environment bond-definitions)
                (let ((focus (core:make-cxx-object 'chem:antechamber-focus-atom-match
                                                   :id (incf *next-chem-info-id*)
                                                   :residue-names residue-list
                                                   :atomic-number focus-atomic-number
                                                   :attached-atoms focus-number-of-attached-atoms
                                                   :attached-hs focus-number-of-attached-hydrogen-atoms
                                                   :attached-ewgs focus-number-of-electron-withdrawing-atoms
                                                   :atomic-property bracketed-atomic-property)))
                  (core:make-cxx-object 'chem:antechamber-root
                                        :id (incf *next-chem-info-id*)
                                        :assign-type (intern type-name :keyword)
                                        :first-test focus
                                        :chain chemical-environment
                                        :after-match-tests bond-definitions
                                        :wild-dict *wild-dict*))))

(esrap:defrule antechamber-line...chemical-environment
    (and type-name/s residue-list/s focus-atomic-number focus-number-of-attached-atoms focus-number-of-attached-hydrogen-atoms focus-number-of-electron-withdrawing-atoms bracketed-atomic-property/s chemical-environment/?s)
  (:destructure (type-name residue-list focus-atomic-number focus-number-of-attached-atoms focus-number-of-attached-hydrogen-atoms focus-number-of-electron-withdrawing-atoms bracketed-atomic-property chemical-environment)
                (let ((focus (core:make-cxx-object 'chem:antechamber-focus-atom-match
                                                   :id (incf *next-chem-info-id*)
                                                   :residue-names residue-list
                                                   :atomic-number focus-atomic-number
                                                   :attached-atoms focus-number-of-attached-atoms
                                                   :attached-hs focus-number-of-attached-hydrogen-atoms
                                                   :attached-ewgs focus-number-of-electron-withdrawing-atoms
                                                   :atomic-property bracketed-atomic-property)))
                  (core:make-cxx-object 'chem:antechamber-root
                                        :id (incf *next-chem-info-id*)
                                        :assign-type (intern type-name :keyword)
                                        :first-test focus
                                        :chain chemical-environment
                                        :wild-dict *wild-dict*))))

(esrap:defrule antechamber-line...bracketed-atomic-property
    (and type-name/s residue-list/s focus-atomic-number focus-number-of-attached-atoms focus-number-of-attached-hydrogen-atoms focus-number-of-electron-withdrawing-atoms bracketed-atomic-property/?s)
  (:destructure (type-name residue-list focus-atomic-number focus-number-of-attached-atoms focus-number-of-attached-hydrogen-atoms focus-number-of-electron-withdrawing-atoms bracketed-atomic-property)
                (let ((focus (core:make-cxx-object 'chem:antechamber-focus-atom-match
                                                   :id (incf *next-chem-info-id*)
                                                   :residue-names residue-list
                                                   :atomic-number focus-atomic-number
                                                   :attached-atoms focus-number-of-attached-atoms
                                                   :attached-hs focus-number-of-attached-hydrogen-atoms
                                                   :attached-ewgs focus-number-of-electron-withdrawing-atoms
                                                   :atomic-property bracketed-atomic-property)))
                  (core:make-cxx-object 'chem:antechamber-root
                                        :id (incf *next-chem-info-id*)
                                        :assign-type (intern type-name :keyword)
                                        :first-test focus
                                        :wild-dict *wild-dict*))))

(esrap:defrule antechamber-line...focus-number-of-electron-withdrawing-atoms
    (and type-name/s residue-list/s focus-atomic-number focus-number-of-attached-atoms focus-number-of-attached-hydrogen-atoms focus-number-of-electron-withdrawing-atoms)
  (:destructure (type-name residue-list focus-atomic-number focus-number-of-attached-atoms focus-number-of-attached-hydrogen-atoms focus-number-of-electron-withdrawing-atoms)
                (let ((focus (core:make-cxx-object 'chem:antechamber-focus-atom-match
                                                   :id (incf *next-chem-info-id*)
                                                   :residue-names residue-list
                                                   :atomic-number focus-atomic-number
                                                   :attached-atoms focus-number-of-attached-atoms
                                                   :attached-hs focus-number-of-attached-hydrogen-atoms
                                                   :attached-ewgs focus-number-of-electron-withdrawing-atoms)))
                  (core:make-cxx-object 'chem:antechamber-root
                                        :id (incf *next-chem-info-id*)
                                        :assign-type (intern type-name :keyword)
                                        :first-test focus
                                        :wild-dict *wild-dict*))))

(esrap:defrule antechamber-line...focus-number-of-attached-hydrogen-atoms
    (and type-name/s residue-list/s focus-atomic-number focus-number-of-attached-atoms focus-number-of-attached-hydrogen-atoms)
  (:destructure (type-name residue-list focus-atomic-number focus-number-of-attached-atoms focus-number-of-attached-hydrogen-atoms)
                (let ((focus (core:make-cxx-object 'chem:antechamber-focus-atom-match
                                                   :id (incf *next-chem-info-id*)
                                                   :residue-names residue-list
                                                   :atomic-number focus-atomic-number
                                                   :attached-atoms focus-number-of-attached-atoms
                                                   :attached-hs focus-number-of-attached-hydrogen-atoms
                                                   :attached-ewgs -1)))
                  (core:make-cxx-object 'chem:antechamber-root
                                        :id (incf *next-chem-info-id*)
                                        :assign-type (intern type-name :keyword)
                                        :first-test focus
                                        :wild-dict *wild-dict*))))

(esrap:defrule antechamber-line...focus-number-of-attached-atoms
    (and type-name/s residue-list/s focus-atomic-number focus-number-of-attached-atoms)
  (:destructure (type-name residue-list focus-atomic-number focus-number-of-attached-atoms)
                (let ((focus (core:make-cxx-object 'chem:antechamber-focus-atom-match
                                                   :id (incf *next-chem-info-id*)
                                                   :residue-names residue-list
                                                   :atomic-number focus-atomic-number
                                                   :attached-atoms focus-number-of-attached-atoms
                                                   :attached-hs -1
                                                   :attached-ewgs -1)))
                  (core:make-cxx-object 'chem:antechamber-root
                                        :id (incf *next-chem-info-id*)
                                        :assign-type (intern type-name :keyword)
                                        :first-test focus
                                        :wild-dict *wild-dict*))))

(esrap:defrule antechamber-line...focus-atomic-number
    (and type-name/s residue-list/s focus-atomic-number)
  (:destructure (type-name residue-list focus-atomic-number)
                (let ((focus (core:make-cxx-object 'chem:antechamber-focus-atom-match
                                                   :id (incf *next-chem-info-id*)
                                                   :residue-names residue-list
                                                   :atomic-number focus-atomic-number
                                                   :attached-atoms -1
                                                   :attached-hs -1
                                                   :attached-ewgs -1)))
                  (core:make-cxx-object 'chem:antechamber-root
                                        :id (incf *next-chem-info-id*)
                                        :assign-type (intern type-name :keyword)
                                        :first-test focus
                                        :wild-dict *wild-dict*))))

(esrap:defrule antechamber-line...residue-list
    (and type-name/s residue-list/s focus-atomic-number)
  (:destructure (type-name residue-list focus-atomic-number)
                (let ((focus (core:make-cxx-object 'chem:antechamber-focus-atom-match
                                                   :id (incf *next-chem-info-id*)
                                                   :residue-names residue-list
                                                   :atomic-number focus-atomic-number)))
                  (core:make-cxx-object 'chem:antechamber-root
                                        :id (incf *next-chem-info-id*)
                                        :assign-type (intern type-name :keyword)
                                        :first-test focus
                                        :wild-dict *wild-dict*))))

(esrap:defrule antechamber-line...type-name
    (and type-name/s focus-atomic-number)
  (:destructure (type-name focus-atomic-number)
                (let ((focus (core:make-cxx-object 'chem:antechamber-focus-atom-match
                                                   :id (incf *next-chem-info-id*)
                                                   :atomic-number focus-atomic-number)))
                  (core:make-cxx-object 'chem:antechamber-root
                                        :id (incf *next-chem-info-id*)
                                        :assign-type (intern type-name :keyword)
                                        :first-test focus
                                        :wild-dict *wild-dict*))))


(defmacro defrule-wild-card-or-number (name)
  (let ((option1-name (gensym))
        (option2-name (gensym))
        (num (gensym)))
    `(progn
       (esrap:defrule ,name
           (or ,option1-name ,option2-name))
       (esrap:defrule ,option1-name
         wild-card/s (:constant -1))
       (esrap:defrule ,option2-name
         number/s (:lambda (,num) ,num)))))

(defrule-wild-card-or-number focus-atomic-number)
(defrule-wild-card-or-number focus-number-of-attached-atoms)
(defrule-wild-card-or-number focus-number-of-attached-hydrogen-atoms)
(defrule-wild-card-or-number focus-number-of-electron-withdrawing-atoms)

(parser.common-rules:defrule/s bracketed-atomic-property
    (or bracketed-atomic-property.wild-card
        bracketed-atomic-property.log-op))

(esrap:defrule bracketed-atomic-property.wild-card
    (and wild-card)
  (:lambda (x)
    (declare (ignore x))
    (core:make-cxx-object 'chem:logical
                          :id (incf *next-chem-info-id*)
                          :op :log-always-true)))

(esrap:defrule bracketed-atomic-property.log-op
    (and "[" log-op "]")
  (:lambda (x)
    (second x)))

(parser.common-rules:defrule/s chemical-environment
    (or chemical-environment.wild-card
        chemical-environment.branches))

(esrap:defrule chemical-environment.wild-card
  wild-card
  (:lambda (s)
    (core:make-cxx-object 'chem:bond-list-match-node
                          :id (incf *next-chem-info-id*))))

(esrap:defrule chemical-environment.branches
    (and branches)
  (:destructure (b)
                b))

(parser.common-rules:defrule/s bond-definitions
    (and ap-tag ":" ap-tag ":" ap-bond)
  (:destructure (tag1 tag2 bond)
                (core:make-cxx-object 'chem:after-match-bond-test
                                      :id (incf *next-chem-info-id*)
                                      :tag1 (intern tag1 :keyword)
                                      :tag2 (intern tag2 :keyword)
                                      :bond bond)))

#||
(esrap:defrule log-op.and
    (and log-op "," log-op)
  (:destructure (left-op and right-op)
                (declare (ignore and))
                (core:make-cxx-object 'chem:logical
                                      :op :log-high-precedence-and
                                      :left left-op
                                      :right right-op)))


(esrap:defrule log-op.or
    (and log-op "." log-op)
  (:destructure (left-op or right-op)
                (declare (ignore or))
                (core:make-cxx-object 'chem:logical
                                      :op :log-or
                                      :left left-op
                                      :right right-op)))
||#


(parser.common-rules:defrule/s atd
  "ATD")

(parser.common-rules:defrule/s ap-tag
    (and (esrap:character-ranges (#\a #\z) (#\A #\Z))
         (* (or (esrap:character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9)))))
  (:text t))

(parser.common-rules:defrule/s type-name
    (and (esrap:character-ranges (#\a #\z) (#\A #\Z))
         (* (or (esrap:character-ranges (#\! #\~))))) ; entire ascii printable range
  (:text t))

(parser.common-rules:defrule/s wild-card
  "*")

(parser.common-rules:defrule/s number
    (+ (esrap:character-ranges (#\0 #\9)))
  (:lambda (x)
    (parse-integer (esrap:text x))))

(esrap:defrule ap-element-2-char
    (or "He" "Li" "Be" "Ne" "Na" "Mg" "Al" "Si" "Cl" "Ar" "Ca" "Sc" "Ti" "Cr" "Mn" "Fe" "Co"
        "Ni" "Cu" "Zn" "Ga" "Ge" "As" "Se" "Br" "Kr" "Rb" "Sr" "Zr" "Nb" "Mo" "Tc" "Ru" "Rh"
        "Pd" "Ag" "Cd" "In" "Sn" "Sb" "Te" "Xe" "Cs" "Ba" "La" "Ce" "Pr" "Nd" "Pm" "Sm" "Eu"
        "Gd" "Tb" "Dy" "Ho" "Er" "Tm" "Yb" "Lu" "Hf" "Ta" "Re" "Os" "Ir" "Pt" "Au" "Hg" "Tl"
        "Pb" "Bi" "Po" "At" "Rn" "Fr" "Ra" "Ac" "Th" "Pa" "Np" "Pu" "Am" "Cm" "Bk" "Cf" "Es"
        "Fm" "Md" "No" "Lr"))

(esrap:defrule ap-element-1-char
    (or "H" "W" "U" "B" "C" "N" "O" "F" "P" "S" "K" "V" "I" "Y"))

(esrap:defrule ap-element
    (or ap-element-2-char
        ap-element-1-char))

(parser.common-rules:defrule/s element
    (or ap-element-2-char
        ap-element-1-char))


(esrap:defrule ap-wild-element
    (or "XA" "XB" "XC" "XD" "XE" "XF" "XG" "XH" "XI" "XJ" "XK" "XL"
        "XM" "XN" "XO" "XP" "XQ" "XR" "XS" "XT" "XU" "XV" "XW" "XX" "XY" "XZ"))

(esrap:defrule ap-element-wild
    (or ap-element-2-char
        ap-wild-element
        ap-element-1-char))

(parser.common-rules:defrule/s wild-element
    (or ap-element-2-char
        ap-wild-element
        ap-element-1-char))

(parser.common-rules:defrule/s terminator
  "&")

(esrap:defrule operator-and
  ","
  (:constant :and))

(esrap:defrule operator-or
  "."
  (:constant :or))


(esrap:defrule log-op
    (or log-op.and
        log-op.or
        atomic-test ))

#+(or)
(esrap:defrule log-op.atomic-test
  atomic-test
  (:lambda (test)
    (architecture.builder-protocol:node* (:unary-operator :op :log-identity)
                                         (1 :operand test))))

(parser.common-rules.operators:define-operator-rules
    (:skippable?-expression (* #\Space))
    (2 log-op.and      operator-and)
  (2 log-op.or       operator-or)
  atomic-test)

(esrap:defrule atomic-test
    (or atomic-test.number.bond
        atomic-test.bond.not-bonded-to-previous
        atomic-test.bond.bonded-to-previous
        atomic-test.bond
        atomic-test.antechamber-arlevel
        atomic-test.antechamber-ring-membership
        atomic-test.number.antechamber-ring-membership
        atomic-test.antechamber-no-ring-membership)
  (:lambda (value)
    (core:make-cxx-object 'chem:logical
                          :id (incf *next-chem-info-id*)
                          :op :log-identity :left value)))

(esrap:defrule atomic-test.bond
  ap-bond
  (:lambda (bond) (chem:make-atom-test-in-bond1 bond)))

(esrap:defrule atomic-test.number.bond
    (and number ap-bond)
  (:destructure (num bond) (chem:make-atom-test-in-bond bond num)))

(esrap:defrule atomic-test.bond.not-bonded-to-previous
    (and ap-bond ap-not-bonded-to-previous)
  (:destructure (bond not-bonded-to-previous)
                (declare (ignore not-bonded-to-previous))
                (chem:make-atom-test-not-bonded-to-previous bond)))

(esrap:defrule ap-not-bonded-to-previous
    "''")

(esrap:defrule atomic-test.bond.bonded-to-previous
    (and ap-bond ap-bonded-to-previous)
  (:destructure (bond bonded-to-previous)
                (declare (ignore bonded-to-previous))
                (chem:make-atom-test-bonded-to-previous bond)))

(esrap:defrule ap-bonded-to-previous
    "'")

(esrap:defrule atomic-test.antechamber-arlevel
    (and ap-antechamber-arlevel)
  (:lambda (x) (chem:make-atom-test-ar-level (first x))))

(esrap:defrule atomic-test.antechamber-ring-membership
  ap-antechamber-ring-membership
  (:lambda (size)
    (chem:make-atom-test-ring-size size)))

(esrap:defrule atomic-test.number.antechamber-ring-membership
    (and number ap-antechamber-ring-membership)
  (:destructure (num size) (chem:make-atom-test-ring-size2 num size)))

(esrap:defrule atomic-test.antechamber-no-ring-membership
    (and ap-antechamber-no-ring-membership)
  (:destructure (dummy) (declare (ignore dummy)) (chem:make-atom-test-no-ring-membership)))

(esrap:defrule branches
    (and "(" branches-list ")")
  (:function second))

(esrap:defrule branches-list
    (or branches-list.and.chain
        chain))

(esrap:defrule branches-list.and.chain
    (and branches-list operator-and chain)
  (:lambda (x)
    (chem:make-branch.left.right (first x) (third x))))

(esrap:defrule chain
    (or chain.chemical-environment-atom.branches
        chain.chemical-environment-atom))

(esrap:defrule chain.chemical-environment-atom
  chemical-environment-atom
  (:lambda (x)
    (chem:make-chain.head x)))

(esrap:defrule chain.chemical-environment-atom.branches
    (and chemical-environment-atom branches)
  (:lambda (x)
    (chem:make-chain.head.tail (first x) (second x))))

(esrap:defrule chemical-environment-atom
    (or chemical-environment-atom.element.number.bracketed-atom-property-or-null.tag-name-or-null
        chemical-environment-atom.element.bracketed-atom-property-or-null.tag-name-or-null))

(esrap:defrule chemical-environment-atom.element.number.bracketed-atom-property-or-null.tag-name-or-null
    (and ap-element-wild number bracketed-atom-property-or-null tag-name-or-null)
  (:destructure (element number property tag)
                (chem:make-antechamber-bond-test (intern element :keyword) number property (intern tag :keyword))))

(esrap:defrule chemical-environment-atom.element.bracketed-atom-property-or-null.tag-name-or-null
    (and ap-element-wild bracketed-atom-property-or-null tag-name-or-null)
  (:destructure (element property tag)
                (chem:make-antechamber-bond-test (intern element :keyword) -1 property (intern tag :keyword))))

(esrap:defrule element
    (or ap-element-wild
        ap-element-wild-card)
  (:function first))

(esrap:defrule bracketed-atom-property-or-null
    (or bracketed-atom-property-or-null.log-op
        bracketed-atom-property-or-null.null))

(esrap:defrule bracketed-atom-property-or-null.null
  skippable?
  (:lambda (x) (declare (ignore x)) (core:make-cxx-object 'chem:logical
                                                          :id (incf *next-chem-info-id*)
                                                          :op :log-always-true)))

(esrap:defrule bracketed-atom-property-or-null.log-op
    (and "[" log-op "]")
  (:function second))

(esrap:defrule tag-name-or-null
    (or tag-name-or-null.ap-tag
        tag-name-or-null.null))

(esrap:defrule tag-name-or-null.null
    (and)
  (:constant ""))

(esrap:defrule tag-name-or-null.ap-tag
    (and "<" ap-tag ">")
  (:function second))

;;; add more options for residue-list if we ever encounter a non * residue-list
(parser.common-rules:defrule/s residue-list
  wild-card
  (:lambda (x)
    (declare (ignore x))
    (core:make-cxx-object 'chem:residue-list)))

(defmacro defrules-simple (name &rest body)
  (let ((subrule-names (loop repeat (length body) collect (gensym))))
    `(progn
       (esrap:defrule ,name (or ,@subrule-names))
       ,@(mapcar (lambda (name parts)
                   `(esrap:defrule ,name ,(car parts) ,(list :constant (second parts))))
                 subrule-names body))))

(defrules-simple ap-bond
    ("SB" :sabsingle-bond)
  ("DB" :sabdouble-bond)
  ("DL" :sabdelocalized-bond)
  ("TB" :sabtriple-bond)
  ("AB" :sabaromatic-bond)
  ("sb" :sabsingle-or-aromatic-bond)
  ("db" :sabdouble-or-aromatic-bond)
  ("tb" :sabtriple-or-aromatic-bond))

(defrules-simple ap-antechamber-ring-membership
    ("RG3" 3)
  ("RG4" 4)
  ("RG5" 5)
  ("RG6" 6)
  ("RG7" 7)
  ("RG" 1)
  )

(esrap:defrule ap-antechamber-no-ring-membership
  "NR"
  (:constant "NR"))

(defrules-simple ap-antechamber-arlevel
    ("AR1" 1)
  ("AR2" 2)
  ("AR3" 3)
  ("AR4" 4)
  ("AR5" 5))

;;(esrap:parse 'chemical-environment "C3")


;;; builder protocol

(defclass builder () ())

(defmethod architecture.builder-protocol:make-node ((builder builder) (kind (eql :unary-operator)) &rest args &key value)
  (declare (optimize (debug 3)))
  (core:make-cxx-object 'chem:logical
                        :id (incf *next-chem-info-id*)
                        :op :log-identity))


(defmethod architecture.builder-protocol:make-node ((builder builder) (kind (eql :binary-operator)) &rest args &key operator bounds)
  (declare (optimize (debug 3)))
  (let ((operator-symbol (cond
                           ((eq operator :and) :log-high-precedence-and)
                           ((eq operator :or) :log-or)
                           (t (error "Illegal binary operator ~s" operator)))))
    (core:make-cxx-object 'chem:logical
                          :id (incf *next-chem-info-id*)
                          :op operator-symbol)))


(defmethod architecture.builder-protocol:relate ((builder builder) (relation (eql :operand)) (left t) (right t) &key)
  (declare (optimize (debug 3)))
  (if (chem:get-left left)
      (chem:set-right left right)
      (chem:set-left left right))
  left)


(defun read-antechamber-type-rules (stream)
  "* Arguments
- filename : A pathname
* Description
Read the contents of the filename into memory and return a buffer-stream on it."
  (let ((builder (make-instance 'builder))
        (data (make-string (file-length stream)))
        (wild-dict (core:make-cxx-object 'chem:wild-element-dict))
        wild-atom atoms
        (type-rules (core:make-cxx-object 'chem:fftypes-db))
        (*next-chem-info-id* 0))
    (read-sequence data stream)
    (architecture.builder-protocol:with-builder builder
      (with-input-from-string (sin data)
        (loop for line = (read-line sin nil :eof)
              for lineno from 1
              for line-len = (if (stringp line) (length line) 0)
              until (eq line :eof)
              do (cond
                   ((and (> line-len #.(length "WILDATOM"))
                         (string= "WILDATOM" line :start2 0 :end2 #.(length "WILDATOM")))
                    (let* ((wild-atom-entry
                             (esrap:parse 'wild-atom-line (string-trim '(#\space #\tab) line)))
                           (wild-atom (intern (car wild-atom-entry) :keyword)))
                      (chem:add-wild-name wild-dict wild-atom)
                      (dolist (e (cdr wild-atom-entry))
                        (chem:add-wild-name-map wild-dict wild-atom (intern e :keyword)))))
                   ((and (> line-len #.(length "ATD"))
                         (string= "ATD" line :start2 0 :end2 #.(length "ATD")))
                    (let* ((*wild-dict* wild-dict)
                           (root (esrap:parse 'antechamber-line (string-trim '(#\space #\tab) line)))
                           (chem-info (core:make-cxx-object 'chem:chem-info
                                                            :root root)))
                      (chem:fftypes-db-add type-rules chem-info)))
                   (t
                    #+(or)(format *debug-io* "Unhandled line: ~a~%" line))))))
    type-rules))


(defgeneric chem:compile-antechamber-type-rule (system rule-string &optional wild-element-dict))

(defparameter *antechamber-builder* (make-instance 'leap.antechamber-type-definition-parser::builder))

(defmethod chem:compile-antechamber-type-rule ((system null) rule-string &optional wild-element-dict)
  (architecture.builder-protocol:with-builder *antechamber-builder*
    (esrap:parse 'antechamber-line (string-trim '(#\space #\tab) rule-string))))


