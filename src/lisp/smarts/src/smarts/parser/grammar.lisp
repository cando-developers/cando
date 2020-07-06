;;;; grammar.lisp --- Grammar for the SMARTS language.
;;;;
;;;; Copyright (C) 2018 Jan Moringen
;;;;
;;;; Author: Jan Moringen <jmoringe@techfak.uni-bielefeld.de>

;;;; Esrap parser for the SMARTS language as an extension of the Esrap
;;;; parser for the SMILES language.
;;;;
;;;; See https://en.wikipedia.org/wiki/Smiles_arbitrary_target_specification
;;;;     https://docs.eyesopen.com/toolkits/cpp/oechemtk/SMARTS.html
;;;;     http://www.daylight.com/cheminformatics/index.html

(cl:in-package #:language.smarts.parser)

(defrule smirks (and smarts ">>" smarts)
  (:destructure (reagent dummy product &bounds start end)
                (declare (ignore dummy))
                (cons reagent product)))

(defrule smarts
    (and atom-pattern (? chain-elements))
  (:destructure (atom elements &bounds start end)
    (if elements
        (bp:node* (:chain :bounds (cons start end))
          (* :element (list* atom elements)))
        atom)))

(defrule chain-pattern
    chain-elements
  (:lambda (elements &bounds start end)
    (bp:node* (:chain :bounds (cons start end))
      (* :element elements))))

(defrule chain-elements
    (+ (or bond-atom-pattern branch)))

(defrule branch
    (and #\( chain-pattern #\)) ; TODO ast node for this?
  (:function second))

(defrule bond-atom-pattern
    (and (? bond-pattern) atom-pattern)
  (:destructure (bond atom &bounds start end)
    (bp:node* (:bond :kind (or bond :single-or-aromatic) :bounds (cons start end))
      (1 :atom atom))))

(defrule atom-pattern ; TODO duplicated form SMILES?
    (and acyclic-atom-pattern (? parser.common-rules:integer-literal/decimal))
  (:destructure (atom label &bounds start end)
    (if label
        (progn
;;          (format t "atom-pattern got label |~s| of type |~s|~%" label (class-of label))
          (bp:node* (:labeled :label label :bounds (cons start end))
                                               (1 :atom atom)))
        atom)))

(defrule acyclic-atom-pattern
    (or modified-atom-pattern wildcard-atom-pattern smiles:atom-symbol))

(defrule modified-atom-pattern
    (and #\[ modified-weak-and-expression #\])
  (:function second)
  (:lambda (expression &bounds start end)
    expression ; returning the expression eliminates useless log-identity nodes
    #+(or)(bp:node* (:bracketed-expression :bounds (cons start end))
      (1 :expression expression))))

(defrule modified-weak-and-expression
    (or weak-and-expression-with-map weak-and-expression))

(defrule weak-and-expression-with-map
    (and weak-and-expression smiles:atom-map-class)
  (:destructure (expression atom-map &bounds start end)
                (bp:node* (:binary-operator :operator :weak-and :bounds (cons start end))
                                                     (* :operand (list expression atom-map)))))



(defrule lisp-func
  (and #\< (esrap:character-ranges (#\a #\z) (#\A #\Z))
       (* (or (esrap:character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9))))
       #\>)
  (:lambda (symbol-parts &bounds start end)
    (let* ((label (string-trim "<>" (esrap:text symbol-parts)))
           (symbol (let ((*package* (find-package :keyword)))
                     (read-from-string label))))
      (bp:node* (:atom :lisp-function symbol
                                                  :bounds (cons start end))))))

;;; SMARTS 4.1 Atomic Primitives

;;; There are no modifiers: 4.3 Logical Operators
;;;
;;; All atomic expressions which are not simple primitives must be
;;; enclosed in brackets. The default operation is & (high precedence
;;; "and"), i.e., two adjacent primitives without an intervening
;;; logical operator must both be true for the expression (or
;;; subexpression) to be true.
(defrule modified-atom-pattern-body
    (or smiles:atom-weight ; TODO this is just a number, i.e. (node* ) is missing
        smiles:atom-symbol
        lisp-func
        
        ; smiles:hydrogen-count
        smiles:charge
        smiles:chirality

        atom-pattern/non-literal

;;; ------ smiles:atom-map-class can't be here it will be strong-and and it needs to be weak-and
;;;        smiles:atom-map-class

        recursive))

(defrule wildcard-atom-pattern
  #\*
  (:lambda (value &bounds start end)
    (declare (ignore value))
    (bp:node* (:atom :kind :wildcard :bounds (cons start end)))))

(macrolet
    ((define-rules (&body clauses)
       (let ((rules '()))
         (flet ((process-clause (name expression
                                 &key
                                 parameter
                                 (kind (make-keyword name)))
                  (let ((rule-name (symbolicate '#:atom-rule- name)))
                    (push rule-name rules)
                    (ecase parameter
                      ((t &optional)
                       `(defrule ,rule-name
                            (and ,expression
                                 ,(if (eq parameter '&optional)
                                      `(? parser.common-rules:integer-literal/decimal)
                                      'parser.common-rules:integer-literal/decimal))
                          (:function second)
                          (:lambda (value &bounds start end)
                            (bp:node* (:atom ',kind value :bounds (cons start end))))))
                      ((nil)
                       `(defrule ,rule-name
                            ,expression
                          (:lambda (value &bounds start end)
                            (declare (ignore value))
                            (bp:node* (:atom :kind ',kind :bounds (cons start end))))))))))
           `(progn
              ,@(map 'list (curry #'apply #'process-clause) clauses)
              (defrule atom-pattern/non-literal (or ,@(nreverse rules))))))))
  (define-rules
    (wildcard                #\*)

    (aromatic                #\a)
    (aliphatic               #\A)

    (degree                  #\D :parameter t)
    (total-hydrogen-count    #\H :parameter &optional)
    (implicit-hydrogen-count #\h :parameter &optional)
    (ring-bond-count         #\R :parameter &optional)
    (ring-connectivity       #\x :parameter &optional) ;; TODO allow x?
    (smallest-ring-size      #\r :parameter &optional)
    (valence                 #\v :parameter t)
    (connectivity            #\X :parameter t)
    (atomic-number           #\# :parameter t)))

;;; SMARTS 4.2 Bonds Primitives

;;; original code without logicals
#+(or)
(defrule bond-pattern
    (or bond-pattern/non-literal smiles:bond))

;;; New code with logicals
(defrule bond-pattern
    (or modified-bond-pattern one-bond-pattern))

(defrule modified-bond-pattern
  (and weak-and-bond-expression)
  (:lambda (expression &bounds start end)
    (bp:node* (:logical-bond-expression :bounds (cons start end))
                                         (1 :bond-expression expression))))

;;; SMARTS 4.3 Logical Operators

(macrolet ((define-operator-rule (name expression
                                  &optional (value (make-keyword name)))
             (let ((rule-name (symbolicate '#:operator- name)))
               `(defrule ,rule-name
                    ,expression
                  (:constant ,value)))))
  (define-operator-rule bond-weak-and     #\;   :bond-weak-and)
  (define-operator-rule bond-or           #\,)
  (define-operator-rule bond-strong-and   #\&   :bond-strong-and)
  (define-operator-rule bond-not          #\!)
  (define-operator-rule bond-implicit-and (and) :bond-implicit-and)) ; same as strong-and

(parser.common-rules.operators:define-operator-rules
    (:skippable?-expression nil)
  (2 weak-and-bond-expression     operator-bond-weak-and)
  (2 or-bond-expression           operator-bond-or)
  (2 strong-and-bond-expression   operator-bond-strong-and)
  (2 implicit-and-bond-expression operator-bond-implicit-and)
  (1 not-bond-expression          operator-bond-not)
  one-bond-pattern)


(defrule one-bond-pattern
    (or bond-pattern/non-literal smiles:bond))

;;; Only additions to SMILES 3.2.2 and higher.
(macrolet
    ((define-rules (&body clauses)
       (let ((rules '()))
         (flet ((process-clause (name expression
                                 &key
                                 (kind (make-keyword name)))
                  (let ((rule-name (symbolicate '#:bond-pattern- name)))
                    (push rule-name rules)
                    `(defrule ,rule-name
                         ,expression
                       (:constant ,kind)))))
           `(progn
              ,@(map 'list (curry #'apply #'process-clause) clauses)
              (defrule bond-pattern/non-literal (or ,@(nreverse rules))))))))
  (define-rules
    (wildcard            #\~)
    (up-or-unspecified   "/?")
    (down-or-unspecified "\\?")
    (same-ring           #\@) ; added in smarts 4.6
    ))

;;; SMARTS 4.3 Logical Operators

(macrolet ((define-operator-rule (name expression
                                  &optional (value (make-keyword name)))
             (let ((rule-name (symbolicate '#:operator- name)))
               `(defrule ,rule-name
                    ,expression
                  (:constant ,value)))))
  (define-operator-rule weak-and     #\;   :weak-and)
  (define-operator-rule or           #\,)
  (define-operator-rule strong-and   #\&   :strong-and)
  (define-operator-rule not          #\!)
  (define-operator-rule implicit-and (and) :implicit-and))

(parser.common-rules.operators:define-operator-rules
    (:skippable?-expression nil)
  (2 weak-and-expression     operator-weak-and)
  (2 or-expression           operator-or)
  (2 strong-and-expression   operator-strong-and)
  (2 implicit-and-expression operator-implicit-and)
  (1 not-expression          operator-not)
  modified-atom-pattern-body)

;;; SMARTS 4.4 Recursive SMARTS

(defrule recursive
    (and #\$ #\( smarts #\))
  (:function third)
  (:lambda (pattern &bounds start end)
    (bp:node* (:recursive :bounds (cons start end))
      (1 :pattern pattern))))
