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

;;;(cl:in-package #:language.smarts.parser)
(in-package :language.smarts.parser)

(defrule smarts
    (and atom-pattern (? chain-elements))
  (:destructure (atom elements &bounds start end)
                (if elements
                    (architecture.builder-protocol:node* (:chain :bounds (cons start end))
                                                         (* :element (list* atom elements)))
                    atom)))

(defrule chain-pattern
  chain-elements
  (:lambda (elements &bounds start end)
    (architecture.builder-protocol:node* (:chain :bounds (cons start end))
                                         (* :element elements))))

(defrule chain-elements
    (+ (or bond-atom-pattern branch)))

(defrule branch
    (and #\( chain-pattern #\)) ; TODO ast node for this?
  (:function second))

(defrule bond-atom-pattern
    (and (? bond-pattern) atom-pattern)
  (:destructure (bond atom &bounds start end)
    (architecture.builder-protocol:node* (:bond :kind (or bond :single) :bounds (cons start end))
      (1 :atom atom))))

(defrule atom-pattern ; TODO duplicated form SMILES?
    (and acyclic-atom-pattern (? parser.common-rules:integer-literal/decimal))
  (:destructure (atom label &bounds start end)
    (if label
        (architecture.builder-protocol:node* (:labeled :label label :bounds (cons start end))
          (1 :atom atom))
        atom)))

(defrule acyclic-atom-pattern
    (or modified-atom-pattern language.smiles.parser:atom-symbol))

(defrule modified-atom-pattern
    (and #\[ weak-and-expression #\])
  (:function second)
  (:lambda (expression &bounds start end)
    (architecture.builder-protocol:node* (:bracketed-expression :bounds (cons start end))
      (1 :expression expression))))

;;; SMARTS 4.1 Atomic Primitives

;;; There are no modifiers: 4.3 Logical Operators
;;;
;;; All atomic expressions which are not simple primitives must be
;;; enclosed in brackets. The default operation is & (high precedence
;;; "and"), i.e., two adjacent primitives without an intervening
;;; logical operator must both be true for the expression (or
;;; subexpression) to be true.
(defrule modified-atom-pattern-body
    (or language.smiles.parser:atom-weight ; TODO this is just a number
        language.smiles.parser:atom-symbol

        ; language.smiles.parser:hydrogen-count
        language.smiles.parser:charge
        language.smiles.parser:chirality

        atom-pattern/non-literal

        ;; clrule
        
        language.smiles.parser:atom-map-class

        recursive))

#+(or)
(progn
  (defrule clrule
      (and "<" symbol-tag ">")
    (:lambda (lt sym gt)
      (architecture.builder-protocol:node* (:atom :clrule sym))))

  (parser.common-rules:defrule/s symbol-tag
      (and (esrap:character-ranges (#\a #\z) (#\A #\Z))
           (* (or (esrap:character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9)))))
    (:text t))
  )


     
(macrolet
    ((define-rules (&body clauses)
       (let ((rules '()))               ; TODO unused
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
                            (architecture.builder-protocol:node* (:atom ',kind value :bounds (cons start end))))))
                      ((nil)
                       `(defrule ,rule-name
                            ,expression
                          (:lambda (value &bounds start end)
                            (declare (ignore value))
                            (architecture.builder-protocol:node* (:atom :kind ',kind :bounds (cons start end))))))))))
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
    ;; TODO allow x?
    (smallest-ring-size      #\r :parameter &optional)
    (valence                 #\v :parameter t)
    (connectivity            #\X :parameter t)
    (atomic-number           #\# :parameter t)))

;;; SMARTS 4.2 Bonds Primitives

(defrule bond-pattern
    (or bond-pattern/non-literal language.smiles.parser:bond))

;;; Only additions to SMILES 3.2.2.
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
    (down-or-unspecified "\\?")))

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
  (1 not-expression          operator-not)
  (2 implicit-and-expression operator-implicit-and)
  modified-atom-pattern-body)

;;; SMARTS 4.4 Recursive SMARTS

(defrule recursive
    (and #\$ #\( smarts #\))
  (:function third)
  (:lambda (pattern &bounds start end)
    (architecture.builder-protocol:node* (:recursive :bounds (cons start end))
      (1 :pattern pattern))))
