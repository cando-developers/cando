;;;; Copyright (c) 2016 Jan Moringen <jmoringe@techfak.uni-bielefeld.de>
;;;;
;;;; Permission is hereby granted, free of charge, to any person
;;;; obtaining a copy of this software and associated documentation files
;;;; (the "Software"), to deal in the Software without restriction,
;;;; including without limitation the rights to use, copy, modify, merge,
;;;; publish, distribute, sublicense, and/or sell copies of the Software,
;;;; and to permit persons to whom the Software is furnished to do so,
;;;; subject to the following conditions:
;;;;
;;;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
;;;; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
;;;; MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
;;;; IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
;;;; CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
;;;; TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
;;;; SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

(cl:in-package #:leap.parser)

;;; Lexical stuff

(defvar *list-context* nil)

(defrule end-of-instruction
    (and (? horizontal-whitespace)
         (or token-\; (& comment) newline <end-of-input>))
  (:constant nil)
  (:error-report nil))

(defrule skippable?
    (? skippable)
  (:error-report nil))

(defrule skippable
    (+ (or whitespace-character comment))
  (:constant nil))

(defrule horizontal-whitespace
    (+ spacoid-character)
  (:constant nil)
  (:error-report nil))

(defrule spacoid-character
    (or #\Space #\Tab)
  (:constant nil))

(defrule whitespace
    (+ whitespace-character)
  (:constant nil))

(defrule whitespace-character
    (or spacoid-character maybe-newline)
  (:constant nil))

(defrule maybe-newline
    newline
  (:when *list-context*))

(defrule newline
    #\Newline
  (:constant nil)
  (:error-report nil))

(defrule comment
    (and shell-style-comment/trimmed
         (* (and newline (? whitespace) shell-style-comment/trimmed)))
  (:text t)
  (:lambda (content &bounds start end)
    (architecture.builder-protocol:node* (:comment :content content
                                                   :bounds  (cons start end))))
  (:error-report nil))

(macrolet ((define-tokens (&rest characters)
             (flet ((make-rule (character)
                      (let ((rule-name (symbolicate '#:token- (string character))))
                        `(defrule/s ,rule-name
                             ,character
                           (:constant nil)))))
               `(progn ,@(mapcar #'make-rule characters)))))
  (define-tokens #\; #\= #\* #\{ #\}))

;;; Entry point

(defrule leap
    (* (or instruction
           comment
           horizontal-whitespace newline token-\;)) ; TODO hacky
  (:lambda (things &bounds start end)
    (let ((instructions (remove nil things)))
      (architecture.builder-protocol:node* (:leap :bounds (cons start end))
        (* :instruction instructions)))))

(defrule instruction
    (and (or assignment function s-expr variable-s-expr) end-of-instruction)
  (:function first))

;;; Function Call

(defrule function
    (or (and function-name/?s arguments)
        (and function-name/?s (and)))
  (:destructure (name arguments &bounds start end)
    (architecture.builder-protocol:node* (:function :name   name
                                                    :bounds (cons start end))
      (* :argument arguments))))

;;; leap.parser::*function-names*

(defvar *function-names/alist* nil
  "This is defined in commands.lisp")

(defun function-name-p (thing)
  (assoc thing *function-names/alist* :test #'string-equal))

(defrule function-name-string
    (+ (character-ranges (#\a #\z) #\- (#\A #\Z) (#\0 #\9)))
  (:text t))

(defrule matter-name-string
    (+ (character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9)))
  (:text t))

(defrule/s function-name
    (function-name-p function-name-string))


(defun file-path-unquoted-p (thing)
  (and (loop for chr across thing
             when (alpha-char-p chr)
               return t)
       (or (position #\/ thing)
           (position #\. thing))))

(defrule file-path-string
    (+ (character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9) #\~ #\- #\_ #\/ #\. ))
  (:text t))

(defrule/s file-path-unquoted
    (file-path-unquoted-p file-path-string))



(defrule arguments
    (and raw-expression (* (and skippable? raw-expression)))
  (:destructure (first rest)
    (when first (list* first (mapcar #'second rest)))))

;;; Assignment

(defrule assignment
    (or (and variable-name/?s token-=/?s (or function raw-expression))
        (and variable-name/?s token-=    (and)))
  (:destructure (name operator value &bounds start end)
    (declare (ignore operator))
    (architecture.builder-protocol:node* (:assignment :name   name
                                                      :bounds (cons start end))
      (architecture.builder-protocol:? :value value))))

(defrule package-name
    (+ (character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9) #\_ #\+ #\-))
  (:lambda (chars)
    (let ((str (string-upcase (esrap:text chars))))
      (let ((pkg (find-package str)))
        (if pkg
            pkg
            (error "Could not find package ~a" str))))))

(defrule variable-name-string
    (+ (character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9) #\_  #\+ #\- #\*))
  (:lambda (chars)
    (let ((str (esrap:text chars)))
      (leap.core:convert-leap-case str))))

(defrule package-colon-variable-name
    (and package-name #\: variable-name-string)
  (:destructure (package colon variable-name-string)
                (declare (ignore colon))
                (let ((sym (find-symbol variable-name-string package)))
                  (if sym
                      sym
                      (error "Could not find symbol ~a in package ~a" variable-name-string package)))))

(defrule current-package-variable-name
    (and variable-name-string)
  (:destructure (variable-name-string)
                (intern variable-name-string leap.core:*variable-package*)))


(defrule/s variable-name
    (or package-colon-variable-name
        current-package-variable-name
        ))

#|
    (or package-colon-variable-name
        current-package-variable-name))
|#
  
(defrule/s keyword
    (and #\: (+ (character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9) #\. #\_ #\+ #\-)))
  (:lambda (chars)
    (let ((str (string-upcase (esrap:text (cdr chars)))))
      (intern str :keyword))))


;;; Expression

(defrule expression
    (or raw-expression function))

(defrule/s raw-expression
    (or s-expr list literal s-expr))

(defun parse-list-elements (text position end)
  (let ((*list-context* t))
    (parse '(and token-{
                 (* (or raw-expression
                        comment
                        horizontal-whitespace newline))
                 token-})
           text :start position :end end :raw t)))

(defrule list
    #'parse-list-elements
  (:function second)
  (:lambda (elements &bounds start end)
    (architecture.builder-protocol:node* (:list :bounds (cons start end))
      (* :element (remove nil elements)))))

(defrule literal
    (or
;;     matter
     file-path-unquoted
     float-literal
     integer-literal/decimal
     string-literal/double-quotes
     string-literal/dollar
     keyword
     variable-name
     dummy)
  (:lambda (value &bounds start end)
    (architecture.builder-protocol:node* (:literal :value  value
                                                   :bounds (cons start end)))))

(defrule string-literal/dollar
    (and #\$ (* (not (or whitespace-character #\, #\;))) (? (or #\, #\;)))
  (:function second)
  (:text t))

(defrule dummy
    token-*)

(defun parse-with-read (text position end)
  (handler-case
      ;; When successful, READ-FROM-STRING returns the read object and
      ;; the position up to which TEXT has been consumed.
      (let ((*package* (find-package :cando-user)))
        (read-from-string text t nil :start position :end end))
    ;; When READ-FROM-STRING fails, indicate the parse failure,
    ;; including CONDITION as explanation.
    #+(or)(stream-error (condition)
            ;; For STREAM-ERRORs, we can try to determine and return the
            ;; exact position of the failure.
            (let ((position (ignore-errors
                             (file-position (stream-error-stream condition)))))
              (values nil position condition)))
    #+(or)(error (condition)
      ;; For general ERRORs, we cannot determine the exact position of
      ;; the failure.
      (values nil nil condition))))

(defrule s-expr
    (and (& #\() #'parse-with-read)
  (:function second)
  (:lambda (expression &bounds start end)
    (architecture.builder-protocol:node*
        (:s-expr :value  expression
                 :bounds (cons start end)))))


(defrule variable-s-expr
    variable-name/?s
  (:lambda (expression &bounds start end)
    (architecture.builder-protocol:node*
        (:s-expr :value  expression
                 :bounds (cons start end)))))



;;; matter expressions

(defun aggregate-name-p (thing)
  (let ((symbol (find-symbol thing :cando-user)))
    (when symbol
      (let ((aggp (typep (symbol-value symbol) 'chem:aggregate)))
        aggp))))

(defrule aggregate-name-string
    (and variable-name-string)
  (:destructure (arg)
                arg))

(defrule aggregate
    (aggregate-name-p aggregate-name-string)
  (:lambda (name)
    (leap.core:lookup-variable (intern name leap.core:*variable-package*))))

(defrule/s aggregate.pdb-sequence-number
    (and aggregate #\. integer-literal/decimal)
  (:destructure (aggregate dot pdb-sequence-number)
                (declare (ignore dot))
                (block residue
                  (chem:do-molecules (mol aggregate)
                    (chem:do-residues (res mol)
                      (when (= (chem:get-id res) pdb-sequence-number)
                        (return-from residue res)))))))

(defrule/s aggregate.residue-name
    (and aggregate #\. matter-name-string)
  (:destructure (aggregate dot residue-name)
                (declare (ignore dot))
                (let ((residue-name (intern residue-name :keyword)))
                  (block residue
                    (chem:do-molecules (mol aggregate)
                      (chem:do-residues (res mol)
                        (when (eq (chem:get-name res) residue-name)
                          (return-from residue res))))))))

(defrule/s aggregate/molecule-number
    (and aggregate #\/ integer-literal/decimal)
  (:destructure (aggregate dot molecule-number)
                (declare (ignore dot))
                (chem:content-with-id aggregate molecule-number)
                #+(or)(block molecule
                  (chem:do-molecules (mol aggregate)
                    (when (= (chem:get-id mol) molecule-number)
                        (return-from molecule mol))))))


(defrule/s aggregate.residue-number.atom-number
    (and aggregate #\. integer-literal/decimal #\. integer-literal/decimal)
  (:destructure (aggregate dot1 residue-number dot2 atom-number)
                (declare (ignore dot1 dot2))
                (block atom
                  (chem:do-molecules (mol aggregate)
                    (chem:do-residues (res mol)
                      (when (= (chem:get-id res) residue-number)
                        (chem:do-atoms (atm res)
                          (when (eq (chem:get-id atm) atom-number)
                            (return-from atom atm)))))))))

(defrule/s aggregate.residue-number.atom-name
    (and aggregate #\. integer-literal/decimal #\. matter-name-string)
  (:destructure (aggregate dot1 residue-number dot2 atom-name)
                (declare (ignore dot1 dot2))
                (let ((atom-name (intern atom-name :keyword)))
                  (block atom
                    (chem:do-molecules (mol aggregate)
                      (chem:do-residues (res mol)
                        (when (= (chem:get-id res) residue-number)
                          (chem:do-atoms (atm res)
                            (when (eq (chem:get-name atm) atom-name)
                              (return-from atom atm))))))))))

(defrule/s aggregate/molecule-number.residue-number.atom-name
    (and aggregate #\/ integer-literal/decimal #\. integer-literal/decimal #\. matter-name-string)
  (:destructure (aggregate at molecule-number dot1 residue-number dot2 atom-name)
                (declare (ignore at dot1 dot2))
                (let ((atom-name (intern atom-name :keyword)))
                  (block atom
                    (chem:do-molecules (mol aggregate)
                      (when (= (chem:get-id mol) molecule-number)
                        (chem:do-residues (res mol)
                          (when (= (chem:get-id res) residue-number)
                            (chem:do-atoms (atm res)
                              (when (eq (chem:get-name atm) atom-name)
                                (return-from atom atm)))))))))))

(defrule/s aggregate/molecule-number.residue-number.atom-number
    (and aggregate #\/ integer-literal/decimal #\. integer-literal/decimal #\. integer-literal/decimal)
  (:destructure (aggregate at molecule-number dot1 residue-number dot2 atom-number)
                (declare (ignore at dot1 dot2))
                (block atom
                  (chem:do-molecules (mol aggregate)
                    (when (= (chem:get-id mol) molecule-number)
                        (chem:do-residues (res mol)
                          (when (= (chem:get-id res) residue-number)
                            (chem:do-atoms (atm res)
                              (when (eq (chem:get-id atm) atom-number)
                                (return-from atom atm))))))))))

(defrule aggregate/molecule-number.residue-number
    (and aggregate #\/ integer-literal/decimal #\. integer-literal/decimal)
  (:destructure (aggregate at molecule-number dot1 residue-number)
                (declare (ignore at dot1 dot2))
                (block residue
                  (chem:do-molecules (mol aggregate)
                    (when (= (chem:get-id mol) molecule-number)
                      (chem:do-residues (res mol)
                        (when (= (chem:get-id res) residue-number)
                          (return-from residue res))))))))

(defrule matter
    (or aggregate/molecule-number.residue-number.atom-number
        aggregate/molecule-number.residue-number.atom-name
        aggregate/molecule-number.residue-number
        aggregate/molecule-number
        aggregate.residue-number.atom-number
        aggregate.residue-number.atom-name
        aggregate.pdb-sequence-number
        aggregate.residue-name
        aggregate
;;;        molecule
;;;        residue
;;;        atom
        ))

                


(defun parse-sub-matter (matter-string)
  (let ((result (esrap:parse 'matter matter-string :junk-allowed t)))
    result))
