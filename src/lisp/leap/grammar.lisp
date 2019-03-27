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
  (:constant nil))

(defrule skippable?
    (? skippable))

(defrule skippable
    (+ (or whitespace-character comment))
  (:constant nil))

(defrule horizontal-whitespace
    (+ spacoid-character)
  (:constant nil))

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
  (:constant nil))

(defrule comment
    (and shell-style-comment/trimmed
         (* (and newline (? whitespace) shell-style-comment/trimmed)))
  (:text t)
  (:lambda (content &bounds start end)
    (architecture.builder-protocol:node* (:comment :content content
                                                   :bounds  (cons start end)))))

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
    (and (or assignment function) end-of-instruction)
  (:function first))

;;; Function Call

(defrule function
    (or (and function-name/?s arguments)
        (and function-name/?s (and)))
  (:destructure (name arguments &bounds start end)
    (architecture.builder-protocol:node* (:function :name   name
                                                    :bounds (cons start end))
      (* :argument arguments))))

(defrule/s function-name
    (function-name-p function-name-string))

(defvar *function-names*
  '("add" "addH" "addIons" "addIons2" "addIonSolv" "addIonsRand" "addPath"
    "addPdbAtomMap" "addPdbResMap" "addAtomTypes" "alias" "alignAxes" "bond"
    "bondByDistance" "center" "charge" "check" "clearPdbAtomMap"
    "clearPdbResMap" "clearVariables" "combine" "copy" "createAtom"
    "createParmset" "createResidue" "createUnit" "crossLink" "debugOff"
    "debugOn" "debugStatus" "deleteBond" "deleteOffLibEntry" "deleteRestraint"
    "desc" "deSelect" "displayPdbAtomMap" "displayPdbResMap" "edit" "flip"
    "groupSelectedAtoms" "help" "impose" "list" "listOff" "loadAmberParams"
    "loadAmberPrep" "loadOff" "loadMol2" "loadMol3" "loadPdb" "loadPdbUsingSeq"
    "logFile" "matchVariables" "measureGeom" "memDebug" "quit" "relax" "remove"
    "restrainAngle" "restrainBond" "restrainTorsion" "saveAmberParm"
    "saveAmberParmNetCDF" "saveAmberParmPert" "saveAmberParmPol"
    "saveAmberParmPolPert" "saveAmberPrep" "saveMol2" "saveMol3" "saveOff"
    "saveOffParm" "savePdb" "scaleCharges" "select" "sequence" "set" "setBox" "setupDefaultPaths"
    "showdefault" "showPaths" "solvateBox" "solvateCap" "solvateDontClip" "solvateOct"
    "solvateShell" "source" "translate" "transform" "update" "verbosity"
    "zMatrix"))

(defvar *function-names/hash-table*
  (alist-hash-table (mapcar (rcurry #'cons t) *function-names*) :test #'equalp))

(defun function-name-p (thing)
  (gethash thing *function-names/hash-table*))

(defrule function-name-string
    (+ (character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9)))
  (:text t))

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

(defrule/s variable-name
    (+ (character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9) #\. #\_))
  (:lambda (chars)
    (let ((str (esrap:text chars)))
      (if (position #\. str)
          str
          (intern (esrap:text chars) :keyword)))))

;;; Expression

(defrule expression
    (or raw-expression function))

(defrule/s raw-expression
    (or list literal))

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
    (or integer-literal/decimal
        string-literal/double-quotes
        string-literal/dollar
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
