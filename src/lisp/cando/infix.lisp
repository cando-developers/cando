(in-package :infix)
;; infix.lisp
;; by AmkG

;; Title	:	Infix notation macro
;; Filename	:	infix.lisp
;; Written by	:	Alan Manuel K. Gloria <almkglor@gmail.com>
;; Corrected by :       Bruno Daniel <bruno.daniel@gmx.net>

;; Copyright status:
;; Err.  I don't care?  Okay, I do care.  To be honest, I'll just
;; be plain happy if somebody ELSE uses this.  Drop me a line,
;; or something.  Or don't.  I built this for my own amusement,
;; so if anyone else is amused too, well, drop me a line!

;; tested on:
;;  clisp
;;  gcl
;;  sbcl
;;  ecl

;; infix notation:
;; Lisp has traditionally used prefix notation for all formulae: in this
;; notation, the operation to be performed is specified before the
;; operands it is to be performed on.  This simple notation gives a
;; consistent, regular syntax which greatly facilitates
;; metaprogramming.

;; However, one weakness of prefix notation is that some simple
;; algebraic/arithmetic operations do not and cannot follow the
;; traditional syntax taught to us during elementary - that is, the
;; notation we are used to in arithmetic is not prefix, but infix
;; notation:
;;    1 + 2 - 4 * 5
;; In Lisp this would be:
;;   (- (+ 1 2) (* 4 5)))
;; Despite the inherent advantages of prefix notation, sometimes
;; it's just much easier to write mathematical formulae a little more
;; like the way we actually write them on paper.

;; This file contains a macro, 'nfx, which allows for such an infix
;; notation.  There are limitations: you need to put a LOT of spaces
;; in the formula:
;; 	(nfx 1+2-3) ;would confuse the Lisp parser
;; 	(nfx 1 + 2 - 3) ;would be understood

;; You can't use certain symbols as variable names:
;; 	(setf + 42)
;; 	(nfx (max + 3)) ; would look awfully like +(max,3), not max(+,3)

;; ...but you can embed prefix notation in it!
;; 	(nfx 1 + (- x 100)) ;it's valid!
;; 	(nfx 1 + (- x (3 * 3))) ;it's ALSO valid!
;; 	(nfx 1 + (- x 3 * 3)) ;err... this can give you unexpected behavior

;; ...also, you can define your own infix symbols using definfix:
;; 	(definfix my-infix-symbol
;; 		:precedence 40
;; 		:function-name +)
;; 	(nfx 1 my-infix-symbol 2)
;; 		=> 3

;; send any feedback, bugreports, bugfixes, and cute girls to:
;; 	almkglor@gmail.com

;; This file defines a macro of the form (nfx ...) whose
;; parameters are a stream of LISP "tokens" (either symbols,
;; constants, or lists).  The stream of tokens is interpreted
;; in infix notation.  The macro then expands to the prefix
;; form equivalent to the infix form
;; Ex.
;; (macroexpand-1 '(nfx foo = 32))
;; 	=> (setf foo 32)
;; (macroexpand-1 '(nfx bar == (3 + foo) * quux ))
;; 	=> (equal bar (* (+ 3 foo) quux))

;; NOTE: detection of function calls
;; function call forms should be supported:
;; (macroexpand-1 '(nfx foo = (max (bar + 32) quux niaw) ))
;; 	=> (setf foo (max (+ bar 32) quux niaw))
;; (macroexpand-1 '(nfx (cdr foo) = (cons (qux + 1) nil)))
;; 	=> (setf (cdr foo) (cons (+ qux 1) nil))

;; If any infix notation is in a function call within an
;; nfx form, it should be within a parentheses:
;; (macroexpand-1 '(nfx (max (bar + foo) (quux + quuux)) ))
;; 	=> (max (+ bar foo) (+ quux quuux))
;; (macroexpand-1 '(nfx (max bar + foo quux + quuux) ))
;; 	=> (max bar + foo quux + quuux)

;; function calls are detected in the following manner:
;; if nfx detects a list in the input stream,
;; 	if the second item is not a registered infix,
;; 		function call, for each element recurse and replace the element
;; 		not a function call, recurse on the list as a new stream

;; this allows us to use prefix operators (such as -) as-is:
;; (macroexpand-1 '(nfx (- bar) == (/ (foo + 1)) ))
;; 	=> (equal (- bar) (/ (+ foo 1)))


;; (nfx-operator-base <action> ...)
;; Handles the operator database

(let ((nfx-data (make-hash-table :test 'eq)))
  (defun nfx-operator-base (action &rest params)
    (labels
        ((getval (sym)
           (gethash sym nfx-data) )
         (add (params)
           (setf (gethash (car params) nfx-data)
                 (apply #'vector params))))
      (case action
        (:get (getval (car params)))
        (:add (add params))))))


;;; the nfx macro

(defun nfx-impl (s)
  (let (opstack
        curop tmp
        (bldg (cons nil nil)))
    (labels
        ( ;; tconc = a cons with pointers to the head and the last element of
         ;; a list. This function appends an element to such a tconc.
         (tconc-append (l v)
           (if (car l)
               (setf (cdr l) (setf (cddr l) (cons v nil)))
               (setf (car l) (setf (cdr l) (cons v nil))) ))
         ;; top-of-stack query
         (top-opstack ()
           (caaar opstack))
         ;; handles any sub-lists
         (enlist (l)
           (when l
             (cond ((getop (cadr l))
                    (cons 'nfx l))
                   (t (mapcar #'(lambda (o) (if (listp o) (enlist o) o))
                              l)))))
         ;; pushes an operation onto the stack - used when the currently
         ;; being built operation is of lower precedence than the
         ;; operation being considered
         (push-oper ()
           (push bldg opstack)
           (setf bldg (cons nil nil)))
         ;; pops off an operation from the stack,
         ;; - used before collapsing the call
         (pop-oper ()
           (let ((top (pop opstack)) )
             (tconc-append top (car bldg))
             (setf bldg top)))
         ;; sub-expression handling code
         (expr (o)
           (if (listp o)
               (enlist o)
               o))
         ;; determine if the specified operation has precedence over the
         ;; operation currently being built
         (precedes? (op1 op2 associativity)
           (or (< (precedence op1) (precedence op2))
               (and (= (precedence op1) (precedence op2))
                    (eq (associativity op1) associativity))))
         ;; funges bldg: collapses fungible operations into one form
         (fungebldg (bldg)
           (mapcan
            #'(lambda (o)
                (cond ((listp o)
                       (setf o (fungebldg o))
                       (cond ((and (fungible (car bldg))
                                   (eq (car o) (car bldg)))
                              (cdr o))
                             (t (cons o nil))))
                      (t (cons o nil))))
            bldg))
         ;; fixes bldg: changes operation objects to their functions
         (fixbldg (bldg)
           (when (vectorp (car bldg))
             (setf (car bldg) (function-name (car bldg))) )
           (mapc
            #'(lambda (o)
                (when (listp o)
                  (fixbldg o) ))
            bldg ))
         ;; error-handling function
         (err (l)
           (error "improper nfx expression:~%~s" l))
         (print-expr-list (l)
           (format t "(")
           (loop :for x in l
              :for i from 0 :do
              (when (plusp i)
                (format t " "))
              (cond ((vectorp x)
                     (format t "~s" (aref x 0)))
                    ((listp x)
                     (print-expr-list x))
                    (t
                     (format t "~s" x))))
           (format t ")"))
         ;;          (print-tconc (tc)
         ;;            (print-expr-list (car tc))
         ;;            (format t "~%"))
         ;; accessor functions
         (precedence (op) (if op (aref op 1) 99999))
         (associativity (op) (aref op 2))
         (function-name (op) (aref op 3))
         (fungible (op) (when (vectorp op) (aref op 4)))
         ;; get the data for a (presumed) infix operator
         (getop (op)
           (nfx-operator-base :get op)))
      (case (length s)
        (1 (expr (car s)))
        (2 (err s))
        (t
         (do ( ;; variable and step list
              (oL s (cddr oL)) )
             ( ;; termination condition
              (null (cddr oL))
              (tconc-append bldg (expr (car oL)))
              (when (cdr oL)
                (err oL) )
              (loop :while opstack :do
                 (pop-oper))
              (fixbldg (fungebldg (caar bldg))))
           (setf curop (getop (cadr oL)))
           (cond (curop
                  (cond ((precedes? curop (caar bldg) :right)
                         (push-oper)
                         (tconc-append bldg curop)
                         (tconc-append bldg (expr (car oL))) )
                        (t
                         (tconc-append bldg (expr (car oL)))
                         ;; collapse while the stacktop is precedent over the
                         ;; current op
                         (loop :while (precedes? (top-opstack) curop :left) :do
                            (pop-oper))
                         (setf tmp (car bldg))
                         (setf bldg (cons nil nil))
                         (tconc-append bldg curop)
                         (tconc-append bldg tmp) )))
                 (t
                  (err oL)))))))))


(defmacro nfx (&rest s)
  (nfx-impl s))

(defmacro infix (&rest s)
  (nfx-impl s))




;; (definfix <operator>
;; 	:precedence <integer>
;; 	[:associativity <:left | :right>]
;; 	[:function-name <actual function>]
;; 	[:fungible <t | nil>])
;; - defines an infix operator with the symbol <operator>,
;;   with a precedence of <integer>.  The smaller the precedence
;;   number, the more precedence it has: * and / have smaller
;;   precedence number than + and -.
;; - associativity defaults to :left, which means that if two
;;   operators of the same precedence are encountered, the first
;;   one resolves first: 1 x 2 x 3 becomes ((1 x 2) x 3).  This
;;   is appropriate for most maths.  :right associativity means
;;   that 1 x 2 x 3 becomes (1 x (2 x 3)).  This is appropriate
;;   for assignment.  IMPORTANT: operators with the same
;;   precedence must have the same associativities!
;; - function-name defaults to the same symbol as the operator.
;;   For example, the function-name of = is setf.
;; - fungible means that if the same operator is encountered
;;   several times, then all inputs are funged into one function
;;   call.  For example, 1 + 2 + 3 becomes (+ 1 2 3), 1 < 2 < 3
;;   becomes (< 1 2 3).  This is not true for assignment:
;;   x = y = z should become (setf x (setf y z)), not (setf x y z)
;; - fungible defaults to t because I noticed that nearly every
;;   single darned operator was fungible.  Except assignment.
;;   Saved some dozen lines of code too.

(defmacro definfix (operator
                    &key precedence (associativity :left)
                    function-name (fungible t))
  `(nfx-operator-base :add ',operator
    ,precedence ,associativity
    ',(if function-name function-name operator )
    ,fungible))

;; Predefined operators
;; Note: precedences are divisibles of 10, in case you
;; want to insert precedence levels between levels.

;; since we expect infix notation only (not infix-postfix),
;; our array accessor is a single infix @

(definfix @ :precedence 10 :function-name aref)
(definfix ** :precedence 20 :function-name expt)
(definfix * :precedence 30)
(definfix / :precedence 29)
(definfix % :precedence 28 :function-name mod)
(definfix + :precedence 40)
(definfix - :precedence 40)
(definfix <= :precedence 60)
(definfix < :precedence 60)
(definfix >= :precedence 60)
(definfix > :precedence 60)

(definfix = :precedence 70 :fungible nil)
(definfix /= :precedence 70)

;; 80-100 should be for the bitwise operators (&, |, ^),
;; once I figure out how CL handles bitops.

(definfix && :precedence 110 :function-name and)
(definfix || :precedence 120 :function-name or)


(set-dispatch-macro-character
 #\# #\n
 #'(lambda (stream c1 c2)
     (declare (ignorable c1 c2))
     (let ((rd (read stream t nil t)))
       (if (listp rd)
           `(nfx ,@rd)
           (progn
             (print 'invalid-#n-usage)
             (print rd)
             (error 'error) )))))


;; (nfx 55 - 3 * 4 + 2000000) used to give -1999957
;; instead of the correct result 2000043. This is now corrected. -- Bruno Daniel
