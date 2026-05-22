(in-package :xtal.symop-parser)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Top-level rule: a symmetry operator is three components separated by commas.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defrule symop
    (and component "," component "," component)
  (:destructure (comp1 comma1 comp2 comma2 comp3)
    ;; Return a list of three coordinate expressions.
    (list comp1 comp2 comp3)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; A component is an addition expression of one or more terms.
;;
;; For example, "1/2-x" parses as two terms: 1/2 and -x,
;; and we return the symbolic expression: (+ 1/2 (- x))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defrule component
    (and term (* operator-term))
  (:destructure (first-term opterms)
    (if (null opterms)
        first-term
        ;; Build an expression of the form:
        ;; (+ first-term term2 term3 …)
        (cons '+ (cons first-term (mapcar (lambda (pair)
                                            ;; pair is a list (op term); a "-" operator wraps the term in a unary minus.
                                            (let ((op (first pair))
                                                  (term (second pair)))
                                              (if (string= op "-")
                                                  (list '- term)
                                                  term)))
                                          opterms))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; An operator-term is an operator followed by a term.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defrule operator-term
    (and op term)
  (:destructure (operator term) (list operator term)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; A term consists of an optional unary sign and a factor.
;;
;; For instance, "-x" is parsed as a sign "-" and a factor "x".
;; If a sign is present, we immediately produce ( - factor ).
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defrule term
    (and (? sign) factor)
  (:lambda (parts)
    (destructuring-bind (maybe-sign fac) parts
      (if maybe-sign
          (if (string= maybe-sign "-")
              (list '- fac)
              fac)  ;; a leading "+" can be dropped
          fac))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; A factor is either a fraction (like "1/2") or an axis symbol ("x","y","z")
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defrule factor
    (or fraction axis))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; An axis is one of "x", "y", or "z". (Case‐sensitive; you can tweak this.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defrule axis (or "x" "y" "z")
  (:lambda (s)
    (intern (string-upcase s) :keyword)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; A fraction is something like "1/2", "1/4", or "3/4".
;; We use a regular expression to match one or more digits, a "/" then one or more digits,
;; and we convert the string to a rational number.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defrule fraction
    (and
     (+ (or "0" "1" "2" "3" "4" "5" "6" "7" "8" "9"))
     "/"
     (+ (or "0" "1" "2" "3" "4" "5" "6" "7" "8" "9"))
     )
  (:destructure (numer slash denom)
    (/ (parse-integer (text numer)) (parse-integer (text denom)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; An operator is either "+" or "-".
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defrule op (or "+" "-"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; A sign is also either "+" or "-".  (Used for unary signs before a factor.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defrule sign (or "+" "-"))


(defun symbolic-to-row-major-matrix (symop)
  "Convert a parsed symmetry operator (3 symbolic expressions) into a 16-element row-major 4x4 matrix."
  (labels ((coeffs (expr)
             ;; Return (x y z offset) from expression
             (let ((x 0) (y 0) (z 0) (offset 0))
               (labels ((walk (e sign)
                          (cond
                            ((numberp e) (incf offset (* sign e)))
                            ((symbolp e)
                             (ecase e
                               (:X (incf x sign))
                               (:Y (incf y sign))
                               (:Z (incf z sign))))
                            ((and (consp e) (eq (first e) '+))
                             (dolist (term (rest e)) (walk term +1)))
                            ((and (consp e) (eq (first e) '-))
                             (walk (second e) (- sign)))
                            (t (error "Unknown term in expr: ~a" e)))))
                 (walk expr +1)
                 (list x y z offset)))))
    (let* ((rows (mapcar #'coeffs symop)) ; 3 rows, each (x y z offset)
           (matrix
             ;; build rows 0–2 and final row
             (append
              (first rows)
              (second rows)
              (third rows)
              '(0 0 0 1))))
      (assert (= (length matrix) 16))
      matrix)))


