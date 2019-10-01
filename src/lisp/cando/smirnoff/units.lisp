(in-package #:smirnoff.units)

(defrule constituent (character-ranges (#\a #\z) (#\A #\Z) #\_))

;;; Given a number, return a QUANTITY with no dimension.
;;; The UNITS library does not provide a way to do this easily,
;;; so we do num*kg*kg^-1 instead. FIXME.
(defun dimensionless (num)
  (units:* (units:make-quantity num units:kilograms)
           (units:power (units:make-quantity 1d0 units:kilograms) -1)))

;;; NOTE: Throughout, we use QUANTITYs rather than UNITs, because we
;;; can multiply them and stuff.

(defun name->unit (name)
  (macrolet ((string-case (&rest clauses)
               `(cond
                  ,@(loop for (string . body) in clauses
                          collect `((string= name ,string) ,@body))
                  (t (error "Unknown unit: ~a" name)))))
    (string-case
     ("kilogram" (units:make-quantity 1d0 units:kilograms))
     ("angstrom" (units:make-quantity 1d0 units:angstroms))
     ("kilocalorie" (units:make-quantity 4184d0 units:joules))
     ("mole" (units:make-quantity 1d0 units:moles))
     ("kilocalories_per_mole"
      (units:make-quantity 1d0 units:kilocalories-per-mole))
     ("degree" (units:make-quantity 1d0 units:degrees))
     ("radian" (units:make-quantity 1d0 units:radians)))))

(defrule primitive (+ constituent)
  (:lambda (list)
    (name->unit (esrap:text list))))

(defrule power (and #\* #\*) (:constant nil))
(defrule mult #\* (:constant nil))

(defrule unit (and primitive (? whitespace+) (? (and power integer-literal/decimal)))
  (:destructure (prim ws maybe-exp)
    (declare (ignore ws))
    (let ((exp (if maybe-exp (second maybe-exp) 1)))
      (units:power prim exp))))

(defrule quantity-unit-part
    (and (? whitespace+) mult (? whitespace+) unit)
  (:function fourth))

(defrule quantity (and double-float-literal (* quantity-unit-part))
  (:destructure (num units)
    ;; NUM is a float, UNITS is a list of QUANTITYs.
    (let ((res (dimensionless num)))
      (loop for unit in units
            do (setf res (units:* res unit)))
      res)))

(defun parse-quantity (string)
  (values (parse 'quantity string)))

;;; Examples:
;;; (parse 'quantity "1.41 * angstrom") => 1.41e-10*m
;;; (parse 'quantity "640.0 * angstrom**-2 * mole**-1 * kilocalorie")
;;;   => 2.67776e+26*kg*s^-2*mole^-1
