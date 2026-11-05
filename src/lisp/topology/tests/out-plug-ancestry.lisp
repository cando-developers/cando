(defpackage #:topology-out-plug-ancestry-tests
  (:use #:cl)
  (:export #:run-tests))

(in-package #:topology-out-plug-ancestry-tests)

(defun run-tests ()
  (flet ((validate (expression)
           (topology::validate-out-plug-ancestry
            (topology::interpret 'ancestry-test expression))))
    ;; An attachment on the root or its child lacks reference atoms.
    (dolist (expression '(((C C :plugs (:+amide)) - O)
                          (O = (C C :plugs (:+ring-amide))
                           - CA.2HA - CB.2HB - (N N :plugs (:+amide)) - HN)))
      (let ((message (handler-case (progn (validate expression) nil)
                       (error (condition) (princ-to-string condition)))))
        (assert (and message (search "needs three reference atoms" message)))))
    ;; Reordered RING-BALA provides ancestry at both outgoing plugs.
    (validate '(HA1 - CA.HA2
                (- (C C :plugs (:+ring-amide)) = O)
                - CB.2HB - (N N :plugs (:+amide)) - HN))
    ;; A preceding residue supplies ancestry through the incoming plug.
    (validate '((N N :plugs (:-amide)) - (C C :plugs (:+amide)) = O))
    ;; A topology without outgoing plugs needs no such references.
    (validate '(O = C)))
  (format t "Outgoing plug ancestry tests passed.~%")
  t)
