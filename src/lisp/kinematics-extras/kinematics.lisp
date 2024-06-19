
(in-package :kin)

(define-condition kinematics-error (error)
  ()
  (:report (lambda (condition stream)
             (format stream "~a ~a" (class-name (class-of condition))))))

(define-condition undefined-internal-coordinates (kinematics-error)
  ((joint :initarg :joint :reader joint)))
