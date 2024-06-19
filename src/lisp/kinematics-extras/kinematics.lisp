
(in-package :kin)

(define-condition kin:kinematics-error (error)
  ()
  (:report (lambda (condition stream)
             (format stream "~a ~a" (class-name (class-of condition))))))

(define-condition kin:undefined-internal-coordinates (kin:kinematics-error)
  ((joint :initarg :joint :reader joint)))
