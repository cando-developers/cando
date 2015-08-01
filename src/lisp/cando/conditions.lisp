(in-package :chem)


(define-condition minimizer-error (error)
  ((message :initarg :message
            :reader minimizer-error-message)))
(define-condition minimizer-exceeded-max-steps-error (minimizer-error) ())
(define-condition minimizer-stuck-error (minimizer-error) ())

