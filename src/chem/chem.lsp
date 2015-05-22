
(define-condition minimizer-error (error) (message))

(define-condition minimizer-exceeded-max-steps (minimizer-error))
(define-condition minimizer-stuck (minimizer-error))
