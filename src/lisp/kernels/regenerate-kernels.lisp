;;;; regenerate-kernels.lisp
;;;;
;;;; Lisp driver for regenerate-kernels.sh.  Not part of the cando build --
;;;; this is run by hand (via `make kernels`) when the kernel DSL definitions
;;;; in energy-kernels.lisp or the mathkernel generator itself change.
;;;;
;;;; Loading :cando-kernels loads energy-kernels.lisp, whose final top-level
;;;; form is (write-all kernels).  write-all takes its output directory from
;;;; the KERNEL_PATH environment variable, which the shell script sets.

(require :asdf)

;; mathkernel's expression-tests.lisp defines test-stretch-diff-x1 twice, which
;; SBCL reports as a full WARNING; ASDF would otherwise escalate that to an
;; error and abort the system load.  Remove this once the duplicate is deleted.
(setf asdf:*compile-file-failure-behaviour* :warn)

(let ((mathkernel (uiop:getenv "MATHKERNEL_DIR"))
      (kernels-dir (uiop:getenv "KERNELS_DIR")))
  (unless (and mathkernel (probe-file mathkernel))
    (format *error-output* "~&regenerate-kernels: MATHKERNEL_DIR ~s not found~%" mathkernel)
    (uiop:quit 2))
  (push (truename mathkernel) asdf:*central-registry*)
  (push (truename kernels-dir) asdf:*central-registry*))

(format t "~&regenerate-kernels: writing to ~a~%" (uiop:getenv "KERNEL_PATH"))

(handler-bind ((warning #'muffle-warning))
  (handler-case
      (let ((*standard-output* (make-broadcast-stream))
            (*trace-output* (make-broadcast-stream)))
        (asdf:load-system :cando-kernels :force t))
    (error (e)
      (format *error-output* "~&regenerate-kernels: FAILED: ~a~%" e)
      (uiop:quit 1))))

;; Report what was actually produced.  A kernel group that fails to build will
;; be silently absent from this list rather than aborting the run, so compare
;; this against the build-kernel-group forms in energy-kernels.lisp.
(let* ((ks (symbol-value (find-symbol "KERNELS" :mathkernel-user)))
       (name-fn (find-symbol "KERNEL-GROUP-NAME" :mathkernel)))
  (format t "~&regenerate-kernels: ~d kernel group~:p written:~%" (length ks))
  (dolist (g (sort (mapcar name-fn ks) #'string<))
    (format t "  ~a~%" g)))

(uiop:quit 0)
