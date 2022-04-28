(format t "Entered cando-startup.lisp~%")
(push :dont-start-cando-user *features*)
(load "sys:extensions;cando;src;lisp;start-cando.lisp")
(format t "Done load start-cando.lisp in cando-startup.lisp~%")

(push (lambda ()
        (ql:quickload :cando-user))
      core:*extension-startup-evals*)
