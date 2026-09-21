(load "~/quicklisp/setup.lisp")

(let ((pn (merge-pathnames "cando-kernels.asd" *load-pathname*)))
  (format t "pn = ~s~%" pn)
  (asdf:load-asd pn)
  (ql:quickload :cando-kernels)
  )
