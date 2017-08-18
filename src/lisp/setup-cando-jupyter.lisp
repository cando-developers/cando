(in-package :common-lisp-user)
(defun configure-jupyter ()
  "This function is run just before the jupyter notebooks kernel starts up to configure the jupyter notebook"
  (format t "Starting nglview~%")
  (asdf:load-asd "~/quicklisp/local-projects/cl-jupyter/nglview/nglview.asd")
  (ql:quickload :nglview)
  (use-package :cando :cl-jupyter-user)
  (shadowing-import '(#:atom #:quit #:desc) :leap)
  (use-package :leap :cl-jupyter-user))
