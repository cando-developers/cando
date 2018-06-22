(in-package :common-lisp-user)
(defun configure-jupyter ()
  "This function is run just before the jupyter notebooks kernel starts up to configure the jupyter notebook"
  (format t "Starting cl-nglview~%")
  (asdf:load-asd "~/quicklisp/local-projects/cl-nglview/cl-nglview/cl-nglview.asd")
  (ql:quickload :cl-nglview)
  (use-package :cando :cl-jupyter-user)
  (shadowing-import '(#:atom #:quit #:desc) :leap)
  (use-package :leap :cl-jupyter-user))
