
(load (merge-pathnames #P"build-cando.lisp" *load-pathname*))
(ql:quickload "cl-jupyter")
(ql:quickload "nglview")
(ql:quickload "trivial-http")
(asdf:load-asd (merge-pathnames #P"build-cando-jupyter.asd" *load-pathname*))
