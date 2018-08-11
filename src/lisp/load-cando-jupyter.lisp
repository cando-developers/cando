
(load (merge-pathnames #P"load-cando.lisp" *load-pathname*))
(asdf:load-asd (merge-pathnames #P"cando-jupyter/cando-jupyter.asd" *load-pathname*))
(asdf:load-asd (merge-pathnames #P"load-cando-jupyter.asd" *load-pathname*))
(ql:quickload "load-cando-jupyter" :verbose t)
