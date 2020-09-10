(in-package :cl-user)

(declaim (optimize (safety 3)))

(load (compile-file "sys:regression-tests;framework.lisp"))

(in-package :clasp-tests)

(reset-clasp-tests)
(load-if-compiled-correctly "source-dir:extensions;cando;src;lisp;regression-tests;leap.lisp")
(load-if-compiled-correctly "source-dir:extensions;cando;src;lisp;regression-tests;spanning-tree.lisp")
