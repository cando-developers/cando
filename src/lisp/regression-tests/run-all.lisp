(in-package :cl-user)

(declaim (optimize (safety 3)))

(load (compile-file "sys:src;lisp;regression-tests;framework.lisp"))

(in-package :clasp-tests)

(reset-clasp-tests)
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;geometry.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;leap.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;spanning-tree.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;energy.lisp")
(ext:quit (if (show-test-summary) 0 1))
