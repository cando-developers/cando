(in-package :cl-user)

(declaim (optimize (safety 3)))

(push :tests *features*)

#+swank
(defmethod gray:stream-interactive-p ((stream swank/gray::slime-output-stream)) t)

(load (compile-file "sys:src;lisp;regression-tests;framework.lisp"))

(in-package :clasp-tests)

(reset-clasp-tests)
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;cremer-pople.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;geometry.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;leap.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;spanning-tree.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;rosetta-nonbond.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;rosetta-elec.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;energy.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;smirnoff-cache.lisp")

#-swank(ext:quit (if (show-test-summary) 0 1))
#+swank(show-test-summary)
