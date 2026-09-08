(in-package :cl-user)

(declaim (optimize (safety 3)))

(push :tests *features*)

#+swank
(defmethod gray:stream-interactive-p ((stream swank/gray::slime-output-stream)) t)

(load (compile-file "sys:src;lisp;regression-tests;framework.lisp"))
(let ((fasl (compile-file
             "sys:extensions;cando;src;lisp;regression-tests;suite.lisp")))
  (unless fasl
    (error "Could not compile the Cando regression-suite runner"))
  (load fasl))

(in-package :clasp-tests)

(reset-clasp-tests)
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;cremer-pople.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;geometry.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;leap.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;spanning-tree.lisp")
#+(or)(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;pairwise-derivatives.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;rosetta-nonbond.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;rosetta-elec.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;energy.lisp")
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;smirnoff-cache.lisp")
;; This test installs a minimal :ROSETTA force field; keep it last.
(load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;lksolvation.lisp")

#-swank(ext:quit (if (show-test-summary) 0 1))
#+swank(show-test-summary)
(let ((success (run-cando-regression-tests)))
  #-swank(ext:quit (if success 0 1))
  #+swank success)
