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
(let ((rosetta-host-p (ext:logical-host-p "ROSETTA")))
  (load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;cremer-pople.lisp")
  (load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;geometry.lisp")
  (load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;leap.lisp")
  (load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;spanning-tree.lisp")
  (load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;pairwise-derivatives.lisp")
  (when rosetta-host-p
    (load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;rosetta-nonbond.lisp")
    (load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;rosetta-elec.lisp")
    (load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;energy.lisp")
    (load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;amber-excluded-atoms.lisp")
    (load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;smirnoff-cache.lisp"))
  (unless rosetta-host-p
    (message :warn
             "Skipping rosetta-nonbond.lisp, rosetta-elec.lisp, energy.lisp, amber-excluded-atoms.lisp, and smirnoff-cache.lisp because the ROSETTA logical host is not configured"))
  ;; This test installs a minimal :ROSETTA force field; keep it last.
  (load-if-compiled-correctly "sys:extensions;cando;src;lisp;regression-tests;lksolvation.lisp"))

#-swank(ext:quit (if (show-test-summary) 0 1))
#+swank(show-test-summary)
