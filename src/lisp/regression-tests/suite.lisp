(in-package #:clasp-tests)

(export '(run-cando-regression-tests
          *cando-regression-test-files*)
        '#:clasp-tests)

(defparameter *cando-regression-test-directory*
  #.(make-pathname :name nil :type nil :version nil
                   :defaults (or *compile-file-truename* *load-truename*))
  "Directory containing the Cando regression test source files.")

(defparameter *cando-regression-test-files*
  '("cremer-pople.lisp"
    "geometry.lisp"
    "leap.lisp"
    "spanning-tree.lisp"
    "pairwise-derivatives.lisp"
    "rosetta-nonbond.lisp"
    "rosetta-elec.lisp"
    "energy.lisp"
    "copy-energy-function.lisp"
    "smirnoff-cache.lisp"
    ;; This cross-repository cache regression loads OPEN-FORCE-FIELD and SPIROS.
    "cached-types.lisp"
    ;; This test installs a minimal :ROSETTA force field; keep it last.
    "lksolvation.lisp")
  "Cando regression files in their required execution order.")

(defun run-cando-regression-tests
    (&key (directory *cando-regression-test-directory*))
  "Compile, load, and report every Cando regression test.

Return true only if all assertions pass and every source file compiles and
loads.  The :TESTS reader feature is bound while the files are compiled so
this works both from the command-line runner and from an interactive image."
  (reset-clasp-tests)
  (let ((*features* (adjoin :tests *features*)))
    (dolist (relative-path *cando-regression-test-files*)
      (load-if-compiled-correctly
       (merge-pathnames relative-path directory))))
  (let ((tests-passed (show-test-summary))
        (files-loaded (null *files-failed-to-compile*))
        (tests-ran (plusp (hash-table-count *test-marker-table*))))
    (and tests-passed files-loaded tests-ran)))
