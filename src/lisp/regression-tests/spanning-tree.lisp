
(in-package #:clasp-tests)


(test spanning-tree
      (let* ((filename (translate-logical-pathname "source-dir:extensions;cando;src;lisp;regression-tests;data;A1.mol2"))
             (agg (chem:load-mol2 filename))
             (mol (chem:content-at agg 0))
             (res (chem:content-at mol 0))
             (atm (chem:content-at res 0))
             (span (chem:spanning-loop/make atm))
             (atoms (chem:all-atoms span)))
        (= (length atoms) 95)))
