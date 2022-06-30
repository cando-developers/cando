
(in-package #:clasp-tests)

(leap:setup-default-paths)

(test-type load-frcmod
           (leap:load-amber-params "frcmod.ff99SB")
           chem:force-field)

(test-type load-frcmod-scee-scnb
           (leap:load-amber-params "frcmod.fb15")
           chem:force-field)

(test-type load-frcmod-pho
           (leap:load-amber-params (namestring (translate-logical-pathname "sys:extensions;cando;src;lisp;regression-tests;data;frcmod.pho")))
           chem:force-field)
