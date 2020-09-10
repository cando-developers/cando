
(in-package #:clasp-tests)

(leap:setup-default-paths)

(test load-frcmod
      (leap:load-amber-params "frcmod.ff99SB"))

(test load-frcmod-scee-scnb
      (leap:load-amber-params "frcmod.fb15"))

(test load-frcmod-pho
      (leap:load-amber-params "./data/frcmod.pho"))
