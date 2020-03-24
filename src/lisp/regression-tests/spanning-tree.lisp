
(in-package #:clasp-tests)


(let* ((mol (chem:load-chemdraw "data/sample-molecule.cdxml"))
       (res (chem:content-at mol 0))
       (atom (chem:content-at res 0))
       (span (chem:make-spanning-tree atom)))
  span)

(test spanning-tree
      (let* ((mol (chem:load-chemdraw "data/sample-molecule.cdxml"))
             (res (chem:content-at mol 0))
             (atom (chem:content-at res 0))
             (span (chem:make-spanning-tree atom)))
        (chem:content-at (chem:content-at mol 0) 0)))
