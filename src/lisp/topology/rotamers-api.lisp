(in-package :topology)

(defun write-backbone-rotamers (oligomer-shape rotamers)
  "Mid-level API: Write backbone rotamers to OLIGOMER-SHAPE.
Handles creation of permissible-backbone-rotamers internally.
For lower-level control, use make-permissible-backbone-rotamers and write-rotamers directly.
ROTAMERS can be :random, nil (first allowed), a rotamer-indexes object,
or a function that receives permissible-backbone-rotamers and returns
a value suitable for write-rotamers."
  (let ((permissible (make-permissible-backbone-rotamers oligomer-shape)))
    (write-rotamers oligomer-shape permissible
                    (etypecase rotamers
                      ((eql :random) :random)
                      (null nil)
                      (rotamer-indexes rotamers)
                      (function (funcall rotamers permissible))))))

(defun write-sidechain-rotamers (oligomer-shape rotamers &key monomers)
  "Mid-level API: Write sidechain rotamers to OLIGOMER-SHAPE.
Queries permissible sidechain rotamers from the current backbone state.
For lower-level control, use make-permissible-sidechain-rotamers and write-rotamers directly.
ROTAMERS can be :random, nil (first allowed), a rotamer-indexes object,
or a function that receives permissible-sidechain-rotamers and returns
a value suitable for write-rotamers.
MONOMERS optionally restricts which sidechain monomers are affected."
  (let ((permissible (make-permissible-sidechain-rotamers oligomer-shape
                                                          :monomers monomers)))
    (write-rotamers oligomer-shape permissible
                    (etypecase rotamers
                      ((eql :random) :random)
                      (null nil)
                      (rotamer-indexes rotamers)
                      (function (funcall rotamers permissible))))))

(defun build-oligomer-shape (oligomer-space &key (oligomer-index 0)
                                                 (backbone-rotamers :random)
                                                 (sidechain-rotamers :random))
  "High-level API: Build an oligomer-shape from OLIGOMER-SPACE and set its rotamers.
Handles oligomer creation, shape construction, and rotamer assignment in one call.
For mid-level control, use write-backbone-rotamers and write-sidechain-rotamers.
For low-level control, use make-permissible-backbone-rotamers, make-permissible-sidechain-rotamers,
and write-rotamers directly.
OLIGOMER-INDEX selects which oligomer to build (default 0).
BACKBONE-ROTAMERS and SIDECHAIN-ROTAMERS can each be :random,
nil (first allowed), a rotamer-indexes object, or a function that
receives the corresponding permissible-rotamers and returns a value
suitable for write-rotamers."
  (let* ((oligomer (make-oligomer oligomer-space oligomer-index))
         (oligomer-shape (make-oligomer-shape oligomer)))
    (write-backbone-rotamers oligomer-shape backbone-rotamers)
    (write-sidechain-rotamers oligomer-shape sidechain-rotamers)
    oligomer-shape))
