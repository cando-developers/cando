(in-package #:topology-shape-key-tests)

(defun test-cache-storage ()
  (topology:define-topology skt-rotamer-cached (C (- H1) (- H2) - H3)
    :cache ((:dihedral :bb1 C) (:ring :ring C)))
  (topology:define-topology skt-rotamer-plain (C (- H1) (- H2) - H3))
  (let* ((cache (topology:make-binned-shape-key-cache '((:bb1 35) (:ring 2 13))))
         (fragment (make-instance 'topology:fragment-internals
                                  :monomer-context :test :trainer-index 0
                                  :internals-values (vector 1d0 2d0 3d0)
                                  :coordinates (vector 4d0 5d0 6d0)))
         (cached (topology:make-fragment-internals-with-shape-key-cache-from-fragment-internals
                  fragment cache))
         (declared-topology (chem:find-topology 'skt-rotamer-cached))
         (plain-topology (chem:find-topology 'skt-rotamer-plain))
         (rotamer (topology:make-backbone-rotamer-from-fragment-internals cached declared-topology))
         (plain (topology:make-backbone-rotamer-from-fragment-internals fragment plain-topology)))
    (assert (typep rotamer 'topology:backbone-with-cache-rotamer))
    (assert (typep plain 'topology:backbone-without-cache-rotamer))
    (assert (signals-error-p
             (lambda () (topology:make-backbone-rotamer-from-fragment-internals fragment declared-topology))))
    (assert (signals-error-p
             (lambda () (topology:make-backbone-rotamer-from-fragment-internals cached plain-topology))))
    (assert (signals-error-p
             (lambda () (topology:make-backbone-rotamer declared-topology
                          :cache (topology:make-binned-shape-key-cache '((:bb1 35)))))))
    (let ((container (make-instance 'topology:backbone-rotamers :monomer-context :test
                                    :rotamer-vector (vector rotamer plain))))
      (uiop:with-temporary-file (:pathname pathname)
        (cando.serialize:save-cando container pathname)
        (dolist (native '(t nil))
          (let ((restored (cando.serialize:load-cando pathname :use-native-reader native)))
            (assert (typep restored 'topology:backbone-rotamers))
            (assert (typep (aref (topology:rotamer-vector restored) 0) 'topology:backbone-with-cache-rotamer))
            (assert (typep (aref (topology:rotamer-vector restored) 1) 'topology:backbone-without-cache-rotamer))))))
    (topology:add-binned-shape-key-entry cache :later '(1))
    (assert (not (member :later (topology:binned-shape-key-names
                                (topology:backbone-shape-key-cache cached)))))
    (topology:add-binned-shape-key-entry (topology:backbone-shape-key-cache cached) :fragment-only '(2))
    (assert (not (member :fragment-only (topology:binned-shape-key-names
                                        (topology:backbone-shape-key-cache rotamer)))))
    (setf (aref (topology:internals-values fragment) 0) 99d0)
    (assert (= 1d0 (aref (topology:internals-values cached) 0)))
    (dolist (object (list cached rotamer))
      (let ((original (topology:backbone-shape-key-cache object)))
        (assert (signals-error-p (lambda () (setf (topology:backbone-shape-key-cache object)
                                                  (topology:make-binned-shape-key nil)))))
        (assert (eq original (topology:backbone-shape-key-cache object))))
      (uiop:with-temporary-file (:pathname pathname)
        (cando.serialize:save-cando object pathname)
        (dolist (native '(t nil))
          (let* ((restored (cando.serialize:load-cando pathname :use-native-reader native))
                 (restored-cache (topology:backbone-shape-key-cache restored)))
            (assert (eq (class-of object) (class-of restored)))
            (assert (eq (class-of restored-cache) (find-class 'topology:binned-shape-key-cache)))
            (assert (equal (topology:binned-shape-key-entries (topology:backbone-shape-key-cache object))
                           (topology:binned-shape-key-entries restored-cache))))))))
  (format t "CACHE-STORAGE-TESTS-PASSED~%")
  t)
