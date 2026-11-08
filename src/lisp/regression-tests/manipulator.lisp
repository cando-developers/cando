(in-package #:clasp-tests)

(defun check-ring-dihedral-classification (memberships-b memberships-c)
  ;; Use distinct objects for atoms and rings, just like the production table.
  (let ((table (make-hash-table :test #'eq))
        (atom-b (gensym "B"))
        (atom-c (gensym "C")))
    (setf (gethash atom-b table) memberships-b
          (gethash atom-c table) memberships-c)
    (topology::inter-ring-dihedral-allowed-p atom-b atom-c table)))

(test-true acyclic-dihedral-driver-allowed
  (check-ring-dihedral-classification nil nil))

(test-true exocyclic-dihedral-driver-allowed
  (let ((ring (list :ring)))
    (and (check-ring-dihedral-classification (list ring) nil)
         (check-ring-dihedral-classification nil (list ring)))))

(test-true separate-ring-dihedral-driver-allowed
  (let ((ring-a (list :a)) (ring-b (list :b)))
    (and (check-ring-dihedral-classification (list ring-a) (list ring-b))
         (check-ring-dihedral-classification (list ring-b) (list ring-a)))))

(test-true within-ring-dihedral-driver-excluded
  (let ((ring (list :ring)))
    (not (check-ring-dihedral-classification (list ring) (list ring)))))

(test-true fused-ring-dihedral-driver-excluded
  (let ((ring-a (list :a)) (ring-b (list :b)) (shared (list :shared)))
    (not (check-ring-dihedral-classification
          (list ring-a shared) (list ring-b shared)))))
