(in-package :cando)



(defun do-stereochemistry-from-structure (molecule ht)
  (let ((cip (chem:make-cip-prioritizer))
        (cip-ht (make-hash-table)))
    (chem:assign-cahn-ingold-prelog-priority-to-atoms-relative-priority cip molecule cip-ht)
    (chem:do-atoms (atm molecule)
      (let ((conf (chem:calculate-stereochemical-configuration atm cip-ht)))
        (unless (eq conf :undefined-configuration)
          (setf (gethash atm ht) conf))))
    cip-ht))

(defgeneric stereochemistry-from-structure (matter)
  (:documentation "Return an alist (atom . configuration-symbol) calculated from the structure"))

(defmethod stereochemistry-from-structure ((molecule chem:molecule))
  (let* ((ht (make-hash-table))
         (cip-ht (do-stereochemistry-from-structure molecule ht)))
    (values ht cip-ht)))

