(in-package :kin)

(defclass fold-tree ()
  ((ataggregate :initarg :ataggregate :accessor aggregate)))


(defun make-fold-tree ()
  (make-instance 'fold-tree
                 :ataggregate (make-instance 'ataggregate)))

(defmethod resize-molecules ((fold-tree fold-tree) num-molecules)
  (resize-molecules (aggregate fold-tree) num-molecules))

(defun build-atmolecule-using-oligomer (fold-tree molecule-id oligomer)
  (let ((atmolecule (atmolecule (ataggregate fold-tree) molecule-id)))
    (resize-residues molecule (chem:oligomer/number-of-monomers oligomer))
    (build-using-oligomer molecule molecule-id oligomer)
    molecule
    ))

(eval-when (:compile-toplevel :execute :load-toplevel)
  (export '(fold-tree)
          )
  )


