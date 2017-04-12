
(in-package :leap.cando-database)
            
(defclass leap-cando-database ()
  ())

(eval-when (:load-toplevel :execute)
  (defparameter chem:*cando-database* (make-instance 'leap-cando-database)))

(defmethod chem:constitution-for-name-or-pdb ((db leap-cando-database) name)
  (let* ((topology (leap.core:lookup-variable name))
         (constitution (chem:get-constitution topology)))
    constitution))

(defmethod chem:monomer-name-for-name-or-pdb ((db leap-cando-database) name)
  name)

(defmethod chem:pdb-name-for-name-or-pdb ((db leap-cando-database) name)
  :UNK)
