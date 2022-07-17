
(in-package :leap.cando-database)




(defmethod chem:constitution-for-name-or-pdb ((db chem:cando-database) name)
  (let* ((topology (leap.core:lookup-variable name))
         (constitution (chem:topology/get-constitution topology)))
    constitution))

(defmethod chem:constitution-for-name-or-pdb ((db null) name)
  (let* ((topology (leap.core:lookup-variable name))
         (constitution (chem:topology/get-constitution topology)))
    constitution))

(defmethod chem:monomer-name-for-name-or-pdb ((db chem:cando-database) name)
  name)

(defmethod chem:pdb-name-for-name-or-pdb ((db chem:cando-database) name)
  :UNK)

#|
(defclass leap-cando-database ()
  ())

(defmethod chem:constitution-for-name-or-pdb ((db leap-cando-database) name)
  (let* ((topology (leap.core:lookup-variable name))
         (constitution (chem:get-constitution topology)))
    constitution))

(defmethod chem:constitution-for-name-or-pdb ((db null) name)
  (let* ((topology (leap.core:lookup-variable name))
         (constitution (chem:get-constitution topology)))
    constitution))

(defmethod chem:monomer-name-for-name-or-pdb ((db leap-cando-database) name)
  name)

(defmethod chem:pdb-name-for-name-or-pdb ((db leap-cando-database) name)
  :UNK)
|#
