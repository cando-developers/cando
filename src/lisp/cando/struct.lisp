(in-package :struct)


(defgeneric find-residue-that-contains-atom (matter atom))
(defmethod find-residue-that-contains-atom ((matter chem:aggregate) atom)
  (chem:map-molecules nil
                      (lambda (m)
                        (let ((res (find-residue-that-contains-atom m atom)))
                          (when res (return-from find-residue-that-contains-atom res))))
                      matter)
  nil)

(defmethod find-residue-that-contains-atom ((matter chem:molecule) atom)
  (chem:map-residues nil
                     (lambda (r)
                       (chem:map-atoms nil
                                       (lambda (at)
                                         (when (eq atom at)
                                           (return-from find-residue-that-contains-atom r)))
                                       r))
                     matter))

(defgeneric remove-atom (matter atom))
(defmethod remove-atom ((matter chem:aggregate) atom)
  (let ((res (find-residue-that-contains-atom matter atom)))
    (remove-atom res atom)))
(defmethod remove-atom ((matter chem:molecule) atom)
  (let ((res (find-residue-that-contains matter atom)))
    (remove-atom res atom)))

(defmethod remove-atom ((residue chem:residue) atom)
  (chem:remove-atom-delete-bonds residue atom))


