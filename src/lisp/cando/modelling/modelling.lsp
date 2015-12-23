(in-package :modelling)

;;; Identifies atoms with property property-name and sorts them by the value
;;; Eg: You could have three atoms with property (:fix a) (:fix b) (:fix c)
;;; The atoms will be returned in the order a,b,c
(defun sorted-atoms-with-property (agg property-name)
  (let* ((sorted-atoms
          (sort (cando:atoms-with-property agg property-name) #'string<
                :key (lambda (a)
                       (string (getf (chem:get-properties a) property-name))))))
    sorted-atoms))

(defun ensure-ovector3 (v)
  (typecase v
    (chem:ovector3 v)
    (vector (chem:make-ovector3 (elt v 0) (elt v 1) (elt v 2)))
    (cons (chem:make-ovector3 (first v) (second v) (third v)))
    (t (error "Add support for converting ~a to chem:ovector3" v))))


;;; Create anchor restraints for the atoms in the fix-atoms list
;;; the anchors are a list of ovector3's or lists of three numbers or v
(defun anchor-atoms (fix-atoms anchors &optional (weight 1.0))
  (loop for atom in fix-atoms
     for pos in anchors
     do (let ((restraint-anchor (make-cxx-object 'chem:restraint-anchor :atom atom :position pos :weight weight)))
          (chem:add-restraint atom restraint-anchor))))


