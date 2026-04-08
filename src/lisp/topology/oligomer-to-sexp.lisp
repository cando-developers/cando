(in-package :topology)

;;; Generate an s-expression representation of an oligomer that shows
;;; the monomer names and coupling names in tree form.
;;; The output format matches the input to make-oligomer-space.

(defun find-coupling-name-between (parent-monomer child-monomer)
  "Find the coupling name for the coupling from PARENT-MONOMER to CHILD-MONOMER."
  (maphash (lambda (plug-name coupling)
             (declare (ignore plug-name))
             (when (and (typep coupling 'topology:directional-coupling)
                        (eq (topology:source-monomer coupling) parent-monomer)
                        (eq (topology:target-monomer coupling) child-monomer))
               (return-from find-coupling-name-between (topology:name coupling))))
           (topology:couplings parent-monomer))
  (error "No coupling found from ~s to ~s" parent-monomer child-monomer))

(defun oligomer-to-sexp-impl (monomer oligomer out-monomers)
  "Recursively build the s-expression for MONOMER and its children.
Each child becomes (:coupling-name child-name child-subtrees...).
Leaf nodes with no children become (:coupling-name child-name)."
  (let* ((name (topology:current-stereoisomer-name monomer oligomer))
         (outs (gethash monomer out-monomers))
         (child-sexps nil))
    (dolist (child outs)
      (let* ((coupling-name (find-coupling-name-between monomer child))
             (child-name (topology:current-stereoisomer-name child oligomer))
             (grandchildren (gethash child out-monomers))
             (child-subtrees (loop for grandchild in grandchildren
                                   collect (let ((gc-coupling (find-coupling-name-between child grandchild)))
                                             (oligomer-to-sexp-subtree gc-coupling grandchild oligomer out-monomers)))))
        (push (list* coupling-name child-name child-subtrees) child-sexps)))
    (if child-sexps
        (list* name (nreverse child-sexps))
        name)))

(defun oligomer-to-sexp-subtree (coupling-name monomer oligomer out-monomers)
  "Build (:coupling-name monomer-name subtrees...) for a child monomer."
  (let* ((name (topology:current-stereoisomer-name monomer oligomer))
         (outs (gethash monomer out-monomers))
         (child-subtrees (loop for child in outs
                               collect (let ((c-name (find-coupling-name-between monomer child)))
                                         (oligomer-to-sexp-subtree c-name child oligomer out-monomers)))))
    (list* coupling-name name child-subtrees)))

(defun oligomer-to-sexp (oligomer)
  "Generate an s-expression representation of an oligomer-shape showing
monomer names and coupling names in tree form.
The output format matches the input to make-oligomer-space.

Example output:
  (:cpent :dkp (:pro4sr :side :3ba) :dkp (:pro4ss :side :4pr))"
  (multiple-value-bind (monomer-shape-vector monomer-shape-info-vector monomer-shape-to-index
                        root in-monomers out-monomers monomer-shape-map)
      (make-monomer-shape-parts oligomer)
    (declare (ignore monomer-shape-vector monomer-shape-info-vector monomer-shape-to-index in-monomers monomer-shape-map))
    (oligomer-to-sexp-impl root oligomer out-monomers)))

(defun oligomer-to-string (oligomer)
  "Generate a compact string representation of an oligomer-shape.
Shows monomer names separated by coupling names."
  (let ((sexp (oligomer-to-sexp oligomer)))
    (let ((*print-pretty* nil)
          (*print-case* :downcase))
      (format nil "~s" sexp))))
