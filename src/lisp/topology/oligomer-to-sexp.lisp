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

(defun oligomer-to-sexp-impl (monomer oligomer oligomer-shape)
  "Recursively build the s-expression for MONOMER and its children."
  (let* ((name (topology:current-stereoisomer-name monomer oligomer))
         (out-monomers (gethash monomer (topology:out-monomers oligomer-shape)))
         (children nil))
    (dolist (child out-monomers)
      (let ((coupling-name (find-coupling-name-between monomer child)))
        (push coupling-name children)
        (let ((child-sexp (oligomer-to-sexp-impl child oligomer oligomer-shape)))
          (push child-sexp children))))
    (if children
        ;; Has children: (name coupling child coupling child ...)
        (list* name (nreverse children))
        ;; Leaf: just the name
        name)))

(defun oligomer-to-sexp (oligomer-shape)
  "Generate an s-expression representation of an oligomer-shape showing
monomer names and coupling names in tree form.
The output format matches the input to make-oligomer-space.

Example output:
  (:cpent :dkp (:pro4sr :side :3ba) :dkp (:pro4ss :side :4pr))"
  (let ((root (topology:the-root-monomer oligomer-shape))
        (oligomer (topology:oligomer oligomer-shape)))
    (oligomer-to-sexp-impl root oligomer oligomer-shape)))

(defun oligomer-to-string (oligomer-shape)
  "Generate a compact string representation of an oligomer-shape.
Shows monomer names separated by coupling names."
  (let ((sexp (oligomer-to-sexp oligomer-shape)))
    (let ((*print-pretty* nil)
          (*print-case* :downcase))
      (format nil "~s" sexp))))
