
(in-package :design)


(defun canonical-sequence-monomer (coupling root monomer-out-couplings unique-ring-couplings)
  (let* ((outs (gethash root monomer-out-couplings))
         (result (list* (chem:current-stereoisomer-name root)
                        (when outs
                          (let ((sorted-outs (sort outs #'string< :key (lambda (coup) (string (chem:get-name coup))))))
                            (loop for sorted-out in sorted-outs
                                  for target-monomer = (chem:get-target-monomer sorted-out)
                                  collect (canonical-sequence-monomer
                                           sorted-out
                                           target-monomer
                                           monomer-out-couplings
                                           unique-ring-couplings)))))))
    (if coupling
      (cons (chem:get-name coupling) result)
      result)))

(defun canonical-sequence (oligomer)
  "Return a canonical sequence for the oligomer that can be compared to other oligomers to
   determine if they are equivalent"
  (let ((root-monomer (chem:root-monomer oligomer)))
    (let ((monomer-out-couplings (make-hash-table))
          (unique-ring-couplings nil))
      (loop for coupling in (chem:couplings-as-list oligomer)
            if (typep coupling 'chem:directional-coupling)
              do (push coupling (gethash (chem:get-source-monomer coupling) monomer-out-couplings))
            else
              do (pushnew coupling unique-ring-couplings))
      (canonical-sequence-monomer nil root-monomer monomer-out-couplings unique-ring-couplings))))
