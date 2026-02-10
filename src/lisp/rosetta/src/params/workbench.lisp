(ql:quickload :rosetta-params)

(progn
  (defparameter aas (rosetta.params:build-capped-residues-in-directories
                     :directories
                     (list
                      (list (probe-file #P"/home/meister/common-lisp/tlt/lk-solvation-model/rosetta-data/peptoids/")
                            (list "011" "014" "308" "309" "332" "403" "506" "631"
                                  "703")))))
  (defparameter rule-pn #p"/home/meister/Development/cando/extensions/cando/src/lisp/rosetta/data/generated-peptoid-type-rules.lisp")
  )

(defun assign-lk-types (rosetta-capped-residue params &key verbose)
  (let* ((matter (rosetta.params:rosetta-capped-residue-molecule rosetta-capped-residue))
         (rings (chem:identify-rings matter))
         (chem:*current-rings* rings))
    (let* ((ff-name (or (chem:force-field-name matter) :mdl))
           (aromaticity (chem:identify-aromatic-rings matter ff-name))
           (chem:*current-aromaticity-information* aromaticity))
      (when verbose
        (format t "aromaticity = ~s~%" aromaticity))
      (rosetta.params:assign-lk-types rosetta-capped-residue params))
    matter))

(defun one-mol (name)
  (loop for aasx in aas
        for name = (rosetta.params:rosetta-capped-residue-name aasx)
        when (member name '(:|008|))
          do (return-from one-mol (rosetta.params:rosetta-capped-residue-molecule aasx))))

(progn
  (defparameter tas (rosetta.params:load-lk-type-assignments :pathname rule-pn))
  (loop for aasx in aas
        for name = (rosetta.params:rosetta-capped-residue-name aasx)
        when (member name '(:|008|))
          do (assign-lk-types aasx tas :verbose t)))

(progn
  (defparameter tas (rosetta.params:load-lk-type-assignments :pathname rule-pn))
  (loop for aasx in aas
        for name = (rosetta.params:rosetta-capped-residue-name aasx)
        unless (member name '(:CYV :CYZ))
          do (assign-lk-types aasx tas)))


