(in-package :aromaticity)


(defun test-oxazole ()
  (let* ((smiles-mol  (smiles:parse-smiles-string "O1C=NC=C1")) ; Kekulé oxazole
         (mol         (smiles:chem-molecule smiles-mol))
         (agg         (let ((a (chem:make-aggregate :test)))
                        (chem:add-matter a mol)
                        a))
         (rings       (chem:identify-rings agg))
         (chem:*current-rings* rings)
         (aromaticity (chem:identify-aromatic-rings agg :smirnoff))
         (chem:*current-aromaticity-information* aromaticity))
    (let* ((atoms (chem:map-atoms 'list #'identity mol))
           (o     (find :O atoms :key #'chem:get-element))
           (n     (find :N atoms :key #'chem:get-element))
           (c2    (find-if (lambda (a)
                             (and (eq (chem:get-element a) :C)
                                  (member o (chem:atom/bonded-atoms-as-list a))
                                  (member n (chem:atom/bonded-atoms-as-list a))))
                           atoms)))
      (list :o-aromatic       (aromaticity:is-aromatic o)
            :n-aromatic       (aromaticity:is-aromatic n)
            :c2-aromatic      (aromaticity:is-aromatic c2)
            :c2-is-ar5        (aromaticity:is-aromatic-type c2 :ar5)
            :c2-by-rule5      (aromaticity:is-aromatic-rule c2 :rule5)
            :o-c2-bond        (chem:bond-order-to o c2)
            :c2-info-string   (aromaticity:aromaticity-information-as-string c2))))
  )


(defun test-aromaticity ()
  (let ((oxazole-expect `(:O-AROMATIC
                           ,(list (make-aromatic-info :aromaticity-type :ar5 :rule-name :rule5))
                           :N-AROMATIC
                           ,(list (make-aromatic-info :aromaticity-type :ar5 :rule-name :rule5))
                           :C2-AROMATIC
                           ,(list (make-aromatic-info :aromaticity-type :ar5 :rule-name :rule5))
                           :C2-IS-AR5 T
                           :C2-BY-RULE5 T
                           :O-C2-BOND :SINGLE-BOND
                           :C2-INFO-STRING
                           "(#S(AROMATICITY:AROMATIC-INFO :AROMATICITY-TYPE :AR5 :RULE-NAME :RULE5))")))
    (if (equalp (test-oxazole) oxazole-expect)
        t
        (format t "Could not match ~%~s~%~s~%" (test-oxazole) oxazole-expect))))
