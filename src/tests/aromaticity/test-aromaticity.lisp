;;; Common code

(defun fake-aromatic-info (&key mol atom-names aromaticity-type (fake-test-name :fake-test))
  "Return fake aromaticity-information in a hash-table.  You provide the atom-names and
what aromaticity information you want for them and this function returns a hash-table
with those atoms with that aromaticity-information."
  (let ((aromaticity-info (make-hash-table)))
    (loop with aromatic-atoms = nil
          for atom-name in atom-names
          for atom = (chem:first-atom-with-name mol atom-name)
          do (aromaticity::set-aromaticity-type aromaticity-info
                                                atom
                                                aromaticity-type
                                                fake-test-name
                                                aromatic-atoms))
    aromaticity-info))

(defun fake-test-aromaticity-rule (&key mol (fake-aromatic-info (make-hash-table))
                                     rule (fake-aromaticity-type :fake-type) (fake-test-name :fake-test))
  "Test for aromaticity using fake aromaticity information"
  (let* ((chem:*current-rings* (chem:identify-rings mol))
         (chem:*current-aromaticity-information* fake-aromatic-info)
         (mol-graph (chem:make-molecule-graph mol))
         (chem:*current-aromaticity-information* fake-aromatic-info)
         (new-ai (make-hash-table)))
    (aromaticity::exhaustively-apply-aromatic-rule new-ai
                                                   mol-graph
                                                   rule
                                                   fake-aromaticity-type
                                                   fake-test-name)
    new-ai))

(defun aromaticity-matches (molecule-name molecule aromaticity-pairs)
  (let ((aromaticity-info (chem:identify-aromatic-rings molecule)))
    (loop for type-rule-atom-names in aromaticity-pairs
          for aromaticity-type = (first type-rule-atom-names)
          for aromaticity-rule = (second type-rule-atom-names)
          for atom-names = (cddr type-rule-atom-names)
          do (loop for atom-name in atom-names
                   for atom = (chem:first-atom-with-name molecule atom-name)
                   unless (aromaticity:is-aromatic-type atom aromaticity-type aromaticity-info)
                     do (error "mismatch in aromaticity information for atom ~a in molecule ~a - it should be ~a and it is ~a" atom molecule-name aromaticity-type (aromaticity:is-aromatic atom aromaticity-info))
                   unless (aromaticity:is-aromatic-rule atom aromaticity-rule aromaticity-info)
                     do (error "mismatch in aromaticity information for atom ~a in molecule ~a - it should be ~a and it is ~a" atom molecule-name aromaticity-rule (aromaticity:is-aromatic atom aromaticity-info))))))
                   
;;; Try some contrived examples
(defparameter *agg* (load-chem-draw-aggregate "polyaromatic.cdxml"))
(chem:fill-in-implicit-hydrogens *agg*)

(defmacro with-aromaticity-test (name aggregate expect)
  (let ((mol-name-gs (gensym))
        (mol-gs (gensym))
        (aromatic-info-gs (gensym)))
    `(let* ((,mol-name-gs ,name)
            (,mol-gs (chem:content-with-name ,aggregate ,mol-name-gs))
            (chem:*current-rings* (chem:identify-rings ,mol-gs))
            (,aromatic-info-gs (chem:identify-aromatic-rings ,mol-gs)))
       (aromaticity-matches ,mol-name-gs ,mol-gs ,expect))))


;;; If everything is working correctly then the following tests should each succeed and return NIL.
;;; If here is a problem with the aromaticity pattern recognition then one or more of the tests will signal an error.
(with-aromaticity-test :ring1 *agg* '((:ar6 :rule1 :c1 :c2 :c7 :c8 :c9 :c10)))
(with-aromaticity-test :ring2 *agg* '((:ar6 :rule1 :c1 :c2 :c3 :c4 :c5 :c6)
                                      (:ar6 :rule2 :c7 :c8 :c9 :c10)))
(with-aromaticity-test :ring4 *agg* '((:ar6 :rule1 :c1 :c2 :c3 :c4 :c5 :c6)
                                      (:ar6 :rule2 :c7 :c8 :c9 :c10
                                       :c11 :c12 :c13 :c14
                                       :C15 :c16 :c17 :c18
                                       )))

(with-aromaticity-test :indole *agg* '((:ar6 :rule1 :c4 :c5 :c6 :c7 :c8 :c9)
                                       (nil nil :N1 :c2 :c3)))


#|


(let ((chem:*current-rings* (chem:identify-rings *ring2*)))
  (chem:identify-aromatic-rings *ring2*))

(let ((chem:*current-rings* (chem:identify-rings *ring4*)))
  (chem:identify-aromatic-rings *ring4*))

(let ((chem:*current-rings* (chem:identify-rings *ring-all*)))
  (chem:identify-aromatic-rings *ring-all*))


(let ((chem:*current-rings* (chem:identify-rings *indole*)))
  (chem:identify-aromatic-rings *indole*))

chem:*verbose*

(fake-test-aromaticity-rule :mol *ring-fake*
                            :rule aromaticity::*rule1*)

(fake-aromatic-info :mol *fm* :atom-names '(:c1 :c2) :aromaticity-type :ar6)

(fake-test-aromaticity-rule :mol *ring-fake*
                            :fake-aromatic-info (fake-aromatic-info :mol *ring-fake*
                                                                    :atom-names '(:c1 :c2)
                                                                    :aromaticity-type :ar6)
                            :rule aromaticity::*rule2*)

(defparameter *fake-ai-c1-c2-ar6*  (fake-aromatic-info :mol *ring-fake*
                                                       :atom-names '(:c1 :c2)
                                                       :aromaticity-type :ar6))

(defparameter *ar6rule* (chem:make-chem-info-graph (chem:compile-smarts "[<ar6>:1]:[<ar6>:2]" :tests aromaticity::*artests*)))

(let ((chem:*current-rings* (chem:identify-rings *fm*))
      (chem:*verbose* 1))
  (fake-test-aromaticity-rule :mol *fm*
                              :fake-aromatic-info *fake-ai-c1-c2-ar6*
                              :rule *ar6rule*))

(let ((chem:*current-rings* (chem:identify-rings *fm*))
      (chem:*current-aromaticity-information* *fake-ai-c1-c2-ar6*)
      (chem:*verbose* 1))
  (format t "chem:*current-rings* -> ~s~%" chem:*current-rings*)
  (format t "chem:*current-aromaticity-information* -> ~s~%" chem:*current-aromaticity-information*)
  (chem:in-aromatic-bond (chem:first-atom-with-name *fm* :c1) (chem:first-atom-with-name *fm* :c2)))

(fake-test-aromaticity-rule :mol *fm*
                            :fake-aromatic-info *fake-ai-c1-c2-ar6*
                            :rule *myrule2*)

(aromaticity:is-aromatic-type (chem:first-atom-with-name *fm* :c7) :ar6 *fake-ai-c1-c2-ar6*)
(aromaticity:is-aromatic-type (chem:first-atom-with-name *fm* :c1) :ar6 *fake-ai-c1-c2-ar6*)
(aromaticity:is-aromatic-type (chem:first-atom-with-name *fm* :c2) :ar6 *fake-ai-c1-c2-ar6*)

(trace aromaticity::is-aromatic-type )

                            
                            
;;; Try a compilicated example
(defparameter *agg* (load-chem-draw-aggregate "ligand8a.cdxml"))
(chem:fill-in-implicit-hydrogens *agg*)
(set-stereoisomer-func (cando:gather-stereocenters *agg*) (lambda (c) :R))
(defparameter *mol* (chem:content-at *agg* 0))
(set-force-field *mol* :smirnoff)
(save-amber-parm *agg* "/tmp/test.parm" "/tmp/test.crd")



|#
