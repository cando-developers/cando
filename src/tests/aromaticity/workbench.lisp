;;; Common code
(progn
  (setf *default-pathname-defaults* #P"~/Development/cando/extensions/cando/src/tests/aromaticity/")
  (setup-default-paths)
  (easy-gaff)
  (load-smirnoff-params #P"~/Development/openforcefield/openforcefield/data/forcefield/smirnoff99Frosst.offxml")
  )



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

;;; Try a compilicated example
(defparameter *agg* (load-chem-draw-aggregate "ligand8a.cdxml"))
(chem:fill-in-implicit-hydrogens *agg*)
(set-stereoisomer-func (cando:gather-stereocenters *agg*) (lambda (c) :R))
(defparameter *mol* (chem:content-at *agg* 0))
(set-force-field *mol* :smirnoff)
(save-amber-parm *agg* "/tmp/test.parm" "/tmp/test.crd")



;;; Try some contrived examples
(defparameter *agg* (load-chem-draw-aggregate "polyaromatic.cdxml"))
(chem:fill-in-implicit-hydrogens *agg*)

(defparameter *ring4* (chem:content-with-name *agg* :ring4))
(defparameter *ring2* (chem:content-with-name *agg* :ring2))
(defparameter *ring-fake* (chem:content-with-name *agg* :ring-fake))
(defparameter *ring-all* (chem:content-with-name *agg* :ring-all))
(defparameter *indole* (chem:content-with-name *agg* :indole))


(chem:identify-aromatic-rings *ring-fake*)
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

                            
                            
