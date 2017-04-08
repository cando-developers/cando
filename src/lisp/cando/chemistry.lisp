
(in-package :cando)

(defun compile-smarts (smarts-str)
  "* Arguments
- smarts-str :: A string
* Description
Compile the string into smarts code."
  (let ((ci (core:make-cxx-object 'chem:chem-info)))
    (chem:compile-smarts ci smarts-str)
    ci))

(defparameter *amide-smarts* (compile-smarts "N[C&D3]~O"))
(defparameter *imidazole-n* (compile-smarts "N([C&D3])[C&D3]"))
(defun calculate-hybridization (atom)
  "* Arguments
- atom : chem:atom
* Description
Calculate the hybridization of the atom based on the number of neighbors
and the element."
  (let ((neighbors (chem:number-of-bonds atom)))
    (case (chem:get-element atom)
      (:C (case neighbors
            (2 :sp)
            (3 :sp2)
            (otherwise :sp3)))
      (:N (cond
            ((or (chem:matches *amide-smarts* atom)
                 (chem:matches *imidazole-n* atom))
             :sp2)
            ((= neighbors 2) :sp2)
            (t :sp3)))
      ((:O :S)
       (cond
         ((= neighbors 1) :sp2)
         (t :sp3))))))
