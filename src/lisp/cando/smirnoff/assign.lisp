(in-package :smirnoff)



(defmethod chem:assign-molecular-force-field-parameters (energy-function (force-field combined-force-field) molecule)
  ;;; The force fields need to be searched in reverse order because by smirnoff
  ;;; later rules over-ride earlier ones.
  (let ((molecule-graph (chem:make-molecule-graph-from-molecule molecule))
        (reversed-force-fields (reverse (chem:force-fields-as-list force-field))))
    (let ((stretch-terms (make-hash-table :test #'equal)))
      (loop for force-field in reversed-force-fields
            do (loop for term in (harmonic
            for stretch-compiled-smarts = 
        
    (let ((
  (error "Finish implementing chem:assign-molecular-force-field-parameters for smirnoff"))
