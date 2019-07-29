(in-package :fep)

(defun only-residue (molecule)
  (unless (= (chem:content-size molecule) 1)
    (error "The molecule ~a must have only one residue" molecule))
  (chem:content-at molecule 0))

(defun fep-calculation-from-ligands (ligands)
  (let* ((fep-structures (loop for molecule in ligands
                              for residue = (only-residue molecule)
                              collect (make-instance 'simple-fep-structure
                                                     :name (chem:get-name molecule)
                                                     :drawing molecule
                                                     :core-residue residue
                                                     :core-residue-name (chem:get-name residue)
                                                     :core-atoms (chem:all-atoms-as-list molecule nil)
                                                     :molecule molecule
                                                     :net-charge 0.0 #|Not always!|#)))
         (calc (make-instance 'fep:fep-calculation
                              :ligands fep-structures
                              :mask-method :closest)))
    calc))
                              
(defmethod calculate-masks (source target (method (eql :closest)))
  "Calculate the mask using the McGregor algorithm to find the largest common subgraph
   and what remains of each molecule are the differences used to construct the mask"
  (let ((source-molecule (molecule source))
        (target-molecule (molecule target)))
    (multiple-value-bind (equivs diff1 diff2)
        (molecule-graph:compare-molecules source-molecule target-molecule)
      (unless (= (chem:content-size source-molecule) 1)
        (error "The source molecule ~a has more than one residue"))
      (unless (= (chem:content-size target-molecule) 1)
        (error "The target molecule ~a has more than one residue"))
      (let ((source-atom-names (mapcar 'chem:get-name diff1))
            (target-atom-names (mapcar 'chem:get-name diff2)))
        (make-instance 'ti-mask
                       :equivalent-atom-names (let (pairs)
                                                (maphash (lambda (pair dummy)
                                                           (declare (ignore dummy))
                                                           (let ((source-atom (car pair))
                                                                 (target-atom (cdr pair)))
                                                             (let ((source-name (chem:get-name source-atom))
                                                                   (target-name (chem:get-name target-atom)))
                                                               (push (cons source-name target-name) pairs))))
                                                         equivs)
                                                pairs)
                       :source source
                       :source-timask-residue-index 1
                       :source-scmask-atom-names source-atom-names
                       :target target
                       :target-timask-residue-index 2
                       :target-scmask-atom-names target-atom-names)))))



(defmethod connect-job-nodes (calculation (pairs list) &key)
  (loop for pair in pairs
        for source-name = (first pair)
        for target-name = (second pair)
        for source-node = (find source-name (nodes (jobs calculation)) :test #'eq :key #'name)
        for target-node = (find target-name (nodes (jobs calculation)) :test #'eq :key #'name)
        do (unless source-node
             (error "Could not find find ~a in nodes" source-name))
        do (unless target-node
             (error "Could not find find ~a in nodes" target-name))
        do (add-job-edge source-node target-node calculation)))
