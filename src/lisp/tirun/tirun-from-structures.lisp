(in-package :tirun)

(defun only-residue (molecule)
  (unless (= (chem:content-size molecule) 1)
    (error "The molecule ~a must have only one residue" molecule))
  (chem:content-at molecule 0))

(defun short-residue-name (index)
  (let ((label (format nil "!~36,2,'0r" index)))
    (intern label :keyword)))

(defun bonded-atoms (atm diffs)
  (let (heavy-match heavy-diff hydrogens)
    (loop for bond in (chem:bonds-as-list atm)
          for other-atm = (chem:bond/get-other-atom bond atm)
          for element = (chem:get-element other-atm)
          for atomic-number = (chem:atomic-number-for-element element)
          do (if (= atomic-number 1)
                 (push other-atm hydrogens)
                 (if (member atm diffs)
                     (push other-atm heavy-diff)
                     (push other-atm heavy-match))))
    (values heavy-match heavy-diff hydrogens)))

(defun find-closest-hydrogen (atm hydrogens)
  (let* ((closest-hydrogen (first hydrogens))
         (closest-distance-squared (chem:distance-squared-between-two-atoms atm closest-hydrogen)))
    (loop for hyd in (cdr hydrogens)
          for dist-squared = (chem:distance-squared-between-two-atoms atm hyd)
          when (< dist-squared closest-distance-squared)
            do (setf closest-distance-squared dist-squared
                     closest-hydrogen hyd))
    (let ((remaining-hydrogens (loop for hyd in hydrogens
                                     unless (eq hyd closest-hydrogen)
                                       collect hyd)))
      (values closest-hydrogen remaining-hydrogens))))


(defun match-hydrogens (hydrogens1 hydrogens2)
  (let ((equivs nil)
        hyd2)
    (loop for hyd1 in hydrogens1
          do (multiple-value-setq (hyd2 hydrogens2)
               (find-closest-hydrogen hyd1 hydrogens2))
          do (push (cons hyd1 hyd2) equivs))
    equivs))

(defun integrate-hydrogens-into-comparison (heavy-equivs heavy-diff1 heavy-diff2)
;;; magically add hydrogens to the comparison
  (flet ((push-all (entries place)
           (loop for atm in entries
                 do (push atm place))
           place))
    (let ((equivs heavy-equivs)
          (diff1 heavy-diff1)
          (diff2 heavy-diff2))
      (loop for pair in heavy-equivs
            for heavy1 = (car pair)
            for heavy2 = (cdr pair)
            do (multiple-value-bind (heavy1-match heavy1-diff hydrogen1)
                   (bonded-atoms heavy1 heavy-diff1)
                 (multiple-value-bind (heavy2-match heavy2-diff hydrogen2)
                     (bonded-atoms heavy2 heavy-diff2)
                   ;; If the number of connected hydrogens are different then
                   ;; consider them all different
                   (if (/= (length hydrogen1) (length hydrogen2))
                       (progn
                         (setf diff1 (push-all hydrogen1 diff1)
                               diff2 (push-all hydrogen2 diff2)))
                       (cond
                         ;; If there are no unmatched connected heavy atoms
                         ;; then any connected hydrogens are matched across
                         ;; the two molecules
                         ((= 0 (length heavy1-diff) (length heavy2-diff))
                          (setf equivs (push-all (match-hydrogens hydrogen1 hydrogen2) equivs)))
                         ;; Otherwise the hydrogens are all considered different
                         (t
                          (setf diff1 (push-all hydrogen1 diff1)
                                diff2 (push-all hydrogen2 diff2))))))))
      (loop for heavy1 in heavy-diff1
            do (multiple-value-bind (heavy-match heavy-diff hydrogens)
                   (bonded-atoms heavy1 heavy-diff1)
                 (setf diff1 (push-all hydrogens diff1))))
      (loop for heavy2 in heavy-diff2
            do (multiple-value-bind (heavy-match heavy-diff hydrogens)
                   (bonded-atoms heavy2 heavy-diff2)
                 (setf diff2 (push-all hydrogens diff2))))
      (values equivs diff1 diff2))))

(defun calculate-masks-for-molecules (source-molecule target-molecule)
  "Calculate the mask using the max-clique algorithm to find the largest common subgraph
   and what remains of each molecule are the differences used to construct the mask"
  (multiple-value-bind (heavy-equivs heavy-diff1 heavy-diff2)
      (molecule-graph.max-clique:compare-molecules-return-largest-connected-cluster-of-equivalent-atoms
       source-molecule
       target-molecule
       :atom-match-callback (lambda (atom1 atom2)
                              (chem:atom-within-angstroms atom1 atom2 0.3))
       :remove-partial-rings t
       :exclude-hydrogens t)
    (integrate-hydrogens-into-comparison heavy-equivs heavy-diff1 heavy-diff2)))

(defmethod calculate-masks (source target (method (eql :closest)))
  "Calculate the mask using the max-clique algorithm to find the largest common subgraph
   and what remains of each molecule are the differences used to construct the mask"
  (let ((source-molecule (molecule source))
        (target-molecule (molecule target)))
    (multiple-value-bind (equivs diff1 diff2)
        (calculate-masks-for-molecules source-molecule target-molecule)
      (unless (= (chem:content-size source-molecule) 1)
        (error "The source molecule ~a has more than one residue"))
      (unless (= (chem:content-size target-molecule) 1)
        (error "The target molecule ~a has more than one residue"))
      (let ((source-atom-names (mapcar 'chem:get-name diff1))
            (target-atom-names (mapcar 'chem:get-name diff2)))
        (make-instance 'ti-mask
                       :equivalent-atom-names (let (pairs)
                                                (loop for pair in equivs
                                                      do (let ((source-atom (car pair))
                                                               (target-atom (cdr pair)))
                                                           (let ((source-name (chem:get-name source-atom))
                                                                 (target-name (chem:get-name target-atom)))
                                                             (push (cons source-name target-name) pairs))))
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
