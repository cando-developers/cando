
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
  (let ((root-monomer (chem:|Oligomer_O::rootMonomer| oligomer)))
    (let ((monomer-out-couplings (make-hash-table))
          (unique-ring-couplings nil))
      (loop for coupling in (chem:|Oligomer_O::couplingsAsList| oligomer)
            if (typep coupling 'chem:directional-coupling)
              do (push coupling (gethash (chem:get-source-monomer coupling) monomer-out-couplings))
            else
              do (pushnew coupling unique-ring-couplings))
      (canonical-sequence-monomer nil root-monomer monomer-out-couplings unique-ring-couplings))))



(defun recursively-build-molecule (prev-monomer
                                   prev-topology
                                   prev-residue
                                   monomer-out-couplings
                                   molecule
                                   monomers-to-residues
                                   monomers-to-topologys
                                   monomer-positions)
  (loop for out-coupling in (gethash prev-monomer monomer-out-couplings)
        for next-monomer = (chem:get-target-monomer out-coupling)
        for next-topology = (let ((next-top (chem:current-topology next-monomer)))
                              (unless (typep next-top 'chem:topology)
                                (error "Unexpected object - expected a chem:topology - got ~s" next-top))
                              next-top)
        for next-stereoisomer-name = (chem:current-stereoisomer-name next-monomer)
        for next-residue = (progn
                             (chem:build-residue-for-monomer-name next-topology next-stereoisomer-name))
        do (setf (gethash next-monomer monomers-to-residues) next-residue
                 (gethash next-monomer monomers-to-topologys) next-topology)
           (let* ((monomer-position (gethash next-monomer monomer-positions))
                  (next-residue-index (chem:content-size molecule)))
             (chem:put-matter molecule monomer-position next-residue))
           (chem:connect-residues prev-topology 
                                  prev-residue
                                  (chem:get-source-plug-name out-coupling)
                                  next-topology
                                  next-residue
                                  (chem:get-target-plug-name out-coupling))
           (recursively-build-molecule next-monomer
                                       next-topology
                                       next-residue
                                       monomer-out-couplings
                                       molecule
                                       monomers-to-residues
                                       monomers-to-topologys
                                       monomer-positions)))

(defun build-molecule (oligomer)
  (let ((root-monomer (chem:|Oligomer_O::rootMonomer| oligomer))
        (monomer-out-couplings (make-hash-table))
        (monomers-to-residues (make-hash-table))
        (monomers-to-topologys (make-hash-table))
        (ring-couplings nil)
        (monomer-positions (make-hash-table)))
    (loop for coupling in (chem:|Oligomer_O::couplingsAsList| oligomer)
          if (typep coupling 'chem:directional-coupling)
            do (push coupling (gethash (chem:get-source-monomer coupling) monomer-out-couplings))
          else
            do (pushnew coupling ring-couplings)) ; Only add ring coupling when unique
    (let* ((molecule (chem:make-molecule :mol))
           (root-topology (let ((top (chem:current-topology root-monomer)))
                            (unless (typep top 'chem:topology)
                              (error "Unexpected object - expected a chem:topology - got ~s" top))
                            top))
           (stereoisomer-name (chem:current-stereoisomer-name root-monomer))
           (root-residue (progn
                           (chem:build-residue-for-monomer-name root-topology stereoisomer-name))))
      (setf (gethash root-monomer monomers-to-residues) root-residue
            (gethash root-monomer monomers-to-topologys) root-topology)
      (let* ((monomer-list (chem:|Oligomer_O::monomersAsList| oligomer))
             (number-of-residues (length monomer-list)))
        (loop for count from 0
              for monomer in monomer-list
              do (setf (gethash monomer monomer-positions) count))
        (let* ((monomer-position (gethash root-monomer monomer-positions))
               (next-residue-index (chem:content-size molecule)))
          (chem:resize-contents molecule number-of-residues)
          (chem:put-matter molecule monomer-position root-residue))
        (recursively-build-molecule root-monomer
                                    root-topology
                                    root-residue
                                    monomer-out-couplings
                                    molecule
                                    monomers-to-residues
                                    monomers-to-topologys
                                    monomer-positions))
      ;; Now close the rings
      (loop for ring-coupling in ring-couplings
            for monomer1 = (chem:get-monomer1 ring-coupling)
            for topology1 = (chem:current-topology monomer1)
            for residue1 = (gethash monomer1 monomers-to-residues)
            for monomer2 = (chem:get-monomer2 ring-coupling)
            for topology2 = (chem:current-topology monomer2)
            for residue2 = (gethash monomer2 monomers-to-residues)
            for plug1name = (chem:get-plug1 ring-coupling)
            for plug2name = (chem:get-plug2 ring-coupling)
            do
            #+(or)(progn
                    (format *debug-io* "ring-coupling: ~a~%" ring-coupling)
                    (format *debug-io* "    topology1: ~a~%" topology1)
                    (format *debug-io* "     residue1: ~a~%" residue1)
                    (format *debug-io* "    plug1name: ~a~%" plug1name)
                    (format *debug-io* "    topology2: ~a~%" topology2)
                    (format *debug-io* "     residue2: ~a~%" residue2)
                    (format *debug-io* "    plug2name: ~a~%" plug2name))
                  (chem:connect-residues topology1
                                         residue1
                                         plug1name
                                         topology2
                                         residue2
                                         plug2name))
      (values molecule monomer-positions))))


(defun build-all-molecules (oligomer &optional (number-of-sequences (chem:|Oligomer_O::numberOfSequences| oligomer)))
  "Return a list of all built molecules for the oligomer.
The number of entries can be limited by passing a lower value for the optional argument **number-of-sequences**
than the (chem:|Oligomer_O::numberOfSequences| oligomer)."
  (loop for index from 0 below number-of-sequences
        for aggregate = (chem:make-aggregate (intern (format nil "seq-~a" index) :keyword))
        for molecule = (progn
                         (chem:|Oligomer_O::gotoSequence| oligomer index)
                         (build-molecule oligomer))
        do (chem:add-matter aggregate molecule)
           (cando:jostle aggregate)
        collect aggregate))
