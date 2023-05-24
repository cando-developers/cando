(in-package :topology)

(defun safe-atom-with-name (residue name)
  (chem:atom-with-name residue name))

(defun build-residue-for-monomer-name (topology monomer-name)
  (let* ((residue (chem:make-residue monomer-name))
         (constitution (topology:constitution topology))
         (constitution-atoms (topology:constitution-atoms constitution))
         (stereoisomer-atoms (topology:stereoisomer-atoms (topology:stereoisomer topology)))
         (atoms (make-array (length constitution-atoms)))
         (named-atoms (make-hash-table)))
    ;; Add atoms to the residue
    (loop for index below (length constitution-atoms)
          for constitution-atom = (elt constitution-atoms index)
          for stereoisomer-atom = (elt stereoisomer-atoms index)
          for atm-name = (topology:atom-name constitution-atom)
          for atm-element = (topology:element constitution-atom)
          for atm = (chem:make-atom atm-name atm-element)
          for atm-type = (when (slot-boundp stereoisomer-atom 'topology:atom-type)
                           (topology:atom-type stereoisomer-atom))
          for atm-charge = (if (slot-boundp stereoisomer-atom 'topology:atom-charge)
                               (topology:atom-charge stereoisomer-atom)
                               0.0)
          do (chem:atom/set-atom-type atm atm-type)
          do (chem:atom/set-charge atm atm-charge)
          do (chem:add-matter residue atm)
          do (setf (elt atoms index) atm
                   (gethash atm-name named-atoms) atm))
    ;; Create bonds between atoms
    (loop for index below (length constitution-atoms)
          for constitution-atom = (elt constitution-atoms index)
          for atm1-index = (topology:index constitution-atom)
          for atm1 = (elt atoms atm1-index)
          do (loop for bond in (topology:bonds constitution-atom)
                   for atm2-index = (topology:to-atom-index bond)
                   for atm2 = (elt atoms atm2-index)
                   when (< atm1-index atm2-index)
                     do (chem:bond-to atm1 atm2 (topology:order bond))))
    ;; Define the :chiral atoms
    (let* ((stereoisomer (stereoisomer topology))
           (stereoconfigurations (stereoconfigurations stereoisomer)))
      (loop for stereoconfiguration in stereoconfigurations
            for atom-name = (atom-name stereoconfiguration)
            for configuration = (configuration stereoconfiguration)
            for atm = (gethash atom-name named-atoms)
            do (chem:set-stereochemistry-type atm :chiral)
            do (chem:set-configuration atm configuration))
      ;; assign stereochem for prochiral atoms so they aren't random
      (chem:do-atoms (atm residue)
        (when (= (chem:number-of-bonds atm) 4)
          (chem:set-hybridization atm :sp3)
          (when (not (eq (chem:get-stereochemistry-type atm) :chiral))
            (chem:set-stereochemistry-type atm :prochiral)
            (chem:set-configuration atm :left-handed))))
      ;; Add restraints to the residue
      ;;  dihedral restraints is all we can specify here for the time being
      (loop for restraint in (restraints topology)
            for rr = (etypecase restraint
                       (dihedral-restraint
                        (let* ((atom1 (safe-atom-with-name residue (atom1-name restraint)))
                               (atom2 (safe-atom-with-name residue (atom2-name restraint)))
                               (atom3 (safe-atom-with-name residue (atom3-name restraint)))
                               (atom4 (safe-atom-with-name residue (atom4-name restraint))))
                          (chem:make-restraint-dihedral atom1 atom2 atom3 atom4
                                                        (dihedral-min-degrees restraint)
                                                        (dihedral-max-degrees restraint)
                                                        (weight restraint)))))
            do (chem:add-restraint residue rr)))
    residue))

(defun build-one-molecule-for-topology (topology)
  "Debugging function - return a molecule that contains a single residue built from topology."
  (let ((residue (build-residue-for-monomer-name topology (name topology)))
        (mol (chem:make-molecule (name topology))))
    (chem:add-matter mol residue)
    mol))

(defun sketch-svg (topology &rest args)
  "Generate an svg sketch of the topology - send args to sketch2d"
  (let ((mol (build-one-molecule-for-topology topology)))
    (funcall (find-symbol "SVG" :sketch2d) (apply (find-symbol "SKETCH2D" :sketch2d) mol args))))


(defun root-monomer (oligomer)
  (loop for index below (length (monomers oligomer))
        for monomer = (elt (monomers oligomer) index)
        unless (has-in-coupling-p monomer)
          do (return-from root-monomer monomer)))




(defun canonical-sequence-monomer (oligomer coupling root monomer-out-couplings unique-ring-couplings)
  (let* ((outs (gethash root monomer-out-couplings))
         (result (list* (current-stereoisomer-name root oligomer)
                        (when outs
                          (let ((sorted-outs (sort outs #'string< :key (lambda (coup) (string (name coup))))))
                            (loop for sorted-out in sorted-outs
                                  for target-monomer = (target-monomer sorted-out)
                                  collect (canonical-sequence-monomer
                                           oligomer
                                           sorted-out
                                           target-monomer
                                           monomer-out-couplings
                                           unique-ring-couplings)))))))
    (if coupling
        (cons (name coupling) result)
        result)))



(defun canonical-sequence (oligomer)
  "Return a canonical sequence for the oligomer that can be compared to other oligomers to
   determine if they are equivalent"
  (let ((root-monomer (root-monomer oligomer)))
    (let ((monomer-out-couplings (make-hash-table))
          (unique-ring-couplings nil))
      (loop for index below (length (couplings oligomer))
            for coupling = (elt (couplings oligomer) index)
            if (typep coupling 'directional-coupling)
              do (push coupling (gethash (source-monomer coupling) monomer-out-couplings))
            else
              do (pushnew coupling unique-ring-couplings))
      (canonical-sequence-monomer oligomer nil root-monomer monomer-out-couplings unique-ring-couplings))))


(defun ordered-sequence-monomers (oligomer coupling root monomer-out-couplings unique-ring-couplings)
  (let* ((outs (gethash root monomer-out-couplings))
         (result (list* root
                        (when outs
                          (let ((sorted-outs (sort outs #'string< :key (lambda (coup) (string (name coup))))))
                            (loop for sorted-out in sorted-outs
                                  for target-monomer = (target-monomer sorted-out)
                                  append (ordered-sequence-monomers
                                           oligomer
                                           sorted-out
                                           target-monomer
                                           monomer-out-couplings
                                           unique-ring-couplings)))))))
    result))

(defun ordered-monomers (oligomer)
  "Return a canonical sequence of monomers"
  (let ((root-monomer (root-monomer oligomer)))
    (let ((monomer-out-couplings (make-hash-table))
          (unique-ring-couplings nil))
      (loop for index below (length (couplings oligomer))
            for coupling = (elt (couplings oligomer) index)
            if (typep coupling 'directional-coupling)
              do (push coupling (gethash (source-monomer coupling) monomer-out-couplings))
            else
              do (pushnew coupling unique-ring-couplings))
      (ordered-sequence-monomers oligomer nil root-monomer monomer-out-couplings unique-ring-couplings))))



(defun recursively-build-molecule (oligomer
                                   prev-monomer
                                   prev-topology
                                   prev-residue
                                   monomer-out-couplings
                                   molecule
                                   monomers-to-residues
                                   monomers-to-topologys
                                   molecule-index
                                   monomer-positions)
  (loop for out-coupling in (gethash prev-monomer monomer-out-couplings)
        for next-monomer = (target-monomer out-coupling)
        for next-topology = (let ((next-top (monomer-topology next-monomer oligomer)))
                              (unless (typep next-top 'topology:topology)
                                (error "Unexpected object - expected a topology - got ~s" next-top))
                              next-top)
        for next-stereoisomer-name = (current-stereoisomer-name next-monomer oligomer)
        for next-residue = (progn
                             (build-residue-for-monomer-name next-topology next-stereoisomer-name))
        do (setf (gethash next-monomer monomers-to-residues) next-residue
                 (gethash next-monomer monomers-to-topologys) next-topology)
           (let* ((monomer-position (gethash next-monomer monomer-positions))
                  (next-residue-index (chem:content-size molecule)))
             (declare (ignore next-residue-index))
             (chem:put-matter molecule (monomer-position-residue-index monomer-position) next-residue))
           (connect-residues prev-topology 
                                  prev-residue
                                  (source-plug-name out-coupling)
                                  next-topology
                                  next-residue
                                  (target-plug-name out-coupling))
           (recursively-build-molecule oligomer
                                       next-monomer
                                       next-topology
                                       next-residue
                                       monomer-out-couplings
                                       molecule
                                       monomers-to-residues
                                       monomers-to-topologys
                                       molecule-index
                                       monomer-positions)))

(defgeneric oligomer-force-field-name (foldamer)
  (:documentation "Return the name of the force field used by the foldamer"))

(defun build-molecule (oligomer aggregate molecule-index monomer-positions)
  (let ((root-monomer (root-monomer oligomer))
        (monomer-out-couplings (make-hash-table))
        (monomers-to-residues (make-hash-table))
        (monomers-to-topologys (make-hash-table))
        (ring-couplings nil))
    (loop for index below (length (couplings oligomer))
          for coupling = (elt (couplings oligomer) index)
          if (typep coupling 'directional-coupling)
            do (push coupling (gethash (source-monomer coupling) monomer-out-couplings))
          else
            do (pushnew coupling ring-couplings)) ; Only add ring coupling when unique
    (let* ((molecule (let ((mol (chem:make-molecule :mol)))
                       (chem:setf-force-field-name mol (oligomer-force-field-name (foldamer (oligomer-space oligomer))))
                       mol))
           (root-topology (let ((top (monomer-topology root-monomer oligomer)))
                            (unless (typep top 'topology:topology)
                              (error "Unexpected object - expected a topology - got ~s" top))
                            top))
           (stereoisomer-name (current-stereoisomer-name root-monomer oligomer))
           (root-residue (progn
                           (build-residue-for-monomer-name root-topology stereoisomer-name))))
      (setf (gethash root-monomer monomers-to-residues) root-residue
            (gethash root-monomer monomers-to-topologys) root-topology)
      (let ((number-of-residues (length (monomers oligomer))))
        (loop for count from 0 below (length (monomers oligomer))
              for monomer = (elt (monomers oligomer) count)
              do (setf (gethash monomer monomer-positions) (make-monomer-position
                                                            :molecule-index molecule-index
                                                            :residue-index count)))
        (let* ((monomer-position (gethash root-monomer monomer-positions))
               (next-residue-index (chem:content-size molecule)))
          (declare (ignore next-residue-index))
          (chem:resize-contents molecule number-of-residues)
          (chem:put-matter molecule (monomer-position-residue-index monomer-position) root-residue))
        (recursively-build-molecule oligomer
                                    root-monomer
                                    root-topology
                                    root-residue
                                    monomer-out-couplings
                                    molecule
                                    monomers-to-residues
                                    monomers-to-topologys
                                    molecule-index
                                    monomer-positions))
      ;; Now close the rings
      (loop for ring-coupling in ring-couplings
            for monomer1 = (monomer1 ring-coupling)
            for topology1 = (monomer-topology monomer1 oligomer)
            for residue1 = (gethash monomer1 monomers-to-residues)
            for monomer2 = (monomer2 ring-coupling)
            for topology2 = (monomer-topology monomer2 oligomer)
            for residue2 = (gethash monomer2 monomers-to-residues)
            for plug1name = (plug1 ring-coupling)
            for plug2name = (plug2 ring-coupling)
            do
            #+(or)(progn
                    (format *debug-io* "ring-coupling: ~a~%" ring-coupling)
                    (format *debug-io* "    topology1: ~a~%" topology1)
                    (format *debug-io* "     residue1: ~a~%" residue1)
                    (format *debug-io* "    plug1name: ~a~%" plug1name)
                    (format *debug-io* "    topology2: ~a~%" topology2)
                    (format *debug-io* "     residue2: ~a~%" residue2)
                    (format *debug-io* "    plug2name: ~a~%" plug2name))
                  (connect-residues topology1
                                    residue1
                                    plug1name
                                    topology2
                                    residue2
                                    plug2name))
      (chem:add-matter aggregate molecule)
      molecule)))


(defun build-all-molecules (oligomer-space &optional (number-of-sequences (number-of-sequences oligomer-space)))
  "Return a list of all built molecules for the oligomer-space.
Also return a hash-table that maps monomers to molecule/residue indices.
The number of entries can be limited by passing a lower value for the optional argument **number-of-sequences**
than the (chem:oligomer/number-of-sequences oligomer)."
  (let* ((monomer-positions (make-hash-table))
         (oligomer (topology:make-oligomer oligomer-space 0))
         (aggregates (loop for index from 0 below number-of-sequences
                           for aggregate = (chem:make-aggregate (intern (format nil "seq-~a" index) :keyword))
                           for molecule = (progn
                                            (topology:goto-sequence oligomer index)
                                            (build-molecule oligomer aggregate index monomer-positions))
                           collect aggregate)))
    (values aggregates monomer-positions)))
