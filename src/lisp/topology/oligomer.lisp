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
      ;; Add properties to the residue
      (let* ((residue-properties (residue-properties constitution))
             (residue-net-charge (getf residue-properties :residue-net-charge)))
        (chem:matter/set-properties residue (append residue-properties (chem:matter/properties residue)))
        (when (numberp residue-net-charge)
          (chem:residue/set-net-charge residue residue-net-charge)))
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


(defun root-monomer (oligomer)
  (loop for index below (length (monomers oligomer))
        for monomer = (elt (monomers oligomer) index)
        unless (has-in-coupling-p monomer)
          do (return-from root-monomer monomer)))

(defun canonical-sequence-monomer (oligomer coupling root monomer-out-couplings unique-ring-couplings monomer-labels)
;;;  (break "Check unique-ring-couplings and monomer-labels")
  (let* ((outs (gethash root monomer-out-couplings))
         (monomer-name (current-stereoisomer-name root oligomer))
         (full-monomer-name (let ((monomer-label (gethash root monomer-labels)))
                              (if monomer-label
                                  (list monomer-name :label monomer-label)
                                  monomer-name)))
         (result (let* ((normal-couplings (list* full-monomer-name
                                                (when outs
                                                  #+(or)
                                                  (let ((sorted-outs (sort outs #'string< :key (lambda (coup) (string (name coup))))))
                                                    (loop for sorted-out in sorted-outs
                                                          for target-monomer = (target-monomer sorted-out)
                                                          collect (canonical-sequence-monomer
                                                                   oligomer
                                                                   sorted-out
                                                                   target-monomer
                                                                   monomer-out-couplings
                                                                   unique-ring-couplings
                                                                   monomer-labels)))
                                                  (let ((parts (loop for sorted-out in outs
                                                                     for target-monomer = (target-monomer sorted-out)
                                                                     collect (canonical-sequence-monomer
                                                                              oligomer
                                                                              sorted-out
                                                                              target-monomer
                                                                              monomer-out-couplings
                                                                              unique-ring-couplings
                                                                              monomer-labels))))
                                                    (sort parts #'< :key #'length)))))
                        )
                   normal-couplings))
         (lresult (cond
                    ((let* ((ring-coupling (gethash root unique-ring-couplings)))
                       (when ring-coupling
                         (let* ((other-monomer-label (let ((oml (gethash (monomer2 ring-coupling) monomer-labels)))
                                                       (unless oml
                                                         (error "Could not find label for ring other side ~s" ring-coupling))
                                                       oml))
                                (ring-info (list :ring
                                                 (plug1 ring-coupling)
                                                 (plug2 ring-coupling)
                                                 other-monomer-label)))
                           (append result (list (list ring-info)))))))
                    ((and (consp result) (consp (car (last result))))
                     (append (subseq result 0 (1- (length result))) (car (last result))))
                    (t (progn
                        result)))))
    (if coupling
        (cons (name coupling) lresult)
        lresult)))

(defun canonical-sequence (oligomer)
  "Return a canonical sequence for the oligomer that can be compared to other oligomers to
   determine if they are equivalent"
  (let ((root-monomer (root-monomer oligomer)))
    (let ((monomer-out-couplings (make-hash-table))
          (unique-ring-couplings (make-hash-table))
          (monomer-labels (make-hash-table)))
      (maphash (lambda (label monomer)
                 (setf (gethash monomer monomer-labels) label))
               (labeled-monomers (oligomer-space oligomer)))
      (loop for index below (length (couplings oligomer))
            for coupling = (elt (couplings oligomer) index)
            if (typep coupling 'directional-coupling)
              do (push coupling (gethash (source-monomer coupling) monomer-out-couplings))
            else
              do (setf (gethash (monomer1 coupling) unique-ring-couplings) coupling))
      (canonical-sequence-monomer oligomer nil root-monomer monomer-out-couplings unique-ring-couplings monomer-labels))))


(defun ordered-sequence-monomers (oligomer coupling root monomer-out-couplings unique-ring-couplings)
  (declare (ignore coupling))
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


(defun build-residue (oligomer
                      monomers-to-residues
                      monomers-to-topologys
                      monomer-positions-accumulator
                      molecule-index
                      molecule
                      prev-topology
                      prev-residue
                      next-monomer
                      &optional out-coupling)
  (let* ((topology (let ((next-top (monomer-topology next-monomer oligomer)))
                     (unless (typep next-top 'topology:topology)
                       (error "Unexpected object - expected a topology - got ~s" next-top))
                     next-top))
         (stereoisomer-name (current-stereoisomer-name next-monomer oligomer))
         (residue (build-residue-for-monomer-name topology stereoisomer-name)))
    (setf (gethash next-monomer monomers-to-residues) residue
          (gethash next-monomer monomers-to-topologys) topology)
    (let ((residue-index (chem:content-size molecule)))
      (setf (gethash next-monomer monomer-positions-accumulator)
            (make-instance 'monomer-position
                           :molecule-index molecule-index
                           :residue-index residue-index)))
    (chem:add-matter molecule residue)
    (when prev-residue
      (connect-residues prev-topology
                        prev-residue
                        (source-plug-name out-coupling)
                        topology
                        residue
                        (target-plug-name out-coupling)))
    (values residue topology)))

(defun recursively-build-molecule (oligomer
                                   monomer-subset
                                   prev-monomer
                                   prev-topology
                                   prev-residue
                                   monomer-out-couplings
                                   molecule
                                   monomers-to-residues
                                   monomers-to-topologys
                                   molecule-index
                                   monomer-positions-accumulator)
  (if prev-monomer
      (let ((next-topology nil)
            (next-residue nil))
        (loop for out-coupling in (gethash prev-monomer monomer-out-couplings)
              for next-monomer = (target-monomer out-coupling)
              do (when (in-monomer-subset monomer-subset next-monomer)
                   (multiple-value-setq (next-residue next-topology)
                     (build-residue oligomer
                                    monomers-to-residues
                                    monomers-to-topologys
                                    monomer-positions-accumulator
                                    molecule-index
                                    molecule
                                    prev-topology
                                    prev-residue
                                    next-monomer
                                    out-coupling)))
              do (recursively-build-molecule oligomer
                                             monomer-subset
                                             next-monomer
                                             next-topology
                                             next-residue
                                             monomer-out-couplings
                                             molecule
                                             monomers-to-residues
                                             monomers-to-topologys
                                             molecule-index
                                             monomer-positions-accumulator)))
      (let* ((root-monomer (root-monomer oligomer))
             (monomer-out-couplings (make-hash-table))
             (monomers-to-topologys (make-hash-table))
             (ring-couplings nil)
             (root-topology nil)
             (root-residue nil))
        (loop for index below (length (couplings oligomer))
              for coupling = (elt (couplings oligomer) index)
              if (typep coupling 'directional-coupling)
                do (push coupling (gethash (source-monomer coupling) monomer-out-couplings))
              else
                do (pushnew coupling ring-couplings)) ; Only add ring coupling when unique
        (when (in-monomer-subset monomer-subset root-monomer)
          (multiple-value-setq (root-residue root-topology)
            (build-residue oligomer
                           monomers-to-residues
                           monomers-to-topologys
                           monomer-positions-accumulator
                           molecule-index
                           molecule
                           prev-topology
                           prev-residue
                           root-monomer)))
        (recursively-build-molecule oligomer
                                    monomer-subset
                                    root-monomer
                                    root-topology
                                    root-residue
                                    monomer-out-couplings
                                    molecule
                                    monomers-to-residues
                                    monomers-to-topologys
                                    molecule-index
                                    monomer-positions-accumulator))))

(defgeneric oligomer-force-field-name (foldamer)
  (:documentation "Return the name of the force field used by the foldamer"))

(defun build-molecule (oligomer
                       &key (aggregate (chem:make-aggregate :dummy))
                         monomer-subset
                         (molecule-index 0)
                         (monomers-to-residues (make-hash-table))
                         (monomer-positions-accumulator (make-hash-table)))
  "Build the molecule for the oligomer using the monomer-subset if provided.
If the caller passes a hash-table for monomers-to-residues, then fill that hash-table with a map
of monomers-to-residues."
  (let ((root-monomer (when (in-monomer-subset monomer-subset (root-monomer oligomer))))
        (monomer-out-couplings (make-hash-table))
        (monomers-to-topologys (make-hash-table))
        (ring-couplings nil)
        (name (name (oligomer-space oligomer))))
    (declare (ignore root-monomer))
    (loop for index below (length (couplings oligomer))
          for coupling = (elt (couplings oligomer) index)
          if (typep coupling 'directional-coupling)
            do (push coupling (gethash (source-monomer coupling) monomer-out-couplings))
          else
            do (pushnew coupling ring-couplings)) ; Only add ring coupling when unique
    (let* ((molecule (let ((mol (chem:make-molecule name)))
                       (chem:setf-force-field-name mol (oligomer-force-field-name (foldamer (oligomer-space oligomer))))
                       mol)))
      (recursively-build-molecule oligomer
                                  monomer-subset
                                  nil
                                  nil
                                  nil
                                  monomer-out-couplings
                                  molecule
                                  monomers-to-residues
                                  monomers-to-topologys
                                  molecule-index
                                  monomer-positions-accumulator)
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
      (maphash (lambda (name monomer)
                 (unless (atom monomer)
                   (error "The monomer ~s cannot be a list" monomer))
                 (when (in-monomer-subset monomer-subset monomer)
                   (let* ((monpos (let ((mp (gethash monomer monomer-positions-accumulator)))
                                    (unless mp (error "Could not find the monomer ~s in the monomer-positions-accumulator ~s"
                                                      monomer
                                                      monomer-positions-accumulator))
                                    mp))
                          (res (chem:content-at molecule (topology:residue-index monpos))))
                     (chem:set-property res :label name))))
               (labeled-monomers (oligomer-space oligomer)))
      (values molecule monomer-positions-accumulator))))


(defun build-all-molecules (oligomer-space
                            &key
                              (monomers-to-residues (make-hash-table) monomers-to-residues-p)
                              (number-of-sequences (number-of-sequences oligomer-space)))
  "Build all of the molecules in the oligomer-space into a single aggregate and return it.
Also return the position of each monomer in the aggregate -
A hash-table that maps monomers to molecule/residue indices.
The number of entries can be limited by passing a lower value for the optional argument **number-of-sequences**
than the (chem:oligomer/number-of-sequences oligomer)."
  (declare (ignore monomers-to-residues-p))
  (let* ((monomer-positions-accumulator (make-hash-table))
         (oligomer (topology:make-oligomer oligomer-space 0))
         (aggregates (loop for index from 0 below number-of-sequences
                           for aggregate = (chem:make-aggregate (intern (format nil "seq-~a" index) :keyword))
                           for molecule = (progn
                                            (topology:goto-sequence oligomer index)
                                            (build-molecule oligomer
                                                            :aggregate aggregate
                                                            :molecule-index index
                                                            :monomers-to-residues monomers-to-residues
                                                            :monomer-positions-accumulator monomer-positions-accumulator))
                           collect aggregate)))
    (values aggregates monomer-positions-accumulator)))
