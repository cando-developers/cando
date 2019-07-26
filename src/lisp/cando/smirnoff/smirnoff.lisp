(in-package :smirnoff)

#+(or)
(defmacro smlog (fmt &rest args))

;;#+(or)
(defmacro smlog (fmt &rest args)
  `(format *debug-io* ,fmt ,@args))

(defstruct element-info element atomic-number atomic-mass)
(defvar *element-info*)

(defgeneric walk-nodes (node))
(defmethod walk-nodes (node)
  (when node
    (let ((children (chem:chem-info-node-children node)))
      (loop for child in children
            do (walk-nodes child)))))

(defmethod walk-nodes ((node chem:atom-test))
  (let ((test (chem:atom-test-type node)))
    ;; ignore everything that isn't element, atomic-number or atomic-mass
    (case test
      (:sapelement (setf (element-info-element *element-info*) (chem:get-symbol-arg node)))
      (:sapatomic-number (setf (element-info-atomic-number *element-info*) (chem:get-int-arg node)))
      (:sapatomic-mass (setf (element-info-atomic-mass *element-info*) (chem:get-int-arg node))))))


(defun smirks-head-atomic-mass (node)
  "Figure out the atomic mass for the first atom using a chem-info object"
  (let ((*element-info* (make-element-info))
        (head (if (typep node 'chem:atom-test)
                  node
                  (first (chem:chem-info-node-children node)))))
    (walk-nodes head)
    (cond
      ((element-info-element *element-info*)
       (chem:atomic-weight-for-element (element-info-element *element-info*)))
      ((element-info-atomic-number *element-info*)
       (let ((element (chem:element-for-atomic-number (element-info-atomic-number *element-info*))))
         (chem:atomic-weight-for-element element)))
      ((element-info-atomic-mass *element-info*)
       (warn "Don't rely on atomic-weight until we add integral atomic weight to elements.cc")
       (element-info-atomic-mass *element-info*))
      (t (error "There must be an element, atomic number or atomic weight in the head of the smirks string")))))

(defun rmin-half-to-nanometers-multiplier (source)
  (let ((rmin-half-unit (rmin-half-unit source)))
    (cond
      ((string-equal rmin-half-unit "angstroms")
       (* 2.0 0.1))
      (t (error "Unrecognized rmin-half-unit ~s" rmin-half-unit)))))

(defun epsilon-to-kj-multiplier (source)
  (let ((epsilon-unit (epsilon-unit source)))
    (cond
      ((string-equal epsilon-unit "kilocalories_per_mole")
       4.184)
      (t (error "Unrecognized epsilon-unit ~s" epsilon-unit)))))


(defclass smirnoff-type-rules ()
  ((rules :initarg :rules :accessor rules)))


#+(or)
(defmethod chem:force-field-type-rules-merged ((first-force-field-part smirnoff-force-field) all-force-field-parts)
  "Accumulate the nonbond terms from first to last and bottom to top.
The first rule that matches is used to assign the type."
  (let ((all-rules (make-array 128 :fill-pointer 0 :adjustable t)))
    (loop for part in force-field-parts
          for vdw-force = (vdw-force part)
          for terms = (terms vdw-force)
          do (loop for index from (1- (length terms)) downto 0
                   for term = (aref terms index)
                   when (vector-push-extend term all-rules)))
    (make-instance 'smirnoff-type-rules :rules all-rules)))


(defmethod chem:assign-force-field-types ((combined-force-field combined-force-field) molecule)
  "The first rule that matches is used to assign the types.
The chem:force-field-type-rules-merged generic function was used to organize the rules."
  (cando:do-atoms (atom molecule)
    (let ((type (loop named assign-type
                      for field in (chem:force-fields-as-list combined-force-field)
                      for vdw-force = (vdw-force field)
                      for terms = (terms vdw-force)
                      do (loop for index from (1- (length terms)) downto 0
                               for term = (aref terms index)
                               for compiled-smirks = (compiled-smirks term)
                               for type = (type term)
                               for match = (chem:matches compiled-smirks atom)
                               when match
                                 do (return-from assign-type type)))))
      (if type
          (chem:set-type atom type)
          (error "Could not set type of atom ~s in force-field ~s" atom :smirnoff)))))
      
   
(defmethod chem:force-field-component-merge ((dest chem:ffnonbond-db) (source vdw-force))
  (let ((rmin-half-to-nanometers-multiplier (rmin-half-to-nanometers-multiplier source))
        (epsilon-to-kj-multiplier (epsilon-to-kj-multiplier source)))
    (loop with terms = (terms source)
          for index below (length terms)
          for term = (aref terms index)
          for rmin-half = (/ (rmin-half term) 2)
          for radius-nanometers = (* rmin-half rmin-half-to-nanometers-multiplier)
          for epsilon = (epsilon term)
          for epsilon-kj = (* epsilon epsilon-to-kj-multiplier)
          for type = (type term)
          when (not (chem:has-type dest type))
            do (let* ((atomic-mass (smirks-head-atomic-mass (compiled-smirks term)))
                      (ffnonbond (chem:make-ffnonbond type
                                                      :mass atomic-mass
                                                      :radius-nanometers radius-nanometers
                                                      :epsilon-kj epsilon-kj)))
                 (chem:ffnonbond-db-add dest ffnonbond)))
    dest))

(defun bonds-hash-table (molecule)
  (let ((bonds (make-hash-table :test #'equal)))
    (chem:map-bonds
     nil
     (lambda (a1 a2 bond-order)
       (let ((key (list a1 a2)))
         (setf (gethash key bonds) nil)))
     molecule)
    bonds))
  

(defmethod chem:generate-molecule-energy-function-tables (energy-function molecule (combined-force-field combined-force-field) active-atoms)
  (let* ((molecule-graph (chem:make-molecule-graph-from-molecule molecule))
         (atom-table (chem:atom-table energy-function))
         (bonds (bonds-hash-table molecule))
         (angles (make-hash-table :test #'equal))
         (ptors (make-hash-table :test #'equal))
         (itors (make-hash-table :test #'equal))
         (*print-pretty* nil))
    (when (chem:verbose 2) (format t "Generating stretch terms~%"))
    (let ((bond-energy (chem:get-stretch-component energy-function))
          )
      ;; Work through the force-fields backwards so that more recent terms will shadow older ones
      (loop for force-field in (reverse (chem:force-fields-as-list combined-force-field))
            for bonds-force = (bonds-force force-field)
            do (loop for term across (terms bonds-force)
                     for compiled-smirks = (compiled-smirks term)
                     for smirks-graph = (chem:make-chem-info-graph compiled-smirks)
                     for hits = (chem:boost-graph-vf2 smirks-graph molecule-graph)
                     do (loop for hit in hits
                              for a1 = (aref hit 1)
                              for a2 = (aref hit 2)
                              for key12 = (list a1 a2)
                                 ;; The vf2 algorithm will find the bond in both directions
                                 ;; one of them will match - so we won't worry about the other one
                              do (multiple-value-bind (val found)
                                     (gethash key12 bonds)
                                   (if found
                                       (progn
                                         (push term (gethash key12 bonds)))
                                       (let ((key21 (reverse key12)))
                                         (multiple-value-bind (val found)
                                             (gethash key21 bonds)
                                           (if found
                                               (progn
                                                 (push term (gethash key21 bonds)))
                                               (error "Could not find keys ~s or ~s in bonds" key12 key21)))))))))
      (maphash (lambda (key terms)
                 (unless terms
                   (error "Could not find stretch parameter for ~s" key))
                 (let* ((a1 (first key))
                        (a2 (second key))
                        (term (first terms))
                        (k (k term))
                        (k-amber (/ k 2.0)) ; Amber drops the factor of 2 in Hooke's law
                        (len (len term)))
                   (chem:add-stretch-term bond-energy atom-table a1 a2 k-amber len)))
               bonds))
    (when (chem:verbose 2) (format t "Generating angle terms~%"))
    (let ((angle-energy (chem:get-angle-component energy-function)))
      (chem:map-angles nil (lambda (a1 a2 a3)
                             (setf (gethash (list a1 a2 a3) angles) nil))
                       molecule)
      ;; Work through the force-fields backwards so that more recent terms will shadow older ones
      (loop for force-field in (reverse (chem:force-fields-as-list combined-force-field))
            for angles-force = (angles-force force-field)
            do (loop for term across (terms angles-force)
                     for smirks = (progn
                                    (smirks term))
                     for compiled-smirks = (compiled-smirks term)
                     for smirks-graph = (chem:make-chem-info-graph compiled-smirks)
                     for hits = (progn
                                  (chem:boost-graph-vf2 smirks-graph molecule-graph))
                     do (loop for hit in hits
                              for a1 = (aref hit 1)
                              for a2 = (aref hit 2)
                              for a3 = (aref hit 3)
                              for key123 = (list a1 a2 a3)
                                 ;; The vf2 algorithm will find the angle in both directions
                                 ;; one of them will match - so we won't worry about the other one
                              do (multiple-value-bind (val found)
                                     (gethash key123 angles)
                                   (if found
                                       (progn
                                         (push term (gethash key123 angles)))
                                       (let ((key321 (nreverse key123)))
                                         (multiple-value-bind (val found)
                                             (gethash key321 angles)
                                           (if found
                                               (progn
                                                 (push term (gethash key321 angles)))
                                               (error "Could not find keys ~s or ~s in angles" key123 key321)))))))))
      (maphash (lambda (key terms)
                 (unless terms
                   (error "Could not find angle parameter for ~s" key))
                 (let* ((a1 (first key))
                        (a2 (second key))
                        (a3 (third key))
                        (term (first terms))
                        (k (k term))
                        (k-amber (/ k 2.0)) ; Amber drops the factor of 2 in Hooke's law
                        (angle-rad (angle-rad term)))
                   (chem:add-angle-term angle-energy atom-table a1 a2 a3 k-amber angle-rad)))
               angles))


    
    (when (chem:verbose 2) (format t "Generating dihedral terms~%"))
    (let ((dihedral-energy (chem:get-dihedral-component energy-function)))
      (chem:map-dihedrals nil (lambda (a1 a2 a3 a4)
                                (setf (gethash (list a1 a2 a3 a4) ptors) nil))
                          molecule)
      ;; Work through the force-fields backwards so that more recent terms will shadow older ones
      (loop for force-field in (reverse (chem:force-fields-as-list combined-force-field))
            for proper-torsion-force = (proper-torsion-force force-field)
            do (loop for term across (terms proper-torsion-force)
                     for smirks = (progn
                                    (smirks term))
                     for compiled-smirks = (compiled-smirks term)
                     for smirks-graph = (chem:make-chem-info-graph compiled-smirks)
                     for hits = (chem:boost-graph-vf2 smirks-graph molecule-graph)
                     do (loop for hit in hits
                              for a1 = (aref hit 1)
                              for a2 = (aref hit 2)
                              for a3 = (aref hit 3)
                              for a4 = (aref hit 4)
                              for key1234 = (list a1 a2 a3 a4)
                              ;; The vf2 algorithm will find the angle in both directions
                              ;; one of them will match - so we won't worry about the other one
                              do (multiple-value-bind (val found)
                                     (gethash key1234 ptors)
                                   (if found
                                       (progn
                                         (push term (gethash key1234 ptors)))
                                       (let ((key4321 (nreverse key1234)))
                                         (multiple-value-bind (val found)
                                             (gethash key4321 ptors)
                                           (if found
                                               (progn
                                                 (push term (gethash key4321 ptors)))
                                               (error "Could not find keys ~s or ~s in ptors" key1234 key4321)))))))))
      (maphash (lambda (key terms)
                 (let* ((a1 (first key))
                        (a2 (second key))
                        (a3 (third key))
                        (a4 (fourth key))
                        (term (first terms)))
                   (if term
                       (loop for part in (parts term)
                             for phase = (phase-rad part)
                             for periodicity = (periodicity part)
                             for k = (k part)
                             for idivf = (idivf part)
                             for v = (/ k idivf)
                             do (chem:add-dihedral-term dihedral-energy atom-table a1 a2 a3 a4 phase t v periodicity))
                       (warn "Could not find proper torsion parameters for ~s" key))))
               ptors)
      (loop for force-field in (reverse (chem:force-fields-as-list combined-force-field))
            for improper-torsion-force = (improper-torsion-force force-field)
            do (loop for term across (terms improper-torsion-force)
                     for smirks = (progn
                                    (smirks term))
                     for compiled-smirks = (compiled-smirks term)
                     for smirks-graph = (chem:make-chem-info-graph compiled-smirks)
                     for hits = (chem:boost-graph-vf2 smirks-graph molecule-graph)
                     do (loop for hit in hits
                              for a1 = (aref hit 1)
                              for a2 = (aref hit 2)
                              for a3 = (aref hit 3)
                              for a4 = (aref hit 4)
                              for key1234 = (list a1 a2 a3 a4)
                              ;; The vf2 algorithm will find the angle in both directions
                              ;; one of them will match - so we won't worry about the other one
                               do (push term (gethash key1234 itors)))))
      (maphash (lambda (key terms)
                 (let* ((a1 (first key))
                        (a2 (second key))
                        (a3 (third key))
                        (a4 (fourth key))
                        (term (first terms)))
                    (loop for part in (parts term)
                         for phase = (phase-rad part)
                         for periodicity = (periodicity part)
                         for k = (k part)
                         for idivf = (idivf part)
                         for v = (/ k idivf)
                         do (chem:add-dihedral-term dihedral-energy atom-table a1 a2 a3 a4 phase nil v periodicity))))
               itors))))
