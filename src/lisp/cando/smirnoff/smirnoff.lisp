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
  (let ((children (chem:chem-info-node-children node)))
    (loop for child in children
          do (walk-nodes child))))

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
  (let ((sigma-unit (sigma-unit source)))
    (cond
      ((string-equal sigma-unit "angstroms")
       (* 2.0 0.1))
      (t (error "Unrecognized sigma-unit ~s" sigma-unit)))))

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
          for nonbonded-force = (nonbonded-force part)
          for terms = (terms nonbonded-force)
          do (loop for index from (1- (length terms)) downto 0
                   for term = (aref terms index)
                   when (vector-push-extend term all-rules)))
    (make-instance 'smirnoff-type-rules :rules all-rules)))


(defmethod chem:assign-force-field-types ((atom-type-rules smirnoff-type-rules) molecule)
  "The first rule that matches is used to assign the types.
The chem:force-field-type-rules-merged generic function was used to organize the rules."
  (cando:do-atoms (atom molecule)
    (let ((type (loop named assign-type
                      for index from 0 below (length atom-type-rules)
                      for term = (aref atom-type-rules index)
                      for compiled-smirks = (compiled-smirks term)
                      for type = (type term)
                      for match = (chem:matches compiled-smirks atom)
                      when match
                        do (return-from assign-type type))))
      (if type
          (chem:set-type atom type)
          (error "Could not set type of atom ~s in force-field ~s" atom :smirnoff)))))
      
   
(defmethod chem:force-field-component-merge ((dest chem:ffnonbond-db) (source nonbonded-force))
  (format t "In chem:nonbond-force-field-component-merge with nonbonded-force: ~s~%" source)
  (let ((rmin-half-to-nanometers-multiplier (rmin-half-to-nanometers-multiplier source))
        (epsilon-to-kj-multiplier (epsilon-to-kj-multiplier source)))
    (loop with terms = (terms source)
          for index below (length terms)
          for term = (aref terms index)
          for rmin-half = (rmin-half term)
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


(defmethod chem:assign-molecular-force-field-parameters (energy-function (force-field combined-force-field) molecule)
  (error "Finish implementing chem:assign-molecular-force-field-parameters for smirnoff"))
