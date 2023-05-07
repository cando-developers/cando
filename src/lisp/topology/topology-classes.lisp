(in-package :topology)

(defclass constitution-atom ()
  ((atom-name :initarg :atom-name :accessor atom-name)
   (index :initarg :index :accessor index)
   (element :initarg :element :accessor element)
   (properties :initform nil :initarg :properties :accessor properties)
   (bonds :initform nil :initarg :bonds :accessor bonds)))

(cando.serialize:make-class-save-load constitution-atom
                            :print-unreadably
                            (lambda (obj stream)
                              (print-unreadable-object (obj stream :type t)
                                (format stream "~a" (atom-name obj)))))

(defclass constitution-bond ()
  ((to-atom-index :initarg :to-atom-index :accessor to-atom-index)
   (order :initarg :order :accessor order)))


(cando.serialize:make-class-save-load constitution-bond)

(defclass constitution-atoms ()
  ((atoms :initarg :atoms :accessor atoms)))

(cando.serialize:make-class-save-load constitution-atoms)

(defclass constitution ()
  ((name :initarg :name :accessor name)
   (constitution-atoms :initarg :constitution-atoms :accessor constitution-atoms)))

(cando.serialize:make-class-save-load
 constitution
  :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~s" (name obj)))))

(defun constitution-atom-named (constitution name)
  (loop for ca across (constitution-atoms constitution)
        when (eq name (atom-name ca))
          do (return ca)
        finally (error "Could not find constitution-atom with name ~a in ~a" name constitution)))

(defclass stereoisomer-atom ()
  ((atom-name :initarg :atom-name :accessor atom-name)
   (constitution-atom-index :initarg :constitution-atom-index :accessor constitution-atom-index)
   (atom-charge :initarg :atom-charge :accessor atom-charge)
   (atom-type :initarg :atom-type :accessor atom-type)))

(cando.serialize:make-class-save-load
 stereoisomer-atom
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~s ~a" (atom-name obj) (constitution-atom-index obj)))))

(defclass stereoisomer-virtual-atom (stereoisomer-atom)
  ())


(defclass stereoconfiguration ()
  ((atom-name :initarg :atom-name :accessor atom-name)
   (configuration :initarg :configuration :accessor configuration)
   (stereochemistry-type :initarg :stereochemistry-type :accessor stereochemistry-type)))

(cando.serialize:make-class-save-load stereoconfiguration
                            :print-unreadably
                            (lambda (obj stream)
                              (print-unreadable-object (obj stream :type t)
                                (format stream "~a ~a" (atom-name obj) (configuration obj)))))

(defclass stereoisomer ()
  ((name :initarg :name :accessor name)
   (pdb :initarg :pdb :accessor pdb)
   (stereoisomer-atoms :initarg :stereoisomer-atoms :accessor stereoisomer-atoms)
   (stereoconfigurations :initarg :stereoconfigurations :accessor stereoconfigurations)
   (complex-restraints :initform nil :initarg :complex-restraints :accessor complex-restraints)))

(cando.serialize:make-class-save-load
 stereoisomer
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~s" (name obj) (stereoconfigurations obj)))))

(defclass restraint ()
  ())

(defclass dihedral-restraint ()
  ((atom1-name :initarg :atom1-name :accessor atom1-name)
   (atom2-name :initarg :atom2-name :accessor atom2-name)
   (atom3-name :initarg :atom3-name :accessor atom3-name)
   (atom4-name :initarg :atom4-name :accessor atom4-name)
   (dihedral-min-degrees :initarg :dihedral-min-degrees
                         :accessor dihedral-min-degrees)
   (dihedral-max-degrees :initarg :dihedral-max-degrees
                         :accessor dihedral-max-degrees)
   (weight :initarg :weight :accessor weight)))

(cando.serialize:make-class-save-load dihedral-restraint
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a ~a ~a ~a ~a ~a"
             (atom1-name obj)
             (atom2-name obj)
             (atom3-name obj)
             (atom4-name obj)
             (dihedral-min-degrees obj)
             (dihedral-max-degrees obj)
             (weight obj)))))

(defclass topology ()
  ((name :initarg :name :accessor name)
   (constitution :initarg :constitution :accessor constitution)
   (stereoisomer :initarg :stereoisomer :accessor stereoisomer)
   (property-list :initform nil :initarg :property-list :accessor property-list)
   (plugs :initform (make-hash-table) :type hash-table :initarg :plugs :accessor plugs)
   (joint-template :initarg :joint-template :accessor joint-template)
   (restraints :initform nil :initarg :restraints :accessor restraints)
   ))


(defgeneric add-plug (topology plug-name plug))

(defmethod add-plug ((topology topology) plug-name plug)
  (setf (gethash plug-name (plugs topology)) plug))

(defgeneric topologyp (obj)
  (:documentation "Return T if the object is a topology"))

(defmethod topologyp ((obj t)) nil)
(defmethod topologyp ((obj topology)) t)

(cando.serialize:make-class-save-load topology
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~s" (name obj)))))

(defun has-plug-named (topology plug-name)
  (gethash plug-name (plugs topology)))

(defun plug-named (topology plug-name)
  (gethash plug-name (plugs topology)))

(defun find-in-plug (topology)
  (maphash (lambda (name plug)
             (declare (ignore name))
             (when (typep plug 'in-plug)
               (return-from find-in-plug plug)))
           (plugs topology)))

(defun out-plugs-as-list (topology)
  (let (out-plugs)
    (maphash (lambda (name plug)
               (declare (ignore name))
               (when (typep plug 'out-plug)
                 (push plug out-plugs)))
             (plugs topology))
    out-plugs))

(defun all-out-plug-names-that-match-in-plug-name (topology in-plug-name)
  (let (out-plugs
        (in-coupling-name (coupling-name in-plug-name)))
    (maphash (lambda (name plug)
               (declare (ignore name))
               (when (and (typep plug 'out-plug)
                          (eq (coupling-name (name plug))
                              in-coupling-name))
                 (push (name plug) out-plugs)))
             (plugs topology))
    out-plugs))

(defun plugs-as-list (topology)
  (let (plugs)
    (maphash (lambda (name plug)
               (declare (ignore name))
               (push plug plugs))
             (plugs topology))
    plugs))

(defclass plug-bond ()
  ((atom-name :initarg :atom-name :accessor atom-name)
   (bond-order :initarg :bond-order :accessor bond-order)))

(cando.serialize:make-class-save-load
 plug-bond
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a" (atom-name obj)))))

(defclass plug ()
  ((name :initarg :name :accessor name)
   (plug-bonds :initform (make-array 2) :initarg :plug-bonds :accessor plug-bonds)))

(cando.serialize:make-class-save-load plug
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a" (name obj) (plug-bonds obj)))))

(defclass in-plug (plug)
  ())

(cando.serialize:make-class-save-load
 in-plug
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a" (name obj) (plug-bonds obj)))))

(defun make-in-plug (name bond0 bond-order0 &optional bond1 bond-order1)
  (make-instance 'in-plug
                 :name name
                 :plug-bonds (if bond1
                                 (vector (make-instance 'plug-bond
                                                        :atom-name bond0
                                                        :bond-order bond-order0)
                                         (make-instance 'plug-bond
                                                        :atom-name bond1
                                                        :bond-order bond-order1))
                                 (vector (make-instance 'plug-bond
                                                        :atom-name bond0
                                                        :bond-order bond-order0)))))

(defclass out-plug (plug)
  ())

(cando.serialize:make-class-save-load
 out-plug
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~a ~a" (name obj) (plug-bonds obj)))))

(defun make-out-plug (name bond0 bond-order0 &optional bond1 bond-order1)
  (make-instance 'out-plug
                 :name name
                 :plug-bonds (if bond1
                                 (vector (make-instance 'plug-bond
                                                        :atom-name bond0
                                                        :bond-order bond-order0)
                                         (make-instance 'plug-bond
                                                        :atom-name bond1
                                                        :bond-order bond-order1))
                                 (vector (make-instance 'plug-bond
                                                        :atom-name bond0
                                                        :bond-order bond-order0)))))

(defun in-plug-name (coupling-name)
  (intern (format nil "~c~a" #\- (string coupling-name)) :keyword))

(defun out-plug-name (coupling-name)
  (intern (format nil "~c~a" #\+ (string coupling-name)) :keyword))

(defun is-in-plug-name (name)
  (let ((first-char (elt (string name) 0)))
    (char= first-char #\-)))

(defun other-plug-name (name)
  (if (is-in-plug-name name)
      (intern (format nil "~c~a" #\+ (subseq (string name) 1)) :keyword)
      (intern (format nil "~c~a" #\- (subseq (string name) 1)) :keyword)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Monomer, coupling, and oligomer
;;;


(defclass monomer ()
  ((id :initform nil
       :initarg :id
       :accessor id)
   (couplings :type hash-table
              :initform (make-hash-table)
              :initarg :couplings
              :accessor couplings)
   (monomers :initarg :monomers :accessor monomers)))

(cando.serialize:make-class-save-load monomer
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream ":id ~a :monomers ~s" (id obj) (monomers obj)))))

(defun number-of-stereoisomers (monomer)
  (length (monomers monomer)))

(defun current-stereoisomer-name (monomer oligomer)
  (let ((monomer-index (position monomer (monomers (oligomer-space oligomer)))))
    (unless monomer-index
      (error "Could not find monomer ~a in ~a" monomer oligomer))
    (elt (monomers monomer) (aref (monomer-indices oligomer) monomer-index))))

(defun monomer-topologys-all-have-plug-named (monomer plug-name)
  (let ((topologys (loop for mon-name in (monomers monomer)
                         collect (chem:find-topology mon-name t))))
    (cond
      ((every (lambda (top) (has-plug-named top plug-name)) topologys)
       t)
      ((every (lambda (top) (not (has-plug-named top plug-name))) topologys)
       nil)
      (t (error "Some topologys of ~s have the plug named ~a and some do not" (monomers monomer) plug-name)))))

(defun monomer-topology (monomer oligomer)
  (let ((monomer-position (position monomer (monomers (oligomer-space oligomer)))))
    (unless monomer-position
      (error "Could not find ~a in ~a" monomer oligomer))
    (let* ((monomer (aref (monomers (oligomer-space oligomer)) monomer-position))
           (monomer-index (aref (monomer-indices oligomer) monomer-position))
           (monomer-name (elt (monomers monomer) monomer-index)))
      (chem:find-topology monomer-name t))))

(defun has-in-coupling-p (monomer)
  (maphash (lambda (plug-name coupling)
             (declare (ignore coupling))
             (when (in-plug-name-p plug-name)
               (return-from has-in-coupling-p t)))
           (couplings monomer)))

(defun monomer-plug-named (monomer plug-name)
  (gethash plug-name (couplings monomer)))

(defgeneric other-monomer (coupling monomer))

(defun monomer-on-other-side (monomer plug-name)
  (let ((coupling (monomer-plug-named monomer plug-name)))
    (unless coupling
      (error "Coupling for plug ~a was not found" plug-name))
    (other-monomer coupling monomer)))


(defclass coupling ()
  ((name :initarg :name :accessor name)))

(defclass directional-coupling (coupling)
  ((source-plug-name :initarg :source-plug-name :accessor source-plug-name)
   (target-plug-name :initarg :target-plug-name :accessor target-plug-name)
   (source-monomer :initarg :source-monomer :accessor source-monomer)
   (target-monomer :initarg :target-monomer :accessor target-monomer)))

(cando.serialize:make-class-save-load directional-coupling
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a ~a ~a" (source-monomer obj) (source-plug-name obj) (target-plug-name obj) (target-monomer obj)))))

(defmethod other-monomer ((coupling directional-coupling) monomer)
  (cond
    ((eq (source-monomer coupling) monomer)
     (target-monomer coupling))
    ((eq (target-monomer coupling) monomer)
     (source-monomer coupling))
    (t (error "Could not find monomer ~a in coupling ~a" monomer coupling))))


(defclass ring-coupling (coupling)
  ((plug1 :initarg :plug1 :accessor plug1)
   (plug2 :initarg :plug2 :accessor plug2)
   (monomer1 :initarg :monomer1 :accessor monomer1)
   (monomer2 :initarg :monomer2 :accessor monomer2)))

(defmethod other-monomer ((coupling ring-coupling) monomer)
  (cond
    ((eq (monomer1 coupling) monomer)
     (monomer2 coupling))
    ((eq (monomer2 coupling) monomer)
     (monomer1 coupling))
    (t (error "Could not find monomer ~a in coupling ~a" monomer coupling))))


(defun has-ring-closing-coupling (monomer)
  (maphash (lambda (plug-name coupling)
             (declare (ignore plug-name))
             (when (typep coupling 'ring-coupling)
               (return-from has-ring-closing-coupling t)))
           (couplings monomer)))

(defclass oligomer-space ()
  ((foldamer :initarg :foldamer :accessor foldamer)
   (monomers :initform (make-array 16 :adjustable t :fill-pointer 0)
             :initarg :monomers :accessor monomers)
   (couplings :initform (make-array 16 :adjustable t :fill-pointer 0)
              :initarg :couplings :accessor couplings)
   (%number-of-sequences :initform nil
                        :initarg :number-of-sequences
                        :accessor %number-of-sequences)))

(cando.serialize:make-class-save-load oligomer-space
  :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t))))


(defun topologys-in-oligomer-space (oligomer-space)
  (loop for monomer across (monomers oligomer-space)
        for monomer-names = (monomers monomer)
        append (loop for monomer-name in monomer-names
                     for topology = (chem:find-topology monomer-name t)
                     collect topology)))

(defun number-of-monomers (oligomer-space)
  (length (monomers oligomer-space)))

(defun add-monomer (oligomer-space monomer)
  (vector-push-extend monomer (monomers oligomer-space)))

(defun calculate-number-of-sequences (oligomer-space)
  (let ((num 1))
  (loop for monomer across (monomers oligomer-space)
        do (setf num (* num (length (monomers monomer)))))
    num))

(defun number-of-sequences (oligomer-space)
  ;; Lazy calculate number-of-sequences
  (when (null (%number-of-sequences oligomer-space))
    (setf (%number-of-sequences oligomer-space) (calculate-number-of-sequences oligomer-space)))
  (%number-of-sequences oligomer-space))

(defun make-oligomer-space (foldamer tree &key (parts *parts*))
  "Make an oligomer-space from a description in the **tree**.
The tree is a nested list of lists that look like 
(component coupling component coupling component ... ).
It starts with the name of a topology or a design:part.
The component can be a name of a single topology, it can be a design:part,
or it can 

Examples:
(make-oligomer-space '(:ccap :default :ala :default :ala :default :ser :default :ncap))
(make-oligomer-space '((:lego3 :label :first)
                 :default :lego3
                 :default :lego3
                 :default (cycle :first)))
"
  (let* ((oligomer-space (make-instance 'oligomer-space
                                        :foldamer foldamer))
         (labels (make-hash-table)))
    (interpret-subtree oligomer-space tree labels :parts parts)
    (setf (%number-of-sequences oligomer-space)
          (calculate-number-of-sequences oligomer-space))
    oligomer-space))

(defclass oligomer ()
  ((monomer-indices :initarg :monomer-indices :accessor monomer-indices)
   (oligomer-space :initarg :oligomer-space :accessor oligomer-space)))

(cando.serialize:make-class-save-load oligomer)

(defun copy-oligomer (oligomer)
  (make-instance 'oligomer
                 :monomer-indices (copy-seq (monomer-indices oligomer))
                 :oligomer-space (oligomer-space oligomer)))

(defun goto-sequence (oligomer index)
  (let* ((bases (loop for monomer across (monomers (oligomer-space oligomer))
                      collect (length (monomers monomer))))
         (digits (sys:positive-integer-to-mixed-base-digits index bases)))
    (loop for monomer-index below (length (monomers (oligomer-space oligomer)))
          for digit in digits
          do (setf (aref (monomer-indices oligomer) monomer-index) digit))))

(defun oligomer-monomer-name-for-monomer (oligomer monomer)
  "Return the monomer name at index in the oligomer"
  (let ((monomer-index (position monomer (monomers oligomer))))
    (elt (monomers monomer) (elt (monomer-indices oligomer) monomer-index))))

(defun oligomer-monomer-name-at-index (oligomer index)
  "Return the monomer name at index in the oligomer"
  (let ((monomer (elt (monomers oligomer) index)))
    (elt (monomers monomer) (elt (monomer-indices oligomer) index))))

(defun oligomer-monomer-names (oligomer)
  "Return a list of monomer names for this oligomer"
  (loop for index below (length (monomers oligomer))
        for monomer across (monomers oligomer)
        collect (elt (monomers monomer) (elt (monomer-indices oligomer) index))))

(defun make-oligomer (oligomer-space &optional (index 0))
  "Build an oligomer in the oligomer space"
  (let ((olig (make-instance 'oligomer
                             :monomer-indices (make-array (length (monomers oligomer-space)) :element-type 'ext:byte32)
                             :oligomer-space oligomer-space)))
    (goto-sequence olig index)
    olig))

(defmethod monomers ((obj oligomer))
  (monomers (oligomer-space obj)))

(defmethod couplings ((obj oligomer))
  (couplings (oligomer-space obj)))


(defun oligomer-space-directional-coupling-iterator-factory (oligomer-space)
  "Return a function that takes an oligomer and returns an
   iterator that iterates over directional couplings in oligomer-space and return (values coupling source-monomer-name target-monomer-name)"
  (let ((all-couplings (loop for coupling across (couplings oligomer-space)
                                   when (typep coupling 'directional-coupling)
                                     collect coupling))
        (monomer-to-index (let ((ht (make-hash-table)))
                            (loop for monomer across (monomers oligomer-space)
                                  for index from 0
                                  do (setf (gethash monomer ht) index))
                            ht)))
    (lambda (oligomer)
      "This is a factory lambda that takes an oligomer and returns an iterator that iterates over the directional couplings in the oligomer"
      (let ((remaining-couplings all-couplings))
        (lambda ()
          "An iterator that iterates over couplings in the oligomer-space/oligomer pair.
              This iterator is closed over an oligomer and it returns the names of monomers in the coupling of the oligomer.
              Couplings are really part of the oligomer-space and so we can precompute them and use them to iterate over multiple oligomers.
              Calling the iterator returns successive (values coupling source-monomer-name target-monomer-name)"
          (when remaining-couplings
            (let* ((coupling (car remaining-couplings))
                   (source-monomer (source-monomer coupling))
                   (source-monomer-index (gethash source-monomer monomer-to-index))
                   (source-monomer-name (elt (monomers source-monomer)
                                             (elt (monomer-indices oligomer) source-monomer-index)))
                   (target-monomer (target-monomer coupling))
                   (target-monomer-index (gethash target-monomer monomer-to-index))
                   (target-monomer-name (elt (monomers target-monomer)
                                             (elt (monomer-indices oligomer) target-monomer-index))))
              (setf remaining-couplings (cdr remaining-couplings))
              (values coupling source-monomer-name target-monomer-name))))))))

(defun formula-for-topology (topology)
  (let* ((ht (make-hash-table))
         (constitution (constitution topology))
         (const-atoms (constitution-atoms constitution)))
    (loop for const-atom across const-atoms
          for element = (element const-atom)
          do (incf (gethash element ht 0)))
    (alexandria:hash-table-alist ht)))

(defun make-topology-from-residue (residue topology-name constitution)
  (let ((stereoisomer-name topology-name))
    (unless constitution
      (error "The constitution must be defined"))
    (let ((ca (topology:constitution-atoms constitution))
          (topology (make-instance 'topology:topology
                                   :name topology-name
                                   :constitution constitution))
          (stereoisomer-atoms (make-array (chem:content-size residue)))
          (stereoconfigurations nil))
      (chem:do-atoms (atm residue)
        (let* ((constitution-index (loop for idx below (length (constitution-atoms constitution))
                                         for constitution-atom = (elt (constitution-atoms constitution) idx)
                                         when (eq (atom-name constitution-atom) (chem:get-name atm))
                                           do (return idx)))
               (sai (make-instance 'topology:stereoisomer-atom
                                   :atom-name (chem:get-name atm)
                                   :atom-type (chem:atom-type atm)
                                   :atom-charge (chem:get-charge atm)
                                   :constitution-atom-index constitution-index)))
          (setf (elt stereoisomer-atoms constitution-index) sai)
          #+(or)(format t "make-topology-from-residue atm ~s configuration ~s stereochemistry-type ~s~%"
                  atm
                  (chem:get-configuration atm)
                  (chem:get-stereochemistry-type atm))

          ))
      (let* ((stereoisomer (make-instance 'stereoisomer
                                          :name (chem:get-name residue)
                                          :pdb (chem:get-pdb-name residue)
                                          :stereoisomer-atoms stereoisomer-atoms
                                          :stereoconfigurations stereoconfigurations)))
        (setf (stereoisomer topology) stereoisomer)
        (values topology constitution)))))


(defun build-residue-for-topology (topology)
  (let* ((stereoisomer (stereoisomer topology))
         (residue-name (name stereoisomer))
         (residue (chem:make-residue residue-name))
         #|         
         CL_DEFMETHOD Residue_sp Topology_O::buildResidueForIsomer(size_t isomer) const ; ; ;
         {                              ; ; ;
         //  core::write_bf_stream(fmt::sprintf("%s:%d Topology_O::buildResidueForIsomer\n" , __FILE__ , __LINE__ )); ; ; ;
         StereoisomerAtoms_sp info = this->_StereoisomerAtomProperties[isomer]; ; ; ;
         //  printf("%s:%d buildResidueForIsomer isomer = %lu  stereoisomerAtoms = %s\n", __FILE__, __LINE__, isomer, _rep_(info).c_str()); ; ; ;
         LOG("creating residue\n");     ; ; ;
         core::Symbol_sp residueName = info->getName(); ; ; ;
         if (residueName.unboundp()) {  ; ; ;
         SIMPLE_ERROR(("residueName for %s was unbound") , _rep_(info)); ; ; ;
         }                              ; ; ;
         Residue_sp res = Residue_O::make(residueName); ; ; ;
         |#
         (constitution-atoms (constitution-atoms (constitution topology)))
         (num-atoms (length constitution-atoms))
         (atoms (make-array num-atoms)))
    ;; Create the atoms
    (chem:resize-contents residue num-atoms)
    (loop for idx below num-atoms
          for ca = (elt constitution-atoms idx)
          for sia = (elt (stereoisomer-atoms stereoisomer) idx)
          for atm = (chem:make-atom (atom-name ca) (element ca))
          do (chem:set-atom-type atm (atom-type sia))
          do (chem:set-charge atm (atom-charge sia))
          do (chem:set-properties atm (copy-seq (properties ca)))
          do (chem:set-content-at residue idx atm)
          do (chem:setf-needs-build atm t)
          )
    #|
    ConstitutionAtoms_sp constitutionAtoms = this->_Constitution->getConstitutionAtoms(); ; ; ;
    size_t numAtoms = constitutionAtoms->numberOfAtoms(); ; ; ;
    gctools::Vec0<Atom_sp> atoms;       ; ; ;
    atoms.resize(numAtoms);             ; ; ;
    res->resizeContents(numAtoms);      ; ; ;
    size_t idx = 0;                     ; ; ;
    for ( size_t idx=0, idxEnd(numAtoms); idx<idxEnd; ++idx ) { ; ; ;
    StereoisomerAtom_sp ai = (*info)[idx]; ; ; ;
    Atom_sp atom = Atom_O::create();    ; ; ;
    ConstitutionAtom_sp ca = (*constitutionAtoms)[ai->constitutionAtomIndex()]; ; ; ;
    atom->setName(ca->_AtomName);       ; ; ;
    atom->setElement(ca->_Element);     ; ; ;
    atom->setType(ai->_AtomType);       ; ; ;
    atom->setCharge(ai->_AtomCharge);   ; ; ;
    atom->setProperties(cl__copy_seq(ca->_Properties)); ; ; ;
    atom->turnOnFlags(needsBuild);      ; ; ;
    //    printf("%s:%d  Creating atom@%d -> %s\n", __FILE__, __LINE__, ai->_ConstitutionAtomIndex, _rep_(atom).c_str()); ; ; ;
    atoms[ai->_ConstitutionAtomIndex] = atom; ; ; ;
    if (ai->_ConstitutionAtomIndex != idx) { ; ; ;
    SIMPLE_ERROR(("The atom %s _ConstitutionAtomIndex %lu does not match the StereoisomerAtoms idx %lu") , _rep_(atom) , ai->_ConstitutionAtomIndex , idx ); ; ; ;
    }                                   ; ; ;
    res->putMatter(idx,atom); // order atoms as in Topology ; ; ;
    }                                   ; ; ;
    |#

    ;; Create the bonds
    (loop for idx below num-atoms
          for from-atm = (chem:content-at residue idx)
          for ca = (elt (constitution-atoms (constitution topology)) idx)
          do (loop for constitution-bond in (bonds ca)
                   for to-atm = (chem:content-at residue (to-atom-index constitution-bond))
                   when (not (chem:is-bonded-to from-atm to-atm))
                     do (let ((order (order constitution-bond)))
                          (chem:bond-to from-atm to-atm order))))
    #|
    for ( size_t i=0, iEnd(constitutionAtoms->numberOfAtoms()); i<iEnd; ++i ) { ; ; ;
    Atom_sp fromAtom = atoms[i];        ; ; ;
    ConstitutionAtom_sp ca = (*constitutionAtoms)[i]; ; ; ;
    //    printf("%s:%d @%zu fromAtom -> %s\n", __FILE__, __LINE__, i, _rep_(fromAtom).c_str()); ; ; ;
    for ( auto bi=ca->_Bonds.begin(); bi!=ca->_Bonds.end(); ++bi ) ; ; ;
    {                                   ; ; ;
    if ((*bi)->_ToAtomIndex>=atoms.size()) { ; ; ;
    SIMPLE_ERROR(("Atom index %d out of bounds (num-atoms %d)") , (*bi)->_ToAtomIndex , atoms.size()); ; ; ;
    }                                   ; ; ;
    Atom_sp toAtom = atoms[(*bi)->_ToAtomIndex]; ; ; ;
    //      printf("%s:%d     @%d toAtom -> %s\n", __FILE__, __LINE__, (*bi)->_ToAtomIndex, _rep_(toAtom).c_str()); ; ; ;
    if ( !fromAtom->isBondedTo(toAtom) ) { ; ; ;
    BondOrder order = (*bi)->_BondOrder; ; ; ;
    Atom_sp tempFromAtom = fromAtom;    ; ; ;
    Atom_sp tempToAtom = toAtom;        ; ; ;
    Bond_O::canonicalizeBondOrder(tempFromAtom,tempToAtom,order); ; ; ;
    tempFromAtom->bondTo(tempToAtom,order); ; ; ;
    }                                   ; ; ;
    }                                   ; ; ;
    }                                   ; ; ;
    |#
    ;; Now add stereochemical restraints
    (let ((stereoconfigurations (stereoconfigurations stereoisomer)))
      (loop for configuration in stereoconfigurations
            for name = (atom-name configuration)
            for residue-atom = (chem:atom-with-name residue name)
            do (chem:set-stereochemistry-type residue-atom (stereochemistry-type configuration))
            do (chem:set-configuration residue-atom (configuration configuration))))
    #|
    // Now add stereochemical restraints ; ; ;
                                        ; ; ;
    //                                  ; ; ;
    // Set chiral restraints            ; ; ;
    //                                  ; ; ;
    Stereoisomer_sp si = this->_Constitution->getStereoisomerWithName(residueName); ; ; ;
    gctools::Vec0<StereoConfiguration_sp>::iterator	sci; ; ; ;
    for (sci=si->_Configurations_begin();sci!=si->_Configurations_end();sci++){ ; ; ;
    core::T_sp name = (*sci)->getAtomName(); ; ; ;
    core::T_mv aa_mv = res->atomWithName(name); ; ; ;
    Atom_sp aa = gc::As<Atom_sp>(aa_mv); ; ; ;
    LOG("Setting the configuration of atom(%s) to(%s)" , aa->description().c_str() , _rep_((*sci)->getConfiguration())  ); // ; ; ;
    auto trans = translate::from_object<ConfigurationEnum,std::true_type>((*sci)->getConfiguration())._v; ; ; ;
    aa->setStereochemistryType(chiralCenter); ; ; ;
    aa->setConfiguration(trans);        ; ; ;
    }                                   ; ; ;
    |#
    ;; Now add dihedral restraints for E/Z pi bonds
    (let ((molecule (chem:make-molecule)))
      (chem:add-matter molecule residue)
      (let ((cip-priorities (chem:assign-priorities-hash-table molecule)))
        (loop for complex-restraint in (complex-restraints stereoisomer)
              do (chem:fill-restraints complex-restraint residue cip-priorities))))
    residue)
)

(defun build-residue-single-name (topology)
  (let ((residue (build-residue-for-topology topology)))
    residue))

;;;
;;; Functions used by C++ code for assigning types
;;;
(defun chem:stereoisomer-atoms (topology)
  (let* ((stereoisomer (stereoisomer topology))
         (stereoisomer-atoms (stereoisomer-atoms stereoisomer)))
    stereoisomer-atoms))

(defun chem:stereoisomer-atom-with-name (stereoisomer-atoms name)
  (find name stereoisomer-atoms :key #'atom-name))

(defun chem:stereoisomer-atom-type (stereoisomer-atom)
  (atom-type stereoisomer-atom))

