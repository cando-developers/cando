(in-package :topology)

(defclass atmatter ()
  ((name :initarg :name :accessor name)))

(defmethod print-object ((obj atmatter) stream)
  (print-unreadable-object (obj stream :type t :identity t)
    (when (slot-boundp obj 'name)
      (format stream "~s" (name obj)))))

(defclass ataggregate (atmatter)
  ((aggregate :initarg :aggregate :accessor aggregate)
   (atmolecules :initarg :atmolecules :initform (make-array 0 :adjustable t) :accessor atmolecules)))

(defun put-atmolecule (ataggregate atmolecule index)
  (check-type ataggregate ataggregate)
  (check-type atmolecule atmolecule)
  (setf (aref (atmolecules ataggregate) index) atmolecule))

(defmethod resize-atmolecules ((ataggregate ataggregate) num-atmolecules)
  (adjust-array (atmolecules ataggregate) num-atmolecules))


(defun walk-atresidue-joints (atresidue callback)
  "Walk the joints of the atresidue"
  (loop for atom-index below (length (joints atresidue))
        for joint = (aref (joints atresidue) atom-index)
        do (unless joint (error "joint is NIL in atresidue ~s" atresidue))
        do (funcall callback joint (list :no-molecule-index :no-residue-index atom-index))))

(defun walk-atmolecule-joints (atmolecule callback)
  "Walk the joints of the atmolecule"
  (loop for residue-index below (length (atresidues atmolecule))
        for atresidue = (aref (atresidues atmolecule) residue-index)
        do (loop for atom-index below (length (joints atresidue))
                 for joint = (aref (joints atresidue) atom-index)
                 do (unless joint (error "joint is NIL in atresidue ~s" atresidue))
                 do (funcall callback joint (list :no-molecule-index residue-index atom-index)))))

(defun walk-ataggregate-joints (ataggregate callback)
  "Walk the joints of the ataggregate"
  (loop for molecule-index below (length (atmolecules ataggregate))
        for atmolecule = (aref (atmolecules ataggregate) molecule-index)
        do (loop for residue-index below (length (atresidues atmolecule))
                 for atresidue = (aref (atresidues atmolecule) residue-index)
                 do (loop for atom-index below (length (joints atresidue))
                          for joint = (aref (joints atresidue) atom-index)
                          do (unless joint (error "joint is NIL in atresidue ~s" atresidue))
                          do (funcall callback joint (list molecule-index residue-index atom-index))))))

(defclass atmolecule (atmatter)
  ((molecule :initarg :molecule :accessor molecule)
   (root-atresidue :initarg :root-atresidue :accessor root-atresidue)
   (atresidues :initarg :atresidues
               :initform (make-array 16 :adjustable t :fill-pointer 0)
               :accessor atresidues)))

(defmethod resize-atresidues ((atmolecule atmolecule) num-atresidues)
  (adjust-array (atresidues atmolecule) num-atresidues))

(defun put-atresidue (atmolecule atresidue index)
  (check-type atmolecule atmolecule)
  (check-type atresidue atresidue)
  (setf (aref (atresidues atmolecule) index) atresidue))

#+(or)
(defun atresidue-factory (residue ring-closing-monomer-map monomer monomer-subset topology)
  (cond
    ((null (in-monomer-subset monomer-subset monomer))
     ;; It's not in the monomer-subset - so forget about it
     nil)
    ((let ((ring-closing-coupling (has-ring-closing-coupling monomer)))
       (when ring-closing-coupling
         (let ((atresidue (make-instance 'ring-closing-atresidue
                                         :name (chem:get-name residue)
                                         :residue residue
                                         :topology topology
                                         :ring-closing-coupling ring-closing-coupling)))
           (setf (gethash monomer ring-closing-monomer-map) atresidue)
           atresidue))))
    (t (make-instance 'atresidue :name (chem:get-name residue)
                                 :residue residue
                                 :topology topology))))

(defun new-atresidue-factory (residue ring-closing-monomer-map oligomer monomer monomer-subset topology)
  (cond
    ((let ((ring-closing-coupling (has-ring-closing-coupling monomer)))
       (when ring-closing-coupling
         (let ((atresidue (make-instance 'ring-closing-atresidue
                                         :name (chem:get-name residue)
                                         :residue residue
                                         :topology topology
                                         :monomer-index (position monomer (monomers (oligomer-space oligomer)))
                                         :ring-closing-coupling ring-closing-coupling)))
           (setf (gethash monomer ring-closing-monomer-map) atresidue)
           atresidue))))
    (t (make-instance 'atresidue :name (chem:get-name residue)
                                 :residue residue
                                 :monomer-index (position monomer (monomers (oligomer-space oligomer)))
                                 :topology topology))))

(defun lookup-atresidue-id (atmolecule atresidue-id)
  (aref (atresidues atmolecule) atresidue-id))

(defun make-ring-closing-connections (oligomer monomer-positions molecule ring-closing-monomer-map)
  (let ((ring-couplings nil))
    (when (> (hash-table-count ring-closing-monomer-map) 0)
      (maphash (lambda (monomer atresidue)
                 (pushnew (topology::ring-closing-coupling atresidue) ring-couplings))
               ring-closing-monomer-map)
      (loop for ring-coupling in ring-couplings
            for monomer1 = (topology::monomer1 ring-coupling)
            for monomer2 = (topology::monomer2 ring-coupling)
            for monomer1-pos = (gethash monomer1 monomer-positions)
            for monomer2-pos = (gethash monomer2 monomer-positions)
            for plug1 = (topology::plug1 ring-coupling)
            for plug2 = (topology::plug2 ring-coupling)
            for topology1 = (chem:find-topology (topology:oligomer-monomer-name oligomer monomer1))
            for topology2 = (chem:find-topology (topology:oligomer-monomer-name oligomer monomer2))
            for residue1 = (chem:content-at molecule (topology:residue-index monomer1-pos))
            for residue2 = (chem:content-at molecule (topology:residue-index monomer2-pos))
            do (format t "Connecting residue ~a to ~a~%" residue1 residue2)
            do (topology:connect-residues topology1 residue1 plug1 topology2 residue2 plug2)
            ))))

(defun build-atmolecule-using-oligomer (oligomer
                                        monomer-subset
                                        molecule
                                        molecule-index
                                        monomer-positions
                                        joint-tree
                                        atom-table
                                        adjustments)
  (let* ((ring-closing-monomer-map (make-hash-table))
         (monomer-out-couplings (make-hash-table))
         (atmolecule (make-instance 'atmolecule :name (chem:get-name molecule) :molecule molecule)))
    (loop for index below (length (couplings oligomer))
          for coupling = (elt (couplings oligomer) index)
          if (typep coupling 'directional-coupling)
            do (push coupling (gethash (source-monomer coupling) monomer-out-couplings)))
    (recursively-build-atmolecule nil   ; root-monomer
                                  monomer-out-couplings
                                  nil   ; nil-or-outgoing-plug-names-to-joint-map
                                  monomer-subset
                                  oligomer
                                  joint-tree
                                  atmolecule
                                  molecule-index
                                  nil   ; atresidue
                                  nil   ; atresidue-index
                                  molecule
                                  ring-closing-monomer-map
                                  nil   ; parent-atresidue
                                  nil   ; coupling
                                  monomer-positions
                                  nil   ; parent-joint
                                  atom-table
                                  adjustments
                                  )
    #+(or)(make-ring-closing-connections oligomer monomer-positions molecule ring-closing-monomer-map)
    atmolecule))

(defclass atresidue (atmatter)
  ((residue :initarg :residue :accessor residue)
   (parent :initarg :parent :initform nil :accessor parent)
   (parent-plug-name :initarg :parent-plug-name :initform nil :accessor parent-plug-name)
   (id :initarg :id :accessor id)
   (children :type hash-table :initform (make-hash-table) :initarg :children :accessor children)
   (stereoisomer-name :initform nil :initarg :stereoisomer-name :accessor stereoisomer-name)
   (topology :initarg :topology :accessor topology)
   (monomer-index :initarg :monomer-index :reader monomer-index
                  :documentation "Only for debugging - stores the position of the monomer that generated this atresidue in the oligomer")
   (rotamer-index :initarg :rotamer-index :accessor rotamer-index
                  :documentation "Only for debugging - stores the rotamer-index that the atresidue was last set to")
   (joints :initarg :joints
                   :initform (make-array 0 :adjustable t)
                   :accessor joints)))

(defun put-joint (atresidue joint index)
  (check-type atresidue atresidue)
  (check-type joint kin:joint)
  (setf (aref (joints atresidue) index) joint))

(defun joint-with-name (atresidue atom-name &optional (errorp t))
  (loop for joint across (joints atresidue)
        when (eq atom-name (kin:joint/name joint))
          do (return-from joint-with-name joint))
  (when errorp
    (error "Could not find joint with name ~a in ~s" atom-name atresidue))
  nil)


(defmethod resize-atatoms ((atresidue atresidue) num-atatoms)
  (adjust-array (joints atresidue) num-atatoms))


(defclass ring-closing-atresidue (atresidue)
  ((ring-closing-coupling :initarg :ring-closing-coupling :accessor ring-closing-coupling)))

(defun verify-atresidue (atresidue parent-atresidue parent-joint)
  (when parent-atresidue
    (unless (eq parent-atresidue (parent atresidue))
      (error "The parent-atresidue ~s has not been set into the (parent atresidue) of ~s"
             parent-atresidue atresidue)))
  (when parent-joint
    (unless (eq parent-joint (kin:parent (elt (joints atresidue) 0)))
      (error "The parent joint of the first joint ~s atresidue ~s has not been set to ~s"
             (elt (joints atresidue) 0) atresidue parent-joint))))

(defun recursively-build-atmolecule (prev-monomer ; nil or current monomer
                                     monomer-out-couplings
                                     nil-or-outgoing-plug-names-to-joint-map
                                     monomer-subset ; monomer-subset
                                     oligomer       ; oligomer
                                     joint-tree ; joint-tree to set root
                                     atmolecule ; current atmolecule we are filling
                                     atmolecule-index
                                     atresidue ; nil or atresidue for current monomer
                                     atresidue-index
                                     molecule ; current molecule cooresponding to atmolecule
                                     ring-closing-monomer-map
                                     parent-atresidue ; nil or parent-atresidue
                                     coupling ; nil or coupling
                                     monomer-positions ; map of monomers to monomer-positions
                                     parent-joint
                                     atom-table
                                     adjustments)
  "Recursively build a atmolecule from an oligomer by linking together kin:atresidues"
  (if prev-monomer
      (let ((next-atresidue nil)
            (next-atresidue-index nil)
            (next-outgoing-plug-names-to-joint-map)
            (out-couplings (gethash prev-monomer monomer-out-couplings)))
        (loop for out-coupling in out-couplings
              for next-monomer = (target-monomer out-coupling)
              for next-out-couplings = (gethash next-monomer monomer-out-couplings)
              do (when (in-monomer-subset monomer-subset next-monomer)
                   (multiple-value-setq (next-atresidue next-atresidue-index next-outgoing-plug-names-to-joint-map)
                     (let* ((next-topology (monomer-topology next-monomer oligomer))
                            (next-monomer-position (gethash next-monomer monomer-positions))
                            (next-residue-index (residue-index next-monomer-position))
                            (next-residue (chem:content-at molecule next-residue-index))
                            (next-atresidue (new-atresidue-factory next-residue
                                                                   ring-closing-monomer-map
                                                                   oligomer
                                                                   next-monomer
                                                                   monomer-subset
                                                                   next-topology))
                            (next-atresidue-index (and next-atresidue (length (atresidues atmolecule))))
                            (out-plug-name (source-plug-name out-coupling))
                            (parent-joint (and nil-or-outgoing-plug-names-to-joint-map
                                               (gethash out-plug-name nil-or-outgoing-plug-names-to-joint-map)))
                            (next-outgoing-plug-names-to-joint-map (maybe-fill-atresidue joint-tree
                                                                                         oligomer
                                                                                         next-monomer
                                                                                         monomer-subset
                                                                                         next-atresidue
                                                                                         parent-joint
                                                                                         atmolecule-index
                                                                                         next-atresidue-index
                                                                                         atom-table
                                                                                         adjustments)))
                       (setf (stereoisomer-name next-atresidue) (current-stereoisomer-name prev-monomer oligomer)
                             (topology next-atresidue) next-topology)
                       (when atresidue
                         (setf (parent next-atresidue) atresidue))
                       (verify-atresidue next-atresidue atresidue parent-joint)
                       (when (/= next-atresidue-index next-residue-index)
                         (error "There is a mismatch between the index of the atresidue ~s at ~s and the residue ~s at ~s"
                                next-atresidue next-atresidue-index
                                next-residue next-residue-index))
                       (vector-push-extend next-atresidue (atresidues atmolecule))
                       (when atresidue
                         (setf (gethash out-plug-name (children atresidue)) next-atresidue
                               (parent-plug-name atresidue) (target-plug-name out-coupling)))
                       (values next-atresidue next-atresidue-index next-outgoing-plug-names-to-joint-map))))
              do (recursively-build-atmolecule next-monomer
                                               monomer-out-couplings
                                               next-outgoing-plug-names-to-joint-map
                                               monomer-subset
                                               oligomer
                                               joint-tree
                                               atmolecule
                                               atmolecule-index
                                               next-atresidue
                                               next-atresidue-index
                                               molecule
                                               ring-closing-monomer-map
                                               atresidue
                                               out-coupling
                                               monomer-positions
                                               parent-joint
                                               atom-table
                                               adjustments)))
      (let* ((prev-monomer (root-monomer oligomer))
             (next-atresidue nil)
             (next-atresidue-index nil)
             (next-outgoing-plug-names-to-joint-map))
        (when (in-monomer-subset monomer-subset prev-monomer)
          (multiple-value-setq (next-atresidue next-atresidue-index next-outgoing-plug-names-to-joint-map)
            (let* ((next-monomer prev-monomer)
                   (next-topology (monomer-topology next-monomer oligomer))
                   (next-monomer-position (gethash next-monomer monomer-positions))
                   (next-residue-index (residue-index next-monomer-position))
                   (next-residue (chem:content-at molecule next-residue-index))
                   (next-atresidue (new-atresidue-factory next-residue
                                                          ring-closing-monomer-map
                                                          oligomer
                                                          next-monomer
                                                          monomer-subset
                                                          next-topology))
                   (next-atresidue-index (and next-atresidue (length (atresidues atmolecule))))
                   (next-outgoing-plug-names-to-joint-map (maybe-fill-atresidue joint-tree
                                                                                oligomer
                                                                                prev-monomer
                                                                                monomer-subset
                                                                                next-atresidue
                                                                                parent-joint
                                                                                atmolecule-index
                                                                                next-atresidue-index
                                                                                atom-table
                                                                                adjustments)))
              (setf (stereoisomer-name next-atresidue) (current-stereoisomer-name prev-monomer oligomer)
                    (topology next-atresidue) next-topology)
              (when parent-atresidue
                (setf (parent next-atresidue) parent-atresidue))
              (verify-atresidue next-atresidue atresidue parent-joint)
              (when (/= next-atresidue-index next-residue-index)
                (error "There is a mismatch between the index of the atresidue ~s at ~s and the residue ~s at ~s"
                       next-atresidue next-atresidue-index
                       next-residue next-residue-index))
              (vector-push-extend next-atresidue (atresidues atmolecule))
              (values next-atresidue next-atresidue-index next-outgoing-plug-names-to-joint-map))))
        (recursively-build-atmolecule prev-monomer
                                      monomer-out-couplings
                                      next-outgoing-plug-names-to-joint-map
                                      monomer-subset
                                      oligomer
                                      joint-tree
                                      atmolecule
                                      atmolecule-index
                                      next-atresidue
                                      next-atresidue-index
                                      molecule
                                      ring-closing-monomer-map
                                      atresidue
                                      nil ; coupling
                                      monomer-positions
                                      nil ; parent-joint
                                      atom-table
                                      adjustments))))

(defun describe-recursively (atresidue prefix stream)
  (princ prefix stream)
  (when (parent-plug-name atresidue)
    (princ (parent-plug-name atresidue) stream)
    (princ " " stream))
  (format stream "~a[~a]" (class-name (class-of atresidue)) (stereoisomer-name atresidue))
  (loop for child-atresidue across (children atresidue)
        do (describe-recursively child-atresidue (concatenate 'string prefix "  ") stream)))


(defmethod print-object ((object atresidue) stream)
  (print-unreadable-object (object stream :type t)
    (format stream "[~s]" (residue object))))

(defun add-joint (atresidue index joint)
  (check-type joint kin:joint)
  (when (< (length (joints atresidue)) index)
    (adjust-array (joints atresidue) (1+ index)))
  (setf (aref (joints atresidue) index) joint))

(defun walk-joints (atresidue function)
  (loop for joint across (joints atresidue)
        for joint-index from 0
        do (funcall function joint-index joint)))



#+(or)(eval-when (:compile-toplevel :execute :load-toplevel)
  (export
   '(
     atmatter
     ataggregate
     atmolecules
     atmolecule
     root-atresidue
     atresidues
     atresidue
     parent
     parent-plug-name
     id
     children
     stereoisomer-name
     topology
     rotamer-index
     joints
     put-joint
     walk-ataggregate-joints
     update-internal-coordinates
     zero-external-coordinates
     ))
  )
