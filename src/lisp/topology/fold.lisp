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


(defun walk-ataggregate-joints (ataggregate callback)
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
               :initform (make-array 0 :adjustable t)
               :accessor atresidues)))

(defmethod resize-atresidues ((atmolecule atmolecule) num-atresidues)
  (adjust-array (atresidues atmolecule) num-atresidues))

(defun put-atresidue (atmolecule atresidue index)
  (check-type atmolecule atmolecule)
  (check-type atresidue atresidue)
  (setf (aref (atresidues atmolecule) index) atresidue))

(defun atresidue-factory (residue ring-closing-monomer-map monomer topology)
  (cond
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
                                        molecule
                                        molecule-index
                                        monomer-positions
                                        joint-tree
                                        atom-table
                                        adjustments
                                        orientation)
  (let* ((root-monomer (root-monomer oligomer))
         (ring-closing-monomer-map (make-hash-table))
         (atmolecule (make-instance 'atmolecule :name (chem:get-name molecule) :molecule molecule))
         (monomer-position (gethash root-monomer monomer-positions))
         (residue-index (residue-index monomer-position))
         (residue (chem:content-at molecule residue-index))
         (topology (monomer-topology root-monomer oligomer))
         (root-atresidue (atresidue-factory residue ring-closing-monomer-map root-monomer topology)))
    (adjust-array (atresidues atmolecule) (number-of-monomers oligomer))
    (put-atresidue atmolecule root-atresidue residue-index)
    (recursively-build-children joint-tree
                                root-atresidue
                                atmolecule
                                molecule
                                ring-closing-monomer-map
                                nil
                                nil
                                root-monomer
                                oligomer
                                monomer-positions
                                molecule-index
                                residue-index
                                nil
                                atom-table
                                adjustments
                                orientation
                                )
    #+(or)(make-ring-closing-connections oligomer monomer-positions molecule ring-closing-monomer-map)
    (setf (root-atresidue atmolecule) root-atresidue)
    atmolecule))

(defclass atresidue (atmatter)
  ((residue :initarg :residue :accessor residue)
   (parent :initarg :parent :initform nil :accessor parent)
   (parent-plug-name :initarg :parent-plug-name :initform nil :accessor parent-plug-name)
   (id :initarg :id :accessor id)
   (children :type hash-table :initform (make-hash-table) :initarg :children :accessor children)
   (stereoisomer-name :initform nil :initarg :stereoisomer-name :accessor stereoisomer-name)
   (topology :initarg :topology :accessor topology)
   (conformation-index :initarg :conformation-index :accessor conformation-index)
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

(defun recursively-build-children (joint-tree
                                   atresidue
                                   atmolecule
                                   molecule
                                   ring-closing-monomer-map
                                   parent-atresidue
                                   coupling
                                   monomer
                                   oligomer
                                   monomer-positions
                                   atmolecule-index
                                   atresidue-index
                                   parent-joint
                                   atom-table
                                   adjustments
                                   orientation)
  "Recursively build a atmolecule from an oligomer by linking together kin:atresidues"
  (when parent-atresidue
    (setf (parent atresidue) parent-atresidue))
  (when coupling
    (setf (parent-plug-name atresidue) (target-plug-name coupling)))
  (let ((outgoing-plug-names-to-joint-map (fill-atresidue joint-tree
                                                          oligomer
                                                          atresidue
                                                          parent-joint
                                                          atmolecule-index
                                                          atresidue-index
                                                          atom-table
                                                          adjustments
                                                          orientation))
        (current-topology (monomer-topology monomer oligomer)))
    (setf (stereoisomer-name atresidue) (current-stereoisomer-name monomer oligomer)
          (topology atresidue) current-topology
          (conformation-index atresidue) 0)
    (maphash (lambda (plug-name coupling)
               (declare (ignore plug-name))
               (unless (typep coupling 'ring-coupling)
                 (let ((directional-coupling coupling))
                   (when (eq (source-monomer directional-coupling) monomer)
                     (let* ((other-monomer (target-monomer directional-coupling))
                            (other-topology (monomer-topology other-monomer oligomer))
                            (other-monomer-position (gethash other-monomer monomer-positions))
                            (other-residue-index (residue-index other-monomer-position))
                            (other-residue (chem:content-at molecule other-residue-index))
                            (other-atresidue (atresidue-factory other-residue
                                                                ring-closing-monomer-map
                                                                other-monomer
                                                                other-topology))
                            (out-plug-name (source-plug-name coupling)))
                       (put-atresidue atmolecule other-atresidue other-residue-index)
                       (setf (gethash out-plug-name (children atresidue)) other-atresidue)
                       (let ((new-parent-joint (gethash out-plug-name outgoing-plug-names-to-joint-map)))
                         (recursively-build-children joint-tree
                                                     other-atresidue
                                                     atmolecule
                                                     molecule
                                                     ring-closing-monomer-map
                                                     atresidue
                                                     directional-coupling
                                                     other-monomer
                                                     oligomer
                                                     monomer-positions
                                                     atmolecule-index
                                                     other-residue-index
                                                     new-parent-joint
                                                     atom-table
                                                     adjustments
                                                     orientation)))))))
             (couplings monomer))))

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
    (format stream "~a[~a]" (class-name (class-of object)) (stereoisomer-name object))))

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
     conformation-index
     joints
     put-joint
     walk-ataggregate-joints
     update-internal-coordinates
     zero-external-coordinates
     ))
  )
