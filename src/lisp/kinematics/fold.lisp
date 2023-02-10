(in-package :kin)

(defclass atmatter ()
  (()))

(defclass ataggregate (atmatter)
  ((aggregate :initarg :aggregate :accessor aggregate)
   (atmolecules :initarg :atmolecules :initform (make-array 0 :adjustable t) :accessor atmolecules)))

(defun atmolecule-aref  (ataggregate atmolecule-id)
  (aref (atmolecules ataggregate) atmolecule-id))

(cando:make-class-save-load ataggregate)

(defun put-atmolecule (ataggregate atmolecule index)
  (check-type ataggregate ataggregate)
  (check-type atmolecule atmolecule)
  (setf (aref (atmolecules ataggregate) index) atmolecule))

(defun atresidue-aref  (atmolecule atresidue-id)
  (aref (atresidues atmolecule) atresidue-id))


(defmethod resize-atmolecules ((ataggregate ataggregate) num-atmolecules)
  (adjust-array (atmolecules ataggregate) num-atmolecules))


(defun walk-ataggregate-joints (ataggregate callback)
  (loop for molecule-index below (length (atmolecules ataggregate))
        for atmolecule = (aref (atmolecules ataggregate) molecule-index)
        do (loop for residue-index below (length (atresidues atmolecule))
                 for atresidue = (aref (atresidues atmolecule) residue-index)
                 do (loop for atom-index below (length (joints atresidue))
                          for joint = (aref (joints atresidue) atom-index)
                          do (funcall callback joint (list molecule-index residue-index atom-index))))))

(defclass atmolecule (atmatter)
  ((molecule :initarg :molecule :accessor molecule)
   (root-atresidue :initarg :root-atresidue :accessor root-atresidue)
   (atresidues :initarg :atresidues
               :initform (make-array 0 :adjustable t)
               :accessor atresidues)))

(cando:make-class-save-load atmolecule)

(defmethod resize-atresidues ((atmolecule atmolecule) num-atresidues)
  (adjust-array (atresidues atmolecule) num-atresidues))

(defun put-atresidue (atmolecule atresidue index)
  (check-type atmolecule atmolecule)
  (check-type atresidue atresidue)
  (setf (aref (atresidues atmolecule) index) atresidue))

(defun atresidue-factory (residue ring-closing-monomer-map monomer topology)
  (cond
    ((chem:monomer/has-ring-closing-out-plug monomer)
     (let ((atresidue (make-instance 'ring-closing-atresidue :residue residue)))
       (setf (gethash ring-closing-monomer-map monomer atresidue) atresidue)
       atresidue))
    (t (make-instance 'atresidue :residue residue
                                 :topology-name (chem:get-name topology)))))

(defun lookup-atresidue-id (atmolecule atresidue-id)
  (aref (atresidues atmolecule) atresidue-id))

(defun make-ring-closing-connections (ring-closing-monomer-map)
  (when (> (hash-table-count ring-closing-monomer-map) 0)
    (error "Implement me")))
             

(defun build-atmolecule-using-oligomer (oligomer molecule molecule-index monomer-positions joint-tree atom-tree)
  (let* ((root-monomer (chem:oligomer/root-monomer oligomer))
         (ring-closing-monomer-map (make-hash-table))
         (atmolecule (make-instance 'atmolecule :molecule molecule))
         (residue-index (gethash root-monomer monomer-positions))
         (residue (chem:content-at molecule residue-index))
         (topology (chem:current-topology root-monomer))
         (root-atresidue (atresidue-factory residue ring-closing-monomer-map root-monomer topology)))
    (adjust-array (atresidues atmolecule) (chem:oligomer/number-of-monomers oligomer))
    (put-atresidue atmolecule root-atresidue residue-index)
    (recursively-build-children joint-tree
                                root-atresidue
                                atmolecule
                                molecule
                                ring-closing-monomer-map
                                nil
                                nil
                                root-monomer
                                monomer-positions
                                molecule-index
                                residue-index
                                nil
                                atom-tree
                                )
    (make-ring-closing-connections ring-closing-monomer-map)
    (setf (root-atresidue atmolecule) root-atresidue)
    atmolecule))

(defclass atresidue (atmatter)
  ((residue :initarg :residue :accessor residue)
   (parent :initarg :parent :initform nil :accessor parent)
   (parent-plug-name :initarg :parent-plug-name :initform nil :accessor parent-plug-name)
   (id :initarg :id :accessor id)
   (children :type hash-table :initform (make-hash-table) :initarg :children :accessor children)
   (stereoisomer-name :initform nil :initarg :stereoisomer-name :accessor stereoisomer-name)
   (topology-name :initarg :topology-name :accessor topology-name)
   (conformation-index :initarg :conformation-index :accessor conformation-index)
   (joints :initarg :joints
                   :initform (make-array 0 :adjustable t)
                   :accessor joints)))

(cando:make-class-save-load
 atresidue
 :print-unreadably (lambda (object stream)
                     (format stream "~a[~a]" (class-name (class-of object)) (chem:get-name (residue object)))))

(defun put-joint (atresidue joint index)
  (check-type atresidue atresidue)
  (check-type joint joint)
  (setf (aref (joints atresidue) index) joint))


(defun joint-count (atresidue)
  (length (joints atresidue)))

(defun joint-aref (atresidue joint-id)
  (aref (joints atresidue) joint-id))



(defmethod resize-atatoms ((atresidue atresidue) num-atatoms)
  (adjust-array (joints atresidue) num-atatoms))


(defclass ring-closing-atresidue (atresidue)
  ())

(defun recursively-build-children (joint-tree
                                   atresidue
                                   atmolecule
                                   molecule
                                   ring-closing-monomer-map
                                   parent-atresidue
                                   coupling
                                   monomer
                                   monomer-positions
                                   atmolecule-index
                                   atresidue-index
                                   parent-joint
                                   atom-tree)
  "Recursively build a atmolecule from an oligomer by linking together atresidues"
  (when parent-atresidue
    (setf (parent atresidue) parent-atresidue))
  (when coupling
    (setf (parent-plug-name atresidue) (chem:get-target-monomer-plug-name coupling)))
  (let ((outgoing-plug-names-to-joint-map (fill-atresidue joint-tree atresidue parent-joint atmolecule-index atresidue-index atom-tree))
        (current-topology (chem:current-topology monomer)))
    (setf (stereoisomer-name atresidue) (chem:current-stereoisomer-name monomer)
          (topology-name atresidue) (chem:get-name current-topology)
          (conformation-index atresidue) 0)
    (loop for plug-and-coupling in (chem:plug-names-and-couplings-as-list monomer)
          for coupling = (second plug-and-coupling)
          unless (chem:coupling/is-ring-closing coupling)
            do (let ((directional-coupling coupling))
                 (when (eq (chem:get-source-monomer directional-coupling) monomer)
                   (let* ((other-monomer (chem::get-target-monomer directional-coupling))
                          (other-topology (chem:current-topology other-monomer))
                          (other-residue-index (gethash other-monomer monomer-positions))
                          (other-residue (chem:content-at molecule other-residue-index))
                          (other-atresidue (atresidue-factory other-residue
                                                              ring-closing-monomer-map
                                                              other-monomer
                                                              other-topology))
                          (out-plug-name (chem:get-source-monomer-plug-name coupling)))
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
                                                   monomer-positions
                                                   atmolecule-index
                                                   other-residue-index
                                                   new-parent-joint
                                                   atom-tree))))))))

(defun describe-recursively (atresidue prefix stream)
  (princ prefix stream)
  (when (parent-plug-name atresidue)
    (princ (parent-plug-name atresidue) stream)
    (princ " " stream))
  (format stream "~a[~a]" (class-name (class-of atresidue)) (stereoisomer-name atresidue))
  (loop for child-atresidue across (children atresidue)
        do (describe-recursively child-atresidue (concatenate 'string prefix "  ") stream)))


(defun add-joint (atresidue index joint)
  (when (< (length (joints atresidue)) index)
    (adjust-array (joints atresidue) (1+ index)))
  (setf (aref (joints atresidue) index) joint))

(defun walk-joints (atresidue function)
  (loop for joint across (joints atresidue)
        for joint-index from 0
        do (funcall function joint-index joint)))



