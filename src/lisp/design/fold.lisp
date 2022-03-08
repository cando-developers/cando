(in-package :kin)

(defclass atmatter ()
  (()))


(defclass ataggregate (atmatter)
  ((aggregate :initarg :aggregate :accessor aggregate)
   (atmolecules :initarg :atmolecules :initform (make-array 0 :adjustable t) :accessor atmolecules)))

(defun put-atmolecule (ataggregate atmolecule index)
  (check-type ataggregate ataggregate)
  (check-type atmolecule atmolecule)
  (setf (aref (atmolecules ataggregate) index) atmolecule))

(defmethod resize-atmolecules ((ataggregate ataggregate) num-atmolecules)
  (adjust-array (atmolecules ataggregate) num-atmolecules))

(defclass atmolecule (atmatter)
  ((molecule :initarg :molecule :accessor molecule)
   (root-atresidue :initarg :root-atresidue :accessor root-atresidue)
   (atresidues :initarg :atresidues
               :initform (make-array 0 :adjustable t)
               :accessor atresidues)))

(defun put-atresidue (atmolecule atresidue index)
  (check-type atmolecule atmolecule)
  (check-type atresidue atresidue)
  (setf (aref (atresidues atmolecule) index) atresidue))

(defun atresidue-factory (residue ring-closing-monomer-map monomer)
  (cond
    ((chem:|Monomer_O::hasRingClosingOutPlug| monomer)
     (let ((atresidue (make-instance 'ring-closing-atresidue :residue residue)))
       (setf (gethash ring-closing-monomer-map monomer atresidue) atresidue)
       atresidue))
    (t (make-instance 'atresidue :residue residue)
       )))


(defun lookup-atresidue-id (atmolecule atresidue-id)
  (aref (atresidues atmolecule) atresidue-id))

(defun make-ring-closing-connections (ring-closing-monomer-map)
  (when (> (hash-table-count ring-closing-monomer-map) 0)
    (error "Implement me")))
             

(defun build-atmolecule-using-oligomer (oligomer molecule monomer-positions)
  (let* ((root-monomer (chem:|Oligomer_O::rootMonomer| oligomer))
         (ring-closing-monomer-map (make-hash-table))
         (atmolecule (make-instance 'atmolecule :molecule molecule))
         (residue (chem:content-at molecule (gethash root-monomer monomer-positions)))
         (root-atresidue (atresidue-factory residue ring-closing-monomer-map root-monomer)))
    (adjust-array (atresidues atmolecule) (chem:|Oligomer_O::numberOfMonomers| oligomer))
    (put-atresidue atmolecule root-atresidue (gethash root-monomer monomer-positions))
    (recursively-build-children root-atresidue atmolecule molecule ring-closing-monomer-map nil nil root-monomer monomer-positions)
    (make-ring-closing-connections ring-closing-monomer-map)
    (setf (root-atresidue atmolecule) root-atresidue)
    atmolecule))

(defclass atresidue (atmatter)
  ((residue :initarg :residue :accessor residue)
   (parent :initarg :parent :initform nil :accessor parent)
   (parent-plug-name :initarg :parent-plug-name :initform nil :accessor parent-plug-name)
   (id :initarg :id :accessor id)
   (children :type hash-table :initform (make-hash-table) :initarg :children :accessor children)
   (stereoisomer-name :initarg :stereoisomer-name :accessor stereoisomer-name)
   (topology :initarg :topology :accessor topology)
   (conformation-index :initarg :conformation-index :accessor conformation-index)
   (joints :initarg :joints
                   :initform (make-array 0 :adjustable t)
                   :accessor joints)))

(defclass ring-closing-atresidue (atresidue)
  ())

(defun recursively-build-children (atresidue
                                   atmolecule
                                   molecule
                                   ring-closing-monomer-map
                                   parent
                                   coupling
                                   monomer
                                   monomer-positions)
  "Recursively build a atmolecule from an oligomer by linking together kin:atresidues"
  (when parent
    (setf (parent atresidue) parent))
  (when coupling
    (setf (parent-plug-name atresidue) (chem:get-target-monomer-plug-name coupling)))
  (setf (stereoisomer-name atresidue) (chem:current-stereoisomer-name monomer)
        (topology atresidue) (chem:current-topology monomer)
        (conformation-index atresidue) 0)
  (loop for plug-and-coupling in (chem:plug-names-and-couplings-as-list monomer)
        for coupling = (second plug-and-coupling)
        unless (chem:|Coupling_O::isRingClosing| coupling)
          do (let ((directional-coupling coupling))
               (when (eq (chem:get-source-monomer directional-coupling)
                         monomer)
                 (let* ((other-monomer (chem::get-target-monomer directional-coupling))
                        (other-residue (chem:content-at molecule (gethash other-monomer monomer-positions)))
                        (other-atresidue (atresidue-factory other-residue
                                                            ring-closing-monomer-map
                                                            other-monomer))
                        (out-plug-name (chem:get-source-monomer-plug-name coupling)))
                   (put-atresidue atmolecule other-atresidue (gethash other-monomer monomer-positions))
                   (setf (gethash out-plug-name (children atresidue)) other-atresidue)
                   (recursively-build-children other-atresidue
                                               atmolecule
                                               molecule
                                               ring-closing-monomer-map
                                               atresidue
                                               directional-coupling
                                               other-monomer
                                               monomer-positions))))))

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
  (when (< (length (joints atresidue)) index)
    (adjust-array (joints atresidue) (1+ index)))
  (setf (aref (joints atresidue) index) joint))

(defun walk-joints (atresidue function)
  (loop for joint across (joints atresidue)
        for joint-index from 0
        do (funcall function joint-index joint)))

                          
(eval-when (:compile-toplevel :execute :load-toplevel)
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
     ))
  )
