
(in-package :kin)

(defclass conformation ()
  ((aggregate :initarg :aggregate :accessor aggregate)
   (ataggregate :initarg :ataggregate :accessor ataggregate)
   (joint-tree :initarg :joint-tree :accessor joint-tree)))

(defun make-conformation (oligomer)
  (let* ((aggregate (chem:make-aggregate))
         (ataggregate (make-instance 'ataggregate :aggregate aggregate))
         (joint-tree (make-joint-tree)))
    (resize-atmolecules ataggregate 1)
    (multiple-value-bind (molecule monomer-positions)
        (design:build-molecule oligomer)
      (format t "monomer-positions: ~a~%" (alexandria:hash-table-alist monomer-positions))
      (multiple-value-bind (molecule molecule-index)
          (chem:add-matter aggregate molecule)
        ;; This is where I would invoke Conformation_O::buildMoleculeUsingOligomer
        ;; Use the monomers-to-topologys
        (let* ((atmolecule (build-atmolecule-using-oligomer oligomer molecule monomer-positions)))
          (put-atmolecule ataggregate atmolecule molecule-index)
          (let ((conf (make-instance 'conformation
                                     :aggregate aggregate
                                     :ataggregate ataggregate
                                     :joint-tree joint-tree)))
            conf
            ))))))






#||
This part is confusing - how do I rearrange it?

void Conformation_O::buildMoleculeUsingOligomer(int moleculeId, chem::Oligomer_sp oligomer)
{_OF();
  ASSERTF(moleculeId<this->numberOfMolecules(),
          BF("Illegal moleculeId[%d]") % moleculeId );
  ChainNode_sp chainNode =
    this->_FoldTree->buildChainUsingOligomer(moleculeId,oligomer);
  LOG(BF("Built FoldTree--->\n%s")%chainNode->_RootMonomerNode->asString());
  this->_JointTree->buildMoleculeUsingChainNode(moleculeId,chainNode,oligomer);
}

||#


(export
 '(
   conformation
   make-conformation
   build-molecule-using-oligomer
   )
   :kin)
