
(in-package :kin)

(defgeneric children (object)
  (:documentation "Return list of children of object"))

(defmethod children ((joint joint))
  (joint-children joint))

(defclass joint-tree ()
  ((root :initarg :root :accessor root)
   (atom-id-map :initarg :atom-id-map
                :initform (chem:make-atom-id-map)
                :accessor atom-id-map)))

(cando:make-class-save-load joint-tree)

(defun make-joint-tree ()
  (make-instance 'joint-tree))


(defun recursively-write-into-atresidue (joint-template
                                         parent-joint
                                         atresidue
                                         atmolecule-index
                                         atresidue-index)
  (let ((joint (write-into-joint-tree joint-template
                                                              parent-joint
                                                              atresidue
                                                              atmolecule-index
                                                              atresidue-index)))
    (loop for child-joint-template in (children joint-template)
          do (recursively-write-into-atresidue child-joint-template
                                               joint
                                               atresidue
                                               atmolecule-index
                                               atresidue-index))
    joint))

(defun fill-atresidue (joint-tree atresidue parent-joint atmolecule-index atresidue-index)
  "Recursively build an atom-tree into the atresidues of an ataggregate.
This duplicates https://github.com/cando-developers/cando/blob/main/src/kinematics/pointTree.cc#L172
JointTree_O::recursivelyBuildMolecule
If (null parent-joint) then this is the root atresidue
"
  (let* ((topology (chem:find-topology (topology-name atresidue) t))
         (constitution (chem:topology/get-constitution topology))
         (constitution-atoms (chem:constitution/get-constitution-atoms constitution)))
    (resize-atatoms atresidue (chem:constitution-atoms/number-of-atoms constitution-atoms))
    (let* ((joint-template (chem:topology/get-property topology :joint-template))
           ;; bondid is represented by a (cons parent child)
           (root-joint (recursively-write-into-atresidue joint-template
                                                         parent-joint
                                                         atresidue
                                                         atmolecule-index
                                                         atresidue-index)))
      (unless (slot-boundp joint-tree 'root)
        (setf (root joint-tree) root-joint))
      (let ((outgoing-plug-names-to-joint-map (make-hash-table)))
        (loop for plug in (chem:topology/plugs-as-list topology)
            when (typep plug 'chem:out-plug)
              do (let* ((out-plug plug)
                        (out-plug-name (chem:plug/get-name out-plug))
                        (atomb0 (chem:plug/get-b0 out-plug))
                        (constitution-bond0-atomid (chem:constitution-atoms/index constitution-atoms atomb0))
                        (joint-bond0-parent (aref (joints atresidue) constitution-bond0-atomid)))
                   (setf (gethash out-plug-name outgoing-plug-names-to-joint-map) joint-bond0-parent)))
        outgoing-plug-names-to-joint-map))))

#|

DONT_OPTIMIZE_WHEN_DEBUG_RELEASE void JointTree_O::recursivelyBuildMolecule(MonomerId monomerId,
                                          MonomerNode_sp monomerNode,
                                          Joint_sp parent,
                                          bool rootNode)
{_OF();
  LOG(BF("recursivelyBuildMolecule into parent: %s monomerNode: %s\n") % _rep_(parent) % _rep_(monomerNode));
  chem::Constitution_mv constitutionAndTopology = monomerNode->identifyConstitutionAndTopology();
  chem::Constitution_sp constitution = constitutionAndTopology;
  chem::Topology_sp topology = constitutionAndTopology.second().as<chem::Topology_O>();
  chem::ConstitutionAtoms_sp constitutionAtoms;
  {_BLOCK_TRACEF(BF("Building constitution[%s] Topology[%s]")
                 % _rep_(constitution->getName()) % _rep_(topology->getName()) );
    constitutionAtoms = constitution->getConstitutionAtoms();
    this->resizeAtoms(monomerId._Chain,monomerId._Monomer,constitutionAtoms->numberOfAtoms());
    core::T_sp ttemplate = topology->getProperty(INTERN_(kw,jointTemplate));
    if (ttemplate.nilp()) {
      SIMPLE_ERROR(BF("The topology %s is missing an :atom-template parameter") % _rep_(topology));
    }
    JointTemplate_sp jointTemplate = gc::As<JointTemplate_sp>(ttemplate);
    BondId_sp incoming = BondId_O::create(parent);
    JointTemplate_O::PlugNamesToBondIdMap outgoing;
	    //
	    // Write the sub tree described by jointTemplate into the JointTree
	    // recursively
    ASSERTF(jointTemplate.notnilp(),BF("The JointTemplate for Topology[%s] is nil")
            % _rep_(topology->getName()) );
    jointTemplate->writeIntoJointTree(this->sharedThis<JointTree_O>(),
                                     monomerId,
                                     incoming,
                                     outgoing,
                                     monomerNode,
                                     rootNode);
  }
	    // Now loop over all the OutPlugs and
	    // create their children and connect them to our OutPlug atoms
	    //
  chem::Topology_O::Plugs::iterator it;
  for ( it=topology->_Plugs.begin(); it!=topology->_Plugs.end(); it++ ) {
    if ( gc::IsA<chem::OutPlug_sp>(it->second) ) {
      chem::OutPlug_sp outPlug = gc::As_unsafe<chem::OutPlug_sp>(it->second);
      core::Symbol_sp atomB0 = outPlug->getB0();
      int constitutionBond0AtomId = constitutionAtoms->index(atomB0);
      chem::AtomId atomId(monomerId._Chain, monomerId._Monomer,constitutionBond0AtomId);
      Joint_sp bond0Parent = this->_AtomMap[atomId];
      MonomerNode_sp nextMonomerNode = monomerNode->_Children.get(outPlug->getName());
      MonomerId nextResidueId = nextMonomerNode->_Id;
      this->recursivelyBuildMolecule(nextResidueId,nextMonomerNode,bond0Parent);
    }
  }
}

void JointTree_O::buildMoleculeUsingChainNode(int moleculeId, ChainNode_sp chainNode, chem::Oligomer_sp oligomer )
{_OF();
  ASSERTF(moleculeId<this->numberOfMolecules(),BF("Illegal moleculeId[%d]") % moleculeId );
  int numResidues = oligomer->numberOfMonomers();
  this->resizeResidues(moleculeId,numResidues);
  MonomerNode_sp monomerNode = chainNode->_RootMonomerNode;
  ASSERTF(this->_Root.notnilp(),BF("The Root of the JointTree cannot be nil"));
  LOG(BF("Building moleculeId[%d]") % moleculeId);
  MonomerId monomerId(moleculeId,monomerNode->_Id._Monomer);
  this->recursivelyBuildMolecule(monomerId,
                                 monomerNode,
                                 this->_Root,
                                 true);
}
|#
