(in-package :topology)

(defgeneric children (object)
  (:documentation "Return list of children of object"))

(defmethod children ((joint kin:joint))
  (kin:joint-children joint))

(defclass joint-tree ()
  ((root-map :initform (make-hash-table) :initarg :root-map :accessor root-map)
   (atom-id-map :initarg :atom-id-map
                :initform (chem:make-atom-id-map)
                :accessor atom-id-map)))

(defmethod print-object ((obj joint-tree) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream ":root-map ~s" (alexandria:hash-table-keys (root-map obj)))))

(defun make-joint-tree ()
  (make-instance 'joint-tree))

(defun recursively-write-into-atresidue (joint-template
                                         parent-joint
                                         atresidue
                                         atmolecule-index
                                         atresidue-index
                                         atom-table
                                         adjustments
                                         orientation)
  (let ((joint (write-into-joint-tree joint-template
                                      parent-joint
                                      atresidue
                                      atmolecule-index
                                      atresidue-index
                                      atom-table
                                      adjustments
                                      orientation)))
    (loop for child-joint-template in (children joint-template)
          do (recursively-write-into-atresidue child-joint-template
                                               joint
                                               atresidue
                                               atmolecule-index
                                               atresidue-index
                                               atom-table
                                               adjustments
                                               orientation))
    joint))

(defun fill-atresidue (joint-tree
                       oligomer
                       atresidue
                       parent-joint
                       atmolecule-index
                       atresidue-index
                       atom-table
                       adjustments
                       orientation)
  "Recursively build an atom-tree into the atresidues of an ataggregate.
This duplicates https://github.com/cando-developers/cando/blob/main/src/kinematics/pointTree.cc#L172
JointTree_O::recursivelyBuildMolecule
If (null parent-joint) then this is the root atresidue
"
  (let* ((topology (topology atresidue))
         (constitution (constitution topology))
         (constitution-atoms (constitution-atoms constitution)))
    (resize-atatoms atresidue (length constitution-atoms))
    (let* ((joint-template (joint-template topology))
           ;; bondid is represented by a (cons parent child)
           (root-joint (recursively-write-into-atresidue joint-template
                                                         parent-joint
                                                         atresidue
                                                         atmolecule-index
                                                         atresidue-index
                                                         atom-table
                                                         adjustments
                                                         orientation)))
      (when (null (gethash oligomer (root-map joint-tree)))
        (setf (gethash oligomer (root-map joint-tree)) root-joint))
      (let ((outgoing-plug-names-to-joint-map (make-hash-table)))
        (maphash (lambda (plug-name plug)
                   (declare (ignore plug-name))
                   (when (typep plug 'out-plug)
                     (let* ((out-plug plug)
                            (out-plug-name (name out-plug))
                            (plug-bond (elt (plug-bonds out-plug) 0))
                            (atomb0-name (atom-name plug-bond))
                            (constitution-bond0-atomid (position atomb0-name constitution-atoms :key #'atom-name))
                            (joint-bond0-parent (if constitution-bond0-atomid
                                                    (aref (joints atresidue) constitution-bond0-atomid)
                                                    (error "Could not find atom ~a in ~a" atomb0-name constitution-atoms))))
                       (setf (gethash out-plug-name outgoing-plug-names-to-joint-map) joint-bond0-parent))))
                 (plugs topology))
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
