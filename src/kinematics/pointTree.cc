/*
    File: jointTree.cc
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#define DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <clasp/core/wrappers.h>
#include <cando/chem/chemPackage.h>
#include <cando/chem/constitution.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/oligomer.h>
#include <cando/chem/topology.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/candoDatabase.h>
#include <cando/kinematics/jointTemplate.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/monomerNode.h>
#include <cando/kinematics/chainNode.h>
#include <cando/kinematics/joint.h>
#include <cando/kinematics/originJumpJoint.h>
#include <cando/kinematics/rootBondedJoint.h>
#include <cando/kinematics/jointTree.h>

namespace kinematics
{

void JointTree_O::fields(core::Record_sp node) {
  node->field_if_not_unbound(INTERN_(kw,root),this->_Root);
};



string AtomHolder::typeAsString() const
{
  switch (this->_Type)
  {
  case unused:
      return "unused";
  case jumpJoint:
      return "jump";
  case originJumpJoint:
      return "origin";
  case bondedJoint:
      return "bonded";
  case delayedBondedJoint:
      return "delayed";
  default:
      return "unknown";
  };
};


void JointTree_O::initialize()
{
  chem::AtomId id;
  this->_Root = this->newOriginJumpJoint(id,_Nil<core::T_O>(),"-root-");
}


Joint_sp JointTree_O::lookup(const chem::AtomId& atomId) const
{
  Joint_sp atom = this->_AtomMap[atomId];
  return atom;
}

Joint_sp JointTree_O::JointTreeLookupAtomid(const chem::AtomId& atomId) const
{
  return this->lookup(atomId);
}


void JointTree_O::resizeMolecules(int numMolecules)
{_OF();
  LOG(BF("Resizing JointTree to %d molecules") % numMolecules);
//	int originalNumMolecules = this->_AtomMap.numberOfMolecules();
  this->_AtomMap.resize(numMolecules);
};

void JointTree_O::resizeResidues(int mol, int numResidues)
{_OF();
  LOG(BF("Resizing JointTree to %d residues in moleculeId[%d]") % numResidues % mol );
  this->_AtomMap.resize(mol,numResidues);
};

void JointTree_O::resizeAtoms(int mol, int res, int numAtoms)
{_OF();
  LOG(BF("Resizing JointTree to %d atoms in residue[%d] moleculeId[%d]")
      % numAtoms 
      % res % mol );
  this->_AtomMap.resize(mol,res,numAtoms);
};

Joint_sp JointTree_O::_initializeNewAtom(Joint_sp atom, const chem::AtomId& atomId )
{_OF();
  JointTree_sp tree = this->sharedThis<JointTree_O>();
  if ( atomId.isDefined() ) this->updateAtomIdMap(atomId,atom);
  return atom;
}


Joint_sp JointTree_O::newJumpJoint(const chem::AtomId& atomId, core::T_sp name, const string& comment)
{_OF();
  return this->_newAtom<JumpJoint_O>(atomId,name,comment);
};

Joint_sp JointTree_O::newOriginJumpJoint(const chem::AtomId& atomId, core::T_sp name,const string& comment)
{_OF();
  return this->_newAtom<OriginJumpJoint_O>(atomId,name,comment);
};


Joint_sp JointTree_O::newBondedJoint(const chem::AtomId& atomId, core::T_sp name, const string& comment)
{_OF();
  return this->_newAtom<BondedJoint_O>(atomId,name,comment);
};

Joint_sp JointTree_O::newRootBondedJoint(const chem::AtomId& atomId, core::T_sp name,  const string& comment,
                                       core::Symbol_sp constitutionName,
                                       core::Symbol_sp topologyName,
                                       chem::Plug_sp inPlug)
{_OF();
  RootBondedJoint_sp atom = gc::As<RootBondedJoint_sp>(this->_newAtom<RootBondedJoint_O>(atomId,name,comment));
  atom->setup(constitutionName,topologyName,inPlug);
  return atom;
};

Joint_sp JointTree_O::newDelayedBondedJoint(const chem::AtomId& atomId,core::T_sp name, const string& comment)
{_OF();
  return this->_newAtom<DelayedBondedJoint_O>(atomId,name,comment);
};

void JointTree_O::updateAtomIdMap(const chem::AtomId& atomId, Joint_sp atomHandle )
{_OF();
  this->_AtomMap[atomId] = atomHandle;
}



#if 0
void JointTree_O::replaceMonomerSubTree(const BondId_sp& incoming, const map<core::Symbol_sp,const BondId_sp>& outgoing, const JointTemplate_sp& newSubTreeTemplate )
{_OF();
  newSubTreeTemplate->writeIntoJointTree(this->sharedThis<JointTree_O>(),incoming,outgoing);
}
#endif



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



string JointTree_O::asString() const
{_OF();
  stringstream ss;
  Joint_sp root = this->_Root;
  root->recursiveDumpChildrenIntoStringStream("", ss);
  return ss.str();
}



void JointTree_O::walkTree(core::Function_sp callback)
{_OF();
  OriginJumpJoint_sp originJumpJoint = gc::As<OriginJumpJoint_sp>(this->_Root);
  originJumpJoint->walkChildren(callback);
}

CL_DEFMETHOD void JointTree_O::walk(core::Function_sp exec)
{_OF();
  this->walkTree(exec);
}

CL_DEFMETHOD void JointTree_O::updateInternalCoords()
{_OF();
  ASSERTF(this->_Root.notnilp(),BF("The Root atom is nil - this should never happen"));
  this->_Root->updateInternalCoords(true,this->sharedThis<JointTree_O>());
}


CL_DEFMETHOD void JointTree_O::walkResidue( int residueId, Joint_sp const& root, core::Function_sp exec)
{_OF();
  root->walkResidueTree(residueId,exec);
}



};
