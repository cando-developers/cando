/*
    File: jointTemplate.cc
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
#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/evaluator.h>
#include <cando/kinematics/jointTemplate.h>
#include <cando/chem/plug.h>
#include <clasp/core/environment.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/constitution.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/kinematics/bondedJoint.h>
#include <cando/kinematics/jointTree.h>
#include <cando/kinematics/monomerNode.h>
#include <cando/kinematics/delayedBondedJoint.h>
#include <clasp/core/wrappers.h>



namespace kinematics
{

// ----------------------------------------------------------------------
//


void Checkpoint_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,constitution_name),this->_ConstitutionName);
  node->field(INTERN_(kw,topology_name),this->_TopologyName);
}

Checkpoint_sp Checkpoint_O::make(const core::Symbol_sp& constitutionName,
                                 const core::Symbol_sp& topologyName )
{
  auto  me = gctools::GC<Checkpoint_O>::allocate_with_default_constructor();
  me->_ConstitutionName = constitutionName;
  me->_TopologyName = topologyName;
  ASSERTF(me->_ConstitutionName.notnilp(),BF("You must provide constitutionName"));
  ASSERTF(me->_TopologyName.notnilp(),BF("You must provide topologyName"));
  return me;
}

void Checkpoint_O::setupDelayedBondedJoint(DelayedBondedJoint_sp atom) const
{_OF();
  chem::CandoDatabase_sp cdb = chem::getCandoDatabase();
  chem::Constitution_sp constitution  = gc::As<chem::Constitution_sp>(core::eval::funcall(chem::_sym_constitutionForNameOrPdb,cdb,this->_ConstitutionName));
  chem::ConstitutionAtoms_sp constitutionAtoms = constitution->getConstitutionAtoms();
  atom->_DelayAtomId = constitutionAtoms->index(this->atomName());
}




// ----------------------------------------------------------------------
//

void CheckpointJoint_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,atom_name),this->_AtomName);
  this->Base::fields(node);
}


CheckpointJoint_sp CheckpointJoint_O::make(core::Symbol_sp atomName)
{
  auto  me  = gctools::GC<CheckpointJoint_O>::allocate_with_default_constructor();
  me->_AtomName = atomName;
  return me;
};


void CheckpointJoint_O::setupDelayedBondedJoint(DelayedBondedJoint_sp atom) const
{_OF();
  this->Base::setupDelayedBondedJoint(atom);
  atom->_DelayType = kw::_sym_delayForInternalResidueAtom;
}



// ----------------------------------------------------------------------
//



void CheckpointOutPlugJoint_O::fields(core::Record_sp node) {
  node->field_if_not_unbound(INTERN_(kw,out_plug),this->_Plug);
  this->Base::fields(node);
}

CheckpointOutPlugJoint_sp CheckpointOutPlugJoint_O::make(chem::OutPlug_sp outPlug)
{
  auto  me  = gctools::GC<CheckpointOutPlugJoint_O>::allocate_with_default_constructor();
  me->_Plug = outPlug;
  ASSERTF(me->_Plug.notnilp(),BF("You must provide outPlug argument"));
  return me;
};

void CheckpointOutPlugJoint_O::setupDelayedBondedJoint(DelayedBondedJoint_sp atom) const
{_OF();
  this->Base::setupDelayedBondedJoint(atom);
  atom->_DelayType = kw::_sym_delayForFollowingResidueBond1;
}


// ----------------------------------------------------------------------
//


void JointTemplate_O::fields(core::Record_sp node) {
  node->field_if_not_nil(INTERN_(kw,parent),this->_Parent);
  node->field(INTERN_(kw,id),this->_Id);
  node->field(INTERN_(kw,name),this->_Name);
  node->field(INTERN_(kw,comment),this->_Comment);
}

CL_DEFMETHOD size_t JointTemplate_O::id() const {
  return this->_Id;
};

string JointTemplate_O::__repr__() const
{
  stringstream ss;
  ss << "#<" << this->className() << " " << _rep_(this->_Name) << ">";
  return ss.str();
}

JointTemplate_sp JointTemplate_O::make(const int id, core::T_sp name, const string& comment, JointTemplate_sp parent)
{
  auto  me  = gctools::GC<JointTemplate_O>::allocate_with_default_constructor();
  me->_Id = id;
  me->_Name = name;
  me->_Comment = comment;
  me->_Parent = parent;
  return me;
};

CL_DOCSTRING(R"(Return the atom name for this joint-template.
This returns the _Name field of the joint-template.)");
CL_DEFMETHOD core::T_sp JointTemplate_O::name() const
{
  return this->_Name;
}

CL_DOCSTRING(R"(Return the atom name for this joint-template. 
You need to pass the constitution-atoms that this joint-template's atom id indexes in to 
because the joint-template doesn't store the atom name, just an index into a constitution-atoms object.)");
CL_LISPIFY_NAME("joint-template-atom-name");
CL_DEFMETHOD core::Symbol_sp JointTemplate_O::jointTemplateAtomName(chem::ConstitutionAtoms_sp constitutionAtoms) const
{
  if ( this->_Id==-1) return nil<core::Symbol_O>();
  chem::ConstitutionAtom_sp ca = (*constitutionAtoms)[this->_Id];
  return ca->atomName();
}
	


gc::Nilable<JointTemplate_sp> JointTemplate_O::parent() const
{
  return this->_Parent;
}
	




// ----------------------------------------------------------------------
//


void BondedJointTemplate_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,children),this->_Children);
  node->field_if_not_default(INTERN_(kw,distance),this->_Distance,0.0);
  node->field_if_not_default(INTERN_(kw,theta),this->_Theta,0.0);
  node->field_if_not_default(INTERN_(kw,phi),this->_Phi,0.0);
  node->field_if_not_unbound(INTERN_(kw,out_plug),this->_OutPlug);
  this->Base::fields(node);
}

BondedJointTemplate_sp BondedJointTemplate_O::make(chem::OutPlug_sp outPlug)
{
  auto  me  = gctools::GC<BondedJointTemplate_O>::allocate_with_default_constructor();
  me->_OutPlug = outPlug;
  return me;
};



core::List_sp BondedJointTemplate_O::children() const {
  core::List_sp l = nil<core::T_O>();
  for ( int i(0), iEnd(this->_Children.size()); i<iEnd; ++i ) {
    l = core::Cons_O::create(this->_Children[i],l);
  }
  return l;
}


void BondedJointTemplate_O::addChildren(Joint_sp me,
                                        MonomerId monomerId,
                                        const JointTree_sp& JointTree,
                                        const BondId_sp& incoming,
                                        const PlugNamesToBondIdMap& outgoing,
                                        MonomerNode_sp monomerNode)
{_OF();
  for ( ChildList::iterator it=this->_Children.begin(); it!=this->_Children.end(); it++ )
  {
    LOG(BF("About to write child[%s]") % (*it)->comment() );
    Joint_sp newChild = (*it)->writeIntoJointTree(JointTree,
                                                 monomerId,
                                                 incoming,
                                                 outgoing,
                                                 monomerNode);
    LOG(BF("Child returned %s") % _rep_(newChild));
	    /* In case the Array was moved by writeIntoJointTree we get bonded again */
    Joint_sp bonded = me;
    ASSERTF(bonded,BF("The dereferenced Jointer is NULL - this should not happen"));
    if ( bonded == newChild)
    {
      stringstream serr;
      serr << _rep_(JointTree);
      LOG(BF("PROBLEM: %s") % serr.str() );
      SIMPLE_ERROR(BF("You are adding an Atom to itself!"));
    }
    bonded->appendChild(newChild);
  }
}


Joint_sp BondedJointTemplate_O::writeIntoJointTree(const JointTree_sp& JointTree,
                                                  MonomerId monomerId,
                                                  const BondId_sp& incoming,
                                                  const PlugNamesToBondIdMap& outgoing,
                                                  MonomerNode_sp monomerNode,
                                                  bool rootNode)
{_OF();
  LOG(BF("BondedJointTemplate_O::writeIntoJointTree this: %s  BondId_sp: %s\n") % _rep_(this->asSmartPtr()) % _rep_(incoming));
  chem::AtomId atomId(monomerId._Chain,monomerId._Monomer,this->_Id);
  Joint_sp ownedBonded = JointTree->newBondedJoint(atomId,this->_Name,this->_Comment);
  monomerNode->addJoint(this->_Id,ownedBonded);
  this->addChildren(ownedBonded,monomerId,JointTree,incoming,outgoing,monomerNode);
  if ( this->outPlug().boundp() )
  {
    this->setupOutPlugJointTree(ownedBonded,
                               JointTree,
                               incoming,
                               outgoing);
  }
  return ownedBonded;
}


void BondedJointTemplate_O::extractInternalCoords(Joint_sp const& atom)
{
  this->_Distance = atom->dof(DofType::distance);
  this->_Theta = atom->dof(DofType::theta);
  this->_Phi = atom->dof(DofType::phi);
}


void BondedJointTemplate_O::setupOutPlugJointTree(Joint_sp owned,
                                                 const JointTree_sp& JointTree,
                                                 const BondId_sp& incoming,
                                                 const PlugNamesToBondIdMap& outgoing )
{_OF();
  if ( this->_OutPlug.boundp() )
  {
    LOG(BF("There is no outplug defined - returning"));
    return;
  }
  LOG(BF("Setting up OutPlug part of JointTree"));
  core::Symbol_sp plugName = this->_OutPlug->getName();
  PlugNamesToBondIdMap::const_iterator findOutPlug = outgoing.find(plugName);
  if ( findOutPlug != outgoing.end() )
  {
    BondId_sp bondId = findOutPlug->second;
    if ( bondId->_Child.boundp() )
    {
		// Attach the next residues root atom to this Bond0 atom
      owned->insertChild(bondId->_Child);
      if ( bondId->_Parent.boundp() )
      {
        bondId->_Parent->eraseChild(bondId->_Child);
      }
    }
  }
}

// ----------------------------------------------------------------------
//

void DelayedBondedJointTemplate_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,checkpoint),this->_Checkpoint);
  this->Base::fields(node);
}

DelayedBondedJointTemplate_sp DelayedBondedJointTemplate_O::make(const Checkpoint_sp& checkpoint)
{
  auto  me  = gctools::GC<DelayedBondedJointTemplate_O>::allocate_with_default_constructor();
  me->_Checkpoint = checkpoint;
  return me;
};


Joint_sp DelayedBondedJointTemplate_O::writeIntoJointTree(const JointTree_sp& JointTree,
                                                         MonomerId monomerId,
                                                         const BondId_sp& incoming,
                                                         const PlugNamesToBondIdMap& outgoing,
                                                         MonomerNode_sp monomerNode,
                                                         bool rootNode)
{
  LOG(BF("BondedJointTemplate_O::writeIntoJointTree this: %s  BondId_sp: %s\n") % _rep_(this->asSmartPtr()) % _rep_(incoming));
  chem::AtomId atomId(monomerId._Chain,monomerId._Monomer,this->_Id);
  Joint_sp ownedBonded = JointTree->newDelayedBondedJoint(atomId,this->_Name,this->_Comment);
  monomerNode->addJoint(this->_Id,ownedBonded);
  this->_Checkpoint->setupDelayedBondedJoint(gc::As<DelayedBondedJoint_sp>(ownedBonded));
  this->addChildren(ownedBonded,monomerId,JointTree,incoming,outgoing,monomerNode);
  if ( this->outPlug().boundp() ) {
    this->setupOutPlugJointTree(ownedBonded,
                               JointTree,
                               incoming,
                               outgoing);
  }
  return ownedBonded;
}





// ----------------------------------------------------------------------
//


void RootBondedJointTemplate_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,constitution_name),this->_ConstitutionName);
  node->field(INTERN_(kw,topology_name),this->_TopologyName);
  node->field(INTERN_(kw,in_plug),this->_InPlug);
  this->Base::fields(node);
}

RootBondedJointTemplate_sp RootBondedJointTemplate_O::make(core::Symbol_sp constitutionName, const core::Symbol_sp topologyName, chem::Plug_sp inPlug)
{
  auto  me  = gctools::GC<RootBondedJointTemplate_O>::allocate_with_default_constructor();
  me->_ConstitutionName = constitutionName;
  ASSERTF(me->_ConstitutionName.notnilp(),BF("You must provide constitutionName"));
  me->_TopologyName = topologyName;
  ASSERTF(me->_TopologyName.notnilp(),BF("You must provide topologyName"));
  me->_InPlug = inPlug;
  ASSERTF(me->_InPlug.notnilp(),BF("You must provide inPlug"));
  return me;
};


Joint_sp RootBondedJointTemplate_O::writeIntoJointTree(const JointTree_sp& JointTree,
                                                      MonomerId monomerId,
                                                      const BondId_sp& incoming,
                                                      const PlugNamesToBondIdMap& outgoing,
                                                      MonomerNode_sp monomerNode,
                                                      bool rootNode)
{_OF();
  LOG(BF("BondedJointTemplate_O::writeIntoJointTree this: %s  BondId_sp: %s\n") % _rep_(this->asSmartPtr()) % _rep_(incoming));
  LOG(BF("Writing %s[%s] into JointTree") % this->className() % this->_Comment );
  chem::AtomId atomId(monomerId._Chain, monomerId._Monomer,this->_Id);
    // The type of the Atom to write depends on what the incoming atom is
    // If it isJump then make this a JumpJoint otherwise make this a BondedJoint
  Joint_sp incomingParent = incoming->_Parent;
  if ( incomingParent.unboundp())
  {	
    SIMPLE_ERROR(BF("The incoming parent has to be defined - it is not "));
  }
  Joint_sp owned;
  if (gc::IsA<JumpJoint_sp>(incomingParent) && rootNode )
  {
    owned = JointTree->newJumpJoint(atomId,this->_Name,this->_Comment);
  } else
  {
    owned = JointTree->newRootBondedJoint(atomId,
                                        this->_Name,
                                        this->_Comment,
                                        this->_ConstitutionName,
                                        this->_TopologyName,
                                        this->_InPlug);
  }
  monomerNode->addJoint(this->_Id,owned);
  this->addChildren(owned,monomerId,JointTree,incoming,outgoing,monomerNode);
    // Attach the new owned atom to the incomingParent as a child
    //
    // QUESTION: If the parent is a Jump atom should I append and otherwise insert????????
    //	I'll insert for now because that should work for every non root residue
//  incomingParent->insertChild(owned);
  // It turned out that inserting was the wrong thing to do - try appending
  incomingParent->appendChild(owned);

    // Detach the old Atom that we are replacing
    //
  if ( incoming->_Child.boundp() )
  {
    incomingParent->eraseChild(incoming->_Child);
  }
  if ( this->outPlug().boundp() )
  {
    this->setupOutPlugJointTree(owned,
                               JointTree,
                               incoming,
                               outgoing);
  }
  return owned;
}







}; /* kinematics */
