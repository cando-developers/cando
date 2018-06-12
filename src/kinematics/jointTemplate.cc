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
#define	DEBUG_LEVEL_NONE

#include <clasp/core/common.h>
#include <clasp/core/evaluator.h>
#include <cando/kinematics/jointTemplate.h>
#include <cando/chem/plug.h>
#include <clasp/core/environment.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/constitution.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/kinematics/bondedAtom.h>
#include <cando/kinematics/atomTree.h>
#include <cando/kinematics/delayedBondedAtom.h>
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
        GC_ALLOCATE(Checkpoint_O, me);
	me->_ConstitutionName = constitutionName;
	me->_TopologyName = topologyName;
	ASSERTF(me->_ConstitutionName.notnilp(),BF("You must provide constitutionName"));
	ASSERTF(me->_TopologyName.notnilp(),BF("You must provide topologyName"));
	return me;
    }

    void Checkpoint_O::setupDelayedBondedAtom(DelayedBondedJoint_sp atom) const
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
        GC_ALLOCATE(CheckpointJoint_O, me );
	me->_AtomName = atomName;
	return me;
    };


    void CheckpointJoint_O::setupDelayedBondedAtom(DelayedBondedJoint_sp atom) const
    {_OF();
	this->Base::setupDelayedBondedAtom(atom);
	atom->_DelayType = kw::_sym_delayForInternalResidueAtom;
    }



// ----------------------------------------------------------------------
//



void CheckpointOutPlugJoint_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,out_plug),this->_Plug);
  this->Base::fields(node);
}

  CheckpointOutPlugJoint_sp CheckpointOutPlugJoint_O::make(chem::OutPlug_sp outPlug)
    {
        GC_ALLOCATE(CheckpointOutPlugJoint_O, me );
	me->_Plug = outPlug;
	ASSERTF(me->_Plug.notnilp(),BF("You must provide outPlug argument"));
	return me;
    };

    void CheckpointOutPlugJoint_O::setupDelayedBondedAtom(DelayedBondedJoint_sp atom) const
    {_OF();
	this->Base::setupDelayedBondedAtom(atom);
	atom->_DelayType = kw::_sym_delayForFollowingResidueBond1;
    }


// ----------------------------------------------------------------------
//


void JointTemplate_O::fields(core::Record_sp node) {
  node->field_if_not_nil(INTERN_(kw,parent),this->_Parent);
  node->field(INTERN_(kw,id),this->_Id);
  node->field(INTERN_(kw,name),this->_Name);
  node->field(INTERN_(kw,comment),this->_Comment);
  this->Base::fields(node);
}

JointTemplate_sp JointTemplate_O::make(const int id, core::T_sp name, const string& comment, JointTemplate_sp parent)
    {
        GC_ALLOCATE(JointTemplate_O, me );
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
      if ( this->_Id==-1) return _Nil<core::Symbol_O>();
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
  node->field(INTERN_(kw,out_plug),this->_OutPlug);
  this->Base::fields(node);
}

  BondedJointTemplate_sp BondedJointTemplate_O::make(chem::OutPlug_sp outPlug)
  {
      GC_ALLOCATE(BondedJointTemplate_O, me );
    me->_OutPlug = outPlug;
    return me;
  };



core::List_sp BondedJointTemplate_O::children() const {
  core::List_sp l = _Nil<core::T_O>();
  for ( int i(0), iEnd(this->_Children.size()); i<iEnd; ++i ) {
    l = core::Cons_O::create(this->_Children[i],l);
  }
  return l;
}


    void BondedJointTemplate_O::addChildren(Joint_sp me,
					   uint moleculeId,
					   uint residueId,
					   const AtomTree_sp& atomTree,
					   const BondId_sp& incoming,
					   const PlugNamesToBondIdMap& outgoing)
    {_OF();
	for ( ChildList::iterator it=this->_Children.begin(); it!=this->_Children.end(); it++ )
	{
	    LOG(BF("About to write child[%s]") % (*it)->comment() );
	    Joint_sp newChild = (*it)->writeIntoAtomTree(atomTree,
							moleculeId,
							residueId,
							incoming,
							outgoing);
	    LOG(BF("Child returned with handle[%d]") % newChild.holderIndex() );
	    /* In case the Array was moved by writeIntoAtomTree we get bonded again */
	    Joint_sp bonded = me;
	    ASSERTF(bonded,BF("The dereferenced pointer is NULL - this should not happen"));
	    if ( bonded == newChild)
	    {
		stringstream serr;
		serr << _rep_(atomTree);
		LOG(BF("PROBLEM: %s") % serr.str() );
		SIMPLE_ERROR(BF("You are adding an Atom to itself!"));
	    }
	    bonded->appendChild(newChild);
	}
    }


    Joint_sp BondedJointTemplate_O::writeIntoAtomTree(const AtomTree_sp& atomTree,
						    uint moleculeId,
						    uint residueId,
						    const BondId_sp& incoming,
						    const PlugNamesToBondIdMap& outgoing,
						    bool rootNode)
    {
	chem::AtomId atomId(moleculeId,residueId,this->_Id);
	Joint_sp ownedBonded = atomTree->newBondedAtom(atomId,this->_Name,this->_Comment);
	this->addChildren(ownedBonded,moleculeId,residueId,atomTree,incoming,outgoing);
	if ( this->outPlug().boundp() )
	{
	    this->setupOutPlugAtomTree(ownedBonded,
				       atomTree,
				       moleculeId,
				       residueId,
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




    void BondedJointTemplate_O::setupOutPlugAtomTree(Joint_sp owned,
						    const AtomTree_sp& atomTree,
						    uint moleculeId,
						    uint residueId,
						    const BondId_sp& incoming,
						    const PlugNamesToBondIdMap& outgoing )
    {_OF();
	if ( this->_OutPlug.boundp() )
	{
	    LOG(BF("There is no outplug defined - returning"));
	    return;
	}
	LOG(BF("Setting up OutPlug part of atomTree"));
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
        GC_ALLOCATE(DelayedBondedJointTemplate_O, me );
	me->_Checkpoint = checkpoint;
	return me;
    };


    Joint_sp DelayedBondedJointTemplate_O::writeIntoAtomTree(const AtomTree_sp& atomTree,
							   uint moleculeId,
							   uint residueId,
							   const BondId_sp& incoming,
							   const PlugNamesToBondIdMap& outgoing,
							   bool rootNode)
    {
	chem::AtomId atomId(moleculeId,residueId,this->_Id);
	Joint_sp ownedBonded = atomTree->newDelayedBondedAtom(atomId,this->_Name,this->_Comment);
	this->_Checkpoint->setupDelayedBondedAtom(gc::As<DelayedBondedJoint_sp>(ownedBonded));
	this->addChildren(ownedBonded,moleculeId,residueId,atomTree,incoming,outgoing);
	if ( this->outPlug().boundp() )
	{
	    this->setupOutPlugAtomTree(ownedBonded,
				       atomTree,
				       moleculeId,
				       residueId,
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
      GC_ALLOCATE(RootBondedJointTemplate_O, me );
    me->_ConstitutionName = constitutionName;
    ASSERTF(me->_ConstitutionName.notnilp(),BF("You must provide constitutionName"));
    me->_TopologyName = topologyName;
    ASSERTF(me->_TopologyName.notnilp(),BF("You must provide topologyName"));
    me->_InPlug = inPlug;
    ASSERTF(me->_InPlug.notnilp(),BF("You must provide inPlug"));
    return me;
  };


Joint_sp RootBondedJointTemplate_O::writeIntoAtomTree(const AtomTree_sp& atomTree,
						    uint moleculeId,
						    uint residueId,
						    const BondId_sp& incoming,
						    const PlugNamesToBondIdMap& outgoing,
						    bool rootNode)
{
    LOG(BF("Writing %s[%s] into AtomTree") % this->className() % this->_Comment );
    chem::AtomId atomId(moleculeId,residueId,this->_Id);
    // The type of the Atom to write depends on what the incoming atom is
    // If it isJump then make this a JumpAtom otherwise make this a BondedAtom
    Joint_sp incomingParent = incoming->_Parent;
    if ( incomingParent.unboundp())
    {	
	SIMPLE_ERROR(BF("The incoming parent has to be defined - it is not "));
    }
    Joint_sp owned;
    if (gc::IsA<JumpJoint_sp>(incomingParent) && rootNode )
    {
      owned = atomTree->newJumpAtom(atomId,this->_Name,this->_Comment);
    } else
    {
	owned = atomTree->newRootBondedAtom(atomId,
                                            this->_Name,
                                            this->_Comment,
					    this->_ConstitutionName,
					    this->_TopologyName,
					    this->_InPlug);
    }
    this->addChildren(owned,moleculeId,residueId,atomTree,incoming,outgoing);
    // Attach the new owned atom to the incomingParent as a child
    //
    // QUESTION: If the parent is a Jump atom should I append and otherwise insert????????
    //	I'll insert for now because that will should work for every non root residue
    incomingParent->insertChild(owned);


    // Detach the old Atom that we are replacing - this should cause a flurry of handle releases
    // as the children ref-counts are decremented to zero and they are deleted
    //
    if ( incoming->_Child.boundp() )
    {
	incomingParent->eraseChild(incoming->_Child);
    }
    if ( this->outPlug().boundp() )
    {
	this->setupOutPlugAtomTree(owned,
				   atomTree,
				   moleculeId,
				   residueId,
				   incoming,
				   outgoing);
    }
    return owned;
}







}; /* kinematics */
