/*
    File: atomTemplate.cc
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
#include <cando/kinematics/atomTemplate.h>
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




#if INIT_TO_FACTORIES


#define ARGS_Checkpoint_O_make "(constitution-name topology-name)"
#define DECL_Checkpoint_O_make ""
#define DOCS_Checkpoint_O_make "make Checkpoint args: constitution-name topology-name"
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

#endif

#ifdef XML_ARCHIVE
    void Checkpoint_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	node->attribute("constitutionName",this->_ConstitutionName);
	node->attribute("topologyName",this->_TopologyName );
	// Archive other instance variables here
    }
#endif



    void Checkpoint_O::setupDelayedBondedAtom(DelayedBondedJoint_sp atom) const
    {_OF();
	chem::CandoDatabase_sp cdb = chem::getCandoDatabase();
        chem::Constitution_sp constitution  = gc::As<chem::Constitution_sp>(core::eval::funcall(chem::_sym_constitutionForNameOrPdb,cdb,this->_ConstitutionName));
	chem::ConstitutionAtoms_sp constitutionAtoms = constitution->getConstitutionAtoms();
        atom->_DelayAtomId = constitutionAtoms->index(this->atomName());
    }




// ----------------------------------------------------------------------
//




#if INIT_TO_FACTORIES

#define ARGS_CheckpointJoint_O_make "(atom-name)"
#define DECL_CheckpointJoint_O_make ""
#define DOCS_CheckpointJoint_O_make "make CheckpointAtom"
CheckpointJoint_sp CheckpointJoint_O::make(core::Symbol_sp atomName)
    {
        GC_ALLOCATE(CheckpointJoint_O, me );
	me->_AtomName = atomName;
	return me;
    };

#else

    core::T_sp CheckpointJoint_O::__init__(core::Function_sp exec, core::Cons_sp args,
					 core::Environment_sp env, core::Lisp_sp lisp)
    {
	this->Base::__init__(exec,args,env,lisp);
	this->_AtomName = translate::from_object<string>::convert(env->lookup(this->Package(),"atomName"));
	ASSERTF(this->_AtomName!="",BF("You must provide atomName"));
	return _Nil<core::T_O>();
    }

#endif

#ifdef XML_ARCHIVE
    void CheckpointJoint_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
	node->attribute("AtomName",this->_AtomName);
    }
#endif



    void CheckpointJoint_O::setupDelayedBondedAtom(DelayedBondedJoint_sp atom) const
    {_OF();
	this->Base::setupDelayedBondedAtom(atom);
	atom->_DelayType = delayForInternalResidueAtom;
    }



// ----------------------------------------------------------------------
//




#if INIT_TO_FACTORIES

#define ARGS_CheckpointOutPlugJoint_O_make "(out-plug)"
#define DECL_CheckpointOutPlugJoint_O_make ""
#define DOCS_CheckpointOutPlugJoint_O_make "make CheckpointOutPlugAtom"
  CheckpointOutPlugJoint_sp CheckpointOutPlugJoint_O::make(const chem::OutPlug_sp& outPlug)
    {
        GC_ALLOCATE(CheckpointOutPlugJoint_O, me );
	me->_Plug = outPlug;
	ASSERTF(me->_Plug.notnilp(),BF("You must provide outPlug argument"));
	return me;
    };

#else

    core::T_sp CheckpointOutPlugJoint_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
	this->Base::__init__(exec,args,env,lisp);
	this->_Plug = translate::from_object<chem::OutPlug_O>::convert(env->lookup(this->Package(),"outPlug"));
	ASSERTF(this->_Plug.notnilp(),BF("You must provide outPlug argument"));
	return _Nil<core::T_O>();
    }

#endif

#ifdef XML_ARCHIVE
    void CheckpointOutPlugJoint_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
	node->archiveObject("OutPlug",this->_Plug);
    }
#endif


    void CheckpointOutPlugJoint_O::setupDelayedBondedAtom(DelayedBondedJoint_sp atom) const
    {_OF();
	this->Base::setupDelayedBondedAtom(atom);
	atom->_DelayType = delayForFollowingResidueBond1;
    }


// ----------------------------------------------------------------------
//




#if INIT_TO_FACTORIES

#define ARGS_AtomTemplate_O_make "(id comment parent)"
#define DECL_AtomTemplate_O_make ""
#define DOCS_AtomTemplate_O_make "make AtomTemplate"
  AtomTemplate_sp AtomTemplate_O::make(const int id, const string& comment, const AtomTemplate_sp parent)
    {
        GC_ALLOCATE(AtomTemplate_O, me );
	me->_Id = id;
	me->_Comment = comment;
	me->_Parent = parent;
	return me;
    };

#else

    core::T_sp AtomTemplate_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
	this->_Id = translate::from_object<int>::convert(env->lookup(this->Package(),"id"));
	this->_Comment = translate::from_object<string>::convert(env->lookup(this->Package(),"comment"));
	this->_Parent = env->lookup(this->Package(),"parent").as<AtomTemplate_O>();
	return _Nil<core::T_O>();
    }

#endif

#ifdef XML_ARCHIVE
    void AtomTemplate_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
	node->archiveWeakPointer("parent",this->_Parent);
	node->attribute("id",this->_Id);
	node->attribute("comment",this->_Comment);
    }
#endif



core::Symbol_sp AtomTemplate_O::atomName(chem::ConstitutionAtoms_sp constitutionAtoms) const
    {
      if ( this->_Id==-1) return _Nil<core::Symbol_O>();
	chem::ConstitutionAtom_sp ca = (*constitutionAtoms)[this->_Id];
	return ca->atomName();
    }
	


gc::Nilable<AtomTemplate_sp> AtomTemplate_O::parent() const
{
  return this->_Parent;
}
	




// ----------------------------------------------------------------------
//




#if INIT_TO_FACTORIES

#define ARGS_BondedAtomTemplate_O_make "(out-plug)"
#define DECL_BondedAtomTemplate_O_make ""
#define DOCS_BondedAtomTemplate_O_make "make BondedAtomTemplate"
  BondedAtomTemplate_sp BondedAtomTemplate_O::make(const chem::OutPlug_sp& outPlug)
  {
      GC_ALLOCATE(BondedAtomTemplate_O, me );
    me->_OutPlug = outPlug;
    return me;
  };

#else

    core::T_sp BondedAtomTemplate_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
	this->Base::__init__(exec,args,env,lisp);
	this->_OutPlug = translate::from_object<chem::OutPlug_sp>::convert(env->lookup(Pkg(),"outPlug"));
	return _Nil<core::T_O>();
    }

#endif

#ifdef XML_ARCHIVE
    void BondedAtomTemplate_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
	node->attribute("dist",this->_Distance);
	node->attribute("theta",this->_Theta);
	node->attribute("phi",this->_Phi);
	node->archiveVector0("children",this->_Children);
	node->attributeIfNotNil("outPlug",this->_OutPlug);
    }
#endif


core::List_sp BondedAtomTemplate_O::children() {
  core::List_sp l = _Nil<core::T_O>();
  for ( int i(0), iEnd(this->_Children.size()); i<iEnd; ++i ) {
    l = core::Cons_O::create(this->_Children[i],l);
  }
  return l;
}


    void BondedAtomTemplate_O::addChildren(Joint_sp me,
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


    Joint_sp BondedAtomTemplate_O::writeIntoAtomTree(const AtomTree_sp& atomTree,
						    uint moleculeId,
						    uint residueId,
						    const BondId_sp& incoming,
						    const PlugNamesToBondIdMap& outgoing,
						    bool rootNode)
    {
	chem::AtomId atomId(moleculeId,residueId,this->_Id);
	Joint_sp ownedBonded = atomTree->newBondedAtom(atomId,this->_Comment);
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




    void BondedAtomTemplate_O::extractInternalCoords(Joint_sp const& atom)
    {
	this->_Distance = atom->dof(DofType::distance);
	this->_Theta = atom->dof(DofType::theta);
	this->_Phi = atom->dof(DofType::phi);
    }




    void BondedAtomTemplate_O::setupOutPlugAtomTree(Joint_sp owned,
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




#if INIT_TO_FACTORIES

#define ARGS_DelayedBondedAtomTemplate_O_make "(checkpoint)"
#define DECL_DelayedBondedAtomTemplate_O_make ""
#define DOCS_DelayedBondedAtomTemplate_O_make "make DelayedBondedAtomTemplate"
    DelayedBondedAtomTemplate_sp DelayedBondedAtomTemplate_O::make(const Checkpoint_sp& checkpoint)
    {
        GC_ALLOCATE(DelayedBondedAtomTemplate_O, me );
	me->_Checkpoint = checkpoint;
	return me;
    };

#else

    core::T_sp DelayedBondedAtomTemplate_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
	this->Base::__init__(exec,args,env,lisp);
	this->_Checkpoint = translate::from_object<Checkpoint_O>::convert(env->lookup(this->Package(),"checkpoint"));
	return _Nil<core::T_O>();
    }

#endif

#ifdef XML_ARCHIVE
    void DelayedBondedAtomTemplate_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
	node->archiveObject("checkpoint",this->_Checkpoint);
    }
#endif



    Joint_sp DelayedBondedAtomTemplate_O::writeIntoAtomTree(const AtomTree_sp& atomTree,
							   uint moleculeId,
							   uint residueId,
							   const BondId_sp& incoming,
							   const PlugNamesToBondIdMap& outgoing,
							   bool rootNode)
    {
	chem::AtomId atomId(moleculeId,residueId,this->_Id);
	Joint_sp ownedBonded = atomTree->newDelayedBondedAtom(atomId,this->_Comment);
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




#if INIT_TO_FACTORIES

#define ARGS_RootBondedAtomTemplate_O_make "(constitution-name topology-name in-plug)"
#define DECL_RootBondedAtomTemplate_O_make ""
#define DOCS_RootBondedAtomTemplate_O_make "make RootBondedAtomTemplate"
  RootBondedAtomTemplate_sp RootBondedAtomTemplate_O::make(const core::Symbol_sp& constitutionName, const core::Symbol_sp topologyName, const chem::Plug_sp& inPlug)
  {
      GC_ALLOCATE(RootBondedAtomTemplate_O, me );
    me->_ConstitutionName = constitutionName;
    ASSERTF(me->_ConstitutionName.notnilp(),BF("You must provide constitutionName"));
    me->_TopologyName = topologyName;
    ASSERTF(me->_TopologyName.notnilp(),BF("You must provide topologyName"));
    me->_InPlug = inPlug;
    ASSERTF(me->_InPlug.notnilp(),BF("You must provide inPlug"));
    return me;
  };

#else

    core::T_sp RootBondedAtomTemplate_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
      this->Base::__init__(exec,args,env,lisp);
      this->_ConstitutionName = translate::from_object<core::Symbol_O>::convert(env->lookup(Pkg(),"constitutionName"));
      ASSERTF(this->_ConstitutionName.notnilp(),BF("You must provide constitutionName"));
      this->_TopologyName = translate::from_object<core::Symbol_O>::convert(env->lookup(Pkg(),"topologyName"));
      ASSERTF(this->_TopologyName.notnilp(),BF("You must provide topologyName"));
      this->_InPlug = translate::from_object<chem::Plug_O>::convert(env->lookup(this->Package(),"inPlug"));
      ASSERTF(this->_InPlug.notnilp(),BF("You must provide inPlug"));
      return _Nil<core::T_O>();
    }

#endif

#ifdef XML_ARCHIVE
    void RootBondedAtomTemplate_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	node->attribute("constitutionName",this->_ConstitutionName);
	node->attribute("topologyName",this->_TopologyName);
	node->archiveObject("plug",this->_InPlug);
    }
#endif



Joint_sp RootBondedAtomTemplate_O::writeIntoAtomTree(const AtomTree_sp& atomTree,
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
	owned = atomTree->newJumpAtom(atomId,this->_Comment);
    } else
    {
	owned = atomTree->newRootBondedAtom(atomId,this->_Comment,
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
