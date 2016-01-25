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
#include <clasp/core/binder.h>
#include <cando/chem/candoDatabase.h>
#include <cando/kinematics/atomTemplate.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/monomerNode.h>
#include <cando/kinematics/chainNode.h>
#include <cando/kinematics/atomo.h>
#include <cando/kinematics/atom.h>
#include <cando/kinematics/atomTree.h>

namespace kinematics
{


    string AtomHolder::typeAsString() const
    {
	switch (this->_Type)
	{
	case unused:
	    return "unused";
	case jumpAtom:
	    return "jump";
	case originJumpAtom:
	    return "origin";
	case bondedAtom:
	    return "bonded";
	case delayedBondedAtom:
	    return "delayed";
	default:
	    return "unknown";
	};
    };


    void AtomTreeWalkFunctor::operator()(kinematics::Atom_sp atom) const
    {
	THROW_HARD_ERROR(BF("Subclass of AtomTreeWalkFunctor must define operator()"));
    }

    void ExecutableAtomTreeWalkFunctor::operator()(kinematics::Atom_sp atom) const
    {
      IMPLEMENT_MEF(BF("Improve this"));
#if 0
	core::Cons_sp args = core::Cons_O::create(atom);
	core::ValueFrame_sp frame(core::ValueFrame_O::create(args,_Nil<core::T_O>());
	core::eval::apply(this->_Callback,frame);
//	this->_Callback->evaluate(args,_lisp->nil<core::Environment_O>(),this->lisp());
#endif
    }




    



    void AtomTree_O::initialize()
    {
	this->_AtomHolders.resize(16);
	this->_JumpAtoms.resize(16);
	this->_OriginJumpAtoms.resize(4);
	this->_BondedAtoms.resize(16);
	this->_RootBondedAtoms.resize(16);
	this->_DelayedBondedAtoms.resize(32);
	chem::AtomId id;
	this->_Root = this->newOriginJumpAtom(id,"-root-");
    }


    RefCountedAtomHandle AtomTree_O::lookup(const chem::AtomId& atomId) const
    {
	RefCountedAtomHandle handle = this->_AtomMap[atomId];
	return handle;
    }

    Atom_sp AtomTree_O::atomTreeLookupAtomid(const chem::AtomId& atomId) const
    {
	RefCountedAtomHandle handle = this->lookup(atomId);
	return Atom_O::create(handle);
    }


    void AtomTree_O::resizeMolecules(int numMolecules)
    {_OF();
	LOG(BF("Resizing AtomTree to %d molecules") % numMolecules);
//	int originalNumMolecules = this->_AtomMap.numberOfMolecules();
	this->_AtomMap.resize(numMolecules);
    };

    void AtomTree_O::resizeResidues(int mol, int numResidues)
    {_OF();
	LOG(BF("Resizing AtomTree to %d residues in moleculeId[%d]") % numResidues % mol );
	this->_AtomMap.resize(mol,numResidues);
    };

    void AtomTree_O::resizeAtoms(int mol, int res, int numAtoms)
    {_OF();
	LOG(BF("Resizing AtomTree to %d atoms in residue[%d] moleculeId[%d]")
	    % numAtoms 
	    % res % mol );
	this->_AtomMap.resize(mol,res,numAtoms);
    };




    Atom_sp AtomTree_O::root() const
    {
	return this->_Root;
    }



    Atom_sp AtomTree_O::_initializeNewAtom(uint holderIndex, const chem::AtomId& atomId )
    {_OF();
	AtomTree_sp tree = this->sharedThis<AtomTree_O>();
	Atom_sp atom = Atom_O::create(tree,holderIndex);
	if ( atomId.isDefined() ) this->updateAtomIdMap(atomId,atom->handle());
	atom->handle().get()->__setAtomHandle(atom->handle());
	LOG(BF("AFTER_ALLOC AtomTree:%s") % this->__repr__());
	return atom;
    }




    Atom_sp AtomTree_O::newJumpAtom(const chem::AtomId& atomId, const string& comment)
    {_OF();
	return this->_newAtom<JumpAtom>(this->_JumpAtoms,atomId,comment);
    };

    Atom_sp AtomTree_O::newOriginJumpAtom(const chem::AtomId& atomId,const string& comment)
    {_OF();
	return this->_newAtom<OriginJumpAtom>(this->_OriginJumpAtoms,atomId,comment);
    };


    Atom_sp AtomTree_O::newBondedAtom(const chem::AtomId& atomId, const string& comment)
    {_OF();
	return this->_newAtom<BondedAtom>(this->_BondedAtoms,atomId,comment);
    };

    Atom_sp AtomTree_O::newRootBondedAtom(const chem::AtomId& atomId, const string& comment,
					  core::Symbol_sp constitutionName,
					  core::Symbol_sp topologyName,
					  chem::Plug_sp inPlug)
    {_OF();
	Atom_sp atom = this->_newAtom<RootBondedAtom>(this->_RootBondedAtoms,atomId,comment);
	RootBondedAtom* rootAtom = dynamic_cast<RootBondedAtom*>(atom->get());
	rootAtom->setup(constitutionName,topologyName,inPlug);
	return atom;
    };

    Atom_sp AtomTree_O::newDelayedBondedAtom(const chem::AtomId& atomId,const string& comment)
    {_OF();
	return this->_newAtom<DelayedBondedAtom>(this->_DelayedBondedAtoms,atomId,comment);
    };




    void AtomTree_O::releaseHandle(RefCountedAtomHandle& handle)
    {_OF();
	// Delete all of the children of the node pointed to by this handle
	Atom* atom = handle.get();
	atom->_releaseAllChildren();
	if ( handle.holder()->_Type == bondedAtom )
	{
	    this->_BondedAtoms.release(handle.holder()->_NodeIndex);
	} else if ( handle.holder()->_Type == rootBondedAtom )
	{
	    this->_RootBondedAtoms.release(handle.holder()->_NodeIndex);
	} else if ( handle.holder()->_Type == delayedBondedAtom )
	{
	    this->_DelayedBondedAtoms.release(handle.holder()->_NodeIndex);
	} else if ( handle.holder()->_Type == jumpAtom )
	{
	    this->_JumpAtoms.release(handle.holder()->_NodeIndex);
	} else if ( handle.holder()->_Type == originJumpAtom )
	{
	    this->_OriginJumpAtoms.release(handle.holder()->_NodeIndex);
	} else
	{
	    SIMPLE_ERROR(BF("Unknown handle type"));
	}
	AtomHolder* holder = handle.holder();
	holder->_Type = unused;
	handle.holder()->_NodeIndex = UndefinedUnsignedInt;
	this->_AtomHolders.release(handle.holderIndex());
	LOG(BF("RELEASED %s") % this->__repr__() );
    }



    void AtomTree_O::updateAtomIdMap(const chem::AtomId& atomId, const WeakAtomHandle& atomHandle )
    {_OF();
	try
	{
	    this->_AtomMap[atomId] = atomHandle;
	}
	catch (core::HardError& err)
	{
	    SIMPLE_ERROR(BF("%s") % err.message());
	}
	catch (...)
	{
	    throw;
	}
    }



#if 0
    void AtomTree_O::replaceMonomerSubTree(const BondId_sp& incoming, const map<core::Symbol_sp,const BondId_sp>& outgoing, const AtomTemplate_sp& newSubTreeTemplate )
    {_OF();
	newSubTreeTemplate->writeIntoAtomTree(this->sharedThis<AtomTree_O>(),incoming,outgoing);
    }
#endif



    void AtomTree_O::recursivelyBuildMolecule(int moleculeId,
					      int residueId,
					      chem::CandoDatabase_sp db,
					      MonomerNode_sp monomerNode,
					      Atom_sp parent,
					      bool rootNode)
    {_OF();
	core::List_sp constitutionAndTopology = monomerNode->identifyConstitutionAndTopology(db);
	chem::Constitution_sp constitution = oCar(constitutionAndTopology).as<chem::Constitution_O>();
	chem::Topology_sp topology = oCadr(constitutionAndTopology).as<chem::Topology_O>();
	{_BLOCK_TRACEF(BF("Building constitution[%s] Topology[%s]")
		       % constitution->getName()->__repr__() % topology->getName()->__repr__() );
	    chem::ConstitutionAtoms_sp constitutionAtoms = constitution->getConstitutionAtoms();
	    this->resizeAtoms(moleculeId,residueId,constitutionAtoms->numberOfAtoms());
	    AtomTemplate_sp atomTemplate = topology->properties()->gethash(_lisp->intern(":atomTemplate")).as<AtomTemplate_O>();
	    BondId_sp incoming = BondId_O::create(parent,_Nil<core::T_O>());
	    AtomTemplate_O::PlugNamesToBondIdMap outgoing;
	    //
	    // Write the sub tree described by atomTemplate into the AtomTree
	    // recursively
	    ASSERTF(atomTemplate.notnilp(),BF("The AtomTemplate for Topology[%s] is nil")
		    % topology->getName()->__repr__() );
	    atomTemplate->writeIntoAtomTree(this->sharedThis<AtomTree_O>(),
					    moleculeId,
					    residueId,
					    incoming,
					    outgoing,
					    rootNode);
	    // Now loop over all the OutPlugs and
	    // create their children and connect them to our OutPlug atoms
	    //
	    chem::Topology_O::Plugs::iterator it;
	    for ( it=topology->_Plugs.begin(); it!=topology->_Plugs.end(); it++ )
	    {
		if ( (it->second)->isAssignableTo<chem::OutPlug_O>() )
		{
		    chem::OutPlug_O::smart_ptr outPlug = (it->second).as<chem::OutPlug_O>();
                    core::Symbol_sp atomB0 = outPlug->getB0();
		    int  constitutionBond0AtomId = constitutionAtoms->index(atomB0);
		    chem::AtomId atomId(moleculeId,residueId,constitutionBond0AtomId);
		    RefCountedAtomHandle bond0AtomHandle =
			RefCountedAtomHandle(this->_AtomMap[atomId]);
		    Atom_sp bond0Parent = Atom_O::create(bond0AtomHandle);
		    MonomerNode_sp nextMonomerNode = monomerNode->_Children.get(outPlug->getName());
		    int nextResidueId = nextMonomerNode->_MonomerId;
		    ASSERT(nextResidueId<(int)this->_AtomMap[moleculeId].size());
		    this->recursivelyBuildMolecule(moleculeId,nextResidueId,db,nextMonomerNode,bond0Parent);
		}
	    }
	}	
    }
					



    void AtomTree_O::buildMoleculeUsingChainNode(int moleculeId, ChainNode_sp chainNode, chem::Oligomer_sp oligomer )
    {_OF();
	ASSERTF(moleculeId<this->numberOfMolecules(),BF("Illegal moleculeId[%d]") % moleculeId );
	int numResidues = oligomer->numberOfMonomers();
	this->resizeResidues(moleculeId,numResidues);
	MonomerNode_sp monomerNode = chainNode->_RootMonomerNode;
	chem::CandoDatabase_sp db = chem::_sym_candoDatabase->symbolValue().as<chem::CandoDatabase_O>();
	ASSERTF(this->_Root.notnilp(),BF("The Root of the AtomTree cannot be nil"));
	LOG(BF("Building moleculeId[%d]") % moleculeId);
	this->recursivelyBuildMolecule(moleculeId,
				       monomerNode->_MonomerId,
				       db,
				       monomerNode,
				       this->_Root,
				       true);
    }



    string AtomTree_O::asString() const
    {_OF();
	stringstream ss;
	Atom* root = this->_Root->get();
	root->recursiveDumpChildrenIntoStringStream("", ss);
	return ss.str();
    }



    string AtomTree_O::__repr__() const
    {_OF();
	stringstream ss;
	ss << (BF("Dump of AtomTree")).str() << std::endl;
	ss << (BF("    Unused list start=%d")
	       % this->_AtomHolders.firstUnusedMember() ).str() << std::endl;;
	for ( int i=0; i<this->_AtomHolders.size(); i++ )
	{
	    AtomHolder const & holder = this->_AtomHolders[i];
	    ss << (BF("[%d] = refCount[%d] type[%s] nodeIndex[%d]")
		   % i
		   % holder.refCount()
		   % holder.typeAsString()
		   % holder.nodeIndex() ).str() << std::endl;
	}
	ss << this->_JumpAtoms.dump() << std::endl;
	ss << this->_OriginJumpAtoms.dump() << std::endl;
	ss << this->_BondedAtoms.dump() << std::endl;
	ss << this->_RootBondedAtoms.dump() << std::endl;
	ss << this->_DelayedBondedAtoms.dump() << std::endl;
	return ss.str();
    }


    void AtomTree_O::walkTree(AtomTreeWalkFunctor& functor)
    {_OF();
	OriginJumpAtom* originJumpAtom = dynamic_cast<OriginJumpAtom*>(this->_Root->get());
	originJumpAtom->walkChildren(functor);
    }

    void AtomTree_O::walk(core::Function_sp exec)
    {_OF();
	ExecutableAtomTreeWalkFunctor functor(exec,_lisp);
	this->walkTree(functor);
    }

    void AtomTree_O::updateInternalCoords()
    {_OF();
	ASSERTF(this->_Root.notnilp(),BF("The Root atom is nil - this should never happen"));
	Stub defaultStub;
	this->_Root->get()->updateInternalCoords(defaultStub,true,this->sharedThis<AtomTree_O>());
    }


    void AtomTree_O::walkResidue( int residueId, Atom_sp const& root, core::Function_sp exec)
    {_OF();
	ExecutableAtomTreeWalkFunctor functor(exec,_lisp);
	root->get()->walkResidueTree(residueId,functor);
    }



};
