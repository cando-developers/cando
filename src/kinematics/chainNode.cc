#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/chem/monomer.h>
#include <cando/chem/coupling.h>
#include <cando/chem/oligomer.h>
#include <cando/kinematics/monomerNode.h>
#include <cando/kinematics/ringClosingMonomerNode.h>
#include <cando/kinematics/chainNode.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(kinematics,ChainNode_O);

    void ChainNode_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<ChainNode_O>()
//	.def_raw("__init__",&ChainNode_O::__init__,"(self)")
	;
    }

    void ChainNode_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,ChainNode,"","",_lisp)
//	.def_raw("__init__",&"(self)")
	;
#endif
    }


#if INIT_TO_FACTORIES

#define ARGS_ChainNode_O_make "()"
#define DECL_ChainNode_O_make ""
#define DOCS_ChainNode_O_make "make ChainNode"
  ChainNode_sp ChainNode_O::make()
    {_G();
      IMPLEMENT_ME();
    };

#else
    core::T_sp ChainNode_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
//      arg = translate::from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void ChainNode_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
	IMPLEMENT_ME();
    }
#endif


    void ChainNode_O::initialize()
    {_OF();
        this->Base::initialize();
    }


    MonomerNode_sp ChainNode_O::monomerNodeFactory(ChainNode_O* chainNode,
						   RingClosingMonomerMap& ringClosingMonomerMap,
						   chem::Monomer_sp monomer, core::Lisp_sp lisp )
    {_G();
	int seqId = monomer->getSequenceNumber();
	if ( monomer->hasRingClosingOutPlug() )
	{
	    RingClosingMonomerNode_sp monomerNode = RingClosingMonomerNode_O::create(seqId);
	    ringClosingMonomerMap[monomer] = monomerNode;
	    chainNode->_IndexedMonomerNodes[seqId] = monomerNode;
	    return monomerNode;
	}
	MonomerNode_sp monomerNode = MonomerNode_O::create(seqId);
	chainNode->_IndexedMonomerNodes[seqId] = monomerNode;
	return monomerNode;
    }

    void ChainNode_O::makeRingClosingConnections(RingClosingMonomerMap& ringClosings)
    {_OF();
	if ( ringClosings.size() == 0 ) return;
	IMPLEMENT_ME(); // Make connections between RingClosingMonomerNodes
    }




    void ChainNode_O::buildUsingOligomer(chem::Oligomer_sp oligomer, int chainId)
    {_OF();
	chem::Monomer_sp rootMonomer = oligomer->rootMonomer();
	RingClosingMonomerMap	ringClosingMonomerMap;
	MonomerNode_sp rootMonomerNode = this->monomerNodeFactory(this,ringClosingMonomerMap,rootMonomer,_lisp);
	rootMonomerNode->recursivelyBuildChildren(this->sharedThis<ChainNode_O>(),
						  ringClosingMonomerMap,MonomerNode_O::_nil,
						  chem::DirectionalCoupling_O::_nil,rootMonomer);
	this->makeRingClosingConnections(ringClosingMonomerMap);
	this->_RootMonomerNode = rootMonomerNode;
    }

    void ChainNode_O::resizeMonomers(int numMonomers)
    {_OF();
	this->_IndexedMonomerNodes.resize(numMonomers);
    }

    MonomerNode_sp ChainNode_O::lookupMonomerId(int monomerId) const
    {_OF();
	ASSERTF(monomerId < (int)this->_IndexedMonomerNodes.size(),
		BF("Illegal monomerId[%d] - there are only %d monomerNodes")
		% monomerId
		% this->_IndexedMonomerNodes.size());
	return this->_IndexedMonomerNodes[monomerId];
    }
    

}; /* kinematics */
