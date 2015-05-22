#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "aggregateNode.h"
#include "chem/oligomer.h"
#include "chainNode.h"
#include "monomerNode.h"
#include "foldTree.h"
#include "core/wrappers.h"
namespace kinematics
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(kinematics,FoldTree_O);

    void FoldTree_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<FoldTree_O>()
//	.def_raw("__init__",&FoldTree_O::__init__,"(self)")
	    .def("fold-tree-lookup-monomer-id",&FoldTree_O::lookupMonomerId)
	;
    }

    void FoldTree_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,FoldTree,"","",_lisp)
//	.def_raw("__init__",&"(self)")
	    .def("fold-tree-lookup-monomer-id",&FoldTree_O::lookupMonomerId)
	;
#endif
    }


#if INIT_TO_FACTORIES

#define ARGS_FoldTree_O_make "()"
#define DECL_FoldTree_O_make ""
#define DOCS_FoldTree_O_make "make FoldTree"
  FoldTree_sp FoldTree_O::make()
    {_G();
      IMPLEMENT_ME();
    };

#else
    core::T_sp FoldTree_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
//      arg = translate::from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void FoldTree_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }
#endif


    void FoldTree_O::initialize()
    {_OF();
        this->Base::initialize();
	this->_AggregateNode = AggregateNode_O::create();
    }


    int FoldTree_O::numberOfChains() const
    {
	return this->_AggregateNode->_Chains.size();
    }

    void FoldTree_O::resizeChains(int numChains)
    {
	this->_AggregateNode->resizeChains(numChains);
    }

    void FoldTree_O::resizeMonomers(int chainId, int numMonomers)
    {_OF();
	ASSERTF(chainId<(int)this->_AggregateNode->_Chains.size(),
		BF("Illegal chainId[%d] - there are only %d chains")
		% chainId
		% this->_AggregateNode->_Chains.size());
	this->_AggregateNode->_Chains[chainId]->resizeMonomers(numMonomers);
    }


    ChainNode_sp FoldTree_O::buildChainUsingOligomer(int chainId, chem::Oligomer_sp oligomer)
    {_OF();
	ASSERTF(chainId<this->numberOfChains(),BF("Illegal chainId[%d]") % chainId );
	this->resizeMonomers(chainId,oligomer->numberOfMonomers());
	this->_AggregateNode->_Chains[chainId]->buildUsingOligomer(oligomer,chainId);
	return this->_AggregateNode->_Chains[chainId];
    }
	
    ChainNode_sp FoldTree_O::getChainNode(int chainId) const
    {_OF();
	ASSERTF(chainId>0 && chainId<this->numberOfChains(),
		BF("Illegal chainId[%d]") % chainId );
	return this->_AggregateNode->_Chains[chainId];
    }


    MonomerNode_sp FoldTree_O::lookupMonomerId(MonomerId const& monomerId) const
    {_OF();
	ASSERTF(monomerId.chainId()<(int)this->_AggregateNode->_Chains.size(),
		BF("Illegal chainId[%d] - there are only %d chains") % monomerId.chainId()
		% this->_AggregateNode->_Chains.size() );
	ChainNode_sp chainNode = this->_AggregateNode->_Chains[monomerId.chainId()];
	MonomerNode_sp monomer = chainNode->lookupMonomerId(monomerId.monomerId());
	return monomer;
    }


}; /* kinematics */
