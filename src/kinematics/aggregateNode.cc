#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/kinematics/chainNode.h>
#include <cando/kinematics/aggregateNode.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(kinematics,AggregateNode_O);

    void AggregateNode_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<AggregateNode_O>()
//	.def_raw("__init__",&AggregateNode_O::__init__,"(self)")
	;
    }

    void AggregateNode_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,AggregateNode,"","",_lisp)
//	.def_raw("__init__",&"(self)")
	;
#endif
    }

#ifdef XML_ARCHIVE
    void AggregateNode_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }
#endif


    void AggregateNode_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    void AggregateNode_O::resizeChains(int numChains)
    {_OF();
	int numMolsOld = this->_Chains.size();
	this->_Chains.resize(numChains);
	if ( numChains > numMolsOld )
	{
	    for ( int i=numMolsOld; i<numChains; i++ )
	    {
		this->_Chains[i] = ChainNode_O::create();
	    }
	}
    }

    

}; /* kinematics */
