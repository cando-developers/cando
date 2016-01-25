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
