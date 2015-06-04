       

#define	DEBUG_LEVEL_NONE

#include <cando/chem/ffVdwDb.h>
//#include "core/archiveNode.h"
#include <clasp/core/wrappers.h>



namespace chem {

    REGISTER_CLASS(chem, FFVdw_O );
    REGISTER_CLASS(chem, FFVdwDb_O );

//
#ifdef XML_ARCHIVE
void	FFVdw_O::archiveBase(core::ArchiveP node) {
	node->attribute("type1",this->_T1);
	node->attribute("type2",this->_T2);
	node->attribute("R",this->_R);
	node->attribute("Eps",this->_Eps);
	node->attribute("M",this->_M);
	node->attribute("N",this->_N);
	node->attribute("Edep",this->_Edep);
    };
#endif

#ifdef XML_ARCHIVE
    void	FFVdwDb_O::archiveBase(core::ArchiveP node)
    {
	this->FFBaseDb_O::archiveBase(node);
	node->attributeVector0(KW("Vdw"),this->_Terms);
    };
#endif


};
