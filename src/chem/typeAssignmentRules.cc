       
#define	DEBUG_LEVEL_FULL

#include <cando/chem/typeAssignmentRules.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/wrappers.h>



namespace chem {


void	OneTypeRule_O::initialize()
{
    this->Base::initialize();
    this->_Match = _Nil<ChemInfo_O>();
}

#ifdef XML_ARCHIVE
    void	OneTypeRule_O::archive(core::ArchiveP node)
{
    node->attribute("assignTypeName",this->_AssignTypeName);
    node->attribute("match",this->_Match);
}
#endif


bool	OneTypeRule_O::testMatch(Atom_sp atom)
{
    return this->_Match->matches(atom);
}




void	TypeAssignmentRules_O::initialize()
{
    this->Base::initialize();
    this->_WildElementDict = _Nil<WildElementDict_O>();
}

#ifdef XML_ARCHIVE
    void	TypeAssignmentRules_O::archive(core::ArchiveP node)
{
    node->attribute("WildElementDict",this->_WildElementDict);
    node->archiveVector0("Rules",this->_Rules);
}
#endif


core::Symbol_sp TypeAssignmentRules_O::calculateType(Atom_sp a)
{
    for ( gctools::Vec0<OneTypeRule_sp>::iterator it=this->_Rules.begin();
    		it!=this->_Rules.end(); it++ )
    {
	if ( (*it)->testMatch(a) )
	{
	    return (*it)->getAssignTypeName();
	}
    }
    return chemkw_intern("UNASSIGNED-TYPE");
}





};


