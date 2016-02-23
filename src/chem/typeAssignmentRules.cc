/*
    File: typeAssignmentRules.cc
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
