/*
    File: virtualAtom.cc
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

#include <cando/chem/virtualAtom.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/cons.h>
#include <cando/chem/residue.h>
#include <cando/chem/calculatePosition.h>
#include <clasp/core/wrappers.h>




namespace chem {



    VirtualAtom_sp VirtualAtom_O::create(MatterName name, CalculatePosition_sp proc)
{
    LOG("status" );
    auto  va  = gctools::GC<VirtualAtom_O>::allocate_with_default_constructor();
    LOG("status" );
    va->setName(name);
    LOG("status" );
    va->_CalculatePosition = proc;
    LOG("status" );
    return va;
}

void VirtualAtom_O::fields(core::Record_sp node) {
  node->field_if_not_unbound(INTERN_(kw,calculate_position),this->_CalculatePosition);
  this->Base::fields(node);
}


//
// Copy constructor
//
VirtualAtom_O::VirtualAtom_O(const VirtualAtom_O& ss) :Atom_O(ss)
{
    this->_CalculatePosition = ss._CalculatePosition;
//    this->_Environment = ss._Environment;
}

//
//	copy
//
//	Allocate a copy of this atom and return it.
//	Keep track of the new atom so that we can redirect pointers to the copy.
//
Matter_sp VirtualAtom_O::copy()
{
    LOG("Copying atom @%X" , this );
    auto aNew = gctools::GC<VirtualAtom_O>::copy( *this); // VirtualAtom_sp aNew = RP_Copy<VirtualAtom_O>(this);
    this->copyAtom = aNew;
    LOG("    copy atom== %s" , aNew->description().c_str() );
    return(aNew);
}



void	VirtualAtom_O::initialize()
{
    this->Base::initialize();
    this->_CalculatePosition = nil<CalculatePosition_O>();
}

#if 0
CL_DEFMETHOD bool VirtualAtom_O::VirtualAtom_equal(core::T_sp obj) const
{
    if ( this->eq(obj) ) return true;
    if ( !obj.isA<VirtualAtom_O>() )
    {
	VirtualAtom_sp other = obj.as<VirtualAtom_O>();
	if ( other->getName() != this->getName() ) return false;
	if ( other->numberOfBonds() != this->numberOfBonds() ) return false;
	return true;
    }
    return false;
}
#endif



#ifdef XML_ARCHIVE
    void	VirtualAtom_O::archive(core::ArchiveP node)
{
    this->Base::archive(node);
    node->attributeIfNotNil("calculatePosition",this->_CalculatePosition);
//    node->attribute("environment",this->_Environment);
}
#endif



void	VirtualAtom_O::calculatePosition(Residue_sp res)
{
    Vector3 vpos = this->_CalculatePosition->calculatePosition(res);
    this->setPosition(vpos);
}










};
