/*
    File: candoDatabaseReference.cc
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

#include <cando/chem/candoDatabaseReference.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/candoDatabase.h>
#include <clasp/core/objRef.h>
#include <clasp/core/wrappers.h>



namespace chem {




CandoDatabaseReference_sp CandoDatabaseReference_O::create(core::LispPtr e,CandoDatabase_sp bdb, const string& link)
{
  auto ref = gctools::GC<CandoDatabaseReference_O>::allocate_with_default_constructor();
    ASSERTNOTNULLP(bdb,"The alchemist database is NULL");
    ref->_Reference = core::ObjRef_O::create(e,link);
    return ref;
}


void	CandoDatabaseReference_O::initialize()
{
    this->Base::initialize();
    this->_Reference = _Nil<core::ObjRef_O>();
}

#ifdef XML_ARCHIVE
    void	CandoDatabaseReference_O::archive(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    if ( node->loading() )
    {
        string refStr;
	node->attribute<string>("ref",refStr);
	this->_Reference = core::ObjRef_O::create(_liO,,refStr);
    } else
    {
        string refStr;
	refStr );
	node->attribute<string>("ref",refStr);
    }
}
#endif


    core::T_sp CandoDatabaseReference_O::oget()
{
    core::T_sp	obj;
    obj = getCandoDatabase();
    ASSERTNOTNULL(obj);
    core::T_sp res = this->_Reference->relativeTo(obj);
    ASSERTNOTNULL(res);
    return res;
}






};
