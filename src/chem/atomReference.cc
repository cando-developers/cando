/*
    File: atomReference.cc
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

#include <string>
#include <clasp/core/common.h>
#include <clasp/core/fileSystem.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/environment.h>
#include <cando/geom/vector3.h>
#include <cando/chem/atomReference.h>
#include <clasp/core/wrappers.h>



namespace chem 
{





    void AtomReferenceBase_O::initialize()
    {
	this->Base::initialize();
    }

#ifdef XML_ARCHIVE
    void AtomReferenceBase_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
    }
#endif






    void AtomReference_O::initialize()
    {
	this->Base::initialize();
    }

#ifdef XML_ARCHIVE
    void AtomReference_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attribute("atom",this->_AtomAlias);
    }
#endif

#if INIT_TO_FACTORIES

#define ARGS_AtomReference_O_make "(alias)"
#define DECL_AtomReference_O_make ""
#define DOCS_AtomReference_O_make "make AtomReference args: alias"
AtomReference_sp AtomReference_O::make(core::T_sp alias)
{
  auto  me  = gctools::GC<AtomReference_O>::allocate_with_default_constructor();
  me->_AtomAlias = alias;
  return me;
};

#else

    core::T_sp AtomReference_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp bargs, core::LispPtr lisp)
    {
	this->_AtomAlias = translate::from_object<string>::convert(bargs->lookup(ChemPkg,"alias"));
	return _Nil<core::T_O>();
    }

#endif


#if 0
bool AtomReference_O::equal(core::T_sp o) const
    {
	if ( this->eq(o) ) return true;
	if ( this->_AtomAlias == o.as<AtomReference_O>()->_AtomAlias ) return true;
	return false;
    }
#endif
};
