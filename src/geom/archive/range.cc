/*
    File: range.cc
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

#include <clasp/core/common.h>
#include <cando/geom/range.h>
//#include "core/archiveNode.h"
#include <clasp/core/symbolTable.h>
#include <clasp/core/multipleValues.h>
//#include "core/archive.h"
#include <clasp/core/lisp.h>
#include <clasp/core/wrappers.h>


namespace geom {

#if 0
    core::T_sp Range_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp environ, core::LispPtr lisp)
    {
      LOG("Argument types %s %s" , core::cl__class_of(core::oCar(args))->classNameAsString() , args->ocadr()->className() );
	SYMBOL_SC_(GeomPkg,begin);
	SYMBOL_SC_(GeomPkg,end);
	SYMBOL_SC_(GeomPkg,step);
	this->_Begin = environ->lookup(_sym_begin)->object().as<Integer_O>()->as_core::LongLongInt();
	this->_End = environ->lookup(_sym_end)->object().as<Integer_O>()->as_core::LongLongInt();
	this->_Step = environ->lookup(_sym_step)->object().as<Integer_O>()->as_core::LongLongInt();
	return _Nil<core::T_O>();
    }
#endif



    Range_sp Range_O::create(core::LispPtr lisp,core::LongLongInt b, core::LongLongInt e)
{
    LOG("Range_O::create with begin=%d, end=%d" , b , e  );
    auto r = gctools::GC<Range_O>::allocate_with_default_constructor();
    r->_Begin = b;
    r->_End = e;
    r->_Step = 1;
    return r;
};
    Range_sp Range_O::create3(core::LispPtr lisp, core::LongLongInt b, core::LongLongInt e, core::LongLongInt step)
{
    LOG("Range_O::create with begin=%d, end=%d, step=%d" , b , e , step  );
    auto r = gctools::GC<Range_O>::allocate_with_default_constructor();
    r->_Begin = b;
    r->_End = e;
    r->_Step = step;
    return r;
};



    core::Cons_sp Range_O::asCons()
{
    core::Cons_sp dummy = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::Cons_O>());
    core::Cons_sp cur = dummy;
    for (int i=this->_Begin; i<this->_End; i+=this->_Step )
    {
	core::Cons_sp one = core::Cons_O::create(core::Integer_O::create(i),_Nil<core::Cons_O>());
	cur->setCdr(one);
	cur = one;
    }
    return dummy->cdr();
}



string	Range_O::__repr__() const
{
    stringstream ss;
    ss << "Range(";
    ss << this->_Begin;
    ss << ",";
    ss << this->_End;
    ss << ",";
    ss << this->_Step;
    ss << ")";
    return ss.str();
}

void	Range_O::initialize()
{
    this->Base::initialize();
    this->_Begin = 0;
    this->_End = 10;
    this->_Step = 1;
}

#ifdef XML_ARCHIVE
    void	Range_O::archive(core::ArchiveP node)
{
    IMPLEMENT_ME();
}
#endif





};
