/*
    File: highlightObject.cc
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

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/symbolTable.h>
#include <cando/geom/highlightObject.h>
#include <clasp/core/multipleValues.h>
#include <clasp/core/executables.h>
#include <clasp/core/environment.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/lisp.h>

// last include is wrappers.h
#include <clasp/core/wrappers.h>


namespace geom {





    void PickableObject_O::initialize()
    {
	this->Base::initialize();
    }
#ifdef XML_ARCHIVE
    void PickableObject_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
    }
#endif





    void ControlObject_O::initialize()
    {
	this->Base::initialize();
    }

#ifdef XML_ARCHIVE
    void ControlObject_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
    }
#endif










    void DragObject_O::initialize()
    {
	this->Base::initialize();
    }

#ifdef XML_ARCHIVE
    void DragObject_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
    }
#endif

    DragObject_sp DragObject_O::make(core::Function_sp start, core::Function_sp drag, core::Function_sp release)
    {
	DragObject_sp me(DragObject_O::create());
	me->_StartCode = start;
	me->_DragCode = drag;
	me->_ReleaseCode = release;
	return me;
    }













    HighlightedObject_O::HighlightedObject_O() : _HighlightedId(0) {};

    HighlightedObject_O::~HighlightedObject_O() {};



#if 0 // old way to initialize
    void	HighlightedObject_O::lispInitialize(core::Cons_sp kargs, core::Lisp_sp env)
{
    this->Base::lispInitialize(kargs,env);
    	// your stuff here
}
#endif
#
#
void	HighlightedObject_O::initialize()
{
    this->Base::initialize();
}

#ifdef XML_ARCHIVE
    void	HighlightedObject_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
}
#endif



    void HighlightedObject_O::setHighlightedId(uint id)
    {
	this->_HighlightedId = id;
    }

    bool HighlightedObject_O::isPicked(uint id)
    {
	return this->_HighlightedId == id;
    }







#if 0
#define ARGS_af_makeHighlightedHolder "(obj)"
#define DECL_af_makeHighlightedHolder ""
#define DOCS_af_makeHighlightedHolder "make HighlightedHolder args: path"
    HighlightedHolder_sp af_makeHighlightedHolder(core::T_sp obj)
    {
	HighlightedHolder_sp me(HighlightedHolder_O::create());
	GlueEnvironment_sp env(GlueEnvironment_O::create()<< _sym_obj <<  obj   ).cons()));
	me->make_init__(core::Function_O::_nil,env->args(),env,_lisp);
	return me;
    }


    T_sp HighlightedHolder_O::make_init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp environ, core::Lisp_sp lisp)
{_OF();
    SYMBOL_SC_(GeomPkg,obj);
    core::T_sp obj = translate::from_object<T_sp>::convert(environ->lookup(_sym_obj));
    this->_HeldObject = obj;
    return _Nil<core::T_O>();
}
#endif

void HighlightedHolder_O::initialize()
{
    this->Base::initialize();
    this->_HeldObject = _Nil<T_O>();
}

#ifdef XML_ARCHIVE
void HighlightedHolder_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->archiveObject("obj",this->_HeldObject);
}
#endif


    string HighlightedHolder_O::__repr__() const
    {
	stringstream ss;
	ss << "HighlightedHolder[" << _rep_(this->_HeldObject) << "]";
	return ss.str();
    }





};
