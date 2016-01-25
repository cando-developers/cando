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
	ss << "HighlightedHolder[" << this->_HeldObject->__repr__() << "]";
	return ss.str();
    }





};




