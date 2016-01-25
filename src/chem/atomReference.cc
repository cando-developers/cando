#define	DEBUG_LEVEL_FULL

#include <string>
#include <clasp/core/common.h>
#include <clasp/core/useBoostRegex.h>
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
    AtomReference_sp AtomReference_O::make(const string& alias)
    {
        GC_ALLOCATE(AtomReference_O, me );
	me->_AtomAlias = alias;
	return me;
    };

#else

    core::T_sp AtomReference_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp bargs, core::Lisp_sp lisp)
    {
	this->_AtomAlias = translate::from_object<string>::convert(bargs->lookup(ChemPkg,"alias"));
	return _Nil<core::T_O>();
    }

#endif


    bool AtomReference_O::equal(core::T_sp o) const
    {
	if ( this->eq(o) ) return true;
	if ( this->_AtomAlias == o.as<AtomReference_O>()->_AtomAlias ) return true;
	return false;
    }

};
