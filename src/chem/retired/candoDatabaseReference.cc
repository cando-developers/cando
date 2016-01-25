#define	DEBUG_LEVEL_NONE

#include <cando/chem/candoDatabaseReference.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/candoDatabase.h>
#include <clasp/core/objRef.h>
#include <clasp/core/wrappers.h>



namespace chem {




CandoDatabaseReference_sp CandoDatabaseReference_O::create(core::Lisp_sp e,CandoDatabase_sp bdb, const string& link)
{
    GC_ALLOCATE(CandoDatabaseReference_O,ref);
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




