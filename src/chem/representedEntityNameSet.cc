#define	DEBUG_LEVEL_FULL
#include <clasp/core/common.h>
#include <cando/adapt/stringSet.h>
#include <cando/adapt/adapters.h>
#include <clasp/core/environment.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/representedEntityNameSet.h>
#include <cando/chem/atomIndexer.h>
#include <cando/chem/candoDatabase.h>
#include <clasp/core/wrappers.h>


namespace chem
{



    void	RepresentedEntityNameSet_O::initialize()
    {_OF();
	LOG(BF("RepresentedEntityNameSet_O::initialize()") );
	this->Base::initialize();
	this->_Representative = _Nil<core::Symbol_O>();
	LOG(BF("Creating %s") % this->description()  );
    }


//
// Copy Constructor
//
    RepresentedEntityNameSet_O::RepresentedEntityNameSet_O(const RepresentedEntityNameSet_O& ss) : EntityNameSet_O(ss)
    {
//    this->_IsConstitutionGroup = ss._IsConstitutionGroup;
	this->_Representative = ss._Representative;
    }


    string RepresentedEntityNameSet_O::description() const
    {
	stringstream ss;
	ss << "( " << this->className() ;
	ss << " :name " << this->_Name->__repr__();
	ss << " :representative " <<this->_Representative->__repr__();
	ss << " :entityNames " << this->_EntityNames->__repr__() << " )";
	return ss.str();
    }



    EntityNameSetBase_sp RepresentedEntityNameSet_O::copy() const
    {
      GC_COPY(RepresentedEntityNameSet_O, newCopy, *this ); // = RP_Copy<RepresentedEntityNameSet_O>(this->const_sharedThis<RepresentedEntityNameSet_O>());
	return newCopy;
    }

CL_LISPIFY_NAME("getRepresentative");
CL_DEFMETHOD     core::Symbol_sp	RepresentedEntityNameSet_O::getRepresentative()
    {
	return this->_Representative;
    }

CL_LISPIFY_NAME("setRepresentative");
CL_DEFMETHOD     void	RepresentedEntityNameSet_O::setRepresentative(core::Symbol_sp cc )
    {
	this->_Representative = cc;
    }

CL_LISPIFY_NAME("hasRepresentative");
CL_DEFMETHOD     bool	RepresentedEntityNameSet_O::hasRepresentative()
    {
	return this->_Representative.notnilp();
    }

#if 0
    RepresentedEntityNameSet_sp RepresentedEntityNameSet_O::asUnexpandedRepresentedEntityNameSet()
    {
	GC_COPY(RepresentedEntityNameSet_O, newGroup, this->sharedThis<RepresentedEntityNameSet_O>() ); // = RP_Copy<RepresentedEntityNameSet_O>(this->sharedThis<RepresentedEntityNameSet_O>());
	newGroup->_ContainedNames = adapt::StringSet_O::create();
	newGroup->_ContainedNames->insert(this->_Name);
	return newGroup;
    }
#endif



#ifdef XML_ARCHIVE
    void	RepresentedEntityNameSet_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attributeIfNotNil("representative", this->_Representative);
    }
#endif




    RepresentativeList_sp RepresentedEntityNameSet_O::expandedRepresentativeList() const
    {_OF();
	RepresentativeList_sp		objList;
	Entity_sp			obj;
	CandoDatabase_sp		bdb;
	bdb = getCandoDatabase();
	if ( this->_Representative.nilp() )
	{
	    LOG(BF("this->_Representative.nilp -- expanding my EntryNames"));
	    return this->Base::expandedRepresentativeList();
	}
	ASSERTP(bdb->recognizesEntityName(this->_Representative),
		"candoDatabase does not recognize entity name("+this->_Representative->__repr__()+")");
	obj = bdb->getEntity(this->_Representative);
	ASSERTF(obj->isTerminalName(),
		BF("Representatives must always be terminal names, this representative(%s) is NOT")% obj->getName() );
	objList = RepresentativeList_O::create();
	GC_COPY(RepresentedEntityNameSet_O, expanded , *this); // = RP_Copy<RepresentedEntityNameSet_O>(this);
	expanded->expandToTerminalEntityNames();
	LOG(BF("I have a representative[%s] so I'm created a RepresentedEntityNameSet: %s") % expanded->__repr__() );
	objList->vectorPushExtend(expanded);
	return objList;
    }

#if INIT_TO_FACTORIES

#define ARGS_RepresentedEntityNameSet_O_make "(representative_name)"
#define DECL_RepresentedEntityNameSet_O_make ""
#define DOCS_RepresentedEntityNameSet_O_make "make RepresentedEntityNameSet"
  RepresentedEntityNameSet_sp RepresentedEntityNameSet_O::make(core::Symbol_sp representativeName)
  {
      GC_ALLOCATE(RepresentedEntityNameSet_O, me );
    me->_Representative = representativeName;
    return me;
  };

#else

    core::T_sp RepresentedEntityNameSet_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
	this->Base::__init__(exec,args,env,lisp);
	this->_Representative = translate::from_object<core::Symbol_O>::convert(env->lookup(Pkg(),"representativeName"));
	return _Nil<core::T_O>();
    }

#endif







}; // namespace chem





