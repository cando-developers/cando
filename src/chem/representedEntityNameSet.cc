/*
    File: representedEntityNameSet.cc
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
	this->_Representative = nil<core::Symbol_O>();
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
	ss << " :name " << _rep_(this->_Name);
	ss << " :representative " <<_rep_(this->_Representative);
	ss << " :entityNames " << _rep_(this->_EntityNames) << " )";
	return ss.str();
    }



    EntityNameSetBase_sp RepresentedEntityNameSet_O::copy() const
    {
      auto  newCopy = gctools::GC<RepresentedEntityNameSet_O>::copy( *this ); // = RP_Copy<RepresentedEntityNameSet_O>(this->const_sharedThis<RepresentedEntityNameSet_O>());
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
	bdb = gc::As<CandoDatabase_sp>(getCandoDatabase());
	if ( this->_Representative.nilp() )
	{
	    LOG(BF("this->_Representative.nilp -- expanding my EntryNames"));
	    return this->Base::expandedRepresentativeList();
	}
	ASSERT(bdb->recognizesEntityName(this->_Representative)); // "candoDatabase does not recognize entity name("+_rep_(this->_Representative)+")");
	obj = bdb->getEntity(this->_Representative);
	ASSERT(obj->isTerminalName()); // BF("Representatives must always be terminal names, this representative(%s) is NOT")% obj->getName() );
	objList = RepresentativeList_O::create();
	auto  expanded  = gctools::GC<RepresentedEntityNameSet_O>::copy( *this); // = RP_Copy<RepresentedEntityNameSet_O>(this);
	expanded->expandToTerminalEntityNames();
	LOG(BF("I have a representative[%s] so I'm created a RepresentedEntityNameSet: %s") % _rep_(expanded) );
	objList->vectorPushExtend(expanded);
	return objList;
    }

#if INIT_TO_FACTORIES

#define ARGS_RepresentedEntityNameSet_O_make "(representative_name)"
#define DECL_RepresentedEntityNameSet_O_make ""
#define DOCS_RepresentedEntityNameSet_O_make "make RepresentedEntityNameSet"
  RepresentedEntityNameSet_sp RepresentedEntityNameSet_O::make(core::Symbol_sp representativeName)
  {
    auto  me  = gctools::GC<RepresentedEntityNameSet_O>::allocate_with_default_constructor();
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
