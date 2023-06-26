#ifdef USE_TOPOLOGY
/*
    File: specificContext.cc
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
#include <cando/chem/specificContext.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/candoDatabase.h>
#include <clasp/core/wrappers.h>


namespace chem {


    SpecificContext_O::SpecificContext_O( const SpecificContext_O& ss ) : SpecificContext_O::Base(ss)
{
    this->_SelfName = ss._SelfName;
    this->_Neighbors = ss._Neighbors;
}


#define ARGS_SpecificContext_O_make "()"
#define DECL_SpecificContext_O_make ""
#define DOCS_SpecificContext_O_make "make SpecificContext"
  SpecificContext_sp SpecificContext_O::make()
  {
    IMPLEMENT_ME();
  };

void	SpecificContext_O::initialize()
{
    this->Base::initialize();
    this->_Neighbors.clear();
    this->_SelfName = nil<core::Symbol_O>();
}


void SpecificContext_O::clear()
{
    this->_Neighbors.clear();
    this->_SelfName = nil<core::Symbol_O>();
}

#ifdef XML_ARCHIVE
    void	SpecificContext_O::archive(core::ArchiveP node)
{
    IMPLEMENT_ME();
}
#endif


CL_LISPIFY_NAME("setSelfName");
CL_DEFMETHOD     void SpecificContext_O::setSelfName(core::Symbol_sp selfName)
{
   this->_SelfName = selfName; 
}

CL_LISPIFY_NAME("pushNeighbor");
CL_DEFMETHOD     void	SpecificContext_O::pushNeighbor(core::Symbol_sp coupling, core::Symbol_sp neighborName)
{
    KeyEntry one;
    one._Coupling = coupling;
    one._Monomer = neighborName;
    this->_Neighbors.push_back(one);
}

    void	SpecificContext_O::setLastNeighborName(core::Symbol_sp neighbor)
{
    ASSERT_gt(this->_Neighbors.size(),0);
    KeyEntry& one = this->_Neighbors.back();
    one._Monomer = neighbor;
}

CL_LISPIFY_NAME("popNeighbor");
CL_DEFMETHOD void	SpecificContext_O::popNeighbor()
{
    this->_Neighbors.pop_back();
}



bool SpecificContext_O::allMonomersInDatabase(CandoDatabase_sp bdb)
{
    if ( !bdb->recognizesMonomerName(this->_SelfName) ) return false;
    for ( gctools::Vec0<KeyEntry>::iterator it = this->_Neighbors.begin();
    		it!=this->_Neighbors.end(); it++ )
    {
	if ( !bdb->recognizesMonomerName(it->_Monomer) ) return false;
    }
    return true;
}

CL_LISPIFY_NAME("asString");
CL_DEFMETHOD string SpecificContext_O::asString()
{
stringstream ss;
    ss << "![";
    ss << this->_SelfName->symbolNameAsString();
    ss << "]";
    for ( gctools::Vec0<KeyEntry>::iterator it = this->_Neighbors.begin();
    		it!=this->_Neighbors.end(); it++ )
    {
	ss << it->_Coupling->symbolNameAsString();
	ss << ".";
	ss << "[";
	ss << it->_Monomer->symbolNameAsString();
	ss << "]";
    }
    return ss.str();
}

    core::Symbol_sp SpecificContext_O::asSymbol()
    {
        string str = this->asString();
        core::Symbol_sp sym = chemkw_intern(str);
        return sym;
    }


CL_DEFMETHOD SpecificContext_sp SpecificContext_O::copySpecificContext()
{
  auto  cpy  = gctools::GC<SpecificContext_O>::copy( *this); // = RP_Copy<SpecificContext_O>(this);
    return cpy;
}





#if INIT_TO_FACTORIES

#define ARGS_SpecificContextSet_O_make "()"
#define DECL_SpecificContextSet_O_make ""
#define DOCS_SpecificContextSet_O_make "make SpecificContextSet"
  SpecificContextSet_sp SpecificContextSet_O::make()
  {
    IMPLEMENT_ME();
    LOG("Creating SpecificContextSet" );
  };

#else

    core::T_sp 	SpecificContextSet_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::LispPtr lisp)
{
    IMPLEMENT_ME();
    LOG("Creating SpecificContextSet" );
    	// your stuff here
}
#endif
void	SpecificContextSet_O::initialize()
{
    this->Base::initialize();
}

#ifdef XML_ARCHIVE
    void	SpecificContextSet_O::archive(core::ArchiveP node)
{
    IMPLEMENT_ME();
}
#endif


    void	SpecificContextSet_O::addWithKey(core::Symbol_sp keyKey, SpecificContext_sp key)
{
    LOG("adding SpecificContext |{}|" , keyKey.c_str()  );
    this->_KeySets.set(keyKey, key );
}


CL_LISPIFY_NAME("addSpecificContext");
CL_DEFMETHOD void	SpecificContextSet_O::add(SpecificContext_sp key)
{
    LOG("adding SpecificContext |{}|" , key->asString().c_str()  );
    this->addWithKey(key->asSymbol(),key);
}


bool	SpecificContextSet_O::contains(SpecificContext_sp key)
{
    return this->_KeySets.count(key->asSymbol())>0;
}

void	SpecificContextSet_O::merge(SpecificContextSet_sp s)
{
    SpecificContextSet_O::iterator it;
    for ( it=s->begin(); it!=s->end(); it++ )
    {
	this->addWithKey(it->first,it->second);
    }
}


void	SpecificContextSet_O::remove(SpecificContextSet_sp s)
{
    SpecificContextSet_O::iterator it;
    for ( it=s->begin(); it!=s->end(); it++ )
    {
	core::Symbol_sp key = it->second->asSymbol();
	if ( this->_KeySets.count(key)>0 )
	{
	    this->_KeySets.remove(key);
	}
    }
}

CL_LISPIFY_NAME("asSymbolSetOfKeys");
CL_DEFMETHOD adapt::SymbolSet_sp SpecificContextSet_O::asSymbolSetOfKeys()
{
    adapt::SymbolSet_sp keys = adapt::SymbolSet_O::create();
    SpecificContextSet_O::iterator it;
    for ( it=this->begin(); it!=this->end(); it++ )
    {
	keys->insert(it->first);
    }
    return keys;
}



CL_LISPIFY_NAME("asCons");
CL_DEFMETHOD core::List_sp SpecificContextSet_O::asCons()
{
    core::Cons_sp first = core::Cons_O::create(nil<core::T_O>(),nil<core::T_O>());
    core::Cons_sp cur = first;
    SpecificContextSet_O::iterator it;
    for ( it=this->begin(); it!=this->end(); it++ )
    {
	core::Cons_sp one = core::Cons_O::create(it->second,nil<core::T_O>());
	cur->setCdr(one);
	cur = one;
    }
    return first->cdr();
}

















};
#endif
