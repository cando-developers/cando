#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/stringSet.h"
#include "core/symbolSet.h"
#include "specificContext.h"
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include "candoDatabase.h"
#include "core/wrappers.h"


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
  {_G();
    IMPLEMENT_ME();
  };

void	SpecificContext_O::initialize()
{_OF();
    this->Base::initialize();
    this->_Neighbors.clear();
    this->_SelfName = _Nil<core::Symbol_O>();
}


void SpecificContext_O::clear()
{_OF();
    this->_Neighbors.clear();
    this->_SelfName = _Nil<core::Symbol_O>();
}

#ifdef XML_ARCHIVE
    void	SpecificContext_O::archive(core::ArchiveP node)
{
    IMPLEMENT_ME();
}
#endif


    void SpecificContext_O::setSelfName(core::Symbol_sp selfName)
{_G();
   this->_SelfName = selfName; 
}

    void	SpecificContext_O::pushNeighbor(core::Symbol_sp coupling, core::Symbol_sp neighborName)
{_G();
    KeyEntry one;
    one._Coupling = coupling;
    one._Monomer = neighborName;
    this->_Neighbors.push_back(one);
}

    void	SpecificContext_O::setLastNeighborName(core::Symbol_sp neighbor)
{_G();
    ASSERT_gt(this->_Neighbors.size(),0);
    KeyEntry& one = this->_Neighbors.back();
    one._Monomer = neighbor;
}

void	SpecificContext_O::popNeighbor()
{_G();
    this->_Neighbors.pop_back();
}



bool SpecificContext_O::allMonomersInDatabase(CandoDatabase_sp bdb)
{_G();
    if ( !bdb->recognizesMonomerName(this->_SelfName) ) return false;
    for ( gctools::Vec0<KeyEntry>::iterator it = this->_Neighbors.begin();
    		it!=this->_Neighbors.end(); it++ )
    {
	if ( !bdb->recognizesMonomerName(it->_Monomer) ) return false;
    }
    return true;
}

string SpecificContext_O::asString()
{_G();
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


SpecificContext_sp SpecificContext_O::copy()
{_G();
    GC_COPY(SpecificContext_O, cpy , *this); // = RP_Copy<SpecificContext_O>(this);
    return cpy;
}





#if INIT_TO_FACTORIES

#define ARGS_SpecificContextSet_O_make "()"
#define DECL_SpecificContextSet_O_make ""
#define DOCS_SpecificContextSet_O_make "make SpecificContextSet"
  SpecificContextSet_sp SpecificContextSet_O::make()
  {_G();
    IMPLEMENT_ME();
    LOG(BF("Creating SpecificContextSet") );
  };

#else

    core::T_sp 	SpecificContextSet_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
{_OF();
    IMPLEMENT_ME();
    LOG(BF("Creating SpecificContextSet") );
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
{_G();
    LOG(BF("adding SpecificContext |%s|") % keyKey.c_str()  );
    this->_KeySets.set(keyKey, key );
}


void	SpecificContextSet_O::add(SpecificContext_sp key)
{_G();
    LOG(BF("adding SpecificContext |%s|") % key->asString().c_str()  );
    this->addWithKey(key->asSymbol(),key);
}


bool	SpecificContextSet_O::contains(SpecificContext_sp key)
{_G();
    return this->_KeySets.count(key->asSymbol())>0;
}

void	SpecificContextSet_O::merge(SpecificContextSet_sp s)
{_G();
    SpecificContextSet_O::iterator it;
    for ( it=s->begin(); it!=s->end(); it++ )
    {
	this->addWithKey(it->first,it->second);
    }
}


void	SpecificContextSet_O::remove(SpecificContextSet_sp s)
{_G();
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

core::SymbolSet_sp SpecificContextSet_O::asSymbolSetOfKeys()
{_G();
    core::SymbolSet_sp keys = core::SymbolSet_O::create();
    SpecificContextSet_O::iterator it;
    for ( it=this->begin(); it!=this->end(); it++ )
    {
	keys->insert(it->first);
    }
    return keys;
}



core::Cons_sp SpecificContextSet_O::asCons()
{_G();
    core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::Cons_O>());
    core::Cons_sp cur = first;
    SpecificContextSet_O::iterator it;
    for ( it=this->begin(); it!=this->end(); it++ )
    {
	core::Cons_sp one = core::Cons_O::create(it->second,_Nil<core::Cons_O>());
	cur->setCdr(one);
	cur = one;
    }
    return first->cdr();
}







    void SpecificContext_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<SpecificContext_O>()
	.def("setSelfName",&SpecificContext_O::setSelfName)
	.def("pushNeighbor",&SpecificContext_O::pushNeighbor)
	.def("popNeighbor",&SpecificContext_O::popNeighbor)
	.def("asString",&SpecificContext_O::asString,"","","",false)
	.def("copy",&SpecificContext_O::copy,"","","",false)
    ;
}

    void SpecificContext_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,SpecificContext,"","",_lisp)
	.def("setSelfName",&SpecificContext_O::setSelfName)
	.def("pushNeighbor",&SpecificContext_O::pushNeighbor)
	.def("popNeighbor",&SpecificContext_O::popNeighbor)
//	.def("asString",&SpecificContext_O::asString,"","","",false)
//	.def("copy",&SpecificContext_O::copy,"","","",false)
    ;
#endif //]
}





    void SpecificContextSet_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<SpecificContextSet_O>()
	.def("asSymbolSetOfKeys",&SpecificContextSet_O::asSymbolSetOfKeys)
	.def("size",&SpecificContextSet_O::size,"","","",false)
	.def("addSpecificContext",&SpecificContextSet_O::add)
	.def("asCons",&SpecificContextSet_O::asCons,"","","",false)
	;
}

    void SpecificContextSet_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,SpecificContextSet,"","",_lisp)
	.def("asSymbolSetOfKeys",&SpecificContextSet_O::asSymbolSetOfKeys)
//	.def("size",&SpecificContextSet_O::size,"","","",false)
	.def("addSpecificContext",&SpecificContextSet_O::add)
//	.def("asCons",&SpecificContextSet_O::asCons,"","","",false)
    ;
#endif //]
}



    EXPOSE_CLASS(chem,SpecificContext_O);
    EXPOSE_CLASS(chem,SpecificContextSet_O);
};




