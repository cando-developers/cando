#define	DEBUG_LEVEL_FULL

#include "alias.h"
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include "core/environment.h"
#include "core/wrappers.h"

namespace chem 
{


    void	Alias_O::initialize()
    {
	this->Base::initialize();
//	this->_WeakFingerprint = _Nil<Fingerprint_O>();
    }

#ifdef XML_ARCHIVE
    void	Alias_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attribute("MonomerAlias",this->_MonomerAlias);
	node->attribute("AtomAlias",this->_AtomAlias);
    }
#endif


/*! \subpage script-functions
 * alias 'monAlias 'atomAlias  -  Define an monomer/atom alias combo
 */


    
    
#define ARGS_af_alias "(monAlias atomAlias)"
#define DECL_af_alias ""
#define DOCS_af_alias "alias"
    core::T_sp af_alias(core::Symbol_sp monAlias, core::Symbol_sp atomAlias)
    {_G();
	return Alias_O::create(_lisp,monAlias,atomAlias);
    }



/*
  __BEGIN_DOC(classes.Alias.!class.Alias)
  \positional{Text::monomerAlias Text::atomAlias}

  Create an Alias object that maintains a \sa{monomerAlias} name and an \sa{atomAlias} name.
  __END_DOC
*/

#if INIT_TO_FACTORIES

#define ARGS_Alias_O_make "(monomer_sym atom_sym)"
#define DECL_Alias_O_make ""
#define DOCS_Alias_O_make "make AtomId args: monomer_sym atom_sym"
    Alias_sp Alias_O::make(core::Symbol_sp monomerSym, core::Symbol_sp atomSym)
    {_G();
        GC_ALLOCATE(Alias_O, me );
	me->_MonomerAlias = monomerSym;
	me->_AtomAlias = atomSym;
	return me;
    };

#else

    core::T_sp Alias_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp bargs, core::Lisp_sp lisp)
    {_G();
	this->_MonomerAlias = translate::from_object<core::Symbol_sp>::convert(bargs->lookup(ChemPkg,"monomerSym"));
	this->_AtomAlias = translate::from_object<core::Symbol_sp>::convert(bargs->lookup(ChemPkg,"atomSym"));
	return _Nil<core::T_O>();
    }

#endif

#if 0
    Fingerprint_sp Alias_O::fingerprint()
    {_OF();
	ASSERTF(this->_WeakFingerprint.pointerp(), BF("The _WeakFingerprint should never be NULL"));
	return this->_WeakFingerprint.lock();
    }
#endif


    string Alias_O::description() const
    {
	stringstream ss;
	ss << "[Alias " << this->_MonomerAlias << "@" << this->_AtomAlias;
	ss << " monomerId("<<this->_MonomerId<<") atomId("<<this->_AtomId<<")";
	return ss.str();
    }


#if 0
    void Alias_O::oligomerChanged(PointProvider_sp alchemist)
    {_G();
	this->_MonomerId = alchemist->getMonomerIdForUniqueAlias(this->sharedThis<Alias_O>());
    }

    void Alias_O::sequenceChanged(PointProvider_sp alchemist)
    {_G();
	this->_AtomId = alchemist->indexOfInterestingAtomWithAlias(this->_MonomerId,this->sharedThis<Alias_O>());
	LOG(BF("[Alias %s@%s] _MonomerId(%d) _AtomId(%d)") % this->_MonomerAlias->__repr__() % this->_AtomAlias->__repr__() % this->_MonomerId % this->_AtomId  );
    }


    Vector3& Alias_O::getUntransformedAtomPosition(PointProvider_sp alchemist)
    {
	return alchemist->getUntransformedAtomPosition(this->_AtomId);
    }


#endif



    bool Alias_O::equal(core::T_sp o) const
    {_G();
	if (!o.isA<Alias_O>()) return false;
	Alias_sp a = o.as<Alias_O>();
	if ( this->_MonomerAlias != a->_MonomerAlias ) return false;
	if ( this->_AtomAlias != a->_AtomAlias ) return false;
	return true;
    }





    void Alias_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<Alias_O>()
//	    .def_raw("core:__init__",&Alias_O::__init__,"(self chem::monomerSym chem::atomSym)")
	    ;
	Defun(alias);
    }

    void Alias_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,Alias,"","",_lisp)
	    ;
#endif
    }


    EXPOSE_CLASS(chem,Alias_O);

};


