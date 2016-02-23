/*
    File: alias.cc
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

#include <cando/chem/alias.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/environment.h>
#include <clasp/core/wrappers.h>
#include <cando/main/foundation.h>
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


    
    
#define ARGS_chem__alias "(monAlias atomAlias)"
#define DECL_chem__alias ""
#define DOCS_chem__alias "alias"
CL_DEFUN core::T_sp chem__alias(core::Symbol_sp monAlias, core::Symbol_sp atomAlias)
    {
	return Alias_O::create(monAlias,atomAlias);
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
    {
        GC_ALLOCATE(Alias_O, me );
	me->_MonomerAlias = monomerSym;
	me->_AtomAlias = atomSym;
	return me;
    };

#else

    core::T_sp Alias_O::__init__(core::Function_sp exec, core::List_sp args, core::Environment_sp bargs, core::Lisp_sp lisp)
    {
	this->_MonomerAlias = translate::from_object<core::Symbol_sp>::convert(bargs->lookup(ChemPkg,"monomerSym"));
	this->_AtomAlias = translate::from_object<core::Symbol_sp>::convert(bargs->lookup(ChemPkg,"atomSym"));
	return _Nil<core::T_O>();
    }

#endif

#if 0
    Fingerprint_sp Alias_O::fingerprint()
    {_OF();
	ASSERTF(this->_WeakFingerprint.pointerp(), BF("The _WeakFingerprint should never be NULL"));
	return this->_WeakFingerprint;
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
    {
	this->_MonomerId = alchemist->getMonomerIdForUniqueAlias(this->sharedThis<Alias_O>());
    }

    void Alias_O::sequenceChanged(PointProvider_sp alchemist)
    {
	this->_AtomId = alchemist->indexOfInterestingAtomWithAlias(this->_MonomerId,this->sharedThis<Alias_O>());
	LOG(BF("[Alias %s@%s] _MonomerId(%d) _AtomId(%d)") % this->_MonomerAlias->__repr__() % this->_AtomAlias->__repr__() % this->_MonomerId % this->_AtomId  );
    }


    Vector3& Alias_O::getUntransformedAtomPosition(PointProvider_sp alchemist)
    {
	return alchemist->getUntransformedAtomPosition(this->_AtomId);
    }


#endif



    bool Alias_O::equal(core::T_sp o) const
    {
	if (!o.isA<Alias_O>()) return false;
	Alias_sp a = o.as<Alias_O>();
	if ( this->_MonomerAlias != a->_MonomerAlias ) return false;
	if ( this->_AtomAlias != a->_AtomAlias ) return false;
	return true;
    }










};
