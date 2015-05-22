#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "atomId.h"
#include "atomIdMap.h"
#include "atom.h"
#include "core/wrappers.h"
namespace chem
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(chem,AtomIdToAtomMap_O);

    void AtomIdToAtomMap_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<AtomIdToAtomMap_O>()
	    .def("lookup-atom",&AtomIdToAtomMap_O::lookupAtom)
	;
    }

    Atom_sp AtomIdToAtomMap_O::lookupAtom(AtomId_sp atomId) const
    {_OF();
	return this->_AtomIdMap[atomId->_AtomId];
    }

    void AtomIdToAtomMap_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,AtomIdToAtomMap,"","",_lisp)
	    ;
#endif
    }

#ifdef XML_ARCHIVE
    void AtomIdToAtomMap_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
	IMPLEMENT_ME();
    }
#endif


    void AtomIdToAtomMap_O::initialize()
    {_OF();
        this->Base::initialize();
    }


    void AtomIdToAtomMap_O::resize(int numMols)
    {_OF();
	this->_AtomIdMap.resize(numMols);
    }

    void AtomIdToAtomMap_O::resize(int mol, int numRes)
    {_OF();
	this->_AtomIdMap.resize(mol,numRes);
    }

    void AtomIdToAtomMap_O::resize(int mol, int res, int numAtoms)
    {_OF();
	this->_AtomIdMap.resize(mol,res,numAtoms);
    }
    void AtomIdToAtomMap_O::set(AtomId const& atomId, Atom_sp atom)
    {_OF();
	this->_AtomIdMap[atomId] = atom;
    }

    

}; /* chem */
