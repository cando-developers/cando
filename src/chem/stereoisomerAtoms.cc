#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/stringSet.h"
#include "core/symbolSet.h"
#include "stereoisomerAtoms.h"
#include "constitutionAtoms.h"
#include "matter.h"
#include "residue.h"
#include "calculatePosition.h"
#include "atom.h"
#include "bond.h"
#include "core/wrappers.h"

namespace chem
{

    EXPOSE_CLASS(chem,StereoisomerAtom_O);

    void StereoisomerAtom_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<StereoisomerAtom_O>()
	    .def("atomName",&StereoisomerAtom_O::atomName)
	    .def("setCharge",&StereoisomerAtom_O::setCharge)
	    .def("getCharge",&StereoisomerAtom_O::getCharge)
	    .def("setAtomType",&StereoisomerAtom_O::setAtomType)
	    .def("getAtomType",&StereoisomerAtom_O::getAtomType)
	    ;
    }

    void StereoisomerAtom_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,StereoisomerAtom,"","",_lisp)
	    .def("atomName",&StereoisomerAtom_O::atomName)
	    .def("setCharge",&StereoisomerAtom_O::setCharge)
	    .def("getCharge",&StereoisomerAtom_O::getCharge)
	    .def("setAtomTypeString",&StereoisomerAtom_O::setAtomTypeString)
	    .def("getAtomTypeString",&StereoisomerAtom_O::getAtomTypeString)
	    ;
#endif
    }

;


    StereoisomerAtom_sp StereoisomerAtom_O::create(ConstitutionAtom_sp constitutionAtom )
    {_G();
	if ( constitutionAtom->isAssignableTo<ConstitutionVirtualAtom_O>() )
	{
	    ConstitutionVirtualAtom_sp cva = constitutionAtom.as<ConstitutionVirtualAtom_O>();
	    StereoisomerVirtualAtom_sp tva = StereoisomerVirtualAtom_O::create(cva);
	    return tva;
	}
	GC_ALLOCATE(StereoisomerAtom_O, ta );
	ta->_AtomName = constitutionAtom->_AtomName;
	ta->_AtomIndex = constitutionAtom->_AtomIndex;
	ta->_Charge = 0.0;
	ta->_AtomType = _Nil<core::Symbol_O>();
	return ta;
    }

    void StereoisomerAtom_O::initialize()
    {
	this->_AtomName = _Nil<core::Symbol_O>();
	this->_AtomIndex = UndefinedUnsignedInt;
	this->_AtomType = _Nil<core::Symbol_O>();
	this->_Charge = 0.0;
    }


#ifdef XML_ARCHIVE
    void StereoisomerAtom_O::archiveBase(core::ArchiveP node)
    {_OF();
	this->Base::archiveBase(node);
	node->attribute("atomName",this->_AtomName);
	node->attribute("ci",this->_AtomIndex);
	node->attribute("charge",this->_Charge);
	node->attribute("type",this->_AtomType);
    }
#endif

    string StereoisomerAtom_O::__repr__() const
    {_G();
	stringstream ss;
	ss << "( " << this->className() << " ";
	ss << " :atomName \"" << this->_AtomName << "\"";
	ss << " :atomIndex " << this->_AtomIndex << " )";
	return ss.str();
    }



    EXPOSE_CLASS(chem,StereoisomerVirtualAtom_O);

    void StereoisomerVirtualAtom_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<StereoisomerVirtualAtom_O>()
	    ;
    }

    void StereoisomerVirtualAtom_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,StereoisomerVirtualAtom,"","",_lisp)
	    ;
#endif
    }

;



    StereoisomerVirtualAtom_sp StereoisomerVirtualAtom_O::create(ConstitutionVirtualAtom_sp atom )
    {_G();
	StereoisomerVirtualAtom_sp ta = StereoisomerVirtualAtom_O::create();
	ta->_AtomName = atom->_AtomName;
	ta->_AtomIndex = atom->_AtomIndex;
	return ta;
    }

    void StereoisomerVirtualAtom_O::initialize()
    {
	this->Base::initialize();
    }

#ifdef XML_ARCHIVE
    void StereoisomerVirtualAtom_O::archiveBase(core::ArchiveP node)
    {_OF();
	this->Base::archiveBase(node);
    }
#endif







    EXPOSE_CLASS(chem,StereoisomerAtoms_O);

    void StereoisomerAtoms_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<StereoisomerAtoms_O>()
	    .def("atomNamesAsSymbolSet",&StereoisomerAtoms_O::atomNamesAsSymbolSet)
	    .def("atomWithName",&StereoisomerAtoms_O::atomWithName)
	    .def("atomWithId",&StereoisomerAtoms_O::atomWithId)
	    .def("stereoisomer-atoms-index",&StereoisomerAtoms_O::index)
	    ;
    }

    void StereoisomerAtoms_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,StereoisomerAtoms,"","",_lisp)
	    .def("atomNamesAsStringSet",&StereoisomerAtoms_O::atomNamesAsStringSet)
	    .def("atomWithName",&StereoisomerAtoms_O::atomWithName)
	    .def("atomWithId",&StereoisomerAtoms_O::atomWithId)
	    .def("stereoisomer-atoms-index",&StereoisomerAtoms_O::index)
	    ;
#endif
    }

;

    StereoisomerAtoms_sp StereoisomerAtoms_O::create(ConstitutionAtoms_sp constitutionAtoms)
    {_G();
	GC_ALLOCATE(StereoisomerAtoms_O, topologyAtoms );
	for ( ConstitutionAtoms_O::iterator it=constitutionAtoms->_Atoms.begin();
	      it!=constitutionAtoms->_Atoms.end(); it++)
	{
	    ConstitutionAtom_sp ca = *it;
	    StereoisomerAtom_sp ta = StereoisomerAtom_O::create(ca);
	    topologyAtoms->_Atoms.push_back(ta);
	}
	return topologyAtoms;
    }


    void StereoisomerAtoms_O::addStereoisomerVirtualAtom(StereoisomerVirtualAtom_sp ca)
    {_OF();
	ConstitutionAtomIndex0N nextIndex = this->_Atoms.size();
	ca->_AtomIndex = nextIndex;
	this->_Atoms.push_back(ca);
    };


    StereoisomerAtom_sp& StereoisomerAtoms_O::operator[](ConstitutionAtomIndex0N idx)
    {_OF();
	ASSERTF(idx>=0&&idx<(int)(this->_Atoms.size()),BF("index[%d] out of range") % idx );
	StereoisomerAtom_sp& ca = this->_Atoms[(uint)idx];
	ASSERTF(idx==ca->_AtomIndex,BF("A mismatch has occured between the index[%d] in the StereoisomerAtoms array and the StereoisomerAtom index[%d] - these have to match or we can't quickly look up atoms by their index") % idx % ca->_AtomIndex );
	return ca;
    }

    core::SymbolSet_sp StereoisomerAtoms_O::atomNamesAsSymbolSet()
    {_OF();
	core::SymbolSet_sp ss = core::SymbolSet_O::create();
	for ( gctools::Vec0<StereoisomerAtom_sp>::const_iterator ci=this->_Atoms.begin();
	      ci!=this->_Atoms.end(); ci++ )
	{
	    ss->insert((*ci)->_AtomName);
	}
	return ss;
    };

    void StereoisomerAtoms_O::initialize()
    {_OF();
    }

#ifdef XML_ARCHIVE
    void StereoisomerAtoms_O::archiveBase(core::ArchiveP n)
    {
	this->Base::archiveBase(n);
	n->archiveVector0("atoms",this->_Atoms);
    }
#endif


    StereoisomerAtom_sp StereoisomerAtoms_O::atomWithName(MatterName nm)
    {_OF();
	for ( gctools::Vec0<StereoisomerAtom_sp>::const_iterator ci=this->_Atoms.begin();
	      ci!=this->_Atoms.end(); ci++ )
	{
	    if ( (*ci)->_AtomName == nm ) return *ci;
	}
	SIMPLE_ERROR(BF("Could not find StereoisomerAtom with name[%s]") % nm );
    }


    StereoisomerAtom_sp StereoisomerAtoms_O::atomWithId(ConstitutionAtomIndex0N id) const
    {_G();
	return this->_Atoms[id];
    }

    int StereoisomerAtoms_O::index(MatterName name) const
    {_OF();
	for ( int idx=0; idx<(int)this->_Atoms.size(); idx++ )
	{
	    const StereoisomerAtom_sp& atom = this->_Atoms[idx];
	    if ( atom->_AtomName == name ) return idx;
	}
	SIMPLE_ERROR(BF("Unknown atom[%s]") % name );
    }


};
