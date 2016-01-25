#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <cando/adapt/stringSet.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/matter.h>
#include <cando/chem/residue.h>
#include <cando/chem/calculatePosition.h>
#include <cando/chem/atom.h>
#include <cando/chem/bond.h>
#include <clasp/core/wrappers.h>

namespace chem
{




;

    ConstitutionBond_sp ConstitutionBond_O::create(BondOrder order, ConstitutionAtomIndex0N index)
    {
	GC_ALLOCATE(ConstitutionBond_O, me );
	me->_BondOrder = order;
	me->_AtomIndex = index;
	return me;
    };


    void ConstitutionBond_O::initialize()
    {_OF();
	this->_AtomIndex = UndefinedUnsignedInt;
	this->_BondOrder = noBond;
    }

#ifdef XML_ARCHIVE
    void ConstitutionBond_O::archiveBase(core::ArchiveP node)
    {_OF();
	this->Base::archiveBase(node);
	node->attribute("ai",this->_AtomIndex);
	node->attributeSymbolEnumHiddenConverter("bo",this->_BondOrder,_sym__PLUS_bondOrderToSymbolConverter_PLUS_);
    }
#endif






;


    ConstitutionAtom_sp ConstitutionAtom_O::create(MatterName name, Element element,
						   ConstitutionAtomIndex0N index )
    {
	ConstitutionAtom_sp ta = ConstitutionAtom_O::create();
	ta->_AtomName = name;
	ta->_Element = element;
	ta->_AtomIndex = index;
	return ta;
    }

    void ConstitutionAtom_O::initialize()
    {
	this->_AtomName = _Nil<MatterName::Type>();
	this->_AtomIndex = UndefinedUnsignedInt;
    }


#ifdef XML_ARCHIVE
    void ConstitutionAtom_O::archiveBase(core::ArchiveP node)
    {_OF();
	this->Base::archiveBase(node);
	node->attribute("atomName",this->_AtomName);
	node->attribute("ci",this->_AtomIndex);
	node->attributeSymbolEnumHiddenConverter("element",this->_Element,_sym__PLUS_elementToSymbolConverter_PLUS_);
	node->archiveVector0("bonds",this->_Bonds);
    }
#endif

    string ConstitutionAtom_O::__repr__() const
    {
	stringstream ss;
	ss << "( " << this->className() << " ";
	ss << " :atomName \"" << this->_AtomName << "\"";
	ss << " :element " << symbolFromElement(this->_Element);
	ss << " :atomIndex " << this->_AtomIndex << " )";
	return ss.str();
    }




;



    ConstitutionVirtualAtom_sp ConstitutionVirtualAtom_O::create(MatterName name
                                                                 , ConstitutionAtomIndex0N index
								 , CalculatePosition_sp code )
    {
	ConstitutionVirtualAtom_sp ta = ConstitutionVirtualAtom_O::create();
	ta->_AtomName = name;
	ta->_AtomIndex = index;
	ta->_CalculatePositionCode = code;
	return ta;
    }

    ConstitutionVirtualAtom_sp ConstitutionVirtualAtom_O::create(MatterName name, CalculatePosition_sp code)
    {
	ConstitutionVirtualAtom_sp ta = ConstitutionVirtualAtom_O::create();
	ta->_AtomName = name;
	ta->_AtomIndex = UndefinedUnsignedInt;
	ta->_CalculatePositionCode = code;
	return ta;
    }

    void ConstitutionVirtualAtom_O::initialize()
    {
	this->Base::initialize();
    }

#ifdef XML_ARCHIVE
    void ConstitutionVirtualAtom_O::archiveBase(core::ArchiveP node)
    {_OF();
	this->Base::archiveBase(node);
    }
#endif











;

    ConstitutionAtoms_sp ConstitutionAtoms_O::create(Residue_sp residue)
    {
	GC_ALLOCATE(ConstitutionAtoms_O, catoms );
	Matter_O::contentIterator it;
	ConstitutionAtomIndex0N index = 0;
	MapAtomsToConstitutionAtomIndex0N atomToIndexMap;
	{_BLOCK_TRACEF(BF("First assign each atom a unique ConstitutionAtomIndex0N"));
	    for ( it=residue->begin_atoms(); it!=residue->end_atoms(); it++, index++ )
	    {
		Atom_sp atom = (*it).as<Atom_O>();
		ConstitutionAtom_sp catom = atom->asConstitutionAtom(index);
		catoms->_Atoms.push_back(catom);
		atomToIndexMap[atom] = index;
	    }
	}
	{_BLOCK_TRACE("Then define the bonds for each atom using the indices");
	    index = 0; // start indexing from 0 again
	    for ( it=residue->begin_atoms(); it!=residue->end_atoms(); it++, index++ )
	    {
		Atom_sp atom = (*it).as<Atom_O>();
		ConstitutionAtom_sp catom = catoms->_Atoms[index];
		atom->defineConstitutionAtomBonding(catom,atomToIndexMap);
	    }
	}
	return catoms;
    }


    void ConstitutionAtoms_O::addConstitutionVirtualAtom(ConstitutionVirtualAtom_sp ca)
    {_OF();
	ConstitutionAtomIndex0N nextIndex = this->_Atoms.size();
	ca->_AtomIndex = nextIndex;
	this->_Atoms.push_back(ca);
    };


    ConstitutionAtom_sp& ConstitutionAtoms_O::operator[](ConstitutionAtomIndex0N idx)
    {_OF();
	ASSERTF(idx>=0&&idx<(int)(this->_Atoms.size()),BF("index[%d] out of range") % idx );
	ConstitutionAtom_sp& ca = this->_Atoms[(uint)idx];
	ASSERTF(idx==ca->_AtomIndex,BF("A mismatch has occured between the index[%d] in the ConstitutionAtoms array and the ConstitutionAtom index[%d] - these have to match or we can't quickly look up atoms by their index") % idx % ca->_AtomIndex );
	return ca;
    }

CL_LISPIFY_NAME("atomNamesAsSymbolSet");
CL_DEFMETHOD     adapt::SymbolSet_sp ConstitutionAtoms_O::atomNamesAsSymbolSet()
    {_OF();
	adapt::SymbolSet_sp ss = adapt::SymbolSet_O::create();
	for ( gctools::Vec0<ConstitutionAtom_sp>::const_iterator ci=this->_Atoms.begin();
	      ci!=this->_Atoms.end(); ci++ )
	{
	    ss->insert((*ci)->_AtomName);
	}
	return ss;
    };

    void ConstitutionAtoms_O::initialize()
    {_OF();
	this->_NetCharge = 0.0;
    }

#ifdef XML_ARCHIVE
    void ConstitutionAtoms_O::archiveBase(core::ArchiveP n)
    {
	this->Base::archiveBase(n);
	n->attribute("netCharge",this->_NetCharge);
	n->archiveVector0("atoms",this->_Atoms);
    }
#endif


CL_LISPIFY_NAME("atomWithName");
CL_DEFMETHOD     ConstitutionAtom_sp ConstitutionAtoms_O::atomWithName(MatterName nm)
    {_OF();
	for ( gctools::Vec0<ConstitutionAtom_sp>::const_iterator ci=this->_Atoms.begin();
	      ci!=this->_Atoms.end(); ci++ )
	{
	    if ( (*ci)->_AtomName == nm ) return *ci;
	}
	SIMPLE_ERROR(BF("Could not find ConstitutionAtom with name[%s]") % nm );
    }


CL_LISPIFY_NAME("atomWithId");
CL_DEFMETHOD     ConstitutionAtom_sp ConstitutionAtoms_O::atomWithId(ConstitutionAtomIndex0N id) const
    {
	return this->_Atoms[id];
    }

CL_LISPIFY_NAME("constitution-atoms-index");
CL_DEFMETHOD     int ConstitutionAtoms_O::index(MatterName name) const
    {_OF();
	for ( int idx=0; idx<(int)this->_Atoms.size(); idx++ )
	{
	    const ConstitutionAtom_sp& atom = this->_Atoms[idx];
	    if ( atom->_AtomName == name ) return idx;
	}
	SIMPLE_ERROR(BF("Unknown atom[%s]") % name );
    }

CL_LISPIFY_NAME("buildResidue");
CL_DEFMETHOD     Residue_sp ConstitutionAtoms_O::buildResidue()
    {_OF();
	Residue_sp res = Residue_O::create();
	res->setNetCharge(this->_NetCharge);
        gctools::Vec0<Atom_sp> atoms;
	atoms.resize(this->_Atoms.size());
	res->resizeContents(this->_Atoms.size());
	for ( gctools::Vec0<ConstitutionAtom_sp>::const_iterator ci=this->_Atoms.begin();
	      ci!=this->_Atoms.end(); ci++ )
	{
	    Atom_sp atom = Atom_O::create();
	    atom->setName((*ci)->_AtomName);
	    atom->setElement((*ci)->_Element);
	    atom->setId((*ci)->_AtomIndex);
	    atoms[(*ci)->_AtomIndex] = atom;
	    res->putMatter((*ci)->_AtomIndex,atom);
	}
	for ( gctools::Vec0<ConstitutionAtom_sp>::const_iterator ci=this->_Atoms.begin();
	      ci!=this->_Atoms.end(); ci++ )
	{
	    ASSERT_lt((*ci)->_AtomIndex,(int)atoms.size());
	    Atom_sp fromAtom = atoms[(*ci)->_AtomIndex];
	    for ( gctools::Vec0<ConstitutionBond_sp>::const_iterator bi=(*ci)->_Bonds.begin(); bi!=(*ci)->_Bonds.end(); bi++ )
	    {
		// Avoid making bonds twice
		if ( (*bi)->_AtomIndex < (*ci)->_AtomIndex ) continue;
		ASSERT_lt((*bi)->_AtomIndex,(int)atoms.size());
		Atom_sp toAtom = atoms[(*bi)->_AtomIndex];
		BondOrder order = (*bi)->_BondOrder;
		fromAtom->bondTo(toAtom,order);
	    }
	}
	return res;
    }

};
