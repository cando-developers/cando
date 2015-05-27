#define DEBUG_LEVEL_FULL
       
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/foundation.h>
#include <clasp/core/hashTableEq.h>
#include <cando/units/symbolTable.h>
#include <cando/chem/atom.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/residue.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyComponent.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/forceField.h>
#include <cando/units/unitsPackage.h>
#include <cando/units/quantity.h>
#include <cando/chem/loop.h>
#include <cando/chem/spanningLoop.h>
#include <cando/chem/ffNonbondDb.h>
#include <clasp/core/profiler.h>
#include <clasp/core/wrappers.h>




namespace chem 
{





    EnergyAtom::EnergyAtom() {};
    EnergyAtom::EnergyAtom(Atom_sp atom,uint coordinateIndex) : Base(atom,coordinateIndex) {} ;
    EnergyAtom::EnergyAtom(ForceField_sp forceField, Atom_sp atom, uint coordinateIndex ) : Base(atom,coordinateIndex)
    {
	this->defineForAtom(forceField,atom,coordinateIndex);
    }


void EnergyAtom::defineForAtom(ForceField_sp forceField, Atom_sp a1, uint coordinateIndex)
{
    this->setupBase(a1,coordinateIndex);
    this->_TypeIndex = forceField->getNonbondDb()->findTypeIndex(a1->getType());
    this->_Charge = a1->getCharge();
}

#ifdef XML_ARCHIVE
    void	EnergyAtom::archive(core::ArchiveP node)
{
    this->Base::archive(node);
}
#endif


#if 0 //[
    adapt::QDomNode_sp	EnergyAtom::asXmlRelativeToContainer(chem::Matter_sp parent)
{
    adapt::QDomNode_sp	node;
    string	path;
    node = chem::new_adapt::QDomNode_sp("EnergyAtom");
    path = this->_Atom->getIdPath(parent);
    node->addAttributeString("storagePath",path);
    node->addAttributeString("atomName",this->_Atom->getName());
    node->addAttributeInt("coordinateIndex",this->_CoordinateIndex);
    node->addAttributeDoubleScientific("x",this->_Atom->getPosition().getX());
    node->addAttributeDoubleScientific("y",this->_Atom->getPosition().getY());
    node->addAttributeDoubleScientific("z",this->_Atom->getPosition().getZ());
    node->addAttributeDoubleScientific("fx",this->_Atom->getForce().getX());
    node->addAttributeDoubleScientific("fy",this->_Atom->getForce().getY());
    node->addAttributeDoubleScientific("fz",this->_Atom->getForce().getZ());
    return node;
}


    void	EnergyAtom::parseFromXmlRelativeToContainer(adapt::QDomNode_sp xml,
							chem::Matter_sp parent)
{_G();
string	path;
double	x,y,z,fx,fy,fz;
Vector3	v;
    LOG(BF("parsing xml with localName=%s") % xml->getLocalName().c_str() );
    path = xml->getAttributeString("storagePath");
    this->_Atom = downcast<chem::Atom_O>(parent->getContentForIdPath(path));
    ASSERTNOTNULL(this->_Atom);
    this->_CoordinateIndex = xml->getAttributeInt("coordinateIndex");
    x = xml->getAttributeDouble("x");
    y = xml->getAttributeDouble("y");
    z = xml->getAttributeDouble("z");
    v.set(x,y,z);
    this->_Atom->setPosition(v);
    fx = xml->getAttributeDouble("fx");
    fy = xml->getAttributeDouble("fy");
    fz = xml->getAttributeDouble("fz");
    v.set(fx,fy,fz);
    this->_Atom->setForce(v);
}
#endif //]

string		EnergyAtom::getResidueAndName()
{
    chem::Residue_sp	res;
    stringstream	ss;
    res = this->atom()->getResidueContainedBy().as<Residue_O>();
    ss.str("");
    ss << ":"<<res->getId()<<"@"<<this->atom()->getName();
    return ss.str();
}


    void AtomTable_O::initialize()
    {
        this->_AtomTableIndices = core::HashTableEq_O::create_default();
    }


    EnergyAtom* AtomTable_O::getEnergyAtomPointer(Atom_sp a)
    {_OF();
        core::T_mv it = this->_AtomTableIndices->gethash(a);
	if ( it.second().nilp() ) // it == this->_AtomTableIndices.end() )
	{
	    SIMPLE_ERROR(BF("Could not find atom[%s] in AtomTable") % a->__repr__() );
	}
	return &this->_Atoms[core::clasp_to_fixnum(it.as<core::Fixnum_O>())];
    }

    int AtomTable_O::addAtomInfo(Atom_sp atom, units::Quantity_sp charge, units::Quantity_sp mass, int typeIndex )
    {_OF();
	int coordinateIndex = this->getNumberOfAtoms()*3;
	EnergyAtom ea(atom,coordinateIndex);
	ea._Charge = charge->value_in_unit_asReal(units::_sym_elementaryCharge->symbolValue().as<units::Unit_O>());
	ea._Mass = mass->value_in_unit_asReal(units::_sym_daltons->symbolValue().as<units::Unit_O>());
	ea._TypeIndex = typeIndex;
	this->add(ea);
	return coordinateIndex;
    }


#ifdef XML_ARCHIVE
    void	AtomTable_O::archiveBase(core::ArchiveP node)
{_G();
    if ( node->loading() )
    {
        this->_Atoms.clear();
        core::VectorNodes::iterator	vi;
        for (vi=node->begin_Children(); vi!=node->end_Children(); vi++ )
        {
	    EnergyAtom ea(_lisp);
	    ea.archive(*vi);
	    this->_Atoms.push_back(ea);
	    (*vi)->setRecognized(true);
	}
    } else
    {
	core::ArchiveP child;
        gctools::Vec0<EnergyAtom>::iterator	eai;
        for ( eai=this->_Atoms.begin(); eai!=this->_Atoms.end(); eai++ )
	{
	    child = node->createChildNode("EnergyAtom");
	    eai->archive(child);
	}
    }
}
#endif


    EXPOSE_CLASS(chem,AtomTable_O);
    void AtomTable_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<AtomTable_O>()
	    .def("addAtomInfo",&AtomTable_O::addAtomInfo)
	    ;
    }

    void AtomTable_O::exposePython(core::Lisp_sp lisp)
    {_G()
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,AtomTable,"","",_lisp)
	    .def("addAtomInfo",&AtomTable_O::addAtomInfo)
	    ;
#endif
    }




void	AtomTable_O::dumpTerms()
{
    gctools::Vec0<EnergyAtom>::iterator	eai;
string				as1,as2,as3,as4;
core::Symbol_sp				str1, str2, str3, str4;
    for ( eai=this->_Atoms.begin(); eai!=this->_Atoms.end(); eai++ ) {
	as1 = atomLabel(eai->atom());
	str1 = eai->atom()->getType();
	_lisp->print(BF("(TERM 0ATOM %-9s %-9s :charge %8.5lf :mass %8.5lf :typeIndex %d)")
                     % as1
                     % _rep_(str1)
                     % eai->_Charge
                     % eai->_Mass
                     % eai->_TypeIndex );
    }
}


    void AtomTable_O::constructFromMatter(Matter_sp matter, ForceField_sp forceField)
    {_OF();
	this->_Atoms.clear();
	this->_AtomTableIndices->clrhash();
	uint idx = 0;
	uint coordinateIndex = 0;
	{_BLOCK_TRACE("Defining ATOMS");
	    Loop loop;
	    Atom_sp a1;
	    loop.loopTopGoal(matter,ATOMS);
	    while ( loop.advanceLoopAndProcess() ) 
	    {
		a1 = loop.getAtom();
		if ( a1->isAssignableTo<VirtualAtom_O>() ) 
		{
		    LOG(BF("Skipping virtual atom[%s]") % a1->__repr__() );
		    continue; // skip virtuals
		}
		LOG(BF("Setting atom[%s] in AtomTable[%d]") % a1->__repr__() % idx );
		this->_AtomTableIndices->setf_gethash(a1,core::clasp_make_fixnum(idx));
		EnergyAtom ea(forceField,a1,coordinateIndex);
		{_BLOCK_TRACE("Building spanning tree for atom");
		    SpanningLoop_sp span = SpanningLoop_O::create(a1);
		    Atom_sp bonded;
		    while ( span->advance() )
		    {
			bonded = span->getAtom();
			if ( bonded == a1 ) continue;
			int backCount = bonded->getBackCount();
			LOG(BF("Looking at atom[%s] at remove[%d]") % bonded->__repr__() % backCount );
			// Once we crawl out 4 bonds we have gone as far as we need
			if ( backCount >= 4 ) 
			{
			    LOG(BF("Hit remove of 4 - terminating spanning loop"));
			    break;
			}
			ASSERT(backCount>0 && backCount<=3);
			ea._AtomsAtRemove[backCount-1].insert(bonded);
		    }
		}
		this->_Atoms.push_back(ea);
		idx++;
		coordinateIndex += 3;
#ifdef DEBUG_ON
		stringstream ss;
                gctools::SmallOrderedSet<Atom_sp>::iterator si;
		for ( int zr = 0; zr<3; zr++ )
		{
		    ss.str("");
		    for ( si = ea._AtomsAtRemove[zr].begin(); si!=ea._AtomsAtRemove[zr].end(); si++ )
		    {
			ss << " " << (*si)->__repr__();
		    }
		    LOG(BF("Atoms at remove %d = %s") % zr % ss.str() );
		}
#endif
	    }
	}
    }




};
