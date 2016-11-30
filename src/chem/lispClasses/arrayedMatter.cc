/*
    File: arrayedMatter.cc
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

#include "arrayedMatter.h"
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include "core/sort.h"
#include "geom/coordinateArray.h"
#include "matter.h"
#include "atom.h"
#include "residue.h"
#include "molecule.h"
#include "core/wrappers.h"

namespace chem {




//
// ArrayedAtom


void	ArrayedAtom_O::initialize()
{
    this->Base::initialize();
    this->_AtomLabel = _Nil<MatterName::Type>();
    this->_Atom.reset();
}

#ifdef XML_ARCHIVE
    void	ArrayedAtom_O::archive(core::ArchiveP node)
{
    node->attribute("name",this->_AtomLabel);
    node->archiveWeakPointer("atom",this->_Atom);
}
#endif


void ArrayedAtom_O::setAtom(Atom_sp a)
{
    this->_AtomLabel = a->getName();
    this->_Atom = a;
}


uint ArrayedAtom_O::calculateHash()
{
    uint hash=0;
    for ( string::iterator it=this->_AtomLabel.begin();
    	it !=this->_AtomLabel.end(); it++ )
    {
	hash += (*it);
	hash = hash << 1;
    }
    if ( hash == UndefinedUnsignedInt ) hash = 1;
    return hash;
}


void ArrayedAtom_O::extractFromCoordinateArray(geom::CoordinateArrayWithHash_sp h, uint idx)
{
    ASSERT(this->_Atom.pointerp());
    this->_Atom->setPosition(h->getElement(idx));
}

void ArrayedAtom_O::writeToCoordinateArray(geom::CoordinateArrayWithHash_sp h, uint idx)
{
    ASSERT(this->_Atom.pointerp());
    h->setElement(idx,this->_Atom->getPosition());
}















//
// ArrayedResidue
//


void	ArrayedResidue_O::initialize()
{
    this->Base::initialize();
    this->_ResidueLabel = "";
    this->_Residue.reset();
    this->_SortedAtoms.clear();
}

#ifdef XML_ARCHIVE
    void	ArrayedResidue_O::archive(core::ArchiveP node)
{
    node->attribute("label",this->_ResidueLabel);
    node->archiveWeakPointer("residue",this->_Residue);
    node->archiveVector0("sortedAtoms",this->_SortedAtoms);
}
#endif


class	OrderByAtomLabel
{
public:
    bool operator()(ArrayedAtom_sp x, ArrayedAtom_sp y )
    {
	return x->_AtomLabel <= y->_AtomLabel;
    }
};
void ArrayedResidue_O::setResidue(Residue_sp res)
{
    if ( res->getUniqueLabel() == "" )
    {
	SIMPLE_ERROR(BF("Every ArrayedResidue must have a label"));
    }
    this->_ResidueLabel = res->getUniqueLabel();
    this->_Residue = res;
    Matter_O::contentIterator it;
    this->_SortedAtoms.clear();
    for ( it=res->begin_atoms(); it!=res->end_atoms(); it++ )
    {
	ArrayedAtom_sp atomRef = ArrayedAtom_O::create();
	Atom_sp a = (*it).as<Atom_O>();
	atomRef->setAtom(a);
	this->_SortedAtoms.push_back(atomRef);
    }
    OrderByAtomLabel orderer;
    core::sort::quickSort(this->_SortedAtoms.begin(),
		    this->_SortedAtoms.end(),orderer);
}

uint ArrayedResidue_O::getNumberOfAtoms()
{
    return this->_SortedAtoms.size();
}


uint ArrayedResidue_O::calculateHash()
{
    uint hash=0;
    for ( gctools::Vec0<ArrayedAtom_sp>::iterator it=this->_SortedAtoms.begin();
    	it !=this->_SortedAtoms.end(); it++ )
    {
	hash += (*it)->calculateHash();
	hash = hash << 1;
    }
    if ( hash == UndefinedUnsignedInt ) hash = 1;
    return hash;
}

uint ArrayedResidue_O::extractFromCoordinateArray(geom::CoordinateArrayWithHash_sp hash, uint start)
{
    for ( gctools::Vec0<ArrayedAtom_sp>::iterator it=this->_SortedAtoms.begin();
    	it !=this->_SortedAtoms.end(); it++ )
    {
	(*it)->extractFromCoordinateArray(hash,start);
	start++;
    }
    return this->_SortedAtoms.size();
}

uint ArrayedResidue_O::writeToCoordinateArray(geom::CoordinateArrayWithHash_sp hash, uint start)
{
    for ( gctools::Vec0<ArrayedAtom_sp>::iterator it=this->_SortedAtoms.begin();
    	it !=this->_SortedAtoms.end(); it++ )
    {
	(*it)->writeToCoordinateArray(hash,start);
	start++;
    }
    return this->_SortedAtoms.size();
}











//
// ArrayedMolecule
//

    core::T_sp 	ArrayedMolecule_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
{
    return _Nil<core::T_O>();
}

void	ArrayedMolecule_O::initialize()
{
    this->Base::initialize();
    this->_SortedResidues.clear();
}

#ifdef XML_ARCHIVE
    void	ArrayedMolecule_O::archive(core::ArchiveP node)
{
    node->archiveVector0("sortedResidues",this->_SortedResidues);
    node->attribute("hash",this->_Hash);
    node->attribute("numberOfAtoms",this->_NumberOfAtoms);
}
#endif

class	OrderByResidueLabel
{
public:
    bool operator()(ArrayedResidue_sp x, ArrayedResidue_sp y )
    {
	return x->_ResidueLabel <= y->_ResidueLabel;
    }
};

void ArrayedMolecule_O::setMolecule(Molecule_sp mol)
{
    this->_Molecule = mol;
    Matter_O::contentIterator it;
    this->_SortedResidues.clear();
    for ( it=mol->begin_residues(); it!=mol->end_residues(); it++ )
    {
	ArrayedResidue_sp linearResidue = ArrayedResidue_O::create();
	Residue_sp res = (*it).as<Residue_O>();
	linearResidue->setResidue(res);
	this->_SortedResidues.push_back(linearResidue);
    }
    OrderByResidueLabel orderer;
    core::sort::quickSort(this->_SortedResidues.begin(),
		    this->_SortedResidues.end(),orderer);
    this->_Hash = this->calculateHash();
    this->_NumberOfAtoms = this->getNumberOfAtoms();
}


uint	ArrayedMolecule_O::getNumberOfAtoms()
{
    uint num = 0;
    for ( gctools::Vec0<ArrayedResidue_sp>::iterator it=this->_SortedResidues.begin();
    		it!=this->_SortedResidues.end(); it++ )
    {
	num += (*it)->getNumberOfAtoms();
    }
    return num;
}

uint	ArrayedMolecule_O::calculateHash()
{
    uint hash = 0;
    for ( gctools::Vec0<ArrayedResidue_sp>::iterator it=this->_SortedResidues.begin();
    		it!=this->_SortedResidues.end(); it++ )
    {
	hash += (*it)->calculateHash();
    }
    if ( hash == UndefinedUnsignedInt ) hash = 1;
    return hash;
}



geom::CoordinateArrayWithHash_sp ArrayedMolecule_O::createCoordinateArray()
{
    uint numberOfAtoms = this->getNumberOfAtoms();
    geom::CoordinateArrayWithHash_sp ca;
    ca = geom::CoordinateArrayWithHash_O::create(_lisp,numberOfAtoms);
    ca->setHash(this->_Hash);
    return ca;
}


void ArrayedMolecule_O::extractFromCoordinateArray(geom::CoordinateArrayWithHash_sp h)
{
    if ( h->getHash() != this->_Hash )
    {
	SIMPLE_ERROR(BF("The hashes don't match"));
    }
    if ( h->size() != this->getNumberOfAtoms() )
    {
	SIMPLE_ERROR(BF("The number of atoms don't match"));
    }
    uint start=0;
    for ( gctools::Vec0<ArrayedResidue_sp>::iterator it=this->_SortedResidues.begin();
    		it!=this->_SortedResidues.end(); it++ )
    {
	start += (*it)->extractFromCoordinateArray(h,start);
    }
}

void ArrayedMolecule_O::writeToCoordinateArray(geom::CoordinateArrayWithHash_sp h)
{
    if ( h->getHash() != this->_Hash )
    {
	SIMPLE_ERROR(BF("The hashes don't match"));
    }
    if ( h->size() != this->getNumberOfAtoms() )
    {
	SIMPLE_ERROR(BF("The number of atoms don't match"));
    }
    uint start=0;
    for ( gctools::Vec0<ArrayedResidue_sp>::iterator it=this->_SortedResidues.begin();
    		it!=this->_SortedResidues.end(); it++ )
    {
	start += (*it)->writeToCoordinateArray(h,start);
    }
}








};
