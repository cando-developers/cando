/*
    File: energyAtomTable.h
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
       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	energyAtomTable.h
 *
 *	Maintains a database of atoms
 */

#ifndef EnergyAtomTable_H
#define	EnergyAtomTable_H //[
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/nativeVector.fwd.h>
#include <cando/units/quantity.fwd.h>
#include <cando/chem/indirectAtomCoordinateReference.h>
#include <cando/geom/vector3.h>


namespace chem
{

SMART(Atom);
SMART(ForceField);
SMART(Matter);

/*! Store a pointer to an Atom and an index into the coordinate vector array
 * which stores coordinates in a 1D array (x1,y1,z1,x2,y2,z2,x3,...,xN,yN,zN)
 */
class EnergyAtom : public IndirectAtomCoordinateReference {
    typedef IndirectAtomCoordinateReference Base;
public:
//        Atom_sp  	_Atom;			//!<Atom
//        uint     	_CoordinateIndex;	//!<index into the coordinate vector (atomNumber)*3
    MatterName		_AtomName;		//!<Atom name
    double		_Charge;		//!<Atom charge in electrons
    double		_Mass;			//!<Atom mass in Daltons
    uint		_TypeIndex;		//!<Type index
public:
	// Temporary variables, not necessary to store
	/*! Sets of all atoms that are bonded to this one at remove 0(bonded),
	  1(ends of angle) and 2(ends of dihedral) (indexed at 0, 1, 2 respectively */
    static const int max_remove = 2;
    gctools::SmallOrderedSet<Atom_sp> _AtomsAtRemoveBondAngle14[max_remove+1]; // s e t<Atom_sp>	_AtomsAtRemoveBondAngle14[3];
public:
//	void	archive( core::ArchiveP node);

	bool inBondOrAngle(Atom_sp a) { return (this->_AtomsAtRemoveBondAngle14[0].contains(a)) || (this->_AtomsAtRemoveBondAngle14[1].count(a)>0);};
	bool relatedBy14(Atom_sp a) { return (this->_AtomsAtRemoveBondAngle14[2].contains(a)); };

    EnergyAtom();
	EnergyAtom(Atom_sp atom, uint coordinateIndex);
	EnergyAtom(ForceField_sp forceField, Atom_sp atom, uint coordinateIndex);

		// methods
	void		defineForAtom(ForceField_sp forceField, Atom_sp atom, uint coordinateIndex);
	string		getResidueAndName();
//	adapt::QDomNode_sp	asXmlRelativeToContainer(chem::Matter_sp parent);
//	void		parseFromXmlRelativeToContainer( adapt::QDomNode_sp xml, chem::Matter_sp parent );
};


EnergyAtom*	_findEnergyAtom(gctools::Vec0<EnergyAtom>& atoms, uint	idx3 );

SMART(AtomTable);
class AtomTable_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,AtomTable_O,"AtomTable",core::CxxObject_O);

 public:
//    void	archiveBase(core::ArchiveP node);

 public:
  void initialize();
 private:
  typedef core::HashTableEq_sp        AtomTable;
  gctools::Vec0<EnergyAtom>	_Atoms;
  AtomTable        _AtomTableIndices; // m a p<Atom_sp,uint>	_AtomTableIndices;
 public:
  typedef gctools::Vec0<EnergyAtom>::iterator iterator;
 public:
  gctools::Vec0<EnergyAtom>&	getVectorEnergyAtoms() { return this->_Atoms;};
  CL_DEFMETHOD uint	getNumberOfAtoms()	{ return this->_Atoms.size();};
  uint	getNVectorSize()	{ return this->_Atoms.size()*3;};
  EnergyAtom*	getEnergyAtomPointer(Atom_sp a);
  EnergyAtom&	energyAtomEntry(uint i) { return this->_Atoms[i]; };
  EnergyAtom*	findEnergyAtomWithCoordinateIndex(uint i) { return &this->_Atoms[i/3]; };

  void	add(EnergyAtom& aa) { this->_Atoms.push_back(aa);};

	/*! Add the info for a single atom and return the coordinate index (index*3) of the atom */
  int addAtomInfo(Atom_sp atom, units::Quantity_sp charge, units::Quantity_sp mass, int typeIndex );

  iterator	begin() { return this->_Atoms.begin(); };
  iterator	end() { return this->_Atoms.end(); };

  void	dumpTerms();
  void constructFromMatter(Matter_sp matter, ForceField_sp forceField, core::T_sp activeAtoms );

  size_t push_back_excluded_atom_indices_and_sort(core::NativeVector_int_sp excludedAtomIndices, size_t atomIndex);
  /*! Calculate the excluded atom list with atom indices starting at 0.
      Atoms with no excluded atoms get a -1 in the excluded_atom_list.
      This is different from AMBER, where it starts counting at 1.
      To generate an AMBER excluded atom list, add 1 to each entry.
      Return Values(number_excluded_atoms,excluded_atom_list).
  */
  core::NativeVector_int_mv calculate_excluded_atom_list();

  // ------------------------------------------------------------
  //
  // Access fields in AtomTable elements directly
  //
  CL_DEFMETHOD MatterName elt_atom_name(int index) { return this->_Atoms[index]._AtomName; };
  CL_DEFMETHOD double elt_charge(int index) { return this->_Atoms[index]._Charge; };
  CL_DEFMETHOD double elt_mass(int index) { return this->_Atoms[index]._Mass; };
  CL_DEFMETHOD int elt_type_index(int index) { return this->_Atoms[index]._TypeIndex; };

  DEFAULT_CTOR_DTOR(AtomTable_O);
};


};


TRANSLATE(chem::AtomTable_O);
#endif //]
