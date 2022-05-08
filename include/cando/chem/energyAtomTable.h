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
#include <cando/chem/molecule.h>
#include <cando/geom/vector3.h>


namespace chem
{

SMART(Atom);
SMART(ForceField);
SMART(Matter);
SMART(FFNonbondDb);

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
  uint		        _TypeIndex;		//!<Type index
  uint                  _AtomicNumber;           //!<Atomic number
  size_t                _Flag;
public:
	// Temporary variables, not necessary to store
	/*! Sets of all atoms that are bonded to this one at remove 0(bonded),
	  1(ends of angle) and 2(ends of dihedral) (indexed at 0, 1, 2 respectively */
  static const int max_remove = 2;
  gctools::SmallOrderedSet<Atom_sp> _AtomsAtRemoveBondAngle14[max_remove+1]; // s e t<Atom_sp>	_AtomsAtRemoveBondAngle14[3];
public:

  bool inBondOrAngle(Atom_sp a) { return (this->_AtomsAtRemoveBondAngle14[0].contains(a)) || (this->_AtomsAtRemoveBondAngle14[1].count(a)>0);};
  bool relatedBy14(Atom_sp a) { return (this->_AtomsAtRemoveBondAngle14[2].contains(a)); };

  EnergyAtom();
  EnergyAtom(Atom_sp atom, uint coordinateIndex);
  EnergyAtom(core::T_sp forceField, Atom_sp atom, uint coordinateIndex);

		// methods
  core::List_sp encode() const;
  void		defineForAtom(core::T_sp forceField, Atom_sp atom, uint coordinateIndex);
  string		getResidueAndName(Residue_sp res);
//	adapt::QDomNode_sp	asXmlRelativeToContainer(chem::Matter_sp parent);
//	void		parseFromXmlRelativeToContainer( adapt::QDomNode_sp xml, chem::Matter_sp parent );
};
};

namespace translate {

template <>
struct	to_object<chem::EnergyAtom >
{
  typedef	core::Cons_sp ExpectedType;
  typedef	core::Cons_sp DeclareType;
  static core::T_sp convert(const chem::EnergyAtom& ea)
  {
    return ea.encode();
  }
};

template <>
struct	from_object<chem::EnergyAtom>
{
  typedef	chem::EnergyAtom	ExpectedType;
  typedef	ExpectedType 		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
          SIMPLE_ERROR(("Implement me"));
        }
};
};

namespace chem {

EnergyAtom*	_findEnergyAtom(gctools::Vec0<EnergyAtom>& atoms, uint	idx3 );

SMART(AtomTable);
class AtomTable_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,AtomTable_O,"AtomTable",core::CxxObject_O);

 public:
  bool fieldsp() const { return true; }
  void fields(core::Record_sp node);
 public:
  void initialize();
public:
  gctools::Vec0<EnergyAtom>	      _Atoms;
  core::HashTableEq_sp                _AtomTableIndices; // m a p<Atom_sp,uint>	_AtomTableIndices;
  core::ComplexVector_int32_t_sp      _ResiduePointers;
  core::ComplexVector_T_sp            _ResidueNames;
  //! Store count of atoms in each molecule. The length of this vector is the number of molecules
  core::ComplexVector_int32_t_sp      _AtomsPerMolecule;
  //! Store the index of the first "solvent" molecule
  core::T_sp                     _firstSolventMoleculeNSPSOL;
  core::T_sp                     _finalSoluteResidueIPTRES;
  core::T_sp                     _totalNumberOfMoleculesNSPM;
  //! Stores actual residues from aggregate
  core::Vector_sp                _Residues;
  core::Vector_sp                _Molecules;
  core::T_sp                     _Matter;
  core::T_sp                     _AggregateName;
  core::T_sp                     _BoundingBox;
  core::T_sp                     _NonbondForceFieldForAggregate;
 public:
  typedef gctools::Vec0<EnergyAtom>::iterator iterator;
public:
  static AtomTable_sp make();
  
 public:
  gctools::Vec0<EnergyAtom>&	getVectorEnergyAtoms() { return this->_Atoms;};
  CL_DEFMETHOD size_t	getNumberOfAtoms()	{ return this->_Atoms.size();};
  CL_DEFMETHOD size_t   getNumberOfMolecules()  { return this->_AtomsPerMolecule->length();};
  core::T_mv getAtomIndexOrNil(Atom_sp atom);

  core::List_sp atomsAtRemoveBondAngle14(Atom_sp atom) const;
  core::T_sp   firstSolventMoleculeNSPSOL() const;
  bool firstSolventMoleculeNSPSOL_boundP() const;
  void set_firstSolventMoleculeNSPSOL(size_t num);
  void makUnbound_firstSolventMoleculeNSPSOL();

  core::T_sp   finalSoluteResidueIPTRES() const;
  bool finalSoluteResidueIPTRES_boundP() const;
  void set_finalSoluteResidueIPTRES(size_t num);
  void makUnbound_finalSoluteResidueIPTRES();

  core::T_sp   totalNumberOfMoleculesNSPM() const;
  bool totalNumberOfMoleculesNSPM_boundP() const;
  void set_totalNumberOfMoleculesNSPM(size_t num);
  void makUnbound_totalNumberOfMoleculesNSPM();

  core::T_sp   aggregateName() const;
  bool aggregateNameBoundP() const;
  void setAggregateName(core::T_sp name);
  void makUnboundAggregateName();

  core::T_sp   boundingBox() const;
  bool boundingBoxBoundP() const;
  void setBoundingBox(core::T_sp name);
  void makUnboundBoundingBox();

  core::T_sp   nonbondForceFieldForAggregate() const;
  bool nonbondForceFieldForAggregateBoundP() const;
  void setNonbondForceFieldForAggregate(core::T_sp forceField);
  void makUnboundNonbondForceFieldForAggregate();

  uint	getNVectorSize()	{ return this->_Atoms.size()*3;};
  EnergyAtom*	getEnergyAtomPointer(Atom_sp a);
  EnergyAtom&	energyAtomEntry(uint i) { return this->_Atoms[i]; };
  EnergyAtom*	findEnergyAtomWithCoordinateIndex(uint i) { return &this->_Atoms[i/3]; };

  void	add(EnergyAtom& aa) { SIMPLE_ERROR(("AtomTable->add(EnergyAtom may cause recursive allocation")); this->_Atoms.push_back(aa);};

	/*! Add the info for a single atom and return the coordinate index (index*3) of the atom */
  int addAtomInfo(Atom_sp atom, units::Quantity_sp charge, units::Quantity_sp mass, int typeIndex, uint atomicNumber );

  iterator	begin() { return this->_Atoms.begin(); };
  iterator	end() { return this->_Atoms.end(); };

  void setAtomFlag(size_t index, size_t flag);
  size_t getAtomFlag(size_t index);

  void	dumpTerms();
  void constructFromMolecule(Molecule_sp mol, core::T_sp nonbondForceField, core::T_sp activeAtoms );

  size_t push_back_excluded_atom_indices_and_sort(core::ComplexVector_int32_t_sp excludedAtomIndices, size_t atomIndex);
  /*! Calculate the excluded atom list with atom indices starting at 0.
      Atoms with no excluded atoms get a -1 in the excluded_atom_list.
      This is different from AMBER, where it starts counting at 1.
      To generate an AMBER excluded atom list, add 1 to each entry.
      Return Values(number_excluded_atoms,excluded_atom_list).
  */
  CL_DEFMETHOD core::T_mv calculate_excluded_atom_list();

  // ------------------------------------------------------------
  //
  // Access fields in AtomTable elements directly
  //
  core::HashTableEq_sp getAtomTableIndices();
  CL_DEFMETHOD MatterName elt_atom_name(int index) { return this->_Atoms[index]._AtomName; };
  CL_DEFMETHOD core::Symbol_sp elt_atom_type(int index);
  CL_DEFMETHOD Atom_sp elt_atom(int index) { return this->_Atoms[index].atom(); };
  CL_DEFMETHOD int elt_atom_coordinate_index_times3(int index) { return this->_Atoms[index].coordinateIndexTimes3(); };
  CL_DEFMETHOD double elt_charge(int index) { return this->_Atoms[index]._Charge; };
  CL_DEFMETHOD double elt_mass(int index) { return this->_Atoms[index]._Mass; };
  CL_DEFMETHOD int elt_type_index(int index) { return this->_Atoms[index]._TypeIndex; };
  CL_DEFMETHOD int elt_atomic_number(int index) { return this->_Atoms[index]._AtomicNumber; };
  core::T_sp atom_table_residue_pointers() const;
  core::T_sp atom_table_residue_names() const;
  core::T_sp atom_table_atoms_per_molecule() const;
  core::T_sp atom_table_residues() const;
  core::T_sp atom_table_molecules() const;
  CL_DEFMETHOD void setf_atom_table_residue_pointers(core::ComplexVector_int32_t_sp val) {this->_ResiduePointers = val; };
  CL_DEFMETHOD void setf_atom_table_residue_names(core::ComplexVector_T_sp val) {this->_ResidueNames = val; };
  CL_DEFMETHOD void setf_atom_table_atoms_per_molecule(core::ComplexVector_int32_t_sp val) {this->_AtomsPerMolecule = val; };
  CL_DEFMETHOD void setf_atom_table_residues(core::Vector_sp val) {this->_Residues = val; };
  CL_DEFMETHOD void setf_atom_table_molecules(core::Vector_sp val) {this->_Molecules = val; };
  size_t getCoordinateIndex(Atom_sp atom);
  size_t getCoordinateIndexForAtomAtIndex(size_t index);
  
 AtomTable_O() : _ResiduePointers(unbound<core::ComplexVector_int32_t_O>()),
                 _ResidueNames(unbound<core::ComplexVector_T_O>()),
                 _AtomsPerMolecule(unbound<core::ComplexVector_int32_t_O>()),
                 _firstSolventMoleculeNSPSOL(unbound<core::T_O>()),
                 _finalSoluteResidueIPTRES(unbound<core::T_O>()),
                 _totalNumberOfMoleculesNSPM(unbound<core::T_O>()),
                 _Residues(unbound<core::Vector_O>()),
                 _Molecules(unbound<core::Vector_O>()),
                 _AggregateName(unbound<core::T_O>()),
                 _BoundingBox(unbound<core::T_O>()),
                 _NonbondForceFieldForAggregate(unbound<core::T_O>()) {};

  virtual void fill_atom_table_from_vectors(core::List_sp values);
//  int residue_index(int atom_index);
};


};


TRANSLATE(chem::AtomTable_O);
#endif //]
