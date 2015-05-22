#ifndef	ArrayedMolecule_H //[
#define ArrayedMolecule_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "core/holder.h"

#include "chemPackage.h"

#include "candoBase/coordinateArray.fwd.h"

namespace chem {

SMART(Atom);
SMART(Residue);
SMART(Molecule);
SMART(Matter);
SMART(ArrayedAtom);
SMART(CoordinateArrayWithHash);


SMART(ArrayedAtom);
class ArrayedAtom_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,ArrayedAtom_O,"ArrayedAtom");
friend class OrderByAtomLabel;
public: // virtual functions inherited from Object
	void	initialize();
//	void	archive(core::ArchiveP node);
//	string	__repr__() const;

protected: // instance variables
    MatterName		_AtomLabel;
    Atom_sp		_Atom;

public:	// Creation class functions

public:
	uint calculateHash();

public:
	ArrayedAtom_O( const ArrayedAtom_O& ss ); //!< Copy constructor
	void extractFromCoordinateArray(candoBase::CoordinateArrayWithHash_sp h, uint start);
	void writeToCoordinateArray(candoBase::CoordinateArrayWithHash_sp h, uint start);
public:
	void setAtom(Atom_sp a);
	DEFAULT_CTOR_DTOR(ArrayedAtom_O);
};



SMART(ArrayedResidue);
class ArrayedResidue_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,ArrayedResidue_O,"ArrayedResidue");
friend class OrderByResidueLabel;

public: // virtual functions inherited from Object
	void	initialize();
//	void	archive(core::ArchiveP node);
//	string	__repr__() const;
protected: // instance variables
	string			_ResidueLabel;
	Residue_wp		_Residue;
    gctools::Vec0<ArrayedAtom_sp>	_SortedAtoms;

public:	// Creation class functions

public:
	uint calculateHash();
	uint extractFromCoordinateArray(candoBase::CoordinateArrayWithHash_sp h, uint start);
	uint writeToCoordinateArray(candoBase::CoordinateArrayWithHash_sp h, uint start);
public:
	ArrayedResidue_O( const ArrayedResidue_O& ss ); //!< Copy constructor
public:
	void setResidue(Residue_sp res);
	uint	getNumberOfAtoms();

	DEFAULT_CTOR_DTOR(ArrayedResidue_O);
};





SMART(ArrayedMolecule);
class ArrayedMolecule_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,ArrayedMolecule_O,"ArrayedMolecule");
    DECLARE_INIT();
public: // virtual functions inherited from Object
	void	initialize();
//	void	archive(core::ArchiveP node);
//	string	__repr__() const;

private: // instance variables
	Molecule_wp		_Molecule;
	uint			_NumberOfAtoms;
	uint			_Hash;
    gctools::Vec0<ArrayedResidue_sp>	_SortedResidues;

public:	// Creation class functions

public:
	ArrayedMolecule_O( const ArrayedMolecule_O& ss ); //!< Copy constructor
public:
	void setMolecule(Molecule_sp mol);

	uint calculateHash();

	uint	getNumberOfAtoms();

	candoBase::CoordinateArrayWithHash_sp	createCoordinateArray();
	void extractFromCoordinateArray(candoBase::CoordinateArrayWithHash_sp h);
	void writeToCoordinateArray(candoBase::CoordinateArrayWithHash_sp h);

	DEFAULT_CTOR_DTOR(ArrayedMolecule_O);
};


};
TRANSLATE(chem::ArrayedAtom_O);
TRANSLATE(chem::ArrayedResidue_O);
TRANSLATE(chem::ArrayedMolecule_O);
#endif //]
