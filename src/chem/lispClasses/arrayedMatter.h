/*
    File: arrayedMatter.h
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
