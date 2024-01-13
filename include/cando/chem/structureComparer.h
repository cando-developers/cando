/*
    File: structureComparer.h
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
       
       

/*
 */

#ifndef StructureComparer_H
#define	StructureComparer_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
//#include	"conformation.h"
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/nVector.h>
#include <cando/geom/coordinateArray.h>



namespace       chem {

extern	string	StructureComparerEntry_Xml;

DEFINE_O_SMART_POINTERS(StructureComparerEntry);
DEFINE_O_SMART_POINTERS(StructureComparer);




/*! Compares one structure (the fixed one) to others (moveable)
 * Calculates the RMS and the matrix to transform the moveable ones
 * to the fixed one.
 */

class StructureComparer_O;
    typedef	gctools::smart_ptr<StructureComparer_O>	StructureComparer_sp;
SMART(StructureComparer);
class StructureComparer_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,StructureComparer_O,"StructureComparer",core::CxxObject_O);

public:
	void initialize();
public:
private:
	Matter_sp				_Matter;
    gctools::SmallOrderedSet<Atom_sp>				_AllAtoms;
    gctools::SmallOrderedSet<Atom_sp>				_SuperposeAtoms;
	geom::SimpleVectorCoordinate_sp			_FixedCoordinates;
private:	// do not archive
	Matrix					_TransformToFixed;



protected:
	void	initializeFixedCoordinates();

public:

		/* Set the Matter object that this StructureComparer
		 * compares conformations of
		 */
	void	setFixedMatter(Matter_sp agg);
	Matter_sp	getMatter() { return this->_Matter;};

	double	calculateRmsWithMatter(Matter_sp matter);
	Matrix	getTransformToFixed() { return this->_TransformToFixed; };

		//! Add a superposition atom
    void	superposeAtomsFromSet(gctools::SmallOrderedSet<Atom_sp> atoms);
		//! Set up to superpose all heavy atoms
	void	superposeAllHeavyAtoms();

    gctools::SmallOrderedSet<Atom_sp>::iterator begin_AllAtoms()
		{ return this->_AllAtoms.begin();};
    gctools::SmallOrderedSet<Atom_sp>::iterator end_AllAtoms()
		{ return this->_AllAtoms.end();};

    gctools::SmallOrderedSet<Atom_sp>::iterator begin_SuperposeAtoms()
		{ return this->_SuperposeAtoms.begin();};
    gctools::SmallOrderedSet<Atom_sp>::iterator end_SuperposeAtoms()
		{ return this->_SuperposeAtoms.end();};

	int	numberOfSuperposeAtoms()	{return this->_SuperposeAtoms.size();};







	DEFAULT_CTOR_DTOR(StructureComparer_O);
};



};


#endif
