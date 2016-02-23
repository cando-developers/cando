/*
    File: superposableConformationCollection.h
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
       
       
#ifndef	SuperposableConformationCollection_H //[
#define SuperposableConformationCollection_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/conformationCollection.h>
#include <cando/adapt/stringList.h>
#include <clasp/core/intArray.h>

#include <cando/chem/chemPackage.h>


namespace chem {



SMART(ZMatrix);




class SuperposableConformationCollection_O;
    typedef	gctools::smart_ptr<SuperposableConformationCollection_O>	SuperposableConformationCollection_sp;
SMART(SuperposableConformationCollection);
class SuperposableConformationCollection_O : public ConformationCollection_O
{
    LISP_CLASS(chem,ChemPkg,SuperposableConformationCollection_O,"SuperposableConformationCollection",ConformationCollection_O);

public:
	void initialize();
public:
//	void	archive(core::ArchiveP node);

private:
	double		_RmsCutOff;
	core::IntArray_sp	_SuperposeAtomIndices;
//	string		_SuperposableFileName;
public:
	typedef	core::IntArray_O::iterator	superposeAtomIndexIterator;
public:

	void	setRmsCutOff(double co) { this->_RmsCutOff = co; };

		//! Clear the superposition atoms
	void	clearSuperposeAtoms();
		//! Add a superposition atom
    void	superposeAtomsFromSet(gctools::SmallOrderedSet<Atom_sp> atoms);
		//! Set up to superpose all heavy atoms
	void	superposeAllHeavyAtoms();
	void	addSuperposeAtom(Atom_sp a);

#if 0 // Figure out how to create an iterator that will let us do this

    gctools::SmallOrderedSet<Atom_sp>::iterator begin_SuperposeAtoms()
		{ return this->_SuperposeAtoms.begin();};
    gctools::SmallOrderedSet<Atom_sp>::iterator end_SuperposeAtoms()
		{ return this->_SuperposeAtoms.end();};
#endif


    gctools::SmallOrderedSet<Atom_sp>	getSuperposeAtoms();

	int	numberOfSuperposeAtoms();

	ConformationCollectionEntry_sp createEntryIfConformationIsNew(Matter_sp matter);

	void	setEntryCoordinatesAsFixedWithinSuperposeEngine(ConformationCollectionEntry_sp entry, SuperposeEngine_sp sup);
	void	setEntryCoordinatesAsMoveableWithinSuperposeEngine(ConformationCollectionEntry_sp entry, SuperposeEngine_sp sup);

	bool canRender() { return true; }
#ifdef RENDER
	geom::Render_sp rendered(core::List_sp options);
#endif


	SuperposableConformationCollection_O( const SuperposableConformationCollection_O& ss ); //!< Copy constructor


	DEFAULT_CTOR_DTOR(SuperposableConformationCollection_O);
};



// Constructor function


};
TRANSLATE(chem::SuperposableConformationCollection_O);
#endif //]
