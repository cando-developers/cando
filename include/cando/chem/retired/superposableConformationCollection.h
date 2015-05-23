       
       
#ifndef	SuperposableConformationCollection_H //[
#define SuperposableConformationCollection_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/conformationCollection.h>
#include <clasp/core/stringList.h>
#include <clasp/core/intArray.h>

#include <cando/chem/chemPackage.h>


namespace chem {



SMART(ZMatrix);




class SuperposableConformationCollection_O;
    typedef	gctools::smart_ptr<SuperposableConformationCollection_O>	SuperposableConformationCollection_sp;
SMART(SuperposableConformationCollection);
class SuperposableConformationCollection_O : public ConformationCollection_O
{
    LISP_BASE1(ConformationCollection_O);
    LISP_CLASS(chem,ChemPkg,SuperposableConformationCollection_O,"SuperposableConformationCollection");

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
	geom::Render_sp rendered(core::Cons_sp options);
#endif


	SuperposableConformationCollection_O( const SuperposableConformationCollection_O& ss ); //!< Copy constructor


	DEFAULT_CTOR_DTOR(SuperposableConformationCollection_O);
};



// Constructor function


};
TRANSLATE(chem::SuperposableConformationCollection_O);
#endif //]
