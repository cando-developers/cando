       
       

/*
 */

#ifndef StructureComparer_H
#define	StructureComparer_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "bond.h"
#include "addon/vector3.h"
//#include	"conformation.h"
#include "atom.h"
#include "residue.h"
#include "nVector.h"
#include "addon/coordinateArray.h"



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
class StructureComparer_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,StructureComparer_O,"StructureComparer");

public:
	void initialize();
public:
//	void	archive(core::ArchiveP node);
private:
	Matter_sp				_Matter;
    gctools::SmallOrderedSet<Atom_sp>				_AllAtoms;
    gctools::SmallOrderedSet<Atom_sp>				_SuperposeAtoms;
	addon::CoordinateArray_sp			_FixedCoordinates;
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


TRANSLATE(chem::StructureComparer_O);
#endif
