#ifndef	indirectAtomCoordinate_H
#define indirectAtomCoordinate_H


#include "chemPackage.h"


namespace chem
{

    SMART(Atom);
/*! Often we want to put a large number of atomic coordinates into a contiguous array
  and access their coordinates indirectly. This base class references an Atom and in index into a
  array of real numbers. The index is * 3 - so the first atom coordinates are at 0,1,2
  the second atom coordinates are at 3,4,5 etc.
*/

    class IndirectAtomCoordinateReference
    {
    private:
	Atom_sp	_SharedAtom;
	uint	_IndexTimes3;
    public:
//	void archive(core::ArchiveP node);
    protected:
	void setupBase(const Atom_sp& atom, uint indexTimes3)
	{
	    this->_SharedAtom = atom;
	    this->_IndexTimes3 = indexTimes3;
	};

    public:
	IndirectAtomCoordinateReference();
    IndirectAtomCoordinateReference(const Atom_sp& atom, uint indexTimes3 ) : _SharedAtom(atom), _IndexTimes3(indexTimes3) {};

	const Atom_sp&	atom() const { return this->_SharedAtom;};
	const uint& coordinateIndexTimes3() const {return this->_IndexTimes3;};

    };

};
#endif // indirectAtomCoordinate
