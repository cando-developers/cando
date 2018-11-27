/*
    File: indirectAtomCoordinateReference.h
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
#ifndef	indirectAtomCoordinate_H
#define indirectAtomCoordinate_H


#include <cando/chem/chemPackage.h>


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
    public:
	Atom_sp	_SharedAtom;
	uint	_IndexTimes3;
    public:
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
