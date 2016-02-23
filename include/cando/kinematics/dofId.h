/*
    File: dofId.h
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
#ifndef	kinematics_dofId_H
#define	kinematics_dofId_H

#include <clasp/core/foundation.h>
#include <cando/geom/vector3.h>

namespace kinematics
{

class DofId
{
private:
    DofType		_Type;
    chem::AtomId	_AtomId;

public:
DofId() : _Type(undefined), _AtomId() {};
DofId(chem::AtomId const& atomId, DofType const& aType) : _Type(aType), _AtomId(atomId) {};

    inline chem::AtomId const& atomId() const { return _AtomId;};
    inline int moleculeId() const { return _AtomId.moleculeId();};
    inline int residueId() const { return _AtomId.residueId();};
    inline chem::ConstitutionAtomIndex0N moleculeId() const { return _AtomId.atomId();};
    inline DofType type() const { return _Type;};

public:
    inline bool operator==(DofId const& other)
    {
	return this->_AtomId == b._AtomId && this->_Type == b._Type;
    };

    inline bool operator!=(DofId const& other)
    {
	return this->_AtomId != b._AtomId || this->_Type != b._Type;
    }

    inline bool operator <(DofId const& other)
    {
	return ( (this._AtomId < other.AtomId) ||
		 (this._AtomId == other._AtomId && this._Type == b._Type));
    }
};

};


#endif // kinematics_dofId_H
