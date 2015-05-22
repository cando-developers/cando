#ifndef	kinematics_dofId_H
#define	kinematics_dofId_H

#include "core/foundation.h"
#include "addon/vector3.h"

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
