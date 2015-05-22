#ifndef _kinematics_rootAtomInfo_H
#define _kinematics_rootAtomInfo_H

#include "core/common.h"
#include "chem/plug.fwd.h"
#include "chem/constitutionAtoms.fwd.h"

namespace kinematics
{

    /*! @class Stores information on the Residue that the containing Atom is the root of */

    class RootAtomInfo
    {
	friend class RootBondedAtom;
	friend class DelayedBondedAtom;
    private:
	core::Symbol_sp			_ConstitutionName;
	core::Symbol_sp			_TopologyName;
	chem::ConstitutionAtomIndex0N	_Bond1Id;
    public:
	void setup(core::Symbol_sp constitutionName,
		   core::Symbol_sp topologyName,
		   chem::Plug_sp inPlug);
    };


};

#endif // _kinematics_rootAtomInfo_H
