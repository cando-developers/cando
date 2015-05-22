#define DEBUG_LEVEL_FULL



#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <cando/kinematics/symbolTable.h>
#include <cando/kinematics/rootBondedAtom.h>

namespace kinematics
{




    void RootBondedAtom::setup(core::Symbol_sp constitutionName,
			       core::Symbol_sp topologyName,
			       chem::Plug_sp inPlug)
    {
	this->_RootInfo.setup(constitutionName,topologyName,inPlug);
    }



    core::Symbol_sp RootBondedAtom::typeSymbol() const
    {_OF();
	return _sym_rootBonded;
    };
};
