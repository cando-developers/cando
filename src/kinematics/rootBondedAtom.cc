#define DEBUG_LEVEL_FULL



#include "core/foundation.h"
#include "core/object.h"
#include "core/lisp.h"
#include "symbolTable.h"
#include "rootBondedAtom.h"

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
