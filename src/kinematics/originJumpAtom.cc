#define DEBUG_LEVEL_FULL



#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <cando/kinematics/symbolTable.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/originJumpAtom.h>

namespace kinematics
{

    void OriginJumpAtom::updateInternalCoords(Stub& stub,
					      bool const recursive,
					      AtomTree_sp at)
    {_OF();
	Stub start;
	if ( recursive )
	{
	    for (int childIdx=0; childIdx<this->_numberOfChildren(); childIdx++ )
	    {
		this->_child(childIdx).get()->updateInternalCoords(start,true,at);
	    }
	}
    }



    core::Symbol_sp OriginJumpAtom::typeSymbol() const
    {_OF();
	return _sym_origin;
    };




};
