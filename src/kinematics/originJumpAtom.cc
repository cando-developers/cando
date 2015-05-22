#define DEBUG_LEVEL_FULL



#include "core/foundation.h"
#include "core/object.h"
#include "core/lisp.h"
#include "symbolTable.h"
#include "stub.h"
#include "originJumpAtom.h"

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
