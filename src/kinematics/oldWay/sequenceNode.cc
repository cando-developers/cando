#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "sequenceNode.h"
#include "core/wrappers.h"
namespace kinematics
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(SequenceNode_O);

    void SequenceNode_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<SequenceNode_O>()
//	.initArgs("(self)")
	;
    }


    Object_sp SequenceNode_O::__init__(core::Executable_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
//      arg = from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _lisp->onil();
    }

    void SequenceNode_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }


    void SequenceNode_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    

}; /* kinematics */
