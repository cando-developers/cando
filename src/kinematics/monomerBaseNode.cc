#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/kinematics/monomerBaseNode.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//





#if INIT_TO_FACTORIES

#define ARGS_MonomerBaseNode_O_make "()"
#define DECL_MonomerBaseNode_O_make ""
#define DOCS_MonomerBaseNode_O_make "make MonomerBaseNode"
  MonomerBaseNode_sp MonomerBaseNode_O::make()
    {
      IMPLEMENT_ME();
    };

#else
    core::T_sp MonomerBaseNode_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
//      this->Base::__init__(exec,args,env,lisp);
//      arg = translate::from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void MonomerBaseNode_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }
#endif


    void MonomerBaseNode_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    

}; /* kinematics */
