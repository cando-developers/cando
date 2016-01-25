#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/kinematics/ringClosingMonomerNode.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//






    RingClosingMonomerNode_sp RingClosingMonomerNode_O::create(int monomerId)
    {
	GC_ALLOCATE(RingClosingMonomerNode_O, monomerNode );
	monomerNode->_MonomerId = monomerId;
	return monomerNode;
    }

#if INIT_TO_FACTORIES

#define ARGS_RingClosingMonomerNode_O_make "()"
#define DECL_RingClosingMonomerNode_O_make ""
#define DOCS_RingClosingMonomerNode_O_make "make RingClosingMonomerNode"
  RingClosingMonomerNode_sp RingClosingMonomerNode_O::make()
    {
      IMPLEMENT_ME();
    };

#else
    core::T_sp RingClosingMonomerNode_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
//      this->Base::__init__(exec,args,env,lisp);
//      arg = translate::from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void RingClosingMonomerNode_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }
#endif


    void RingClosingMonomerNode_O::initialize()
    {_OF();
        this->Base::initialize();
    }



    

}; /* kinematics */
