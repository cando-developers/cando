#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/kinematics/bondId.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(kinematics,BondId_O);

    void BondId_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<BondId_O>()
//	.def_raw("__init__",&BondId_O::__init__,"(self)")
	;
    }

    void BondId_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,BondId,"","",_lisp)
//	.def_raw("__init__",&"(self)")
	;
#endif
    }


    BondId_sp BondId_O::create(Atom_sp parent, Atom_sp child)
    {_G();
	GC_ALLOCATE(BondId_O, bondId );
	bondId->_Parent = parent;
	bondId->_Child = child;
	return bondId;
    }

#if INIT_TO_FACTORIES

  #define ARGS_BondId_O_make "()"
#define DECL_BondId_O_make ""
#define DOCS_BondId_O_make "make BondId_O"
  BondId_sp BondId_O::make()
  {_G();
    IMPLEMENT_ME();
  };

#else
    core::T_sp BondId_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
//      arg = translate::from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void BondId_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }
#endif


    void BondId_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    

}; /* kinematics */
