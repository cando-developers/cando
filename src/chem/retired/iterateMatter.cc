#define	DEBUG_LEVEL_FULL
#include <clasp/core/common.h>
#include <cando/chem/iterateMatter.h>
#include <cando/chem/symbolTable.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/wrappers.h>


namespace chem
{


    SYMBOL_SC_(ChemPkg,iterateMatterSymbolConverter);

    void IterateMatter_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<IterateMatter_O>()
	    .def("initTopAndGoal",&IterateMatter_O::initTopAndGoal)
	    .def("advance",&IterateMatter_O::advance)
	    ;
    }

    void IterateMatter_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,IterateMatter,"","",_lisp)
	    .def("initTopAndGoal",&IterateMatter_O::initTopAndGoal)
	    .def("advance",&IterateMatter_O::advance)
	    ;
#endif
    }




#if INIT_TO_FACTORIES

#define ARGS_IterateMatter_O_make "()"
#define DECL_IterateMatter_O_make ""
#define DOCS_IterateMatter_O_make "make IterateMatter"
  IterateMatter_sp IterateMatter_O::make()
  {_G();
    IMPLEMENT_ME();
  };

#else

    core::T_sp  IterateMatter_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
	IMPLEMENT_ME();
//	this->Base::oldLispInitialize(kargs,env);
	// nothing
    }
#endif
    IterateMatter_sp IterateMatter_O::create(core::Lisp_sp e,Matter_sp top, int goal)
    {
	GC_ALLOCATE(IterateMatter_O, m );
	m->_Top = top;
	m->_Goal = goal;
	return m;
    };



    void	IterateMatter_O::initialize()
    {
	this->Base::initialize();
    }

#ifdef XML_ARCHIVE
    void	IterateMatter_O::archive(core::ArchiveP node)
    {
	IMPLEMENT_ME();
    }
#endif



    void IterateMatter_O::initTopAndGoal(Matter_sp top, int goal)
    {_G();
	ANN(top);
	this->_Top = top;
	if ( goal == ALL_MATTERS )
	{
	    goal = MOLECULES|RESIDUES|ATOMS;
	}
	this->_Goal = goal;

    }

    void	IterateMatter_O::first()
    {
	this->_IsDone = false;
	this->_Loop.loopTopGoal(this->_Top,this->_Goal);
	this->next();
    }

    void	IterateMatter_O::next()
    {
	this->advance();
    };

    bool	IterateMatter_O::advance()
    {_G();
	if ( !this->_IsDone )
	{
	    bool d = this->_Loop.advance();
	    if ( !d ) 
	    {
		LOG(BF("Advanced loop and hit end") );
		this->_IsDone = true;
	    }
	}
	return this->_IsDone;
    }

    core::T_sp	IterateMatter_O::currentObject()
    {_OF();
	if ( this->_IsDone ) return _Nil<Matter_O>();
	Matter_sp res = this->_Loop.getMatter();
	ASSERTNOTNULL(res);
	return res;
    }

    bool	IterateMatter_O::isDone()
    {
	return this->_IsDone;
    }








    IterateAtoms_sp IterateAtoms_O::create(Matter_sp top, core::Lisp_sp e)
    {_G();
	GC_ALLOCATE(IterateAtoms_O, ia );
	ia->initTopAndGoal(top,ATOMS);
	return ia;
    }



    
    
#define ARGS_af_createForMatter "(matter)"
#define DECL_af_createForMatter ""
#define DOCS_af_createForMatter "createForMatter"
    core::T_sp af_createForMatter(Matter_sp matter)
    {_G();
	GC_ALLOCATE(IterateAtoms_O, ia );
	ia->initTopAndGoal(matter,ATOMS);
	return ia;
    }

#if INIT_TO_FACTORIES

#define ARGS_IterateAtoms_O_make "(matter)"
#define DECL_IterateAtoms_O_make ""
#define DOCS_IterateAtoms_O_make "make IterateAtoms"
  IterateAtoms_sp IterateAtoms_O::make(Matter_sp matter)
  {_G();
      GC_ALLOCATE(IterateAtoms_O, me );
    me->initTopAndGoal(matter,ATOMS);
    return me;
  };

#else

    core::T_sp 	IterateAtoms_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
	this->Base::__init__(exec,args,env,lisp);
	Matter_sp matter = env->lookup(lisp->internWithPackageName(ChemPkg,"matter")).as<Matter_O>();
	this->initTopAndGoal(matter,ATOMS);
	return _Nil<core::T_O>();
    }

#endif

    core::T_sp	IterateAtoms_O::currentObject()
    {_G();
	if ( this->isDone() ) return _Nil<core::T_O>();
	return this->_Loop.getAtom();
    }





    IterateResidues_sp IterateResidues_O::create(core::Lisp_sp e,Matter_sp top)
    {_G();
	GC_ALLOCATE(IterateResidues_O, ia );
	ia->initTopAndGoal(top,RESIDUES);
	return ia;
    }


    
    
    
#define ARGS_IterateResidues_O_createForMatter "(matter)"
#define DECL_IterateResidues_O_createForMatter ""
#define DOCS_IterateResidues_O_createForMatter "IterateResidues_O_createForMatter"
    core::T_sp IterateResidues_O::createForMatter(Matter_sp matter)
    {_G();
	return IterateResidues_O::create(_lisp,matter);
    }


    core::T_sp	IterateResidues_O::currentObject()
    {_G();
	if ( this->isDone() ) return _Nil<core::T_O>();
	return this->_Loop.getResidue();
    }




    IterateBonds_sp IterateBonds_O::create(core::Lisp_sp e,Matter_sp top)
    {_G();
	GC_ALLOCATE(IterateBonds_O, ia );
	ia->initTopAndGoal(top,BONDS);
	return ia;
    }




    
    
#define ARGS_IterateBonds_O_createForMatter "(matter)"
#define DECL_IterateBonds_O_createForMatter ""
#define DOCS_IterateBonds_O_createForMatter "IterateBonds_O_createForMatter"
    core::T_sp IterateBonds_O::createForMatter(Matter_sp matter)
    {_G();
	return IterateBonds_O::create(_lisp,matter);
    }



    core::T_sp	IterateBonds_O::currentObject()
    {_G();
	return this->sharedThis<IterateBonds_O>();
    }


    Atom_sp	IterateBonds_O::getAtom1()
    {_G();
	return this->_Loop.getAtom1();
    }

    Atom_sp	IterateBonds_O::getAtom2()
    {_G();
	return this->_Loop.getAtom2();
    }

    int IterateBonds_O::getBondOrder()
    {_G();
	return this->_Loop.getBondOrder();
    }



    void IterateMatter_O::lisp_initGlobals(core::Lisp_sp lisp)
    {_G();
	SYMBOL_EXPORT_SC_(ChemPkg,allMatter);
	SYMBOL_EXPORT_SC_(ChemPkg,atoms);
	SYMBOL_EXPORT_SC_(ChemPkg,residues);
	SYMBOL_EXPORT_SC_(ChemPkg,molecules);
	SYMBOL_EXPORT_SC_(ChemPkg,bonds);
	core::enum_<int>(_sym_iterateMatterSymbolConverter,"IterateMatter")
	    .value(_sym_allMatter,ALL_MATTERS)	// MOLECULES/RESIDUES/ATOMS
	    .value(_sym_atoms,ATOMS)
	    .value(_sym_residues,RESIDUES)
	    .value(_sym_molecules,MOLECULES)
	    .value(_sym_bonds,BONDS)
	    ;
    }





    
    void IterateAtoms_O::exposeCando(core::Lisp_sp lisp)
	{
	    core::class_<IterateAtoms_O>()
//		.def_raw("core:__init__",&IterateAtoms_O::__init__,"(self matter)")
		;
	    Defun(createForMatter);
	}

    void IterateAtoms_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,IterateAtoms,"","",_lisp)
	    ;
#endif //]
    }







    void IterateResidues_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<IterateResidues_O>()
	    ;
	af_def(ChemPkg,"residues",&IterateResidues_O::createForMatter,
	       ARGS_IterateResidues_O_createForMatter,
	       DECL_IterateResidues_O_createForMatter,
	       DOCS_IterateResidues_O_createForMatter );
    }

    void IterateResidues_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,IterateResidues,"","",_lisp)
	    ;
#endif //]
    }





    void IterateBonds_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<IterateBonds_O>()
	    .def("getAtom1",&IterateBonds_O::getAtom1)
	    .def("getAtom2",&IterateBonds_O::getAtom2)
	    .def("getBondOrder",&IterateBonds_O::getBondOrder)
	    ;
	core::af_def(ChemPkg,"bonds", &IterateBonds_O::createForMatter,
		     ARGS_IterateBonds_O_createForMatter,
		     DECL_IterateBonds_O_createForMatter,
		     DOCS_IterateBonds_O_createForMatter );
    }

    void IterateBonds_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,IterateBonds,"","",_lisp)
	    ;
#endif //]
    }






    EXPOSE_CLASS_AND_GLOBALS(chem,IterateMatter_O);
    EXPOSE_CLASS(chem,IterateAtoms_O);
    EXPOSE_CLASS(chem,IterateResidues_O);
    EXPOSE_CLASS(chem,IterateBonds_O);


};




