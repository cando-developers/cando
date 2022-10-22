/*
    File: iterateMatter.cc
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#define	DEBUG_LEVEL_NONE
#include <clasp/core/common.h>
#include <cando/chem/iterateMatter.h>
#include <clasp/core/symbolTable.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/wrappers.h>


namespace chem
{


    SYMBOL_SC_(ChemPkg,iterateMatterSymbolConverter);






#if INIT_TO_FACTORIES

#define ARGS_IterateMatter_O_make "()"
#define DECL_IterateMatter_O_make ""
#define DOCS_IterateMatter_O_make "make IterateMatter"
  IterateMatter_sp IterateMatter_O::make()
  {
    IMPLEMENT_ME();
  };

#else

    core::T_sp  IterateMatter_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::LispPtr lisp)
    {
	IMPLEMENT_ME();
//	this->Base::oldLispInitialize(kargs,env);
	// nothing
    }
#endif
    IterateMatter_sp IterateMatter_O::create(Matter_sp top, int goal)
    {
      auto  m  = gctools::GC<IterateMatter_O>::allocate_with_default_constructor();
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



CL_LISPIFY_NAME("initTopAndGoal");
CL_DEFMETHOD     void IterateMatter_O::initTopAndGoal(Matter_sp top, int goal)
    {
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

CL_LISPIFY_NAME("advance");
CL_DEFMETHOD     bool	IterateMatter_O::advance()
    {
	if ( !this->_IsDone )
	{
	    bool d = this->_Loop.advance();
	    if ( !d ) 
	    {
		LOG("Advanced loop and hit end" );
		this->_IsDone = true;
	    }
	}
	return this->_IsDone;
    }

    core::T_sp	IterateMatter_O::currentObject()
    {
	if ( this->_IsDone ) return nil<Matter_O>();
	Matter_sp res = this->_Loop.getMatter();
	ASSERTNOTNULL(res);
	return res;
    }

    bool	IterateMatter_O::isDone()
    {
	return this->_IsDone;
    }








    IterateAtoms_sp IterateAtoms_O::create(Matter_sp top)
    {
      auto  ia  = gctools::GC<IterateAtoms_O>::allocate_with_default_constructor();
	ia->initTopAndGoal(top,ATOMS);
	return ia;
    }



    
    
#define ARGS_chem__create_for_matter "(matter)"
#define DECL_chem__create_for_matter ""
#define DOCS_chem__create_for_matter "createForMatter"
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__create_for_matter(Matter_sp matter)
    {
      auto  ia  = gctools::GC<IterateAtoms_O>::allocate_with_default_constructor();
	ia->initTopAndGoal(matter,ATOMS);
	return ia;
    }

#if INIT_TO_FACTORIES

#define ARGS_IterateAtoms_O_make "(matter)"
#define DECL_IterateAtoms_O_make ""
#define DOCS_IterateAtoms_O_make "make IterateAtoms"
  IterateAtoms_sp IterateAtoms_O::make(Matter_sp matter)
  {
    auto  me  = gctools::GC<IterateAtoms_O>::allocate_with_default_constructor();
    me->initTopAndGoal(matter,ATOMS);
    return me;
  };

#else

    core::T_sp 	IterateAtoms_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::LispPtr lisp)
    {
	this->Base::__init__(exec,args,env,lisp);
	Matter_sp matter = env->lookup(lisp->internWithPackageName(ChemPkg,"matter")).as<Matter_O>();
	this->initTopAndGoal(matter,ATOMS);
	return _Nil<core::T_O>();
    }

#endif

    core::T_sp	IterateAtoms_O::currentObject()
    {
	if ( this->isDone() ) return nil<core::T_O>();
	return this->_Loop.getAtom();
    }





    IterateResidues_sp IterateResidues_O::create(Matter_sp top)
    {
      auto  ia  = gctools::GC<IterateResidues_O>::allocate_with_default_constructor();
	ia->initTopAndGoal(top,RESIDUES);
	return ia;
    }


    
    
    
#define ARGS_IterateResidues_O_createForMatter "(matter)"
#define DECL_IterateResidues_O_createForMatter ""
#define DOCS_IterateResidues_O_createForMatter "IterateResidues_O_createForMatter"
CL_LISPIFY_NAME(residues);
DOCGROUP(cando);
CL_DEFUN core::T_sp IterateResidues_O::createForMatter(Matter_sp matter)
    {
	return IterateResidues_O::create(matter);
    }


    core::T_sp	IterateResidues_O::currentObject()
    {
	if ( this->isDone() ) return nil<core::T_O>();
	return this->_Loop.getResidue();
    }



CL_LISPIFY_NAME(bonds);
DOCGROUP(cando);
CL_DEFUN IterateBonds_sp IterateBonds_O::make(Matter_sp top)
{
  auto  ia  = gctools::GC<IterateBonds_O>::allocate_with_default_constructor();
  ia->initTopAndGoal(top,BONDS);
  return ia;
}




    
    
#define ARGS_IterateBonds_O_createForMatter "(matter)"
#define DECL_IterateBonds_O_createForMatter ""
#define DOCS_IterateBonds_O_createForMatter "IterateBonds_O_createForMatter"
    core::T_sp IterateBonds_O::createForMatter(Matter_sp matter)
    {
	return IterateBonds_O::make(matter);
    }



    core::T_sp	IterateBonds_O::currentObject()
    {
	return this->sharedThis<IterateBonds_O>();
    }


CL_LISPIFY_NAME("getAtom1");
CL_DEFMETHOD     Atom_sp	IterateBonds_O::getAtom1()
    {
	return this->_Loop.getAtom1();
    }

CL_LISPIFY_NAME("getAtom2");
CL_DEFMETHOD     Atom_sp	IterateBonds_O::getAtom2()
    {
	return this->_Loop.getAtom2();
    }

CL_LISPIFY_NAME("getBondOrder");
CL_DEFMETHOD     int IterateBonds_O::getBondOrder()
    {
	return this->_Loop.getBondOrder();
    }


SYMBOL_EXPORT_SC_(ChemPkg,allMatter);
SYMBOL_EXPORT_SC_(ChemPkg,atoms);
SYMBOL_EXPORT_SC_(ChemPkg,residues);
SYMBOL_EXPORT_SC_(ChemPkg,molecules);
SYMBOL_EXPORT_SC_(ChemPkg,bonds);
CL_BEGIN_ENUM(LoopEnum,_sym_iterateMatterSymbolConverter,"IterateMatter");
CL_VALUE_ENUM(_sym_atoms,atoms);
CL_VALUE_ENUM(_sym_residues,residues);
CL_VALUE_ENUM(_sym_molecules,molecules);
CL_VALUE_ENUM(_sym_bonds,bonds);
CL_END_ENUM(_sym_iterateMatterSymbolConverter);






    























};
