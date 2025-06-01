/*
    File: complexRestraints.cc
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

#include <cando/chem/complexRestraints.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/array.h>
#include <cando/chem/restraint.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>
#include <cando/chem/chemInfo.h>
#
// last include is wrappers.h
#include <clasp/core/wrappers.h>


namespace chem {












#define ARGS_ComplexRestraint_O_make "()"
#define DECL_ComplexRestraint_O_make ""
#define DOCS_ComplexRestraint_O_make "make ComplexRestraint"
  ComplexRestraint_sp ComplexRestraint_O::make()
  {
    IMPLEMENT_ME();
  };


void	ComplexRestraint_O::initialize()
{
    this->Base::initialize();
}

#ifdef XML_ARCHIVE
    void	ComplexRestraint_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
}
#endif


;

#if INIT_TO_FACTORIES

SYMBOL_EXPORT_SC_(ChemKwPkg,E);
SYMBOL_EXPORT_SC_(ChemKwPkg,Z);

#define ARGS_RestrainedPiBond_O_make "(config pi_atoms)"
#define DECL_RestrainedPiBond_O_make ""
#define DOCS_RestrainedPiBond_O_make "make RestrainedPiBond"
RestrainedPiBond_sp RestrainedPiBond_O::make(core::Symbol_sp config, core::List_sp piAtoms)
  {
    auto  me  = gctools::GC<RestrainedPiBond_O>::allocate_with_default_constructor();
    me->_Configuration = config;
    if ( me->_Configuration != chemkw::_sym_E && me->_Configuration != chemkw::_sym_Z )
    {
	SIMPLE_ERROR(("config:  must be ""E"" or ""Z""" ));
    }
    if ( core::cl__length(piAtoms) != 2 )
    {
	SIMPLE_ERROR("You must supply two atom names in piAtoms: argument");
    }
    me->_PiAtomX = oCar(piAtoms).as<MatterName::Type>();
    me->_PiAtomY = oCadr(piAtoms).as<MatterName::Type>();
    return me;
  };

#else

    core::T_sp 	RestrainedPiBond_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::LispPtr lisp)
{
    this->_Configuration = env->lookup(Pkg(),"config").as<core::Str_O>()->get();
    core::Cons_sp piAtoms = env->lookup(Pkg(),"piAtoms").as<core::Cons_O>();
    if ( this->_Configuration != "E" && this->_Configuration != "Z" )
    {
	SIMPLE_ERROR(("config:  must be ""E"" or ""Z""" ));
    }
    if ( core::cl__length(piAtoms) != 2 )
    {
	SIMPLE_ERROR("You must supply two atom names in piAtoms: argument");
    }
    this->_PiAtomX = piAtoms->car<core::Str_O>()->get();
    this->_PiAtomY = piAtoms->cadr<core::Str_O>()->get();
    return _Nil<core::T_O>();
}

#endif

void	RestrainedPiBond_O::fillOneDihedralRestraint(Residue_sp residue,
				Atom_sp p1X, Atom_sp pX, Atom_sp pY, Atom_sp p1Y,
				double deg, double weight)
{
    RestraintDihedral_sp dih = RestraintDihedral_O::create();
    dih->setAtomA(p1X);
    dih->setAtomB(pX);
    dih->setAtomC(pY);
    dih->setAtomD(p1Y);
    dih->setParameters(deg,weight);
    residue->addRestraint(dih);
}



/*! Create four dihedral restraints for each RestrainedPiBond entry
 * Only one trans bond was defined, the other three interactions are inferred
 *
 *	x1   y1
 *	  X=Y
 *	x2   y2
 *
 * Depending on if _Configuration is "E" or "Z"
 * four dihedral restraints are set up with:
 * for "E"  x1&y1 trans
 * for "Z"  x1&y1 cis
 */
size_t	RestrainedPiBond_O::fillRestraints(Residue_sp residue, core::HashTable_sp cip)
{
  size_t terms = 0;
  core::Symbol_sp config = this->_Configuration;
  Atom_sp pX  = gc::As_unsafe<Atom_sp>(residue->atomWithName(this->getPiAtomX()));
  Atom_sp pY  = gc::As_unsafe<Atom_sp>(residue->atomWithName(this->getPiAtomY()));

  Atom_sp p1X = pX->highestPriorityNeighborThatIsnt(pY,cip);
  Atom_sp p2X = pX->lowestPriorityNeighborThatIsnt(pY,cip);
  if ( p2X == p1X )
    {
      p2X = nil<Atom_O>();
    }
  Atom_sp p1Y = pY->highestPriorityNeighborThatIsnt(pX,cip);
  Atom_sp p2Y = pY->lowestPriorityNeighborThatIsnt(pX,cip);
  if ( p2Y == p1Y )
    {
      p2Y = nil<Atom_O>();
    }
  if ( pX->numberOfBonds() < 2 || pX->numberOfBonds() > 3  )
    {
      stringstream ss;
      ss << "Dihedral restraints X atom: " << this->_PiAtomX;
      ss << " must have between 2 and 3 bonds - it has " << pX->numberOfBonds();
      SIMPLE_ERROR("Dihedral restraints Y atom: {}" , core::_rep_(this->_PiAtomY) );
      ss << " must have between 2 and 3 bonds - it has " << pY->numberOfBonds();
      SIMPLE_ERROR("{}" , ss.str() );
    }
  double transDeg = -180.0;
  double cisDeg = 0.0;
  double weight = 1.0;
  if ( config == chemkw::_sym_Z )
    {
      this->fillOneDihedralRestraint(residue,p1X,pX,pY,p1Y,cisDeg,weight);
      if ( p2Y.notnilp() ) 
	{
          this->fillOneDihedralRestraint(residue,p1X,pX,pY,p2Y,transDeg,weight);
          terms++;
	}
      if ( p2X.notnilp() ) 
	{
          this->fillOneDihedralRestraint(residue,p2X,pX,pY,p1Y,transDeg,weight);
          terms++;
	}
      if ( p2X.notnilp() && p2Y.notnilp() )
	{
          this->fillOneDihedralRestraint(residue,p2X,pX,pY,p2Y,cisDeg,weight);
          terms++;
	}
    } else if ( config == chemkw::_sym_E )
    {
      this->fillOneDihedralRestraint(residue,p1X,pX,pY,p1Y,transDeg,weight);
      if ( p2Y.notnilp() )
	{
          this->fillOneDihedralRestraint(residue,p1X,pX,pY,p2Y,cisDeg,weight);
          terms++;
	}
      if ( p2X.notnilp() )
	{
          this->fillOneDihedralRestraint(residue,p2X,pX,pY,p1Y,cisDeg,weight);
          terms++;
	}
      if ( p2X.notnilp() && p2Y.notnilp() )
	{
          this->fillOneDihedralRestraint(residue,p2X,pX,pY,p2Y,transDeg,weight);
          terms++;
	}
    } else
    {
      SIMPLE_ERROR("Illegal Pi bond configuration: {} must be E or Z" , core::_rep_(config));
    }
  return terms;
}





bool RestrainedExoCyclicAtom_O::_LazyInitializedSmarts = false;
SYMBOL_EXPORT_SC_(ChemPkg,STARAtomExoToSixMemberedRingSTAR);

void	RestrainedExoCyclicAtom_O::lazyInitializeSmarts()
{
  SIMPLE_ERROR("Fix me - my smarts is wrong");
    if ( !RestrainedExoCyclicAtom_O::_LazyInitializedSmarts )
    {
	RestrainedExoCyclicAtom_O::_LazyInitializedSmarts = true;
        SmartsRoot_sp atomExoToSixMemberedRing = chem__compile_smarts("[*]1[*:1][*:2][*:3][*:4][*:5][*:6]9",nil<core::T_O>());
        chem::_sym_STARAtomExoToSixMemberedRingSTAR->defparameter(atomExoToSixMemberedRing);
    }
}

void	RestrainedExoCyclicAtom_O::initialize()
{
    this->_ExoCyclicAtomName = nil<MatterName::Type>();
    this->_Axial = true;
}

#if INIT_TO_FACTORIES

#define ARGS_RestrainedExoCyclicAtom_O_make "()"
#define DECL_RestrainedExoCyclicAtom_O_make ""
#define DOCS_RestrainedExoCyclicAtom_O_make "make RestrainedExoCyclicAtom"
  RestrainedExoCyclicAtom_sp RestrainedExoCyclicAtom_O::make()
  {
    IMPLEMENT_ME();
  };

#else

    core::T_sp 	RestrainedExoCyclicAtom_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::LispPtr lisp)
{
    IMPLEMENT_ME();
//    this->Base::oldLispInitialize(kargs,env);
}

#endif

#ifdef XML_ARCHIVE
void	RestrainedExoCyclicAtom_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("exoCyclicAtomName",this->_ExoCyclicAtomName);
    node->attribute("axial",this->_Axial);
}
#endif



size_t RestrainedExoCyclicAtom_O::fillRestraints(Residue_sp residue, core::HashTable_sp cip )
{
  SIMPLE_ERROR("Fix the implementation of {}" , __FUNCTION__);
  this->lazyInitializeSmarts();
  if ( !residue->hasAtomWithName(this->_ExoCyclicAtomName) )
    {
      SIMPLE_ERROR("Residue({}) doesn't have atom with name({})" , residue->description() , core::_rep_(this->_ExoCyclicAtomName) );
    }
  Atom_sp exoCyclicAtom = gc::As_unsafe<Atom_sp>(residue->atomWithName(this->_ExoCyclicAtomName));
  Root_sp atomExoToSixMemberedRing = gctools::As<Root_sp>(chem::_sym_STARAtomExoToSixMemberedRingSTAR->symbolValue());
  core::T_mv match_mv = chem__chem_info_match(atomExoToSixMemberedRing,exoCyclicAtom);
  if ( match_mv.nilp() )
    {
      SIMPLE_ERROR("In residue({}) the atom with name({}) is not exo-cyclic to a six-membered ring" , residue->description() , _rep_(this->_ExoCyclicAtomName) );
    }
  //ChemInfoMatch_sp match = gc::As<ChemInfoMatch_sp>(values.second(match_mv.number_of_values()));
  FIX_ME(); // is the above correct?
}


};
