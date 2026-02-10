/*
    File: energyNonbond14.h
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

//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	energyNonbond14.h
 *
 *	Maintains a database of stretch types
 */

#pragma once

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>
#include "clasp/core/ql.h"
//#include "geom/render.fwd.h"// energyNonbond14.h wants DisplayList needs render.fwd.h

namespace       chem {
  FORWARD(EnergyFunction); // Declares class EnergyFunction_O {} and EnergyFunction_sp
  FORWARD(EnergyNonbond14);
  FORWARD(FFNonbondDb);
  FORWARD(AtomTable);

// Bring in the AmberNonbond14Term kernel term
#include "cando/chem/energyKernels/amber_nonbond14.h"


/*! A Nonbond term, this is subclassed in AmberFunction to deal with Amber specific stuff.
 * I should combine energyFunction.cc and amberFunction.cc together.
 */
class EnergyNonbond14 : public EnergyTerm {
public:
  Atom_sp      	_Atom1_enb;
  Atom_sp      	_Atom2_enb;
  AmberNonbond14Term	term;
public:
  string	className()	{ return "EnergyNonbond"; };
  Atom_sp	getAtom1() { return this->_Atom1_enb; };
  Atom_sp	getAtom2() { return this->_Atom2_enb; };
  bool	defineForAtomPair(core::T_sp	forceField,
                          Atom_sp        a1,
                          Atom_sp        a2,
                          size_t         a1CoordinateIndexTimes3,
                          size_t         a2CoordinateIndexTimes3,
                          EnergyNonbond14_sp energyNonbond14,
                          core::HashTable_sp atomTypes,
                          core::T_sp     keepInteraction
                          );

public:
  core::List_sp encode() const;
  void decode(core::List_sp alist);
};
};

namespace translate {

  template <>
struct	to_object<chem::EnergyNonbond14>
  {
    typedef	core::Cons_sp ExpectedType;
    typedef	core::Cons_sp DeclareType;
    static core::T_sp convert(const chem::EnergyNonbond14& nonbond)
    {
      return nonbond.encode();
    }
  };

  template <>
    struct	from_object<chem::EnergyNonbond14>
  {
    typedef	chem::EnergyNonbond14	ExpectedType;
    typedef	ExpectedType 		DeclareType;
    DeclareType _v;
    from_object(core::T_sp o)
    {
      SIMPLE_ERROR("Implement me");
    }
  };
};

namespace chem {


class EnergyNonbond14_O : public EnergyComponent_O
{
  LISP_CLASS(chem,ChemPkg,EnergyNonbond14_O,"EnergyNonbond14",EnergyComponent_O);

 public:
  virtual bool restraintp() const override {return false;};
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);

 public: // virtual functions inherited from Object
  void	initialize();
//	string	__repr__() const;
 public:
  typedef EnergyNonbond14 TermType;
 public: // instance variables
  gctools::Vec0<TermType>         _Terms14;
  // Pairlist
  AtomTable_sp                    _AtomTable;
  core::T_sp                      _NonbondForceField;
  core::HashTable_sp              _AtomTypes;
  core::T_sp                      _KeepInteractionFactory;
  //  Tables for nonbonded calculation using excluded atoms and the Amber way
  size_t                          _ntypes;          // ntypes
  core::SimpleVector_sp           _atom_name_vector;  // atom-name-vector
  NVector_sp                      _charge_vector;          // charge-vector
  NVector_sp                      _mass_vector;            // masses
  core::SimpleVector_int32_t_sp   _atomic_number_vector;    // vec

  // Excluded atom table
  core::SimpleVector_int32_t_sp   _NumberOfExcludedAtomIndexes;
  core::SimpleVector_int32_t_sp   _ExcludedAtomIndexes;
  size_t                          _InteractionsKept;
  size_t                          _InteractionsDiscarded;

 public:
  virtual std::string implementation_details() const;
  virtual std::string descriptionOfContents() const;
  typedef gctools::Vec0<TermType>::iterator iterator;
   static EnergyNonbond14_sp make(SetupAccumulator& setupAcc);
 public:
  virtual size_t numberOfTerms() { return this->_Terms14.size(); };
  void callForEachTerm(core::Function_sp callback);
 public:

  size_t runTestCalls(core::T_sp stream, chem::NVector_sp coords) const;

  CL_DEFMETHOD core::T_sp keepInteractionFactory() const { return this->_KeepInteractionFactory; };
 public:
  void addTerm14(const TermType& term);
  virtual void dumpTerms(core::HashTable_sp atomTypes);

  virtual double evaluateAllComponent( ScoringFunction_sp scorer,
                                       NVector_sp 	pos,
                                       core::T_sp energyScale,
                                       core::T_sp componentEnergy,
                                       bool 		calcForce,
                                       gc::Nilable<NVector_sp> 	force,
                                       bool		calcDiagonalHessian,
                                       bool		calcOffDiagonalHessian,
                                       gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                       gc::Nilable<NVector_sp>	hdvec,
                                       gc::Nilable<NVector_sp> dvec,
                                       core::T_sp activeAtomMask,
                                       core::T_sp debugInteractions );

  double debugAllComponent( ScoringFunction_sp scorer,
                            NVector_sp 	pos,
                            core::T_sp energyScale,
                            core::T_sp componentEnergy,
                            bool 		calcForce,
                            gc::Nilable<NVector_sp> 	force,
                            bool		calcDiagonalHessian,
                            bool		calcOffDiagonalHessian,
                            gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                            gc::Nilable<NVector_sp>	hdvec,
                            gc::Nilable<NVector_sp> dvec,
                            core::T_sp activeAtomMask,
                            core::T_sp debugInteractions );

  core::T_sp getFFNonbondDb();

  void construct14InteractionTerms(AtomTable_sp atomTable, Matter_sp matter, core::T_sp nbforceField, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes );

  EnergyComponent_sp copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc);

 public:
  EnergyNonbond14_O( const EnergyNonbond14_O& ss ); //!< Copy constructor

  EnergyNonbond14_O() :
      _KeepInteractionFactory(nil<core::T_O>())
  {};
};

};

