/*
    File: energyComponentSubclassTemplate.h
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
 *	energyComponentSubclassTemplate.h
 *
 *	TEMPLATE for a new EnergyComponent_O subclass.
 *
 *	HOW TO USE
 *	  1. Copy energyComponentSubclassTemplate.{fwd.h,h,cc}
 *	  2. Rename EnergyComponentSubclassTemplate -> YourClass, and the include guards
 *	  3. Add src/chem/yourClass.cc to src/chem/cscript.lisp
 *	  4. Replace TermType with your term struct, or delete the term machinery entirely if
 *	     your component holds something else (a composite over other components, say)
 *
 *	This file compiles as-is, so it stays honest.  Every method here is one the base class
 *	either requires (pure virtual) or defaults in a way that will bite you if left alone -
 *	NUMBEROFTERMS and DESCRIPTIONOFCONTENTS both SUBCLASS_MUST_IMPLEMENT on the base.
 */

#ifndef EnergyComponentSubclassTemplate_H  //[
#define	EnergyComponentSubclassTemplate_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/energyComponent.h>

namespace       chem {

  FORWARD(EnergyFunction);

/*! One term.  Replace with whatever this component stores.  Terms are POD held in a
    gctools::Vec0, so no smart pointers in here - keep atom references as coordinate indices
    (the I1,I2,... convention every other component uses). */
struct TermEnergyComponentSubclassTemplate
{
  REAL		k;	// force constant, or whatever parameterizes one term
  INT		I1;	// coordinate index x3 of the first atom
  INT		I2;	// coordinate index x3 of the second atom
};

SMART(EnergyComponentSubclassTemplate);
class EnergyComponentSubclassTemplate_O : public EnergyComponent_O
{
  LISP_CLASS(chem,ChemPkg,EnergyComponentSubclassTemplate_O,"EnergyComponentSubclassTemplate",EnergyComponent_O);

public:
  typedef TermEnergyComponentSubclassTemplate	TermType;

public: // instance variables
  gctools::Vec0<TermType>	_Terms;

public:
  typedef gctools::Vec0<TermType>::iterator iterator;
  iterator begin() { return this->_Terms.begin(); };
  iterator end() { return this->_Terms.end(); };

public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);

public:
  /*! TRUE if this component is a restraint rather than a physical energy.  Restraints are
      disabled when computing conformational energies - see EnergyComponent_O. */
  virtual bool restraintp() const { return false; };
  virtual bool is_restraint() const { return false; };

  virtual std::string descriptionOfContents() const;
  virtual std::string implementation_details() const;
  virtual size_t numberOfTerms() { return this->_Terms.size(); };

  void addTerm(const TermType& term);

  virtual void dumpTerms(core::HashTable_sp atomTypes);
  virtual void atomsForEachTerm(core::Function_sp callback);
  virtual core::List_sp extract_vectors_as_alist() const;

  virtual size_t runTestCalls(core::T_sp stream, chem::NVector_sp pos) const;
  virtual void emitTestCalls(core::T_sp stream, chem::NVector_sp pos) const;

  virtual void setupHessianPreconditioner(NVector_sp nvPosition,
                                          AbstractLargeSquareMatrix_sp m,
                                          core::T_sp activeAtomMask );

  /*! The one that matters.  Accumulate this component's energy into the return value, and when
      CALCFORCE/CALC*HESSIAN are set, into FORCE/HESSIAN.  Respect ACTIVEATOMMASK.  A composite
      component forwards this to its children instead of walking _Terms. */
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

  virtual void compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos );

  /*! Return a copy keeping only the interactions KEEPINTERACTIONFACTORY selects.  See
      specializeKeepInteractionFactory / skipInteraction in energyFunction.cc - the factory is
      called with this component's CLASS and returns T, NIL, or a per-term predicate. */
  EnergyComponent_sp copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc);

  virtual string __repr__() const;

public:
  EnergyComponentSubclassTemplate_O( const EnergyComponentSubclassTemplate_O& ss ); //!< Copy constructor

  EnergyComponentSubclassTemplate_O() {};
  virtual ~EnergyComponentSubclassTemplate_O() {};
};

};

#endif //]