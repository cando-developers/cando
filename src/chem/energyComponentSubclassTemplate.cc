/*
    File: energyComponentSubclassTemplate.cc
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
 *	energyComponentSubclassTemplate.cc
 *
 *	TEMPLATE for a new EnergyComponent_O subclass - see the header for how to use it.
 */
#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/energyComponentSubclassTemplate.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/nVector.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>

namespace chem {

// DO NOT define the copy constructor.  It is declared in the header, and defining it emits a
// call to EnergyComponent_O's copy constructor, which is declared but never defined anywhere -
// the link fails with an undefined reference.  Every existing component declares one and leaves
// it undefined; build copies through make() and the add* methods instead.

void EnergyComponentSubclassTemplate_O::fields(core::Record_sp node)
{
  // Terms are POD in a Vec0 - serialize whatever your component needs to reconstruct itself.
  // Delete this method entirely if the component is rebuilt rather than loaded.
  this->Base::fields(node);
}

void EnergyComponentSubclassTemplate_O::addTerm(const TermType& term)
{
  this->_Terms.push_back(term);
}

std::string EnergyComponentSubclassTemplate_O::descriptionOfContents() const
{
  return fmt::format("{} terms", this->_Terms.size());
}

std::string EnergyComponentSubclassTemplate_O::implementation_details() const
{
  return "energyComponentSubclassTemplate.cc - replace me";
}

string EnergyComponentSubclassTemplate_O::__repr__() const
{
  // Not className() - it is virtual but NOT const on EnergyComponent_O (energyComponent.h:146),
  // so it cannot be called from a const method.
  stringstream ss;
  ss << "#<EnergyComponentSubclassTemplate " << this->descriptionOfContents() << ">";
  return ss.str();
}

void EnergyComponentSubclassTemplate_O::dumpTerms(core::HashTable_sp atomTypes)
{
  // Print one line per term.  Used by the energy-function dump routines.
  for ( auto& term : this->_Terms ) {
    core::clasp_write_string(fmt::format("term k({}) I1({}) I2({})\n", term.k, term.I1, term.I2));
  }
}

void EnergyComponentSubclassTemplate_O::atomsForEachTerm(core::Function_sp callback)
{
  // Hand the callback the atoms of each term.  The energy-function walkers use this to relate
  // terms back to atoms - a component that does not implement it is invisible to them.
  IMPLEMENT_ME();
}

core::List_sp EnergyComponentSubclassTemplate_O::extract_vectors_as_alist() const
{
  IMPLEMENT_ME();
}

size_t EnergyComponentSubclassTemplate_O::runTestCalls(core::T_sp stream, chem::NVector_sp pos) const
{
  IMPLEMENT_ME();
}

void EnergyComponentSubclassTemplate_O::emitTestCalls(core::T_sp stream, chem::NVector_sp pos) const
{
  IMPLEMENT_ME();
}

void EnergyComponentSubclassTemplate_O::setupHessianPreconditioner(NVector_sp nvPosition,
                                                                   AbstractLargeSquareMatrix_sp m,
                                                                   core::T_sp activeAtomMask )
{
  // Only needed if this component contributes to the preconditioner.  Leaving it empty is
  // legitimate; the base does not require it.
}

double EnergyComponentSubclassTemplate_O::evaluateAllComponent( ScoringFunction_sp scorer,
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
                                                                core::T_sp debugInteractions )
{
  double totalEnergy = 0.0;
  if ( !this->isEnabled() ) return totalEnergy;
  this->_Evaluations++;
  // Walk _Terms and accumulate.  A COMPOSITE component ignores _Terms entirely and instead
  // forwards this whole call to each enabled child, summing the results - that is the only
  // difference between a leaf component and a group.
  for ( auto& term : this->_Terms ) {
    (void)term;
    // energy += ...;  and when calcForce, write into (*force)[term.I1 + 0..2] etc.
  }
  return totalEnergy;
}

void EnergyComponentSubclassTemplate_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos )
{
  // Finite-difference check of the analytical force/hessian.  Worth implementing while the
  // derivatives are fresh - every existing component has one.
  IMPLEMENT_ME();
}

EnergyComponent_sp EnergyComponentSubclassTemplate_O::copyFilter(core::T_sp keepInteractionFactory,
                                                                 SetupAccumulator& setupAcc)
{
  // Return a copy holding only the terms KEEPINTERACTIONFACTORY selects.  Use
  // specializeKeepInteractionFactory(keepInteractionFactory, this->__class()) to get the
  // per-term predicate, then skipInteraction(...) per term - see energyFunction.cc:996-1030.
  IMPLEMENT_ME();
}

}; // namespace chem
