/*
    File: energyComponentGroup.cc
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
 *	energyComponentGroup.cc
 *
 *	A composite EnergyComponent_O - see the header for why it exists.
 */
#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/energyComponentGroup.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/nVector.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>

namespace chem {

// The copy constructor is DECLARED in the header and deliberately not defined here.  Defining
// it would emit a call to EnergyComponent_O's copy constructor, which is declared but has no
// definition anywhere - the link fails with an undefined reference to it.  Every other
// component does the same: EnergyDihedralRestraint_O declares one and never defines it.
// copyFilter builds copies through make() + addComponent instead.

CL_LISPIFY_NAME(make-energy-component-group);
CL_DEF_CLASS_METHOD
EnergyComponentGroup_sp EnergyComponentGroup_O::make()
{
  return gctools::GC<EnergyComponentGroup_O>::allocate();
}

void EnergyComponentGroup_O::fields(core::Record_sp node)
{
  node->field_if_not_empty(INTERN_(kw,components),this->_Components);
  node->field_if_not_nil(INTERN_(kw,name),this->_GroupName);
  this->Base::fields(node);
}

void EnergyComponentGroup_O::addComponent(EnergyComponent_sp c)
{
  // Push-then-assign, the same dance Joint_O::_appendChild uses: Vec0 wants an element
  // constructed in place before a GC-managed pointer is written into it.
  EnergyComponent_sp empty;
  size_t index = this->_Components.size();
  this->_Components.push_back(empty);
  this->_Components[index] = c;
}

core::List_sp EnergyComponentGroup_O::components() const
{
  core::List_sp result = nil<core::T_O>();
  for ( int i = (int)this->_Components.size()-1; i >= 0; --i ) {
    result = core::Cons_O::create(this->_Components[i],result);
  }
  return result;
}

bool EnergyComponentGroup_O::restraintp() const
{
  for ( auto& c : this->_Components ) if ( c->restraintp() ) return true;
  return false;
}

bool EnergyComponentGroup_O::is_restraint() const
{
  for ( auto& c : this->_Components ) if ( c->is_restraint() ) return true;
  return false;
}

size_t EnergyComponentGroup_O::numberOfTerms()
{
  size_t total = 0;
  for ( auto& c : this->_Components ) total += c->numberOfTerms();
  return total;
}

std::string EnergyComponentGroup_O::descriptionOfContents() const
{
  // The NAME first when there is one - a blueprint holds ~940 of these and they are otherwise
  // indistinguishable in a backtrace.
  if (this->_GroupName.notnilp())
    return fmt::format("{} {} components", _rep_(this->_GroupName), this->_Components.size());
  return fmt::format("{} components", this->_Components.size());
}

std::string EnergyComponentGroup_O::implementation_details() const
{
  return "energyComponentGroup.cc - composite over other energy components";
}

string EnergyComponentGroup_O::__repr__() const
{
  // Not className() - it is virtual but NOT const on EnergyComponent_O (energyComponent.h:146),
  // so it cannot be called from a const method.
  stringstream ss;
  ss << "#<EnergyComponentGroup " << this->descriptionOfContents() << ">";
  return ss.str();
}

void EnergyComponentGroup_O::dumpTerms(core::HashTable_sp atomTypes)
{
  for ( auto& c : this->_Components ) c->dumpTerms(atomTypes);
}

void EnergyComponentGroup_O::atomsForEachTerm(core::Function_sp callback)
{
  for ( auto& c : this->_Components ) c->atomsForEachTerm(callback);
}

core::List_sp EnergyComponentGroup_O::extract_vectors_as_alist() const
{
  core::List_sp result = nil<core::T_O>();
  for ( int i = (int)this->_Components.size()-1; i >= 0; --i ) {
    result = core::Cons_O::create(this->_Components[i]->extract_vectors_as_alist(),result);
  }
  return result;
}

size_t EnergyComponentGroup_O::runTestCalls(core::T_sp stream, chem::NVector_sp pos) const
{
  size_t total = 0;
  for ( auto& c : this->_Components ) total += c->runTestCalls(stream,pos);
  return total;
}

void EnergyComponentGroup_O::emitTestCalls(core::T_sp stream, chem::NVector_sp pos) const
{
  for ( auto& c : this->_Components ) c->emitTestCalls(stream,pos);
}

void EnergyComponentGroup_O::setupHessianPreconditioner(NVector_sp nvPosition,
                                                        AbstractLargeSquareMatrix_sp m,
                                                        core::T_sp activeAtomMask )
{
  for ( auto& c : this->_Components ) c->setupHessianPreconditioner(nvPosition,m,activeAtomMask);
}

double EnergyComponentGroup_O::evaluateAllComponent( ScoringFunction_sp scorer,
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
  // A disabled CHILD contributes nothing, which is how one slot is evaluated without the
  // others.  The group's own _Enabled gates the whole set.
  for ( auto& c : this->_Components ) {
    if ( !c->isEnabled() ) continue;
    totalEnergy += c->evaluateAllComponent(scorer,pos,energyScale,componentEnergy,
                                           calcForce,force,
                                           calcDiagonalHessian,calcOffDiagonalHessian,
                                           hessian,hdvec,dvec,
                                           activeAtomMask,debugInteractions);
  }
  return totalEnergy;
}

void EnergyComponentGroup_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos )
{
  for ( auto& c : this->_Components ) c->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
}

EnergyComponent_sp EnergyComponentGroup_O::copyFilter(core::T_sp keepInteractionFactory,
                                                      SetupAccumulator& setupAcc)
{
  // Filter each child and keep the group's shape.  A child whose copyFilter yields NIL - every
  // interaction rejected - is dropped rather than kept empty.
  EnergyComponentGroup_sp copy = EnergyComponentGroup_O::make();
  copy->_Enabled = this->_Enabled;
  copy->_Scale = this->_Scale;
  copy->_GroupName = this->_GroupName;   // a filtered copy is still the same group
  for ( auto& c : this->_Components ) {
    EnergyComponent_sp filtered = c->copyFilter(keepInteractionFactory,setupAcc);
    if ( filtered.notnilp() ) copy->addComponent(filtered);
  }
  return copy;
}

}; // namespace chem
