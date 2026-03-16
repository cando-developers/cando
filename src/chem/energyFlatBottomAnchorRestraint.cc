/*
    File: energyFlatBottomAnchorRestraint.cc
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

This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University
at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#define DEBUG_LEVEL_NONE

#include <clasp/core/foundation.h>
#include <clasp/core/ql.h>
#include <cando/chem/energyFlatBottomAnchorRestraint.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/nVector.h>
#include <clasp/core/wrappers.h>

namespace chem {

// Pull in the generated kernel (defines Flat_Bottom_Anchor_Restraint<HESSIAN>)
#include <cando/chem/energyKernels/flat_bottom_anchor_restraint.c>

#define FBAR_APPLY_ATOM_MASK(I1) \
  if (hasActiveAtomMask && !bitvectorActiveAtomMask->testBit(I1/3)) continue;

core::List_sp EnergyFlatBottomAnchorRestraint::encode() const {
  return core::Cons_O::createList(
    core::Cons_O::create(INTERN_(kw, ka), core::clasp_make_double_float(this->term.ka)),
    core::Cons_O::create(INTERN_(kw, R0), core::clasp_make_double_float(this->term.r0)),
    core::Cons_O::create(INTERN_(kw, i1), core::make_fixnum(this->term.i3x1)),
    core::Cons_O::create(INTERN_(kw, xa), core::clasp_make_double_float(this->term.xa)),
    core::Cons_O::create(INTERN_(kw, ya), core::clasp_make_double_float(this->term.ya)),
    core::Cons_O::create(INTERN_(kw, za), core::clasp_make_double_float(this->term.za)),
    core::Cons_O::create(INTERN_(kw, atom1), this->_Atom1));
}

void EnergyFlatBottomAnchorRestraint_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, terms), this->_Terms);
  this->Base::fields(node);
}

void EnergyFlatBottomAnchorRestraint_O::initialize() {
  this->Base::initialize();
}

std::string EnergyFlatBottomAnchorRestraint_O::descriptionOfContents() const {
  stringstream ss;
  ss << ":enabled " << ((this->_Enabled) ? "T" : "NIL");
  ss << " number-of-terms " << this->_Terms.size();
  return ss.str();
}

void EnergyFlatBottomAnchorRestraint_O::addTerm(const EnergyFlatBottomAnchorRestraint& r) {
  this->_Terms.push_back(r);
}

size_t EnergyFlatBottomAnchorRestraint_O::addFlatBottomAnchorRestraintTerm(
    EnergyFunction_sp energyFunction, Atom_sp a1,
    Vector3 center, double R0, double forceConstant) {
  AtomTable_sp atomTable = energyFunction->atomTable();
  EnergyAtom* ea = atomTable->getEnergyAtomPointer(a1);
  this->_Terms.emplace_back(a1, forceConstant, R0, center,
                            ea->coordinateIndexTimes3());
  return this->_Terms.size() - 1;
}

void EnergyFlatBottomAnchorRestraint_O::dumpTerms(core::HashTable_sp atomTypes) {}

void EnergyFlatBottomAnchorRestraint_O::setupHessianPreconditioner(
    NVector_sp nvPosition, AbstractLargeSquareMatrix_sp m,
    core::T_sp activeAtomMask) {
  return;
}

double EnergyFlatBottomAnchorRestraint_O::evaluateAllComponent(
    ScoringFunction_sp scorer, NVector_sp pos,
    core::T_sp energyScale, core::T_sp energyComponents,
    bool calcForce, gc::Nilable<NVector_sp> force,
    bool calcDiagonalHessian, bool calcOffDiagonalHessian,
    gc::Nilable<AbstractLargeSquareMatrix_sp> hessian,
    gc::Nilable<NVector_sp> hdvec, gc::Nilable<NVector_sp> dvec,
    core::T_sp activeAtomMask, core::T_sp debugInteractions) {

  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  this->_Evaluations++;

  auto evalType = determineEnergyComponentEvalType(force, hdvec, dvec);
  double totalEnergy = 0.0;
  DOUBLE* position = &(*pos)[0];
  DOUBLE* rforce = NULL;
  DOUBLE* rdvec = NULL;
  DOUBLE* rhdvec = NULL;
  DOUBLE Energy = 0.0;
  Flat_Bottom_Anchor_Restraint<NoHessian> kernel;

  if (evalType == energyEval) {
    for (auto cri = this->_Terms.begin(); cri != this->_Terms.end(); cri++) {
      FBAR_APPLY_ATOM_MASK(cri->term.i3x1);
      Energy = kernel.energy(cri->term, position, &totalEnergy);
    }
  } else if (evalType == gradientEval) {
    rforce = &(*force)[0];
    for (auto cri = this->_Terms.begin(); cri != this->_Terms.end(); cri++) {
      FBAR_APPLY_ATOM_MASK(cri->term.i3x1);
      Energy = kernel.gradient(cri->term, position, &totalEnergy, rforce);
    }
  } else {
    rforce = &(*force)[0];
    rdvec = &(*dvec)[0];
    rhdvec = &(*hdvec)[0];
    for (auto cri = this->_Terms.begin(); cri != this->_Terms.end(); cri++) {
      FBAR_APPLY_ATOM_MASK(cri->term.i3x1);
      Energy = kernel.hessian(cri->term, position, &totalEnergy,
                              rforce, NoHessian(), rdvec, rhdvec);
    }
  }
  maybeSetEnergy(energyComponents,
                 EnergyFlatBottomAnchorRestraint_O::static_classSymbol(),
                 totalEnergy);
  return totalEnergy;
}

EnergyComponent_sp EnergyFlatBottomAnchorRestraint_O::copyFilter(
    core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc) {
  EnergyFlatBottomAnchorRestraint_sp copy =
      EnergyFlatBottomAnchorRestraint_O::create();
  copyEnergyComponent(copy, this->asSmartPtr());
  core::T_sp keepInteraction = specializeKeepInteractionFactory(
      keepInteractionFactory, EnergyFlatBottomAnchorRestraint_O::staticClass());
  for (auto edi = this->_Terms.begin(); edi != this->_Terms.end(); edi++) {
    Atom_sp a1 = edi->_Atom1;
    if (skipInteraction(keepInteraction, a1)) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}

} // namespace chem
