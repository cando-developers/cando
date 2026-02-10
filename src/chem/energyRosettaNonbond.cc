/*
    File: energyRosettaNonbond.cc
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
#include <clasp/core/bformat.h>
#include <cando/chem/energyRosettaNonbond.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/array.h>
#include <clasp/core/evaluator.h>
#include <cando/geom/color.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/nativeVector.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/matter.h>
#include <cando/chem/bond.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/loop.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>

namespace chem {

#include <cando/chem/energyKernels/rosetta_nonbond_dd_cutoff.c>

EnergyRosettaNonbond_sp EnergyRosettaNonbond_O::make(AtomTable_sp atomTable, core::T_sp nbForceField,
                                  core::HashTable_sp atomTypes, core::T_sp keepInteractionFactory,
                                  SetupAccumulator& acc ,
                                  core::T_sp tcoordinates) {
  auto obj = EnergyRosettaNonbond_O::create();
  obj->_Parameters.do_apply(acc);
  if (keepInteractionFactory.nilp()) return obj;
  obj->_AtomTable = atomTable;
  obj->_NonbondForceField = nbForceField;
  obj->_AtomTypes = atomTypes;
  obj->_KeepInteractionFactory = keepInteractionFactory;
  if (tcoordinates.notnilp()) {
    obj->rebuildPairList(tcoordinates);
  }
  return obj;
}

std::string EnergyRosettaNonbond_O::implementation_details() const {
  Rosetta_Nonbond_Dd_Cutoff<NoHessian> nonbond;
  std::stringstream ss;
  ss << nonbond.description();
  return ss.str();
}

std::string EnergyRosettaNonbond_O::descriptionOfContents() const {
  stringstream ss;
  ss << ":enabled " << ((this->_Enabled) ? "T" : "NIL");
  ss << " number-of-terms " << this->_Terms.size();
  return ss.str();
}

}

// #define DEBUG_NONBOND_TERM 1
#define LOG_ENERGY(...)
// #define LOG_ENERGY(...) core::clasp_write_string(fmt::format(__VA_ARGS__))

namespace chem {

SYMBOL_EXPORT_SC_(ChemPkg, EnergyRosettaNonbond);
SYMBOL_EXPORT_SC_(ChemPkg, energyRosettaNonbond);

core::List_sp EnergyRosettaNonbond::encode() const {
  ql::list ll;
  this->term.encode(ll);
  ll << INTERN_(kw, atom1) << this->_Atom1_enb
     << INTERN_(kw, atom2) << this->_Atom2_enb;
  return ll.cons();
}

void EnergyRosettaNonbond::decode(core::List_sp alist) { SIMPLE_ERROR("Implement decode of EnergyRosettaNonbond"); }

core::T_sp debug_rosetta_nonbond(double Energy, double x1, double y1, double z1, double x2, double y2, double z2,
                                 const rosetta_nonbond_term& term) {
  double r2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
  ql::list ll;
  ll << INTERN_(kw, Energy) << mk_double_float(Energy) << INTERN_(kw, x1) << mk_double_float(x1) << INTERN_(kw, y1)
     << mk_double_float(y1) << INTERN_(kw, z1) << mk_double_float(z1) << INTERN_(kw, x2) << mk_double_float(x2)
     << INTERN_(kw, y2) << mk_double_float(y2) << INTERN_(kw, z2) << mk_double_float(z2)
     << INTERN_(kw, sigma) << mk_double_float(term.sigma) << INTERN_(kw, epsilon) << mk_double_float(term.epsilon)
     << INTERN_(kw, r2) << mk_double_float(r2);
  return ll.cons();
}

#define NONBOND_DEBUG_INTERACTIONS(term)                                                                                         \
  if (doDebugInteractions) {                                                                                                      \
    core::eval::funcall(debugInteractions, chem::_sym_EnergyRosettaNonbond,                                                       \
                        debug_rosetta_nonbond(Energy,                                                                            \
                                              position[term.i3x1], position[term.i3x1 + 1], position[term.i3x1 + 2],             \
                                              position[term.i3x2], position[term.i3x2 + 1], position[term.i3x2 + 2],             \
                                              term),                                                                             \
                        core::make_fixnum(term.i3x1), core::make_fixnum(term.i3x2));                                              \
  }

struct NoFiniteDifference {
  static void maybeTestFiniteDifference(ScoringFunction_sp score, core::T_sp energyScale, int I1, int I2, core::T_sp activeAtomMask, double x1,
                                        double y1, double z1, double x2, double y2, double z2, double dA, double dC, double dQ1Q2,
                                        double fx1, double fy1, double fz1, double fx2, double fy2, double fz2, int index, size_t &fails,
                                        bool debugForce) {}
};

inline double calculate_dQ1Q2(double electrostaticScale, double electrostaticModifier, double charge1, double charge2) {
  return electrostaticScale * electrostaticModifier * charge1 * charge2;
}

template <class MaybeFiniteDiff>
double template_evaluateUsingTerms(EnergyRosettaNonbond_O* mthis,
                                   const gctools::Vec0<EnergyRosettaNonbond>& terms,
                                   core::T_sp termSymbol,
                                   ScoringFunction_sp score, NVector_sp nvposition,
                                   core::T_sp energyScale, core::T_sp energyComponents,
                                   bool calcForce, gc::Nilable<NVector_sp> force, bool calcDiagonalHessian,
                                   bool calcOffDiagonalHessian, gc::Nilable<AbstractLargeSquareMatrix_sp> hessian,
                                   gc::Nilable<NVector_sp> hdvec, gc::Nilable<NVector_sp> dvec, core::T_sp activeAtomMask,
                                   core::T_sp debugInteractions, size_t& fails, size_t& index, bool debugForce = false) {
  const rosetta_nonbond_parameters& params = mthis->_Parameters;
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  auto evalType = determineEnergyComponentEvalType(force, hdvec, dvec);
  double totalNonbondEnergy = 0.0;
  DOUBLE* position = &(*nvposition)[0];
  DOUBLE* rforce = NULL;
  DOUBLE* rdvec = NULL;
  DOUBLE* rhdvec = NULL;
  DOUBLE Energy = 0.0;
  Rosetta_Nonbond_Dd_Cutoff<NoHessian> nonbond;

#define KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(I1, I2)                                                            \
  if (hasActiveAtomMask && !(bitvectorActiveAtomMask->testBit(I1 / 3) && bitvectorActiveAtomMask->testBit(I2 / 3))) continue;

  if (evalType == energyEval) {
    for (auto si = terms.begin(); si != terms.end(); si++) {
      KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(si->term.i3x1, si->term.i3x2);
      Energy = nonbond.energy(params, si->term, position, &totalNonbondEnergy);
      NONBOND_DEBUG_INTERACTIONS(si->term);
    }
  } else if (evalType == gradientEval) {
    rforce = &(*force)[0];
    for (auto si = terms.begin(); si != terms.end(); si++) {
      KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(si->term.i3x1, si->term.i3x2);
      Energy = nonbond.gradient(params, si->term, position, &totalNonbondEnergy, rforce);
      NONBOND_DEBUG_INTERACTIONS(si->term);
    }
  } else {
    rforce = &(*force)[0];
    rdvec = &(*dvec)[0];
    rhdvec = &(*hdvec)[0];
    for (auto si = terms.begin(); si != terms.end(); si++) {
      KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(si->term.i3x1, si->term.i3x2);
      Energy = nonbond.hessian(params, si->term, position, &totalNonbondEnergy, rforce, NoHessian(), rdvec, rhdvec);
      NONBOND_DEBUG_INTERACTIONS(si->term);
    }
  }
  maybeSetEnergy(energyComponents, termSymbol, totalNonbondEnergy);
  return totalNonbondEnergy;
}

bool EnergyRosettaNonbond::defineForAtomPair(core::T_sp forceField, Atom_sp a1, Atom_sp a2,
                                             size_t i3x1, size_t i3x2,
                                             EnergyRosettaNonbond_sp energyRosettaNonbond,
                                             core::HashTable_sp atomTypes, core::T_sp keepInteraction,
                                             const rosetta_nonbond_parameters& params) {
  this->_Atom1_enb = a1;
  this->_Atom2_enb = a2;
  core::Symbol_sp t1 = a1->getType(atomTypes);
  core::Symbol_sp t2 = a2->getType(atomTypes);
  ASSERT(forceField && forceField.notnilp());
  core::T_sp tffNonbond1 = core::eval::funcall(_sym_find_type, forceField, t1);
  core::T_sp tffNonbond2 = core::eval::funcall(_sym_find_type, forceField, t2);
  if (tffNonbond1.nilp() || tffNonbond2.nilp()) return false;
  FFNonbond_sp ffNonbond1 = gc::As<FFNonbond_sp>(tffNonbond1);
  FFNonbond_sp ffNonbond2 = gc::As<FFNonbond_sp>(tffNonbond2);

  double rstar = ffNonbond1->getRadius_Angstroms() + ffNonbond2->getRadius_Angstroms();
  double epsilonij = sqrt(ffNonbond1->getEpsilon_kcal() * ffNonbond2->getEpsilon_kcal());
  double rtothe6th = rstar * rstar * rstar * rstar * rstar * rstar;
  double rtothe12th = rtothe6th * rtothe6th;
  double parmA = epsilonij * rtothe12th;
  double parmC = 2.0 * epsilonij * rtothe6th;

  this->term = rosetta_nonbond_term(params, parmA, parmC, i3x1, i3x2);
  return true;
}

void EnergyRosettaNonbond_O::initialize() { this->Base::initialize(); }

void EnergyRosettaNonbond_O::addTerm(const EnergyRosettaNonbond& term) { this->_Terms.push_back(term); }

void EnergyRosettaNonbond_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, terms), this->_Terms);
  node->field(INTERN_(kw, AtomTable), this->_AtomTable);
  node->field(INTERN_(kw, KeepInteractionFactory), this->_KeepInteractionFactory);
  this->_Parameters.fields(node);
  this->Base::fields(node);
}

void EnergyRosettaNonbond_O::dumpTerms(core::HashTable_sp atomTypes) {
  gctools::Vec0<EnergyRosettaNonbond>::iterator eni;
  string as1, as2;
  string str1, str2;
  core::clasp_write_string(fmt::format("Dumping {} terms\n", this->_Terms.size()));
  for (eni = this->_Terms.begin(); eni != this->_Terms.end(); eni++) {
    as1 = _rep_(eni->_Atom1_enb->getName());
    as2 = _rep_(eni->_Atom2_enb->getName());
    if (as1 < as2) {
      str1 = as1;
      str2 = as2;
    } else {
      str2 = as1;
      str1 = as2;
    }
    core::clasp_write_string(fmt::format("TERM ROSETTA_NONBOND {:<9} - {:<9}\n", str1, str2));
  }
}

void EnergyRosettaNonbond_O::callForEachTerm(core::Function_sp callback) {
  for (auto eni = this->_Terms.begin(); eni != this->_Terms.end(); eni++) {
    core::eval::funcall(callback, eni->_Atom1_enb, eni->_Atom2_enb,
                        core::make_fixnum(eni->term.i3x1),
                        core::make_fixnum(eni->term.i3x2));
  }
}


EnergyComponent_sp EnergyRosettaNonbond_O::copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc) {
  EnergyRosettaNonbond_sp copy = EnergyRosettaNonbond_O::create();
  copyEnergyComponent(copy, this->asSmartPtr());

  copy->_Parameters = this->_Parameters;

  copy->_AtomTable = this->_AtomTable;
  copy->_NonbondForceField = this->_NonbondForceField;
  copy->_AtomTypes = this->_AtomTypes;
  copy->_KeepInteractionFactory = keepInteractionFactory;
  if (this->_DisplacementBuffer.notnilp()) {
    copy->_DisplacementBuffer = copy_nvector(gc::As<NVector_sp>(this->_DisplacementBuffer));
  }

  core::T_sp keepInteraction = specializeKeepInteractionFactory(keepInteractionFactory, EnergyRosettaNonbond_O::staticClass());
  if (keepInteraction == _lisp->_true()) {
    for (auto edi = this->_Terms.begin(); edi != this->_Terms.end(); edi++) {
      copy->_Terms.push_back(*edi);
    }
  } else {
    for (auto edi = this->_Terms.begin(); edi != this->_Terms.end(); edi++) {
      Atom_sp a1 = edi->_Atom1_enb;
      Atom_sp a2 = edi->_Atom2_enb;
      size_t ia1 = edi->term.i3x1;
      size_t ia2 = edi->term.i3x2;
      if (skipInteraction_EnergyNonbond(keepInteraction, a1, a2, core::make_fixnum(ia1), core::make_fixnum(ia2))) continue;
      copy->_Terms.push_back(*edi);
    }
  }
  copy->_Parameters.do_apply(setupAcc);
  copy->_DisplacementBuffer = nil<core::T_O>();
  copy->_Terms.clear();
  return copy;
}

core::T_mv EnergyRosettaNonbond_O::maybeRebuildPairList(core::T_sp tcoordinates) {
  auto coords = gc::As<NVector_sp>(tcoordinates);
  if (this->_DisplacementBuffer.nilp()) {
    return this->rebuildPairList(tcoordinates);
  } else if (gc::IsA<NVector_sp>(this->_DisplacementBuffer)) {
    NVector_sp nvDisplacementBuffer = gc::As_unsafe<NVector_sp>(this->_DisplacementBuffer);
    if (nvDisplacementBuffer->size() != coords->size()) {
      SIMPLE_ERROR("The size of the _DispacementBuffer({}) MUST match the size of the coordinatess({})",
                   nvDisplacementBuffer->size(), coords->size());
    }
    double skinThickness = this->_Parameters.rpairlist - this->_Parameters.rcut;
    double movedTrigger = 0.5 * skinThickness;
    double movedTrigger2 = movedTrigger * movedTrigger;
    vecreal* raw_db = &(*nvDisplacementBuffer)[0];
    vecreal* raw_coords = &(*coords)[0];
    for (size_t ci = 0; ci < nvDisplacementBuffer->size(); ci += 3) {
      const vecreal& dx = raw_db[ci];
      const vecreal& dy = raw_db[ci + 1];
      const vecreal& dz = raw_db[ci + 2];
      const vecreal& cx = raw_coords[ci];
      const vecreal& cy = raw_coords[ci + 1];
      const vecreal& cz = raw_coords[ci + 2];
      if (dx == cx && dy == cy && dz == cz) continue;
      vecreal dist2 = (dx - cx) * (dx - cx) + (dy - cy) * (dy - cy) + (dz - cz) * (dz - cz);
      if (dist2 > movedTrigger2) {
        return this->rebuildPairList(tcoordinates);
      }
    }
    return Values0<core::T_O>();
  }
  SIMPLE_ERROR("{}: We should never get here", __FUNCTION__);
}

core::T_mv EnergyRosettaNonbond_O::rebuildPairList(core::T_sp tcoordinates) {
  this->_DisplacementBuffer = copy_nvector(gc::As<NVector_sp>(tcoordinates));
  size_t interactionsKept = 0;
  size_t interactionsDiscarded = 0;
  size_t totalInteractions = 0;
  this->_Terms.clear();
  if (this->_KeepInteractionFactory.notnilp()) {
    core::T_sp keepInteraction = specializeKeepInteractionFactory(this->_KeepInteractionFactory, EnergyRosettaNonbond_O::staticClass());
    bool hasKeepInteractionFunction = gc::IsA<core::Function_sp>(keepInteraction);
    double r_pairlist2 = this->_Parameters.rpairlist * this->_Parameters.rpairlist;
    const rosetta_nonbond_parameters& params = this->_Parameters;
    if (this->_AtomTable->getNumberOfAtoms() >= 2) {
      vecreal* coords = NULL;
      if (gc::IsA<NVector_sp>(tcoordinates)) {
        coords = &(*gc::As_unsafe<NVector_sp>(tcoordinates))[0];
      }
      gctools::Vec0<EnergyAtom>::iterator iea1;
      gctools::Vec0<EnergyAtom>::iterator iea2;
      for (iea1 = this->_AtomTable->begin(); iea1 != this->_AtomTable->end() - 1; iea1++) {
        Vector3 pos1;
        if (coords) {
          pos1.set(coords[iea1->coordinateIndexTimes3()],
                   coords[iea1->coordinateIndexTimes3() + 1],
                   coords[iea1->coordinateIndexTimes3() + 2]);
        } else {
          pos1 = iea1->atom()->getPosition();
        }
        for (iea2 = iea1 + 1; iea2 != this->_AtomTable->end(); iea2++) {
          if (!(iea1->inBondOrAngle(iea2->atom()) || iea1->relatedBy14(iea2->atom()))) {
            totalInteractions++;
            Vector3 pos2;
            if (coords) {
              pos2.set(coords[iea2->coordinateIndexTimes3()],
                       coords[iea2->coordinateIndexTimes3() + 1],
                       coords[iea2->coordinateIndexTimes3() + 2]);
            } else {
              pos2 = iea2->atom()->getPosition();
            }
            Vector3 posDelta = pos1 - pos2;
            double dist2 = posDelta.dotProduct(posDelta);
            if (dist2 < r_pairlist2) {
              bool keep = true;
              if (hasKeepInteractionFunction) {
                core::T_sp result = core::eval::funcall(keepInteraction,
                                                        iea1->atom(), iea2->atom(),
                                                        core::make_fixnum(iea1->coordinateIndexTimes3()),
                                                        core::make_fixnum(iea2->coordinateIndexTimes3()));
                keep = result.notnilp();
              }
              if (keep) {
                EnergyRosettaNonbond term;
                term.defineForAtomPair(this->_NonbondForceField,
                                       iea1->atom(), iea2->atom(),
                                       iea1->coordinateIndexTimes3(),
                                       iea2->coordinateIndexTimes3(),
                                       this->asSmartPtr(),
                                       this->_AtomTypes,
                                       keepInteraction,
                                       params);
                this->addTerm(term);
                ++interactionsKept;
              } else {
                ++interactionsDiscarded;
              }
            }
          }
        }
      }
    }
  }
  return Values(core::clasp_make_fixnum(interactionsKept),
                core::clasp_make_fixnum(interactionsDiscarded),
                core::clasp_make_fixnum(totalInteractions));
}

// Evaluate

double EnergyRosettaNonbond_O::evaluateAllComponent(ScoringFunction_sp score,
                                                    NVector_sp pos,
                                                    core::T_sp energyScale,
                                                    core::T_sp energyComponents,
                                                    bool calcForce,
                                                    gc::Nilable<NVector_sp> force,
                                                    bool calcDiagonalHessian,
                                                    bool calcOffDiagonalHessian,
                                                    gc::Nilable<AbstractLargeSquareMatrix_sp> hessian,
                                                    gc::Nilable<NVector_sp> hdvec,
                                                    gc::Nilable<NVector_sp> dvec,
                                                    core::T_sp activeAtomMask,
                                                    core::T_sp debugInteractions) {
  this->_Evaluations++;
  double energy = 0.0;
  size_t fails = 0;
  size_t index = 0;
  this->maybeRebuildPairList(pos);
  energy += template_evaluateUsingTerms<NoFiniteDifference>(this, this->_Terms, _sym_energyRosettaNonbond,
                                                            score, pos, energyScale, energyComponents, calcForce, force,
                                                            calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec,
                                                            activeAtomMask, debugInteractions, fails, index);
  return energy;
}

double EnergyRosettaNonbond_O::debugAllComponent(ScoringFunction_sp score,
                                                 NVector_sp pos,
                                                 core::T_sp energyScale,
                                                 core::T_sp energyComponents,
                                                 bool calcForce,
                                                 gc::Nilable<NVector_sp> force,
                                                 bool calcDiagonalHessian,
                                                 bool calcOffDiagonalHessian,
                                                 gc::Nilable<AbstractLargeSquareMatrix_sp> hessian,
                                                 gc::Nilable<NVector_sp> hdvec,
                                                 gc::Nilable<NVector_sp> dvec,
                                                 core::T_sp activeAtomMask,
                                                 core::T_sp debugInteractions) {
  return this->evaluateAllComponent(score, pos, energyScale, energyComponents, calcForce, force,
                                    calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec,
                                    activeAtomMask, debugInteractions);
}

}; // namespace chem
