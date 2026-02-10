/*
    File: energyRosettaLKSolvation.cc
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
#include <cando/chem/energyRosettaLKSolvation.h>
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

#include <cando/chem/energyKernels/rosetta_lk_solvation.c>

EnergyRosettaLKSolvation_sp EnergyRosettaLKSolvation_O::make(AtomTable_sp atomTable, core::T_sp nbForceField,
                                                             core::HashTable_sp atomTypes, core::T_sp keepInteractionFactory,
                                                             SetupAccumulator& acc,
                                                             core::T_sp tcoordinates) {
  auto obj = EnergyRosettaLKSolvation_O::create();
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

std::string EnergyRosettaLKSolvation_O::implementation_details() const {
  Rosetta_Lk_Solvation<NoHessian> lk;
  std::stringstream ss;
  ss << lk.description();
  return ss.str();
}

std::string EnergyRosettaLKSolvation_O::descriptionOfContents() const {
  stringstream ss;
  ss << ":enabled " << ((this->_Enabled) ? "T" : "NIL");
  ss << " number-of-terms " << this->_Terms.size();
  return ss.str();
}

}

// #define DEBUG_LK_SOLVATION_TERM 1
#define LOG_ENERGY(...)
// #define LOG_ENERGY(...) core::clasp_write_string(fmt::format(__VA_ARGS__))

namespace chem {

SYMBOL_EXPORT_SC_(ChemPkg, EnergyRosettaLKSolvation);
SYMBOL_EXPORT_SC_(ChemPkg, energyRosettaLKSolvation);

core::List_sp EnergyRosettaLKSolvation::encode() const {
  ql::list ll;
  this->term.encode(ll);
  ll << INTERN_(kw, atom1) << this->_Atom1_enb
     << INTERN_(kw, atom2) << this->_Atom2_enb;
  return ll.cons();
}

void EnergyRosettaLKSolvation::decode(core::List_sp alist) { SIMPLE_ERROR("Implement decode of EnergyRosettaLKSolvation"); }

core::T_sp debug_rosetta_lk_solvation(double Energy, double x1, double y1, double z1, double x2, double y2, double z2,
                                      const rosetta_lk_solvation_term& term) {
  double r2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
  ql::list ll;
  ll << INTERN_(kw, Energy) << mk_double_float(Energy) << INTERN_(kw, x1) << mk_double_float(x1) << INTERN_(kw, y1)
     << mk_double_float(y1) << INTERN_(kw, z1) << mk_double_float(z1) << INTERN_(kw, x2) << mk_double_float(x2)
     << INTERN_(kw, y2) << mk_double_float(y2) << INTERN_(kw, z2) << mk_double_float(z2)
     << INTERN_(kw, sigmaij) << mk_double_float(term.sigmaij)
     << INTERN_(kw, pref_i_j) << mk_double_float(term.pref_i_j)
     << INTERN_(kw, pref_j_i) << mk_double_float(term.pref_j_i)
     << INTERN_(kw, r2) << mk_double_float(r2);
  return ll.cons();
}

#define LK_SOLVATION_DEBUG_INTERACTIONS(term)                                                                                     \
  if (doDebugInteractions) {                                                                                                      \
    core::eval::funcall(debugInteractions, chem::_sym_EnergyRosettaLKSolvation,                                                  \
                        debug_rosetta_lk_solvation(Energy,                                                                       \
                                                   position[term.i3x1], position[term.i3x1 + 1], position[term.i3x1 + 2],       \
                                                   position[term.i3x2], position[term.i3x2 + 1], position[term.i3x2 + 2],       \
                                                   term),                                                                        \
                        core::make_fixnum(term.i3x1), core::make_fixnum(term.i3x2));                                              \
  }

struct NoFiniteDifference {
  static void maybeTestFiniteDifference(ScoringFunction_sp score, core::T_sp energyScale, int I1, int I2, core::T_sp activeAtomMask, double x1,
                                        double y1, double z1, double x2, double y2, double z2, double dA, double dC, double dQ1Q2,
                                        double fx1, double fy1, double fz1, double fx2, double fy2, double fz2, int index, size_t &fails,
                                        bool debugForce) {}
};

template <class MaybeFiniteDiff>
double template_evaluateUsingTerms(EnergyRosettaLKSolvation_O* mthis,
                                   const gctools::Vec0<EnergyRosettaLKSolvation>& terms,
                                   core::T_sp termSymbol,
                                   ScoringFunction_sp score, NVector_sp nvposition,
                                   core::T_sp energyScale, core::T_sp energyComponents,
                                   bool calcForce, gc::Nilable<NVector_sp> force, bool calcDiagonalHessian,
                                   bool calcOffDiagonalHessian, gc::Nilable<AbstractLargeSquareMatrix_sp> hessian,
                                   gc::Nilable<NVector_sp> hdvec, gc::Nilable<NVector_sp> dvec, core::T_sp activeAtomMask,
                                   core::T_sp debugInteractions, size_t& fails, size_t& index, bool debugForce = false) {
  const rosetta_lk_solvation_parameters& params = mthis->_Parameters;
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  auto evalType = determineEnergyComponentEvalType(force, hdvec, dvec);
  double totalEnergy = 0.0;
  DOUBLE* position = &(*nvposition)[0];
  DOUBLE* rforce = NULL;
  DOUBLE* rdvec = NULL;
  DOUBLE* rhdvec = NULL;
  DOUBLE Energy = 0.0;
  Rosetta_Lk_Solvation<NoHessian> lk;

#define KERNEL_TERM_LK_SOLVATION_APPLY_ATOM_MASK(I1, I2)                                                          \
  if (hasActiveAtomMask && !(bitvectorActiveAtomMask->testBit(I1 / 3) && bitvectorActiveAtomMask->testBit(I2 / 3))) continue;

  if (evalType == energyEval) {
    for (auto si = terms.begin(); si != terms.end(); si++) {
      KERNEL_TERM_LK_SOLVATION_APPLY_ATOM_MASK(si->term.i3x1, si->term.i3x2);
      Energy = lk.energy(params, si->term, position, &totalEnergy);
      LK_SOLVATION_DEBUG_INTERACTIONS(si->term);
    }
  } else if (evalType == gradientEval) {
    rforce = &(*force)[0];
    for (auto si = terms.begin(); si != terms.end(); si++) {
      KERNEL_TERM_LK_SOLVATION_APPLY_ATOM_MASK(si->term.i3x1, si->term.i3x2);
      Energy = lk.gradient(params, si->term, position, &totalEnergy, rforce);
      LK_SOLVATION_DEBUG_INTERACTIONS(si->term);
    }
  } else {
    rforce = &(*force)[0];
    rdvec = &(*dvec)[0];
    rhdvec = &(*hdvec)[0];
    for (auto si = terms.begin(); si != terms.end(); si++) {
      KERNEL_TERM_LK_SOLVATION_APPLY_ATOM_MASK(si->term.i3x1, si->term.i3x2);
      Energy = lk.hessian(params, si->term, position, &totalEnergy, rforce, NoHessian(), rdvec, rhdvec);
      LK_SOLVATION_DEBUG_INTERACTIONS(si->term);
    }
  }
  maybeSetEnergy(energyComponents, termSymbol, totalEnergy);
  return totalEnergy;
}

static bool lookup_lk_solvation_parameters(FFLKSolvation_sp ffLKSolvation,
                                           double& lk_dgfree,
                                           double& lk_lambda,
                                           double& lj_radius,
                                           double& lk_volume) {
  // TODO: Wire real LK solvation per-atom parameters when available.
  lk_dgfree = ffLKSolvation->getLK_dgfree();
  lk_lambda = ffLKSolvation->getLK_lambda();
  lk_volume = ffLKSolvation->getLK_volume();
  lj_radius = ffLKSolvation->getLJ_radius();
  return true;
}

bool EnergyRosettaLKSolvation::defineForAtomPair(core::T_sp forceField, Atom_sp a1, Atom_sp a2,
                                                 size_t i3x1, size_t i3x2,
                                                 EnergyRosettaLKSolvation_sp energyRosettaLKSolvation,
                                                 core::HashTable_sp atomTypes, core::T_sp keepInteraction,
                                                 const rosetta_lk_solvation_parameters& params) {
  this->_Atom1_enb = a1;
  this->_Atom2_enb = a2;
  core::Symbol_sp t1 = a1->getPropertyOrDefault(INTERN_(kw,lk_solvation_atom_type),nil<core::Symbol_O>());
  core::Symbol_sp t2 = a2->getPropertyOrDefault(INTERN_(kw,lk_solvation_atom_type),nil<core::Symbol_O>());
  ASSERT(forceField && forceField.notnilp());
  core::T_sp tffLKSolvation1 = core::eval::funcall(_sym_find_lksolvation_type, forceField, t1);
  core::T_sp tffLKSolvation2 = core::eval::funcall(_sym_find_lksolvation_type, forceField, t2);
  if (tffLKSolvation1.nilp() || tffLKSolvation2.nilp()) {
    SIMPLE_ERROR("Could not find one of the LKSolvation parameters");
  }
  auto ffLKSolvation1 = gc::As<FFLKSolvation_sp>(tffLKSolvation1);
  auto ffLKSolvation2 = gc::As<FFLKSolvation_sp>(tffLKSolvation2);

  double lk_dgfree_i = 0.0;
  double lk_lambda_i = 1.0;
  double lj_radius_i = 0.0;
  double lk_volume_i = 0.0;
  bool ok1 = lookup_lk_solvation_parameters(ffLKSolvation1, lk_dgfree_i, lk_lambda_i, lj_radius_i, lk_volume_i);

  double lk_dgfree_j = 0.0;
  double lk_lambda_j = 1.0;
  double lj_radius_j = 0.0;
  double lk_volume_j = 0.0;
  bool ok2 = lookup_lk_solvation_parameters(ffLKSolvation2, lk_dgfree_j, lk_lambda_j, lj_radius_j, lk_volume_j);

  this->term = rosetta_lk_solvation_term(params,
                                         lk_dgfree_i, lk_lambda_i, lj_radius_i, lk_volume_i,
                                         lk_dgfree_j, lk_lambda_j, lj_radius_j, lk_volume_j,
                                         i3x1, i3x2);
  return (ok1 && ok2);
}

void EnergyRosettaLKSolvation_O::initialize() { this->Base::initialize(); }

void EnergyRosettaLKSolvation_O::addTerm(const EnergyRosettaLKSolvation& term) { this->_Terms.push_back(term); }

void EnergyRosettaLKSolvation_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, terms), this->_Terms);
  node->field(INTERN_(kw, AtomTable), this->_AtomTable);
  node->field(INTERN_(kw, KeepInteractionFactory), this->_KeepInteractionFactory);
  this->_Parameters.fields(node);
  this->Base::fields(node);
}

void EnergyRosettaLKSolvation_O::dumpTerms(core::HashTable_sp atomTypes) {
  gctools::Vec0<EnergyRosettaLKSolvation>::iterator eni;
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
    core::clasp_write_string(fmt::format("TERM ROSETTA_LK_SOLVATION {:<9} - {:<9}\n", str1, str2));
  }
}

void EnergyRosettaLKSolvation_O::callForEachTerm(core::Function_sp callback) {
  for (auto eni = this->_Terms.begin(); eni != this->_Terms.end(); eni++) {
    core::eval::funcall(callback, eni->_Atom1_enb, eni->_Atom2_enb,
                        core::make_fixnum(eni->term.i3x1),
                        core::make_fixnum(eni->term.i3x2));
  }
}

EnergyComponent_sp EnergyRosettaLKSolvation_O::copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc) {
  EnergyRosettaLKSolvation_sp copy = EnergyRosettaLKSolvation_O::create();
  copyEnergyComponent(copy, this->asSmartPtr());

  copy->_Parameters = this->_Parameters;

  copy->_AtomTable = this->_AtomTable;
  copy->_NonbondForceField = this->_NonbondForceField;
  copy->_AtomTypes = this->_AtomTypes;
  copy->_KeepInteractionFactory = keepInteractionFactory;
  if (this->_DisplacementBuffer.notnilp()) {
    copy->_DisplacementBuffer = copy_nvector(gc::As<NVector_sp>(this->_DisplacementBuffer));
  }

  core::T_sp keepInteraction = specializeKeepInteractionFactory(keepInteractionFactory, EnergyRosettaLKSolvation_O::staticClass());
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

core::T_mv EnergyRosettaLKSolvation_O::maybeRebuildPairList(core::T_sp tcoordinates) {
  auto coords = gc::As<NVector_sp>(tcoordinates);
  if (this->_DisplacementBuffer.nilp()) {
    return this->rebuildPairList(tcoordinates);
  } else if (gc::IsA<NVector_sp>(this->_DisplacementBuffer)) {
    NVector_sp nvDisplacementBuffer = gc::As_unsafe<NVector_sp>(this->_DisplacementBuffer);
    if (nvDisplacementBuffer->size() != coords->size()) {
      SIMPLE_ERROR("The size of the _DispacementBuffer({}) MUST match the size of the coordinatess({})",
                   nvDisplacementBuffer->size(), coords->size());
    }
    double r_pairlist = this->_Parameters.rpairlist;
    double r_cut = this->_Parameters.r_solv_high;
    double skinThickness = r_pairlist - r_cut;
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

core::T_mv EnergyRosettaLKSolvation_O::rebuildPairList(core::T_sp tcoordinates) {
  this->_DisplacementBuffer = copy_nvector(gc::As<NVector_sp>(tcoordinates));
  size_t interactionsKept = 0;
  size_t interactionsDiscarded = 0;
  size_t totalInteractions = 0;
  this->_Terms.clear();
  if (this->_KeepInteractionFactory.notnilp()) {
    core::T_sp keepInteraction = specializeKeepInteractionFactory(this->_KeepInteractionFactory, EnergyRosettaLKSolvation_O::staticClass());
    bool hasKeepInteractionFunction = gc::IsA<core::Function_sp>(keepInteraction);
    double r_pairlist = this->_Parameters.rpairlist;
    double r_pairlist2 = r_pairlist * r_pairlist;
    const rosetta_lk_solvation_parameters& params = this->_Parameters;
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
                EnergyRosettaLKSolvation term;
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

double EnergyRosettaLKSolvation_O::evaluateAllComponent(ScoringFunction_sp score,
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
  energy += template_evaluateUsingTerms<NoFiniteDifference>(this, this->_Terms, _sym_energyRosettaLKSolvation,
                                                            score, pos, energyScale, energyComponents, calcForce, force,
                                                            calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec,
                                                            activeAtomMask, debugInteractions, fails, index);
  return energy;
}

double EnergyRosettaLKSolvation_O::debugAllComponent(ScoringFunction_sp score,
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
