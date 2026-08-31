/*
    File: energyRosettaElec.cc
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
#include <cando/chem/energyRosettaElec.h>
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
#include <cando/chem/pairList.h>
#include <clasp/core/wrappers.h>

namespace chem {

#include <cando/chem/energyKernels/rosetta_elec_cutoff.c>

EnergyRosettaElec_sp EnergyRosettaElec_O::make(EnergyFunction_sp energyFunction,
                                               core::T_sp keepInteractionFactory,
                                               SetupAccumulator& acc) {
  core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyRosettaElec_O::staticClass() );
  if (keepInteraction.notnilp()) {
    auto obj = EnergyRosettaElec_O::create();
    obj->_Parameters.do_apply(acc);
    obj->_AtomTable = energyFunction->_AtomTable;
    obj->_NonbondForceField = energyFunction->_AtomTable->nonbondForceFieldForAggregate();
    obj->_AtomTypes = energyFunction->_AtomTypes;
    obj->_KeepInteractionFactory = keepInteractionFactory;
    return obj;
  }
  SIMPLE_ERROR("Mismatch between keepInteractionFactory (says don't create EnergyRosettaElec_O) and EnergyRosettaElec_O::make which says make it");
}

std::string EnergyRosettaElec_O::implementation_details() const {
  Rosetta_Elec_Cutoff<NoHessian> elec;
  std::stringstream ss;
  ss << elec.description();
  return ss.str();
}

std::string EnergyRosettaElec_O::descriptionOfContents() const {
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

SYMBOL_EXPORT_SC_(ChemPkg, EnergyRosettaElec);
SYMBOL_EXPORT_SC_(ChemPkg, energyRosettaElec);

core::List_sp EnergyRosettaElec::encode() const
{
  ql::list ll;
  this->term.encode(ll);
  ll << INTERN_(kw, atom1) << this->_Atom1_enb
     << INTERN_(kw, atom2) << this->_Atom2_enb;
  return ll.cons();
}

void EnergyRosettaElec::decode(core::List_sp alist) { SIMPLE_ERROR("Implement decode of EnergyRosettaElec"); }

core::T_sp debug_rosetta_elec(double Energy, double x1, double y1, double z1, double x2, double y2, double z2,
                              const rosetta_elec_term& term) {
  double r2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
  ql::list ll;
  ll << INTERN_(kw, Energy) << mk_double_float(Energy) << INTERN_(kw, x1) << mk_double_float(x1) << INTERN_(kw, y1)
     << mk_double_float(y1) << INTERN_(kw, z1) << mk_double_float(z1) << INTERN_(kw, x2) << mk_double_float(x2)
     << INTERN_(kw, y2) << mk_double_float(y2) << INTERN_(kw, z2) << mk_double_float(z2)
     << INTERN_(kw, kqq) << mk_double_float(term.kqq) << INTERN_(kw, r2) << mk_double_float(r2);
  return ll.cons();
}

#define ELEC_DEBUG_INTERACTIONS(term)                                                                                            \
  if (doDebugInteractions) {                                                                                                      \
    core::eval::funcall(debugInteractions, chem::_sym_EnergyRosettaElec,                                                          \
                        debug_rosetta_elec(Energy,                                                                               \
                                           position[term.i3x1], position[term.i3x1 + 1], position[term.i3x1 + 2],                \
                                           position[term.i3x2], position[term.i3x2 + 1], position[term.i3x2 + 2],                \
                                           term),                                                                                \
                        core::make_fixnum(term.i3x1), core::make_fixnum(term.i3x2));                                              \
  }

struct NoFiniteDifference {
  static void maybeTestFiniteDifference(ScoringFunction_sp score, core::T_sp energyScale, int I1, int I2, core::T_sp activeAtomMask, double x1,
                                        double y1, double z1, double x2, double y2, double z2, double dA, double dC, double dQ1Q2,
                                        double fx1, double fy1, double fz1, double fx2, double fy2, double fz2, int index, size_t &fails,
                                        bool debugForce) {}
};

template <class MaybeFiniteDiff>
double template_evaluateUsingTerms(EnergyRosettaElec_O* mthis,
                                   const gctools::Vec0<EnergyRosettaElec>& terms,
                                   core::T_sp termSymbol,
                                   ScoringFunction_sp score, NVector_sp nvposition,
                                   core::T_sp energyScale, core::T_sp energyComponents,
                                   bool calcForce, gc::Nilable<NVector_sp> force, bool calcDiagonalHessian,
                                   bool calcOffDiagonalHessian, gc::Nilable<AbstractLargeSquareMatrix_sp> hessian,
                                   gc::Nilable<NVector_sp> hdvec, gc::Nilable<NVector_sp> dvec, core::T_sp activeAtomMask,
                                   core::T_sp debugInteractions, size_t& fails, size_t& index, bool debugForce = false) {
  rosetta_elec_parameters& params = mthis->_Parameters;
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  auto evalType = determineEnergyComponentEvalType(force, hdvec, dvec);
  double totalEnergy = 0.0;
  DOUBLE* position = &(*nvposition)[0];
  DOUBLE* rforce = NULL;
  DOUBLE* rdvec = NULL;
  DOUBLE* rhdvec = NULL;
  DOUBLE Energy = 0.0;
  Rosetta_Elec_Cutoff<NoHessian> elec;

#define KERNEL_TERM_ELEC_APPLY_ATOM_MASK(I1, I2)                                                                  \
  if (hasActiveAtomMask && !activeAtomMaskAnyAtomIsActive(bitvectorActiveAtomMask, I1, I2)) continue;

  if (evalType == energyEval) {
    for (auto si = terms.begin(); si != terms.end(); si++) {
      KERNEL_TERM_ELEC_APPLY_ATOM_MASK(si->term.i3x1, si->term.i3x2);
      Energy = elec.energy(params, si->term, position, &totalEnergy);
      ELEC_DEBUG_INTERACTIONS(si->term);
    }
  } else if (evalType == gradientEval) {
    rforce = &(*force)[0];
    for (auto si = terms.begin(); si != terms.end(); si++) {
      KERNEL_TERM_ELEC_APPLY_ATOM_MASK(si->term.i3x1, si->term.i3x2);
      Energy = elec.gradient(params, si->term, position, &totalEnergy, rforce);
      ELEC_DEBUG_INTERACTIONS(si->term);
    }
  } else {
    rforce = &(*force)[0];
    rdvec = &(*dvec)[0];
    rhdvec = &(*hdvec)[0];
    for (auto si = terms.begin(); si != terms.end(); si++) {
      KERNEL_TERM_ELEC_APPLY_ATOM_MASK(si->term.i3x1, si->term.i3x2);
      Energy = elec.hessian(params, si->term, position, &totalEnergy, rforce, NoHessian(), rdvec, rhdvec);
      ELEC_DEBUG_INTERACTIONS(si->term);
    }
  }
  maybeSetEnergy(energyComponents, termSymbol, totalEnergy);
  return totalEnergy;
}

bool EnergyRosettaElec::defineForAtomPair(core::T_sp forceField, Atom_sp a1, Atom_sp a2,
                                          size_t i3x1, size_t i3x2,
                                          EnergyRosettaElec_sp energyRosettaElec,
                                          core::HashTable_sp atomTypes, core::T_sp keepInteraction,
                                          const rosetta_elec_parameters& params) {
  this->_Atom1_enb = a1;
  this->_Atom2_enb = a2;

  double amber_charge_conversion_18dot2223 =
      core::Number_O::as_double_float(gc::As<core::Number_sp>(_sym_STARamber_charge_conversion_18_DOT_2223STAR->symbolValue()));
  double dQ1Q2Scale = amber_charge_conversion_18dot2223 * amber_charge_conversion_18dot2223;
  double electrostaticScale = 1.0;

  double charge1 = a1->getCharge();
  double charge2 = a2->getCharge();
  double kqq = calculate_dQ1Q2(electrostaticScale, dQ1Q2Scale, charge1, charge2);

  this->term = rosetta_elec_term(params, kqq, i3x1, i3x2);
  return true;
}

void EnergyRosettaElec_O::initialize() { this->Base::initialize(); }

void EnergyRosettaElec_O::ensureParameterCache() {
  AtomTable_sp at = this->_AtomTable;
  if (at.nilp()) return;
  size_t n = at->getNumberOfAtoms();
  if (this->_CachedForAtomTable == at && this->_CachedCharge.size() == n) return;  // still valid

  // Hoisted out of the per-pair path: defineForAtomPair did this symbolValue
  // lookup and boxed-number unbox on EVERY atom pair.
  double conv = core::Number_O::as_double_float(
      gc::As<core::Number_sp>(_sym_STARamber_charge_conversion_18_DOT_2223STAR->symbolValue()));
  this->_DQ1Q2Scale = conv * conv;

  this->_CachedCharge.assign(n, 0.0);
  auto& energyAtoms = at->getVectorEnergyAtoms();
  for (size_t i = 0; i < n; i++) {
    this->_CachedCharge[i] = energyAtoms[i].atom()->getCharge();
  }

  // rosetta_elec_term is exactly linear in kqq - every field is a linear
  // combination of e_rmin/e_rlow/de_rlow/e_rhi/de_rhi with params-only
  // coefficients, and each of those is proportional to kqq.  So one prototype
  // at kqq==1 serves every pair and the three exp() calls happen once.
  this->_Prototype = rosetta_elec_term(this->_Parameters, 1.0, 0, 0);
  this->_PrototypeValid = true;
  this->_CachedForAtomTable = at;
}

void EnergyRosettaElec_O::addTerm(const EnergyRosettaElec& term) { this->_Terms.push_back(term); }

void EnergyRosettaElec_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, terms), this->_Terms);
  node->field(INTERN_(kw, AtomTable), this->_AtomTable);
  node->field(INTERN_(kw, NonbondForceField), this->_NonbondForceField);
  node->field(INTERN_(kw, AtomTypes), this->_AtomTypes);
  this->_Parameters.fields(node);
  this->Base::fields(node);
}

/*! ATOMS first, then their I3 values - the convention EnergyComponent_O::atomsForEachTerm
 *  documents.  A caller taking &rest reads any component without knowing which one it has.
 */
void EnergyRosettaElec_O::atomsForEachTerm(core::Function_sp callback) {
  for (auto eni = this->_Terms.begin(); eni != this->_Terms.end(); eni++) {
    core::eval::funcall(callback, eni->_Atom1_enb,
                          eni->_Atom2_enb,
                          core::make_fixnum(eni->term.i3x1),
                          core::make_fixnum(eni->term.i3x2));
  }
}

void EnergyRosettaElec_O::dumpTerms(core::HashTable_sp atomTypes) {
  gctools::Vec0<EnergyRosettaElec>::iterator eni;
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
    core::clasp_write_string(fmt::format("TERM ROSETTA_ELEC {:<9} - {:<9}\n", str1, str2));
  }
}

void EnergyRosettaElec_O::callForEachTerm(core::Function_sp callback) {
  for (auto eni = this->_Terms.begin(); eni != this->_Terms.end(); eni++) {
    core::eval::funcall(callback, eni->_Atom1_enb, eni->_Atom2_enb,
                        core::make_fixnum(eni->term.i3x1),
                        core::make_fixnum(eni->term.i3x2));
  }
}

EnergyComponent_sp EnergyRosettaElec_O::copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc) {
  EnergyRosettaElec_sp copy = EnergyRosettaElec_O::create();
  copyEnergyComponent(copy, this->asSmartPtr());

  copy->_Parameters = this->_Parameters;
  copy->_AtomTable = this->_AtomTable;
  copy->_NonbondForceField = this->_NonbondForceField;
  copy->_AtomTypes = this->_AtomTypes;
  copy->_Matter1 = this->_Matter1;
  copy->_Matter2 = this->_Matter2;
  copy->_KeepInteractionFactory = keepInteractionFactory;
  copy->_Parameters.do_apply(setupAcc);
  copy->invalidatePairList();
  copy->_Terms.clear();
  return copy;
}

CL_LAMBDA((self chem:energy-rosetta-elec) mat1 mat2 energy-function keep-interaction-factory);
CL_DEFMETHOD void EnergyRosettaElec_O::constructNonbondTermsBetweenMatters(Matter_sp mat1, Matter_sp mat2,
                                                                            EnergyFunction_sp energyFunction,
                                                                            core::T_sp keepInteractionFactory) {
  this->_Matter1 = mat1;
  this->_Matter2 = mat2;
  this->_KeepInteractionFactory = keepInteractionFactory;
  this->_Terms.clear();
  this->_AtomTable = energyFunction->_AtomTable;
  this->_AtomTypes = energyFunction->atomTypes();
  this->_NonbondForceField = this->_AtomTable->nonbondForceFieldForAggregate();
  this->invalidatePairList();
}

core::T_mv EnergyRosettaElec_O::rebuildPairListBetweenMatters(core::T_sp tcoordinates) {
  core::T_sp keepInteractionFactory = this->_KeepInteractionFactory;
  if (keepInteractionFactory.nilp()) return Values0<core::T_O>();
  NVector_sp coords = gc::As<NVector_sp>(tcoordinates);
  core::T_sp keepInteraction = specializeKeepInteractionFactory(keepInteractionFactory, EnergyRosettaElec_O::staticClass());
  Matter_sp mat1 = gc::As<Matter_sp>(this->_Matter1);
  Matter_sp mat2 = gc::As<Matter_sp>(this->_Matter2);
  bool hasKeepInteractionFunction = gc::IsA<core::Function_sp>(keepInteraction);
  rosetta_elec_parameters& params = this->_Parameters;
  double rpairlist2 = params.rpairlist * params.rpairlist;
  auto atomTable = this->_AtomTable;
  size_t interactionsKept = 0;
  size_t interactionsDiscarded = 0;
  {
    this->_Terms.clear();
    Loop lMat1(mat1, ATOMS);
    while (lMat1.advanceLoopAndProcess()) {
      Atom_sp a1 = lMat1.getAtom();
      size_t i3x1 = atomTable->getCoordinateIndexTimes3(a1);
      Vector3 v1(coords, i3x1, Safe());
      Loop lMat2(mat2, ATOMS);
      while (lMat2.advanceLoopAndProcess()) {
        Atom_sp a2 = lMat2.getAtom();
        size_t i3x2 = atomTable->getCoordinateIndexTimes3(a2);
        Vector3 v2(coords, i3x2, Safe());
        Vector3 vdiff = v1 - v2;
        double dist2 = vdiff.dotProduct(vdiff);
        if (dist2 < rpairlist2) {
          if (hasKeepInteractionFunction) {
            core::T_sp result = core::eval::funcall(keepInteraction, a1, a2,
                                                    core::make_fixnum(i3x1),
                                                    core::make_fixnum(i3x2));
            if (result.notnilp()) {
              EnergyRosettaElec term;
              term.defineForAtomPair(this->_NonbondForceField, a1, a2,
                                     i3x1, i3x2,
                                     this->asSmartPtr(),
                                     this->_AtomTypes,
                                     keepInteraction,
                                     params);
              this->addTerm(term);
              ++interactionsKept;
            } else {
              ++interactionsDiscarded;
            }
          } else {
            EnergyRosettaElec term;
            term.defineForAtomPair(this->_NonbondForceField, a1, a2,
                                   i3x1, i3x2,
                                   this->asSmartPtr(),
                                   this->_AtomTypes,
                                   keepInteraction,
                                   params);
            this->addTerm(term);
            ++interactionsKept;
          }
        }
      }
    }
  }
  size_t totalInteractions = interactionsKept + interactionsDiscarded;
  return Values(core::clasp_make_fixnum(interactionsKept),
                core::clasp_make_fixnum(interactionsDiscarded),
                core::clasp_make_fixnum(totalInteractions));
}

core::T_mv EnergyRosettaElec_O::maybeRebuildPairList(core::T_sp tcoordinates) {
  // Shared implementation - see maybeRebuildPairListImpl in pairList.h.
  return maybeRebuildPairListImpl(this, tcoordinates);
}

core::T_mv EnergyRosettaElec_O::rebuildPairList(core::T_sp tcoordinates) {
#if 1
  return rebuildPairListImpl(this, tcoordinates);
#else
  this->_DisplacementBuffer = copy_nvector(gc::As<NVector_sp>(tcoordinates));
  if (this->_Matter1.notnilp()) {
    ASSERT(this->_Matter2.notnilp());
    return this->rebuildPairListBetweenMatters(tcoordinates);
  }
  size_t interactionsKept = 0;
  size_t interactionsDiscarded = 0;
  size_t totalInteractions = 0;
  this->_Terms.clear();
  if (this->_KeepInteractionFactory.notnilp()) {
    core::T_sp keepInteraction = specializeKeepInteractionFactory(this->_KeepInteractionFactory, EnergyRosettaElec_O::staticClass());
    bool hasKeepInteractionFunction = gc::IsA<core::Function_sp>(keepInteraction);
    rosetta_elec_parameters& params = this->_Parameters;
    double r_pairlist2 = params.rpairlist*params.rpairlist;
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
                EnergyRosettaElec term;
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
#endif
}

// Evaluate

double EnergyRosettaElec_O::evaluateAllComponent(ScoringFunction_sp score,
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
  energy += template_evaluateUsingTerms<NoFiniteDifference>(this, this->_Terms, _sym_energyRosettaElec,
                                                            score, pos, energyScale, energyComponents, calcForce, force,
                                                            calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec,
                                                            activeAtomMask, debugInteractions, fails, index);
  return energy;
}

void EnergyRosettaElec_O::setupHessianPreconditioner(NVector_sp nvPosition, AbstractLargeSquareMatrix_sp m, core::T_sp activeAtomMask) {
  return; // not used for preconditioner
}

double EnergyRosettaElec_O::debugAllComponent(ScoringFunction_sp score,
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
