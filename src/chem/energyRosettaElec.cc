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
#include <clasp/core/wrappers.h>

namespace chem {

#include <cando/chem/energyKernels/rosetta_elec_cutoff.c>

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
  ll << INTERN_(kw, kqq) << core::clasp_make_double_float(this->term.kqq)
     << INTERN_(kw, e_rmin) << core::clasp_make_double_float(this->term.e_rmin)
     << INTERN_(kw, aa_low) << core::clasp_make_double_float(this->term.aa_low)
     << INTERN_(kw, bb_low) << core::clasp_make_double_float(this->term.bb_low)
     << INTERN_(kw, cc_low) << core::clasp_make_double_float(this->term.cc_low)
     << INTERN_(kw, dd_low) << core::clasp_make_double_float(this->term.dd_low)
     << INTERN_(kw, aa_high) << core::clasp_make_double_float(this->term.aa_high)
     << INTERN_(kw, bb_high) << core::clasp_make_double_float(this->term.bb_high)
     << INTERN_(kw, cc_high) << core::clasp_make_double_float(this->term.cc_high)
     << INTERN_(kw, dd_high) << core::clasp_make_double_float(this->term.dd_high)
     << INTERN_(kw, i1) << core::make_fixnum(this->term.i3x1)
     << INTERN_(kw, i2) << core::make_fixnum(this->term.i3x2)
     << INTERN_(kw, atom1) << this->_Atom1_enb
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

inline double calculate_dQ1Q2(double electrostaticScale, double electrostaticModifier, double charge1, double charge2) {
  return electrostaticScale * electrostaticModifier * charge1 * charge2;
}

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
  if (hasActiveAtomMask && !(bitvectorActiveAtomMask->testBit(I1 / 3) && bitvectorActiveAtomMask->testBit(I2 / 3))) continue;

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

void EnergyRosettaElec_O::addTerm(const EnergyRosettaElec& term) { this->_Terms.push_back(term); }

void EnergyRosettaElec_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, terms), this->_Terms);
  node->field(INTERN_(kw, AtomTable), this->_AtomTable);
  node->field(INTERN_(kw, KeepInteractionFactory), this->_KeepInteractionFactory);
  node->field(INTERN_(kw, elec_weight), this->_Parameters.elec_weight);
  node->field(INTERN_(kw, eps_core), this->_Parameters.eps_core);
  node->field(INTERN_(kw, eps_solvent), this->_Parameters.eps_solvent);
  node->field(INTERN_(kw, rmin), this->_Parameters.rmin);
  node->field(INTERN_(kw, rlow), this->_Parameters.rlow);
  node->field(INTERN_(kw, rhi), this->_Parameters.rhi);
  node->field(INTERN_(kw, rcut), this->_Parameters.rcut);
  node->field(INTERN_(kw, rpairlist), this->_Parameters.rpairlist);
  node->field(INTERN_(kw, InteractionsKept), this->_InteractionsKept);
  node->field(INTERN_(kw, InteractionsDiscarded), this->_InteractionsDiscarded);
  this->Base::fields(node);
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

void EnergyRosettaElec_O::constructNonbondTermsFromAtomTable(AtomTable_sp atomTable, core::T_sp nbForceField,
                                                             core::HashTable_sp atomTypes, core::T_sp keepInteractionFactory,
                                                             core::T_sp tcoordinates) {
  if (keepInteractionFactory.nilp()) return;
  this->_AtomTable = atomTable;
  this->_NonbondForceField = nbForceField;
  this->_AtomTypes = atomTypes;
  this->_KeepInteractionFactory = keepInteractionFactory;
  if (tcoordinates.notnilp()) {
    this->rebuildPairList(tcoordinates);
  }
}

EnergyRosettaElec_sp EnergyRosettaElec_O::copyFilter(core::T_sp keepInteractionFactory) {
  EnergyRosettaElec_sp copy = EnergyRosettaElec_O::create();
  copyEnergyComponent(copy, this->asSmartPtr());

  copy->_Parameters = this->_Parameters;
  copy->_AtomTable = this->_AtomTable;
  copy->_NonbondForceField = this->_NonbondForceField;
  copy->_AtomTypes = this->_AtomTypes;
  copy->_KeepInteractionFactory = keepInteractionFactory;
  copy->_Nonbond_r_pairlist = this->_Nonbond_r_pairlist;
  copy->_InteractionsKept = this->_InteractionsKept;
  copy->_InteractionsDiscarded = this->_InteractionsDiscarded;
  if (this->_DisplacementBuffer.notnilp()) {
    copy->_DisplacementBuffer = copy_nvector(gc::As<NVector_sp>(this->_DisplacementBuffer));
  }

  core::T_sp keepInteraction = specializeKeepInteractionFactory(keepInteractionFactory, EnergyRosettaElec_O::staticClass());
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
  return copy;
}

core::T_mv EnergyRosettaElec_O::maybeRebuildPairList(core::T_sp tcoordinates) {
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

core::T_mv EnergyRosettaElec_O::rebuildPairList(core::T_sp tcoordinates) {
  this->_DisplacementBuffer = copy_nvector(gc::As<NVector_sp>(tcoordinates));
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
  this->_InteractionsKept = interactionsKept;
  this->_InteractionsDiscarded = interactionsDiscarded;
  return Values(core::clasp_make_fixnum(interactionsKept),
                core::clasp_make_fixnum(interactionsDiscarded),
                core::clasp_make_fixnum(totalInteractions));
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
