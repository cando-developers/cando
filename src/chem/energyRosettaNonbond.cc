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

#include <atomic>
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
#include <cando/chem/pairList.h>
#include <clasp/core/wrappers.h>

namespace chem {

#include <cando/chem/energyKernels/rosetta_nonbond_dd_cutoff.c>

EnergyRosettaNonbond_sp EnergyRosettaNonbond_O::make(EnergyFunction_sp energyFunction,
                                                     core::T_sp keepInteractionFactory,
                                                     SetupAccumulator& acc) {
  core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyRosettaNonbond_O::staticClass() );
  if (keepInteraction.notnilp()) {
    auto obj = EnergyRosettaNonbond_O::create();
    obj->_Parameters.do_apply(acc);
    obj->_AtomTable = energyFunction->_AtomTable;
    obj->_NonbondForceField = energyFunction->_AtomTable->nonbondForceFieldForAggregate();
    obj->_AtomTypes = energyFunction->_AtomTypes;
    obj->_KeepInteractionFactory = keepInteractionFactory;
    return obj;
  }
  SIMPLE_ERROR("Mismatch between keepInteractionFactory (says don't create EnergyRosettaNonbond_O) and EnergyRosettaNonbond_O::make which says make it");
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
SYMBOL_EXPORT_SC_(ChemPkg, energyRosettaNonbondFaRep);

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
  typedef Rosetta_Nonbond_Dd_Cutoff<NoHessian> KernelType;
  // acc[0] is the total energy; acc[1] is the fa_rep decomposition.
  double acc[KernelType::EnergyAccumulatorSize] = {0.0};
  DOUBLE* position = &(*nvposition)[0];
  DOUBLE* rforce = NULL;
  DOUBLE* rdvec = NULL;
  DOUBLE* rhdvec = NULL;
  DOUBLE Energy = 0.0;
  KernelType nonbond;

#define KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(I1, I2)                     \
  if (hasActiveAtomMask && !(bitvectorActiveAtomMask->testBit(I1 / 3) && bitvectorActiveAtomMask->testBit(I2 / 3))) continue;

  if (evalType == energyEval) {
    for (auto si = terms.begin(); si != terms.end(); si++) {
      KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(si->term.i3x1, si->term.i3x2);
      Energy = nonbond.energy(params, si->term, position, acc);
      NONBOND_DEBUG_INTERACTIONS(si->term);
    }
  } else if (evalType == gradientEval) {
    rforce = &(*force)[0];
    for (auto si = terms.begin(); si != terms.end(); si++) {
      KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(si->term.i3x1, si->term.i3x2);
      Energy = nonbond.gradient(params, si->term, position, acc, rforce);
      NONBOND_DEBUG_INTERACTIONS(si->term);
    }
  } else {
    rforce = &(*force)[0];
    rdvec = &(*dvec)[0];
    rhdvec = &(*hdvec)[0];
    for (auto si = terms.begin(); si != terms.end(); si++) {
      KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(si->term.i3x1, si->term.i3x2);
      Energy = nonbond.hessian(params, si->term, position, acc, rforce, NoHessian(), rdvec, rhdvec);
      NONBOND_DEBUG_INTERACTIONS(si->term);
    }
  }
  maybeSetEnergy(energyComponents, termSymbol, acc[0]);
  maybeSetEnergy(energyComponents, _sym_energyRosettaNonbondFaRep,
                 acc[KernelType::EnergyAccumulatorIndex_fa_rep]);
  mthis->_LastFaRep = acc[KernelType::EnergyAccumulatorIndex_fa_rep]; 
  return acc[0];
}

bool EnergyRosettaNonbond::defineForAtomPair(core::T_sp forceField, Atom_sp a1, Atom_sp a2,
                                             size_t i3x1, size_t i3x2,
                                             EnergyRosettaNonbond_sp energyRosettaNonbond,
                                             core::HashTable_sp atomTypes, core::T_sp keepInteraction,
                                             const rosetta_nonbond_parameters& params) {
  this->_Atom1_enb = a1;
  this->_Atom2_enb = a2;
  core::T_sp t1 = a1->getType(atomTypes);
  core::T_sp t2 = a2->getType(atomTypes);
  ASSERT(forceField && forceField.notnilp());
  core::T_sp tffNonbond1 = core::eval::funcall(_sym_find_type, forceField, t1);
  core::T_sp tffNonbond2 = core::eval::funcall(_sym_find_type, forceField, t2);
  if (tffNonbond1.nilp() || tffNonbond2.nilp()) return false;
  FFNonbond_sp ffNonbond1 = gc::As<FFNonbond_sp>(tffNonbond1);
  FFNonbond_sp ffNonbond2 = gc::As<FFNonbond_sp>(tffNonbond2);

  double parmA, parmC;
  if (!combineNonbondParams(ffNonbond1->getRadius_Angstroms(), ffNonbond1->getEpsilon_kcal(),
                            ffNonbond2->getRadius_Angstroms(), ffNonbond2->getEpsilon_kcal(),
                            parmA, parmC)) return false;
  this->term = rosetta_nonbond_term(params, parmA, parmC, i3x1, i3x2);
  return true;
}

void EnergyRosettaNonbond_O::initialize() { this->Base::initialize(); }

void EnergyRosettaNonbond_O::addTerm(const EnergyRosettaNonbond& term) { this->_Terms.push_back(term); }

void EnergyRosettaNonbond_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, terms), this->_Terms);
  node->field(INTERN_(kw, AtomTable), this->_AtomTable);
  node->field(INTERN_(kw, NonbondForceField), this->_NonbondForceField);
  node->field(INTERN_(kw, AtomTypes), this->_AtomTypes);
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

void EnergyRosettaNonbond_O::atomsForEachTerm(core::Function_sp callback) {
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
  copy->_Matter1 = this->_Matter1;
  copy->_Matter2 = this->_Matter2;
  copy->_KeepInteractionFactory = keepInteractionFactory;
  copy->_Parameters.do_apply(setupAcc);
  copy->invalidatePairList();
  copy->_Terms.clear();
  copy->invalidateParameterCache();
  return copy;
}

CL_LAMBDA((self chem:energy-rosetta-nonbond) mat1 mat2 energy-function keep-interaction-factory);
CL_DEFMETHOD void EnergyRosettaNonbond_O::constructNonbondTermsBetweenMatters(Matter_sp mat1, Matter_sp mat2,
                                                                               EnergyFunction_sp energyFunction,
                                                                               core::T_sp keepInteractionFactory) {
  this->invalidateParameterCache();
  this->_Matter1 = mat1;
  this->_Matter2 = mat2;
  this->_KeepInteractionFactory = keepInteractionFactory;
  this->_Terms.clear();
  this->_AtomTable = energyFunction->_AtomTable;
  this->_AtomTypes = energyFunction->atomTypes();
  this->_NonbondForceField = this->_AtomTable->nonbondForceFieldForAggregate();
  this->invalidatePairList();
}

core::T_mv EnergyRosettaNonbond_O::maybeRebuildPairList(core::T_sp tcoordinates) {
  // Shared implementation - see maybeRebuildPairListImpl in pairList.h.  The threshold is still
  // 0.5*(rpairlist - rcut) off _Parameters; only the way the drift is measured has changed.
  return maybeRebuildPairListImpl(this, tcoordinates);
}

void EnergyRosettaNonbond_O::ensureParameterCache() {
  AtomTable_sp at = this->_AtomTable;
  if (at.nilp()) return;
  size_t n = at->getNumberOfAtoms();
  // The generation test is what makes the shared table safe: the atom-table pointer can be
  // unchanged while its slots have been rebuilt underneath us.
  if (this->_CachedForAtomTable == at
      && this->_CachedNBGeneration == at->_NBGeneration
      && this->_TypeSlot.size() == n) return;   // still valid

  FFNonbondDb_sp db = gc::As<FFNonbondDb_sp>(this->_NonbondForceField);
  auto& energyAtoms = at->getVectorEnergyAtoms();

  // ---- SHARED across every component over this atom table ----
  //
  // An atom's nonbond type is a property of the ATOM and the force field, not of the component
  // asking, so this mapping is identical for all of a blueprint's slot-group components.  Building
  // it per component cost TWO hash lookups per atom - Atom_O::getType is a gethash (atom.cc:1712)
  // and chem__FFNonbond_findType is a second one keyed on the type STRING - and every slot group
  // owns its own copies of these components, each starting with an invalid cache.  That is ~2820
  // cold starts over ~6700 atoms, ~40 million lookups, measured at 49% of the single-scan fill.
  // Built once here it is ~6700.
  if (!at->nbTypeSlotsValidFor(this->_NonbondForceField)) {
    at->_NBTypeSlot.assign(n, -1);
    at->_NBUniq.clear();
    // One lookup per atom, then collapse identical (radius,epsilon) into slots.  Two distinct types
    // with identical parameters merging is harmless - they produce the same term.
    for (size_t i = 0; i < n; i++) {
      core::T_sp type = energyAtoms[i].atom()->getType(this->_AtomTypes);
      core::T_sp tff  = chem__FFNonbond_findType(db, type);
      if (tff.nilp()) continue;                     // slot stays -1
      FFNonbond_sp ff = gc::As<FFNonbond_sp>(tff);
      double r = ff->getRadius_Angstroms();
      double e = ff->getEpsilon_kcal();
      int slot = -1;
      size_t nslots = at->_NBUniq.size() / 2;       // 2 doubles per slot - see energyAtomTable.h
      for (size_t s = 0; s < nslots; ++s)
        if (at->_NBUniq[2*s+0] == r && at->_NBUniq[2*s+1] == e) { slot = (int)s; break; }
      if (slot < 0) {
        slot = (int)nslots;
        at->_NBUniq.push_back(r);
        at->_NBUniq.push_back(e);
      }
      at->_NBTypeSlot[i] = slot;
    }
    at->_NBCachedForceField = this->_NonbondForceField;
    at->_NBGeneration++;   // every rebuild is a new generation, so stale copies are detectable
  }

  // Copy rather than alias, so the accessors in the header keep indexing this->_TypeSlot
  // unchanged.  n ints per component is nothing next to the lookups just avoided.  Element-wise
  // rather than whole-vector assignment - gctools::Vec0 is not std::vector.
  this->_TypeSlot.assign(n, -1);
  for (size_t i = 0; i < n; ++i) this->_TypeSlot[i] = at->_NBTypeSlot[i];

  // ---- PER COMPONENT: the term cache depends on _Parameters, which is NOT shared ----
  size_t nt = at->_NBUniq.size() / 2;      // 2 doubles per slot
  this->_NTypeSlots = nt;
  this->_TermCache.assign(nt*nt, rosetta_nonbond_term());
  this->_TermCacheValid.assign(nt*nt, 0);
  for (size_t s1 = 0; s1 < nt; ++s1)
    for (size_t s2 = 0; s2 < nt; ++s2) {
      double parmA, parmC;
      if (combineNonbondParams(at->_NBUniq[2*s1+0], at->_NBUniq[2*s1+1],
                               at->_NBUniq[2*s2+0], at->_NBUniq[2*s2+1], parmA, parmC)) {
        this->_TermCache[s1*nt + s2] =
            rosetta_nonbond_term(this->_Parameters, parmA, parmC, 0, 0);
        this->_TermCacheValid[s1*nt + s2] = 1;
      }
    }
  // Some atom types carry no Lennard-Jones parameters at all: AMBER gives HO
  // and HW radius 0 and epsilon 0, putting all of the sterics on the heavy
  // atom they are bonded to.  combineNonbondParams refuses those pairs because
  // the term derives sigma6 = 2A/C and epsilon = C^2/(4A), which is 0/0 when
  // both are zero.  Skipping such a pair yields exactly the zero the term
  // would have evaluated to, so no energy is lost - and electrostatics lives
  // in EnergyRosettaElec_O, so no charge interaction goes with it.
  //
  // Report once per process rather than once per energy function: hundreds of
  // scan energy functions are built, but a *new* offender is still worth
  // seeing.  One zero-parameter type invalidates its whole row and column, so
  // b such types account for nt^2 - (nt-b)^2 pairs; anything above that count
  // means a pair failed for some other reason and deserves a look.
  {
    static std::atomic<bool> s_reported{false};
    size_t invalid = 0;
    for (size_t k = 0; k < nt*nt; ++k) if (!this->_TermCacheValid[k]) ++invalid;
    if (invalid && !s_reported.exchange(true)) {
      fmt::print(stderr, "RosettaNonbond: {} of {} type pairs have zero Lennard-Jones parameters"
                         " - skipped, each contributes exactly zero\n", invalid, nt*nt);
      size_t zeroTypes = 0;
      for (size_t s = 0; s < nt; ++s) {
        size_t bad = 0;
        for (size_t s2 = 0; s2 < nt; ++s2) if (!this->_TermCacheValid[s*nt + s2]) ++bad;
        if (bad != nt) continue;                          // this slot combines with nothing
        ++zeroTypes;
        for (size_t i = 0; i < n; i++) {
          if (this->_TypeSlot[i] == (int)s) {
            core::T_sp type = energyAtoms[i].atom()->getType(this->_AtomTypes);
            fmt::print(stderr, "   slot {}: radius={:.4f} epsilon={:.4f} type={} example-atom={}\n",
                       s, at->_NBUniq[2*s+0], at->_NBUniq[2*s+1],
                       _rep_(type), _rep_(energyAtoms[i].atom()));
            break;
          }
        }
      }
      size_t good = nt - zeroTypes;
      size_t accounted = nt*nt - good*good;
      if (invalid != accounted)
        fmt::print(stderr, "   NOTE: {} zero-parameter type(s) account for only {} of the {} invalid"
                           " pairs - the rest fail to combine for another reason\n",
                   zeroTypes, accounted, invalid);
    }
  }
  this->_CachedForAtomTable = at;
  this->_CachedNBGeneration = at->_NBGeneration;
}



core::T_mv EnergyRosettaNonbond_O::rebuildPairList(core::T_sp tcoordinates) {
  return rebuildPairListImpl(this, tcoordinates);
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

void EnergyRosettaNonbond_O::setupHessianPreconditioner(NVector_sp nvPosition, AbstractLargeSquareMatrix_sp m, core::T_sp activeAtomMask) {
  return; // not used for preconditioner
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
