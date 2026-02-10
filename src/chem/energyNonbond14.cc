/*
    File: energyNonbond.cc
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
#include <cando/chem/energyNonbond14.h>
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

#include <cando/chem/energyKernels/amber_nonbond14.c>


EnergyNonbond14_sp EnergyNonbond14_O::make(SetupAccumulator& acc) {
  auto obj = EnergyNonbond14_O::create();
  return obj;
}

std::string EnergyNonbond14_O::implementation_details() const {
  Amber_Nonbond14<NoHessian> nonbond;

  std::stringstream ss;
  ss << nonbond.description();
  return ss.str();
}

std::string EnergyNonbond14_O::descriptionOfContents() const {
  stringstream ss;
  ss << ":enabled " << ((this->_Enabled) ? "T" : "NIL");
  ss << " number-of-14terms " << this->_Terms14.size();
  return ss.str();
}


size_t EnergyNonbond14_O::runTestCalls(core::T_sp stream, chem::NVector_sp coords) const
{
  #define POS_SIZE 6
  double energy_new;
  double energy_ground;
  double position[POS_SIZE];
  double force_new[POS_SIZE];
  double force_ground[POS_SIZE];
  double hessian_new[POS_SIZE*POS_SIZE];
  double hessian_ground[POS_SIZE*POS_SIZE];
  double dvec_new[POS_SIZE];
  double dvec_ground[POS_SIZE];
  double hdvec_new[POS_SIZE];
  double hdvec_ground[POS_SIZE];
  size_t idx=0;
  size_t errs = 0;
  Amber_Nonbond14<double*> nonbond14;
  for ( auto si=this->_Terms14.begin();
        si!=this->_Terms14.end(); si++ ) {
    position[0] = coords[si->term.i3x1];
    position[1] = coords[si->term.i3x1+1];
    position[2] = coords[si->term.i3x1+2];
    position[3] = coords[si->term.i3x2];
    position[4] = coords[si->term.i3x2+1];
    position[5] = coords[si->term.i3x2+2];
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    nonbond14.gradient( si->term, position, &energy_new, force_new );
    nonbond14.gradient_fd( si->term, position, &energy_ground, force_ground );
    if (!test_match( stream, "nonbond_dd_cutoff_force", POS_SIZE,
                     force_new, force_ground,
                     0, 0,
                     0, 0 )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH nonbond_dd_cutoff_gradient #{} a = {}  c = {}\n",
                              idx, si->term.a, si->term.c ), stream );
    }
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    nonbond14.hessian( si->term, position, &energy_new, force_new, hessian_new, dvec_new, hdvec_new );
    nonbond14.hessian_fd( si->term, position, &energy_ground, force_ground, hessian_ground, dvec_ground, hdvec_ground );
    if (!test_match( stream, "nonbond_dd_cutoff_hessian", POS_SIZE,
                     force_new, force_ground,
                     hessian_new, hessian_ground,
                     hdvec_new, hdvec_ground )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH nonbond_dd_cutoff_hessian #{} a = {}  c = {}\n",
                              idx, si->term.a, si->term.c ), stream );
    }
    idx++;
  }
  core::print(fmt::format("{} errors = {}\n", Amber_Nonbond14<double*>().description(), errs), stream);
  return errs;
}



}



// #define DEBUG_NONBOND_TERM 1
#define LOG_ENERGY(...)
// #define LOG_ENERGY(...) core::clasp_write_string(fmt::format(__VA_ARGS__))

namespace chem {

#if 0
CL_DEFMETHOD
void EnergyNonbond14_O::set_nonbond_pairlist_parameters( double r_switch, double r_cut, double r_pairlist, double distance_dielectric ) {
  this->_Nonbond_r_switch = r_switch;
  this->_Nonbond_r_cut = r_cut;
  this->_Nonbond_r_pairlist = r_pairlist;
  if (distance_dielectric<VERYSMALL) {
    SIMPLE_ERROR("distance-dielectric must be a reasonable value - passed {}", distance_dielectric );
  }
  this->_Nonbond_invdd = 1.0/distance_dielectric;
}
#endif


EnergyFunction_sp energyFunctionNonbond14Parameters(ScoringFunction_sp score,
                                                  core::T_sp energyScale,
                                                  double &dielectricConstant,
                                                  double &dQ1Q2Scale,
                                                  double &cutoff) {
  auto energyFunction = gc::As<EnergyFunction_sp>(score);
  cutoff = energyScaleNonbondCutoff(energyScale);
  dielectricConstant = energyScaleDielectricConstant(energyScale);
  double amber_charge_conversion_18dot2223 =
      core::Number_O::as_double_float(gc::As<core::Number_sp>(_sym_STARamber_charge_conversion_18_DOT_2223STAR->symbolValue()));
  dQ1Q2Scale = amber_charge_conversion_18dot2223 * amber_charge_conversion_18dot2223;
  return energyFunction;
}

core::T_sp debug_nonbond14(double Energy, double x1, double y1, double z1, double x2, double y2, double z2, 
                           double kqq, double dA, double dC ) {
  double r2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
  ql::list ll;
  ll << INTERN_(kw, Energy) << mk_double_float(Energy) << INTERN_(kw, x1) << mk_double_float(x1) << INTERN_(kw, y1)
     << mk_double_float(y1) << INTERN_(kw, z1) << mk_double_float(z1) << INTERN_(kw, x2) << mk_double_float(x2) << INTERN_(kw, y2)
     << mk_double_float(y2) << INTERN_(kw, z2) << mk_double_float(z2)
     << INTERN_(kw, kqq) << mk_double_float(kqq)
     << INTERN_(kw, a) << mk_double_float(dA)
     << INTERN_(kw, c) << mk_double_float(dC)
     << INTERN_(kw, r2) << mk_double_float(r2)
      ;
  return ll.cons();
}

SYMBOL_EXPORT_SC_(ChemPkg, EnergyNonbond14);

core::T_sp nonbond14_type(bool is14) {
  return chem::_sym_EnergyNonbond14;
}
#define NONBOND_APPLY_ATOM_MASK(I1, I2)                                                                                            \
  if (hasActiveAtomMask && !(bitvectorActiveAtomMask->testBit(I1 / 3) && bitvectorActiveAtomMask->testBit(I2 / 3)))                \
    goto SKIP_term;

/* No periodic boundary conditons/bounding-box is used in this code.
   So make these macros do nothing.
   See energyPeriodicBoundaryConditionsNonbond.cc for counter-example.
*/
core::List_sp EnergyNonbond14::encode() const {
  ql::list ll;
  this->term.encode(ll);
  ll << INTERN_(kw, charge1) << core::clasp_make_double_float(this->_Atom1_enb->getCharge())
     << INTERN_(kw, charge2) << core::clasp_make_double_float(this->_Atom2_enb->getCharge())
     << INTERN_(kw,atom1) << this->_Atom1_enb
     << INTERN_(kw,atom2) << this->_Atom2_enb
    ;
  return ll.cons();
}

void EnergyNonbond14::decode(core::List_sp alist) { SIMPLE_ERROR("Implement decode of EnergyNonbond"); }

#define ENERGY_FUNCTION score, energyScale, I1, I2, activeAtomMask

#define NONBOND_DEBUG_INTERACTIONS(term)                                                                                         \
  if (doDebugInteractions) {                                                                                                       \
    core::eval::funcall(debugInteractions, nonbond14_type(interactionsAre14),                                                        \
                        debug_nonbond14(Energy, \
                                        position[term.i3x1],position[term.i3x1+1],position[term.i3x1+2], \
                                        position[term.i3x2],position[term.i3x2+1],position[term.i3x2+2], \
                                        term.kqq, term.a, term.c), \
                        core::make_fixnum(term.i3x1), core::make_fixnum(term.i3x2));                                                             \
  }

double _evaluateEnergyOnly_Nonbond(ScoringFunction_sp score, core::T_sp energyScale, int I1, int I2, core::T_sp activeAtomMask, num_real x1, num_real y1,
                                   num_real z1, num_real x2, num_real y2, num_real z2, num_real dA, num_real dC, num_real dQ1Q2);

struct NoFiniteDifference {
  static void maybeTestFiniteDifference(ScoringFunction_sp score, core::T_sp energyScale, int I1, int I2, core::T_sp activeAtomMask, double x1, double y1,
                                        double z1, double x2, double y2, double z2, double dA, double dC, double dQ1Q2, double fx1,
                                        double fy1, double fz1, double fx2, double fy2, double fz2, int index, size_t &fails,
                                        bool debugForce) {}
};


inline double calculate_dQ1Q2(double electrostaticScale, double electrostaticModifier, double charge1, double charge2) {
  return electrostaticScale * electrostaticModifier * charge1 * charge2;
}

/*! The main nonbond code using excluded atoms.
    It is a template function so that template arguments can inline or elide testing code.
 */

/*! The core nonbond code using pairwise terms.
    It as a template function so that template arguments can inline or elide testing code.
 */
template <class MaybeFiniteDiff>
double template_evaluateUsingTerms(EnergyNonbond14_O *mthis,
                                   bool interactionsAre14,
                                   const gctools::Vec0<EnergyNonbond14>& terms,
                                   core::T_sp termSymbol,
                                   ScoringFunction_sp score, NVector_sp nvposition,
                                   core::T_sp energyScale, core::T_sp energyComponents,
                                   bool calcForce, gc::Nilable<NVector_sp> force, bool calcDiagonalHessian,
                                   bool calcOffDiagonalHessian, gc::Nilable<AbstractLargeSquareMatrix_sp> hessian,
                                   gc::Nilable<NVector_sp> hdvec, gc::Nilable<NVector_sp> dvec, core::T_sp activeAtomMask,
                                   core::T_sp debugInteractions, size_t &fails, size_t &index, bool debugForce = false) {
  double dielectricConstant;
  double dQ1Q2Scale;
  double cutoff;
  EnergyFunction_sp energyFunction = energyFunctionNonbond14Parameters(score, energyScale, dielectricConstant, dQ1Q2Scale, cutoff);
  double nonbondCutoffSquared = cutoff * cutoff;
#define CUTOFF_SQUARED nonbondCutoffSquared
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  double energyElectrostatic = 0.0;
  KahanSummation energyVdw;
  double energyVdw14 = 0.0;
  double energyElectrostatic14 = 0.0;
  double vdwScale = energyScaleVdwScale(energyScale);
  double eelScale = energyScaleElectrostaticScale(energyScale);
  double DIELECTRIC = energyScaleDielectricConstant(energyScale);
  auto evalType = determineEnergyComponentEvalType(force,hdvec,dvec);
  double totalNonbondEnergy = 0.0;
  DOUBLE* position = &(*nvposition)[0];
  DOUBLE* rforce = NULL;
  DOUBLE* rhessian = NULL; // &(*hessian)[0];
  DOUBLE* rdvec = NULL;
  DOUBLE* rhdvec = NULL;
  DOUBLE  Energy = 0.0;
  Amber_Nonbond14<NoHessian> nonbond14;

  // If you are going to use openmp here, you need to control access to the force and hessian
  // arrays so that only one thread updates each element at a time.
  {

#define KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(I1, I2)                     \
    if (hasActiveAtomMask && !(bitvectorActiveAtomMask->testBit(I1 / 3) && bitvectorActiveAtomMask->testBit(I2 / 3)))  continue;

    if (evalType==energyEval) {
      for (auto si = terms.begin(); si != terms.end(); si++ ) {
        KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(si->term.i3x1,si->term.i3x2);
        Energy = nonbond14.energy( si->term,
                                   position,
                                   &totalNonbondEnergy);
        NONBOND_DEBUG_INTERACTIONS(si->term);
      }
    } else if (evalType==gradientEval) {
      rforce = &(*force)[0];
      for (auto si = terms.begin(); si != terms.end(); si++ ) {
        KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(si->term.i3x1,si->term.i3x2);
        Energy = nonbond14.gradient(si->term,
                                    position,
                                    &totalNonbondEnergy,
                                    rforce );
        NONBOND_DEBUG_INTERACTIONS(si->term);
      }
    } else { // if (hasHdAndD) {
      rforce = &(*force)[0];
      rdvec = &(*dvec)[0];
      rhdvec = &(*hdvec)[0];
      for (auto si = terms.begin(); si != terms.end(); si++ ) {
        KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(si->term.i3x1,si->term.i3x2);
        Energy = nonbond14.hessian( si->term,
                                    position,
                                    &totalNonbondEnergy,
                                    rforce,
                                    NoHessian(),
                                    rdvec,
                                    rhdvec);
        NONBOND_DEBUG_INTERACTIONS(si->term);
      }
    }
  }
  maybeSetEnergy(energyComponents, termSymbol, totalNonbondEnergy );
  return totalNonbondEnergy;
}

SYMBOL_EXPORT_SC_(ChemPkg, find_type);

double EnergyNonbond14_O::evaluateAllComponent(ScoringFunction_sp score,
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
  //  printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  double energy = 0.0;
  size_t fails = 0;
  size_t index = 0;
  // Evaluate the 1-4 terms
  energy += template_evaluateUsingTerms<NoFiniteDifference>(this, true, this->_Terms14, _sym_energyNonbond14, 
                                                            score, pos, energyScale, energyComponents, calcForce, force,
                                                            calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec,
                                                            activeAtomMask, debugInteractions, fails, index);
  return energy;
}


#if 0
double _evaluateEnergyOnly_Nonbond(ScoringFunction_sp score, core::T_sp energyScale, int I1, int I2, core::T_sp activeAtomMask, num_real x1, num_real y1,
                                   num_real z1, num_real x2, num_real y2, num_real z2, num_real dA, num_real dC, num_real dQ1Q2) {
  double dielectricConstant;
  double dQ1Q2Scale;
  double cutoff;
  auto energyFunction = energyFunctionNonbondParameters(score, energyScale, dielectricConstant, dQ1Q2Scale, cutoff);
  double vdwScale = energyScaleVdwScale(energyScale);
  double eelScale = energyScaleElectrostaticScale(energyScale);
  double DIELECTRIC = energyScaleDielectricConstant(energyScale);
#define CUTOFF_SQUARED (cutoff * cutoff)
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
#undef NONBOND_SET_PARAMETER
#define NONBOND_SET_PARAMETER(x)                                                                                                   \
  {}
#undef NONBOND_SET_POSITION
#define NONBOND_SET_POSITION(x, ii, of)                                                                                            \
  {}
#undef NONBOND_ENERGY_ACCUMULATE
#define NONBOND_ENERGY_ACCUMULATE(e)                                                                                               \
  {}
#undef NONBOND_FORCE_ACCUMULATE
#define NONBOND_FORCE_ACCUMULATE(i, o, v)                                                                                          \
  {}
#undef NONBOND_DIAGONAL_HESSIAN_ACCUMULATE
#define NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(i1, o1, i2, o2, v)                                                                     \
  {}
#undef NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1, o1, i2, o2, v)                                                                 \
  {}
#undef NONBOND_EVDW_ENERGY_ACCUMULATE
#define NONBOND_EVDW_ENERGY_ACCUMULATE(x)                                                                                          \
  {}
#undef NONBOND_EEEL_ENERGY_ACCUMULATE
#define NONBOND_EEEL_ENERGY_ACCUMULATE(x)                                                                                          \
  {}

#undef NONBOND_CALC_FORCE // Don't calculate FORCE or HESSIAN

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_Nonbond_termDeclares.cc>
#pragma clang diagnostic pop
#define NONBOND_DEBUG_INTERACTIONS(term)
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
#undef NONBOND_DEBUG_INTERACTIONS
#undef CUTOFF_SQUARED
#undef DIELECTRIC
  return Energy;
}
#endif




bool EnergyNonbond14::defineForAtomPair(core::T_sp forceField, Atom_sp a1, Atom_sp a2,
                                        size_t i3x1, size_t i3x2,
                                        EnergyNonbond14_sp energyNonbond,
                                        core::HashTable_sp atomTypes, core::T_sp keepInteraction ) {
  double amber_charge_conversion_18dot2223 =
      core::Number_O::as_double_float(gc::As<core::Number_sp>(_sym_STARamber_charge_conversion_18_DOT_2223STAR->symbolValue()));
  double dQ1Q2Scale = amber_charge_conversion_18dot2223 * amber_charge_conversion_18dot2223;
  num_real epsilonij;
  num_real vdwScale;
  num_real electrostaticScale;
  this->_Atom1_enb = a1;
  this->_Atom2_enb = a2;
  core::Symbol_sp t1 = a1->getType(atomTypes);
  core::Symbol_sp t2 = a2->getType(atomTypes);
  LOG("Defining nonbond between types: {} - {}", _rep_(t1), _rep_(t2));
  ASSERT(forceField && forceField.notnilp());
  LOG("forceField @{}   .notnilp()->{}", (void *)forceField.raw_(), forceField.notnilp());
  core::T_sp tffNonbond1 = core::eval::funcall(_sym_find_type, forceField, t1);
  core::T_sp tffNonbond2 = core::eval::funcall(_sym_find_type, forceField, t2);
  if (tffNonbond1.nilp()) {
    //     	SIMPLE_ERROR(("Unknown force field type(",iea1->_Atom->getType(atomTypes).c_str(),") for non-bonded interaction"));
    return false;
  }
  if (tffNonbond2.nilp()) {
    //     	SIMPLE_ERROR(("Unknown force field type(",iea2->_Atom->getType(atomTypes).c_str(),") for non-bonded interaction"));
    return false;
  }
  // See AMBER manual 1-4 Non-Bonded Interaction Scaling
  // Section 15.1.6
  // vdwScale is 1/2 of vdw
  // electrostaticScale is 1.0/1.2 of electrostatic
  vdwScale = 1.0 / 2.0;
  electrostaticScale = 1.0 / 1.2;
  LOG("vdwScale = {}", (double)(vdwScale));
  LOG("electrostaticScale = {}", (double)(electrostaticScale));
  FFNonbond_sp ffNonbond1 = gc::As<FFNonbond_sp>(tffNonbond1);
  FFNonbond_sp ffNonbond2 = gc::As<FFNonbond_sp>(tffNonbond2);
  double term_a, term_c, term_kqq;
  {
    double rstar = ffNonbond1->getRadius_Angstroms() + ffNonbond2->getRadius_Angstroms();
    epsilonij = sqrt(ffNonbond1->getEpsilon_kcal() * ffNonbond2->getEpsilon_kcal());
    double rtothe6th = rstar*rstar*rstar*rstar*rstar*rstar;
    double rtothe12th = rtothe6th*rtothe6th;
    double parmA = epsilonij * rtothe12th; // pow(this->_RStar, 12.0);
    double parmC = 2.0 * epsilonij * rtothe6th; // pow(this->_RStar, 6.0);
    term_a = parmA * vdwScale;
    term_c = parmC * vdwScale;
  }
  {
    double charge1 = a1->getCharge();
    double charge2 = a2->getCharge();
    term_kqq = calculate_dQ1Q2(electrostaticScale, dQ1Q2Scale, charge1, charge2 );
    LOG("Calc dQ1Q2 electrostaticScale= {}", (double)(electrostaticScale));
    LOG("Calc dQ1Q2 Charge1 = {}", (double)(charge1));
    LOG("Calc dQ1Q2 Charge2 = {}", (double)(charge2));
    LOG("dQ1Q2 = {}", (double)(term_kqq));
  }
  this->term = AmberNonbond14Term( term_a, term_c, term_kqq, i3x1, i3x2 );
  return true;
}


void EnergyNonbond14_O::dumpTerms(core::HashTable_sp atomTypes) {
  IMPLEMENT_ME();
  #if 0
  gctools::Vec0<Amber_Nonbond14>::iterator eni;
  string as1, as2, as3, as4;
  string str1, str2, str3, str4;
  core::clasp_write_string(fmt::format("Dumping {} terms\n", this->_Terms.size()));
  for (eni = this->_Terms14.begin(); eni != this->_Terms14.end(); eni++) {
    as1 = _rep_(eni->_Atom1_enb->getName());
    as2 = _rep_(eni->_Atom2_enb->getName());
    if (as1 < as2) {
      str1 = as1;
      str2 = as2;
    } else {
      str2 = as1;
      str1 = as2;
    }
    core::clasp_write_string(fmt::format("TERM 4CALC14 {:<9} - {:<9}\n", str1, str2));
  }
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
    core::clasp_write_string(fmt::format("TERM 5NONBOND {:<9} - {:<9}\n", str1, str2));
  }
#endif
}


SYMBOL_EXPORT_SC_(KeywordPkg, nonbond);
SYMBOL_EXPORT_SC_(KeywordPkg, vdwRadii);
SYMBOL_EXPORT_SC_(KeywordPkg, distance);
SYMBOL_EXPORT_SC_(KeywordPkg, force);

void EnergyNonbond14_O::initialize() {
  this->Base::initialize();
}

void EnergyNonbond14_O::addTerm14(const EnergyNonbond14 &term) { this->_Terms14.push_back(term); }

void EnergyNonbond14_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, terms14), this->_Terms14);
  node->field(INTERN_(kw, AtomTable), this->_AtomTable);
  node->field(INTERN_(kw, ntypes), this->_ntypes);
  node->field(INTERN_(kw, atom_name_vector), this->_atom_name_vector);         // atom-name-vector
  node->field(INTERN_(kw, charge_vector), this->_charge_vector);               // charge-vector
  node->field(INTERN_(kw, mass_vector), this->_mass_vector);                   // masses
  node->field(INTERN_(kw, atomic_number_vector), this->_atomic_number_vector); // vec
  node->field(INTERN_(kw, KeepInteractionFactory), this->_KeepInteractionFactory);
  this->Base::fields(node);
}

void EnergyNonbond14_O::construct14InteractionTerms(AtomTable_sp atomTable, Matter_sp matter, core::T_sp forceField,
                                                  core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes) {
  if (keepInteractionFactory.nilp()) return;
  core::T_sp keepInteraction = specializeKeepInteractionFactory(keepInteractionFactory,EnergyNonbond14_O::staticClass());
  bool hasKeepInteractionFunction = gc::IsA<core::Function_sp>(keepInteraction);
  EnergyNonbond14 energyNonbond;
  Loop loop;
  loop.loopTopGoal(matter, PROPERS);
  size_t terms(0);
  // Don't add duplicate terms
  std::set<pair<size_t, size_t>> duplicates;
  while (loop.advanceLoopAndProcess()) {
    Atom_sp a1 = loop.getAtom1();
    Atom_sp a4 = loop.getAtom4();
    auto ea1 = atomTable->getEnergyAtomPointer(a1);
    auto ea4 = atomTable->getEnergyAtomPointer(a4);
    size_t ia1 = ea1->_IndexTimes3;
    size_t ia4 = ea4->_IndexTimes3;
    if (skipInteraction_EnergyNonbond(keepInteraction, a1, a4, core::make_fixnum(ia1), core::make_fixnum(ia4))) continue;
    if (ia1 > ia4) {
      size_t ttt = ia1;
      ia1 = ia4;
      ia4 = ttt;
    }
    auto ia1ia4 = std::make_pair(ia1, ia4);
    if (!duplicates.contains(ia1ia4)) {
      duplicates.insert(ia1ia4);
      LOG("About to addTerm");
      energyNonbond.defineForAtomPair(forceField,
                                      ea1->atom(), ea4->atom(),
                                      ea1->coordinateIndexTimes3(),
                                      ea4->coordinateIndexTimes3(),
                                      this->asSmartPtr(),
                                      atomTypes,
                                      keepInteraction );
      this->addTerm14(energyNonbond);
      ++terms;
    }
  }
  if (chem__verbose(0))
    core::clasp_write_string(fmt::format("Built 14 interaction table with {} terms\n", terms));
}




CL_DOCSTRING(R"doc(For each nonbond term call the CALLBACK function.
The CALLBACK function must accept two arguments:
1. The first atom of the nonbond term
2. The second atom of the nonbond term
3. The Ix3 coordinate for the first atom
4. The Ix3 coordinate for the second atom)doc");
CL_DEFMETHOD void EnergyNonbond14_O::callForEachTerm(core::Function_sp callback) {
  for (auto eni = this->_Terms14.begin(); eni != this->_Terms14.end(); eni++) {
    core::eval::funcall(callback,eni->_Atom1_enb,eni->_Atom2_enb,
                        core::make_fixnum(eni->term.i3x1),
                        core::make_fixnum(eni->term.i3x2));
  }
}



EnergyComponent_sp EnergyNonbond14_O::copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc) {
  EnergyNonbond14_sp copy = EnergyNonbond14_O::create();
  copyEnergyComponent( copy, this->asSmartPtr() );
  core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyNonbond14_O::staticClass() );

#define COPY_FIELD(xxx) { copy->xxx = this->xxx; }
  COPY_FIELD(_AtomTable);
  COPY_FIELD(_ntypes);
  COPY_FIELD(_NonbondForceField);
  COPY_FIELD(_AtomTypes);
  COPY_FIELD(_atom_name_vector);
  COPY_FIELD(_charge_vector);
  COPY_FIELD(_mass_vector);
  COPY_FIELD(_atomic_number_vector);
  COPY_FIELD(_NumberOfExcludedAtomIndexes);
  COPY_FIELD(_ExcludedAtomIndexes);
  COPY_FIELD(_InteractionsKept);
  COPY_FIELD(_InteractionsDiscarded);

  // Copy the interaction factory
  copy->_KeepInteractionFactory = keepInteractionFactory;
  if (keepInteraction==_lisp->_true()) {
    // Copy the 1-4 interactions
    for (auto edi = this->_Terms14.begin(); edi != this->_Terms14.end(); edi++) {
      copy->_Terms14.push_back(*edi);
    }
  } else {
    SIMPLE_ERROR("EnergyNonbond uses excluded atoms and we don't support keepInteraction = {}", _rep_(keepInteraction));
  }
  return copy;
}




}; // namespace chem
