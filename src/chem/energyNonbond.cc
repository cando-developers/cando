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
#define DEBUG_LEVEL_FULL

#include <clasp/core/foundation.h>
#include <clasp/core/bformat.h>
#include <cando/chem/energyNonbond.h>
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

#include "cando/chem/energyKernels/nonbond.c"

#if 1
#include "cando/chem/energyKernels/nonbond_dd_cutoff_gpt.c"
#else
#include "cando/chem/energyKernels/nonbond_dd_cutoff.c"
#endif

}




namespace chem {

std::string EnergyNonbond_O::description() const {
  Nonbond<NoHessian> nonbond;
  Nonbond_Dd_Cutoff<NoHessian> nonbond_dd_cutoff;

  std::stringstream ss;
  ss << nonbond.description() << " & " << nonbond_dd_cutoff.description();
  return ss.str();
}


size_t EnergyNonbond_O::runTestCalls(core::T_sp stream, chem::NVector_sp coords) const
{
  if (this->_UsesExcludedAtoms) {
    SIMPLE_ERROR("You cannot runTestCalls with _UseExcludedAtoms = true");
  }
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
  Nonbond_Dd_Cutoff<double*> nonbond_dd_cutoff;
  for ( auto si=this->_Terms.begin();
        si!=this->_Terms.end(); si++ ) {
    position[0] = coords[si->term.I1];
    position[1] = coords[si->term.I1+1];
    position[2] = coords[si->term.I1+2];
    position[3] = coords[si->term.I2];
    position[4] = coords[si->term.I2+1];
    position[5] = coords[si->term.I2+2];
    double invdd = this->_Nonbond_invdd;
    energy_new = 0.0;
    energy_ground = 0.0;
#define POW2(xx) (xx*xx)
    //#define NONBOND_DD_ARGS si->term.dA, si->term.dC, si->term.dQ1Q2, this->_Nonbond_invdd, this->_Nonbond_r_switch, POW2(this->_Nonbond_r_switch), this->_Nonbond_r_cut, POW2(this->_Nonbond_r_cut), this->_Nonbond_r_cut-this->_Nonbond_r_switch, 0, 3
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    nonbond_dd_cutoff.gradient( si->term.dA, si->term.dC, si->term.dQ1Q2, this->_Nonbond_invdd, this->_Nonbond_r_switch, POW2(this->_Nonbond_r_switch), this->_Nonbond_r_cut, POW2(this->_Nonbond_r_cut), this->_Nonbond_r_cut-this->_Nonbond_r_switch, 0, 3, position, &energy_new, force_new, hessian_new, dvec_new, hdvec_new );
    nonbond_dd_cutoff.gradient_fd( si->term.dA, si->term.dC, si->term.dQ1Q2, this->_Nonbond_invdd, this->_Nonbond_r_switch, POW2(this->_Nonbond_r_switch), this->_Nonbond_r_cut, POW2(this->_Nonbond_r_cut), this->_Nonbond_r_cut-this->_Nonbond_r_switch, 0, 3, position, &energy_ground, force_ground, hessian_ground, dvec_ground, hdvec_ground );
    if (!test_match( stream, "nonbond_dd_cutoff_force", POS_SIZE,
                     force_new, force_ground,
                     0, 0,
                     0, 0 )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH nonbond_dd_cutoff_gradient #{} dA = {}  dC = {}\n",
                              idx, si->term.dA, si->term.dC ), stream );
    }
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    nonbond_dd_cutoff.hessian( si->term.dA, si->term.dC, si->term.dQ1Q2, this->_Nonbond_invdd, this->_Nonbond_r_switch, POW2(this->_Nonbond_r_switch), this->_Nonbond_r_cut, POW2(this->_Nonbond_r_cut), this->_Nonbond_r_cut-this->_Nonbond_r_switch, 0, 3, position, &energy_new, force_new, hessian_new, dvec_new, hdvec_new );
    nonbond_dd_cutoff.hessian_fd( si->term.dA, si->term.dC, si->term.dQ1Q2, this->_Nonbond_invdd, this->_Nonbond_r_switch, POW2(this->_Nonbond_r_switch), this->_Nonbond_r_cut, POW2(this->_Nonbond_r_cut), this->_Nonbond_r_cut-this->_Nonbond_r_switch, 0, 3, position, &energy_ground, force_ground, hessian_ground, dvec_ground, hdvec_ground );
    if (!test_match( stream, "nonbond_dd_cutoff_hessian", POS_SIZE,
                     force_new, force_ground,
                     hessian_new, hessian_ground,
                     hdvec_new, hdvec_ground )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH nonbond_dd_cutoff_hessian #{} dA = {}  dC = {}\n",
                              idx, si->term.dA, si->term.dC ), stream );
    }
    idx++;
  }
  core::print(fmt::format("nonbond_dd_cutoff errors = {}\n", errs), stream);
  return errs;
}



}


#if 0
int	_areValuesCloseNB( double numVal, double analVal, const char* funcName, const char* termName, int index )
{
  double	rel = 0.0;
  if ( fabs(numVal)<FINITE_DIFFERENCE_TOLERANCE && fabs(analVal)<FINITE_DIFFERENCE_TOLERANCE ) {
    return true; // both close to zero
  }
  core::lisp_write(fmt::format("NB close? {}[{}] analVal->{}  numVal->{}\n",termName,index,analVal,numVal));
  if ( numVal*analVal > 0.0 ) {
    double avg = (fabs(numVal)+fabs(analVal));
    rel = (fabs(numVal-analVal)/avg);
    if ( rel<FINITE_DIFFERENCE_TOLERANCE ) {
      return true;	// they are within a few percent of each other
    }
  }
  return false;
}

#undef TEST_FORCE
#define TEST_FORCE(func, delta, argLow, argHigh, term, idx)                                                                        \
  {                                                                                                                                \
    num_real eLow = func argLow;                                                                                                   \
    num_real eHigh = func argHigh;                                                                                                 \
    num_real numForce = -(eHigh - eLow) / (delta);                                                                                 \
    LOG(("eHigh = %le eLow = %le delta = %le"), (eHigh), (eLow), (delta));                                                         \
    LOG(("numForce = %le"), numForce);                                                                                             \
    if (!_areValuesCloseNB(numForce, term, #func, #term, idx)) {                                                                   \
      fails++;                                                                                                                     \
    };                                                                                                                             \
  }
#endif

// #define DEBUG_NONBOND_TERM 1
#define LOG_ENERGY(...)
// #define LOG_ENERGY(...) core::clasp_write_string(fmt::format(__VA_ARGS__))

namespace chem {

CL_DEFMETHOD
void EnergyNonbond_O::set_nonbond_pairlist_parameters( double r_switch, double r_cut, double r_pairlist, double distance_dielectric ) {
  this->_Nonbond_r_switch = r_switch;
  this->_Nonbond_r_cut = r_cut;
  this->_Nonbond_r_pairlist = r_pairlist;
  if (distance_dielectric<VERYSMALL) {
    SIMPLE_ERROR("distance-dielectric must be a reasonable value - passed {}", distance_dielectric );
  }
  this->_Nonbond_invdd = 1.0/distance_dielectric;
}

#define PBX(_delta_) (_delta_)
#define PBY(_delta_) (_delta_)
#define PBZ(_delta_) (_delta_)
// FIXME: Disabled BAIL_OUT_IF_CUTOFF
#define BAIL_OUT_IF_CUTOFF(deltaSquared)
// #define BAIL_OUT_IF_CUTOFF(deltaSquared) if (deltaSquared>CUTOFF_SQUARED) goto SKIP_term;

EnergyFunction_sp energyFunctionNonbondParameters(ScoringFunction_sp score,
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

core::T_sp debug_nonbond(double Energy, double x1, double y1, double z1, double x2, double y2, double z2, 
                         double dQ1Q2, double dA, double dC/*, double Eeel, double Evdw, double fx1, double fy1,
                                                             double fz1, double fx2, double fy2, double fz2*/) {
  ql::list ll;
  ll << INTERN_(kw, Energy) << mk_double_float(Energy) << INTERN_(kw, x1) << mk_double_float(x1) << INTERN_(kw, y1)
     << mk_double_float(y1) << INTERN_(kw, z1) << mk_double_float(z1) << INTERN_(kw, x2) << mk_double_float(x2) << INTERN_(kw, y2)
     << mk_double_float(y2) << INTERN_(kw, z2) << mk_double_float(z2)
     << INTERN_(kw, dq1q2) << mk_double_float(dQ1Q2) << INTERN_(kw, da)
     << mk_double_float(dA) << INTERN_(kw, dc) << mk_double_float(dC)
#if 0
     << INTERN_(kw, eel) << mk_double_float(Eeel)
     << INTERN_(kw, evdw) << mk_double_float(Evdw) << INTERN_(kw, fx1) << mk_double_float(fx1) << INTERN_(kw, fy1)
     << mk_double_float(fy1) << INTERN_(kw, fz1) << mk_double_float(fz1) << INTERN_(kw, fx2) << mk_double_float(fx2)
     << INTERN_(kw, fy2) << mk_double_float(fy2) << INTERN_(kw, fz2) << mk_double_float(fz2)
#endif
      ;
  return ll.cons();
}

SYMBOL_EXPORT_SC_(ChemPkg, EnergyNonbond);
SYMBOL_EXPORT_SC_(ChemPkg, EnergyNonbond14);
SYMBOL_EXPORT_SC_(ChemPkg, trapEnergyNonbondDefineForAtomPair );
SYMBOL_EXPORT_SC_(ChemPkg, trapEnergyNonbondNotDefineForAtomPair );

core::T_sp nonbond_type(bool is14) {
  if (is14)
    return chem::_sym_EnergyNonbond14;
  return chem::_sym_EnergyNonbond;
}
#define NONBOND_APPLY_ATOM_MASK(I1, I2)                                                                                            \
  if (hasActiveAtomMask && !(bitvectorActiveAtomMask->testBit(I1 / 3) && bitvectorActiveAtomMask->testBit(I2 / 3)))                \
    goto SKIP_term;

/* No periodic boundary conditons/bounding-box is used in this code.
   So make these macros do nothing.
   See energyPeriodicBoundaryConditionsNonbond.cc for counter-example.
*/
core::List_sp EnergyNonbond::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw, da), core::clasp_make_double_float(this->term.dA)),
                                  core::Cons_O::create(INTERN_(kw, dc), core::clasp_make_double_float(this->term.dC)),
                                  core::Cons_O::create(INTERN_(kw, i1), core::make_fixnum(this->term.I1)),
                                  core::Cons_O::create(INTERN_(kw, i2), core::make_fixnum(this->term.I2)),
                                  core::Cons_O::create(INTERN_(kw, charge1), core::clasp_make_double_float(this->_Atom1_enb->getCharge())),
                                  core::Cons_O::create(INTERN_(kw, charge2), core::clasp_make_double_float(this->_Atom2_enb->getCharge())),
                                  core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1_enb),
                                  core::Cons_O::create(INTERN_(kw,atom2), this->_Atom2_enb)
                                  );
}

void EnergyNonbond::decode(core::List_sp alist) { SIMPLE_ERROR("Implement decode of EnergyNonbond"); }

#define ENERGY_FUNCTION score, energyScale, I1, I2, activeAtomMask

#define NONBOND_DEBUG_INTERACTIONS(I1, I2)                                                                                         \
  if (doDebugInteractions) {                                                                                                       \
    core::eval::funcall(debugInteractions, nonbond_type(interactionsAre14),                                                        \
                        debug_nonbond(Energy, \
                                      position[I1],position[I1+1],position[I1+2], \
                                      position[I2],position[I2+1],position[I2+2], \
                                      dQ1Q2, dA, dC /*, Eeel, Evdw, fx1, fy1, fz1, fx2, fy2, fz2 */),       \
                        core::make_fixnum(I1), core::make_fixnum(I2));                                                             \
  }

double _evaluateEnergyOnly_Nonbond(ScoringFunction_sp score, core::T_sp energyScale, int I1, int I2, core::T_sp activeAtomMask, num_real x1, num_real y1,
                                   num_real z1, num_real x2, num_real y2, num_real z2, num_real dA, num_real dC, num_real dQ1Q2);

void calculateFiniteDifferenceForces(ScoringFunction_sp score, core::T_sp energyScale, core::T_sp activeAtomMask, int I1, int I2, double x1, double y1,
                                     double z1, double x2, double y2, double z2, double dA, double dC, double dQ1Q2, double &fdfx1,
                                     double &fdfy1, double &fdfz1, double &fdfx2, double &fdfy2, double &fdfz2) {
  double delta = 0.00001;
  double deltaTimes2 = delta * 2.0;
  double ehigh;
  double elow;
  ehigh = _evaluateEnergyOnly_Nonbond(score, energyScale, I1, I2, activeAtomMask, x1 + delta, y1, z1, x2, y2, z2, dA, dC, dQ1Q2);
  elow = _evaluateEnergyOnly_Nonbond(score, energyScale, I1, I2, activeAtomMask, x1 - delta, y1, z1, x2, y2, z2, dA, dC, dQ1Q2);
  fdfx1 = -(ehigh - elow) / deltaTimes2;
  ehigh = _evaluateEnergyOnly_Nonbond(score, energyScale, I1, I2, activeAtomMask, x1, y1 + delta, z1, x2, y2, z2, dA, dC, dQ1Q2);
  elow = _evaluateEnergyOnly_Nonbond(score, energyScale, I1, I2, activeAtomMask, x1, y1 - delta, z1, x2, y2, z2, dA, dC, dQ1Q2);
  fdfy1 = -(ehigh - elow) / deltaTimes2;
  ehigh = _evaluateEnergyOnly_Nonbond(score, energyScale, I1, I2, activeAtomMask, x1, y1, z1 + delta, x2, y2, z2, dA, dC, dQ1Q2);
  elow = _evaluateEnergyOnly_Nonbond(score, energyScale, I1, I2, activeAtomMask, x1, y1, z1 - delta, x2, y2, z2, dA, dC, dQ1Q2);
  fdfz1 = -(ehigh - elow) / deltaTimes2;
  ehigh = _evaluateEnergyOnly_Nonbond(score, energyScale, I1, I2, activeAtomMask, x1, y1, z1, x2 + delta, y2, z2, dA, dC, dQ1Q2);
  elow = _evaluateEnergyOnly_Nonbond(score, energyScale, I1, I2, activeAtomMask, x1, y1, z1, x2 - delta, y2, z2, dA, dC, dQ1Q2);
  fdfx2 = -(ehigh - elow) / deltaTimes2;
  ehigh = _evaluateEnergyOnly_Nonbond(score, energyScale, I1, I2, activeAtomMask, x1, y1, z1, x2, y2 + delta, z2, dA, dC, dQ1Q2);
  elow = _evaluateEnergyOnly_Nonbond(score, energyScale, I1, I2, activeAtomMask, x1, y1, z1, x2, y2 - delta, z2, dA, dC, dQ1Q2);
  fdfy2 = -(ehigh - elow) / deltaTimes2;
  ehigh = _evaluateEnergyOnly_Nonbond(score, energyScale, I1, I2, activeAtomMask, x1, y1, z1, x2, y2, z2 + delta, dA, dC, dQ1Q2);
  elow = _evaluateEnergyOnly_Nonbond(score, energyScale, I1, I2, activeAtomMask, x1, y1, z1, x2, y2, z2 - delta, dA, dC, dQ1Q2);
  fdfz2 = -(ehigh - elow) / deltaTimes2;
}

struct NoFiniteDifference {
  static void maybeTestFiniteDifference(ScoringFunction_sp score, core::T_sp energyScale, int I1, int I2, core::T_sp activeAtomMask, double x1, double y1,
                                        double z1, double x2, double y2, double z2, double dA, double dC, double dQ1Q2, double fx1,
                                        double fy1, double fz1, double fx2, double fy2, double fz2, int index, size_t &fails,
                                        bool debugForce) {}
};

struct DebugFiniteDifference {
  static void maybeTestFiniteDifference(ScoringFunction_sp score, core::T_sp energyScale, int I1, int I2, core::T_sp activeAtomMask, double x1, double y1,
                                        double z1, double x2, double y2, double z2, double dA, double dC, double dQ1Q2, double fx1,
                                        double fy1, double fz1, double fx2, double fy2, double fz2, int index, size_t &fails,
                                        bool debugForce) {
#undef TEST_DIAGONAL_HESSIAN
#define TEST_DIAGONAL_HESSIAN(x, y, o1, o2, o3, d, i)
#undef TEST_OFF_DIAGONAL_HESSIAN
#define TEST_OFF_DIAGONAL_HESSIAN(x, y, o1, o2, o3, o4, d, i)
#define ENERGY_FUNCTION score, energyScale, I1, I2, activeAtomMask
#include <cando/chem/energy_functions/_Nonbond_debugFiniteDifference.cc>
#undef ENERGY_FUNCTION
    if (debugForce) {
      double fdfx1;
      double fdfy1;
      double fdfz1;
      double fdfx2;
      double fdfy2;
      double fdfz2;
      calculateFiniteDifferenceForces(score, energyScale, activeAtomMask, I1, I2, x1, y1, z1, x2, y2, z2, dA, dC, dQ1Q2, fdfx1, fdfy1, fdfz1,
                                      fdfx2, fdfy2, fdfz2);
      core::lisp_write(fmt::format("({0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} "
                                   "{20} {21} {22}) ; {23}:{24}:{25} \n",
                                   I1, I2,                                   // 0-1
                                   x1, y1, z1, x2, y2, z2,                   // 2-7
                                   dA, dC, dQ1Q2,                            // 8-10
                                   fx1, fy1, fz1, fx2, fy2, fz2,             // 11
                                   fdfx1, fdfy1, fdfz1, fdfx2, fdfy2, fdfz2, // 17
                                   __FILE__, __LINE__, __FUNCTION__          // 23
                                   ));
    }
  }
};

inline double calculate_dQ1Q2(double electrostaticScale, double electrostaticModifier, double charge1, double charge2) {
  return electrostaticScale * electrostaticModifier * charge1 * charge2;
}

/*! The main nonbond code using excluded atoms.
    It is a template function so that template arguments can inline or elide testing code.
 */
template <class MaybeFiniteDiff>
//__attribute__((optnone))
double template_evaluateUsingExcludedAtoms(EnergyNonbond_O *mthis, ScoringFunction_sp score, NVector_sp pos,
                                           core::T_sp energyScale,
                                           core::T_sp energyComponents, bool calcForce, gc::Nilable<NVector_sp> force,
                                           bool calcDiagonalHessian, bool calcOffDiagonalHessian,
                                           gc::Nilable<AbstractLargeSquareMatrix_sp> hessian, gc::Nilable<NVector_sp> hdvec,
                                           gc::Nilable<NVector_sp> dvec, core::T_sp activeAtomMask, core::T_sp debugInteractions,
                                           size_t &fails, size_t &index, bool debugForce = false) {
  double dielectricConstant;
  double dQ1Q2Scale;
  double cutoff;
  EnergyFunction_sp energyFunction = energyFunctionNonbondParameters(score, energyScale, dielectricConstant, dQ1Q2Scale, cutoff);
#define CUTOFF_SQUARED (cutoff * cutoff)
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  //  printf("%s:%d In evaluateUsingExcludedAtoms starting this->_DebugEnergy -> %d\n", __FILE__, __LINE__, this->_DebugEnergy );
  if (!mthis->_iac_vec) {
    SIMPLE_ERROR("The nonbonded excluded atoms parameters have not been set up");
  }
  core::SimpleVector_int32_t_sp numberOfExcludedAtoms = mthis->_NumberOfExcludedAtomIndexes;
  core::SimpleVector_int32_t_sp excludedAtomIndexes = mthis->_ExcludedAtomIndexes;

  //  printf("%s:%d electrostaticcharge %lf\n", __FILE__, __LINE__, electrostaticScale );

  bool hasForce = force.notnilp();
  bool hasHessian = hessian.notnilp();
  bool hasHdAndD = (hdvec.notnilp()) && (dvec.notnilp());

  double totalEnergy = 0.0;
  DOUBLE* position = &(*pos)[0];
  DOUBLE* rforce = NULL;
  DOUBLE* rhessian = NULL; // &(*hessian)[0];
  DOUBLE* rdvec = NULL;
  DOUBLE* rhdvec = NULL;

  DOUBLE  Energy = 0.0;

  num_real dA, dC, dQ1Q2;
  double vdwScale = energyScaleVdwScale(energyScale);
  double eelScale = energyScaleElectrostaticScale(energyScale);
  double DIELECTRIC = energyScaleDielectricConstant(energyScale);
  int i3x1, i3x2;
  int i = 0;
  int endIndex = pos->length() / 3;
  //  int endIndex = mthis->_AtomTable->getNumberOfAtoms();
  int excludedAtomIndex = 0;
  int nlocaltype = 0;
  // Find the max local type
  for (i = 0; i < mthis->_iac_vec->length(); ++i) {
    if (nlocaltype < (*mthis->_iac_vec)[i]) {
      nlocaltype = (*mthis->_iac_vec)[i];
    }
  }

#define KERNEL_EXCLUDED_NONBOND_APPLY_ATOM_MASK(I1, I2)                                                              \
  if (hasActiveAtomMask && !(bitvectorActiveAtomMask->testBit(I1 / 3) && bitvectorActiveAtomMask->testBit(I2 / 3)))  continue;

  bool interactionsAre14 = false;
  Nonbond<NoHessian> nonbond;
  int index1_end = endIndex - 1;
  if (!hasForce) {
    for (int index1 = 0; index1 < index1_end; ++index1) {
      LOG("{} ====== top of outer loop - index1 = {}\n", __FUNCTION__, index1);
      // Skip 0 in excluded atom list that amber requires
      bool has_excluded_atoms = ((*excludedAtomIndexes)[excludedAtomIndex] >= 0);
      int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
      num_real charge11 = (*mthis->_charge_vector)[index1];
      for (int index2 = index1 + 1, index2_end(endIndex); index2 < index2_end; ++index2) {
        int maybe_excluded_atom = (*excludedAtomIndexes)[excludedAtomIndex];
        // state TOP-INNER
        if (numberOfExcludedAtomsRemaining > 0 && maybe_excluded_atom == index2) {
          ++excludedAtomIndex;
          --numberOfExcludedAtomsRemaining;
          continue;
        }
        int localindex1 = (*mthis->_iac_vec)[index1];
        int localindex2 = (*mthis->_iac_vec)[index2];
        dA = (*mthis->_cn1_vec)[(*mthis->_ico_vec)[nlocaltype * (localindex1 - 1) + localindex2 - 1] - 1];
        dC = (*mthis->_cn2_vec)[(*mthis->_ico_vec)[nlocaltype * (localindex1 - 1) + localindex2 - 1] - 1];
        num_real charge22 = (*mthis->_charge_vector)[index2];
        dQ1Q2 = calculate_dQ1Q2(1.0, dQ1Q2Scale, charge11, charge22);
        i3x1 = index1 * 3;
        i3x2 = index2 * 3;
        KERNEL_EXCLUDED_NONBOND_APPLY_ATOM_MASK(i3x1,i3x2);
        Energy = nonbond.energy( dA, dC, dQ1Q2,
                                 i3x1, i3x2,
                                 position,
                                 &totalEnergy,
                                 NULL, // force
                                 NoHessian(), // hessian
                                 NULL, // dvec
                                 NULL ); // hdvec
        NONBOND_DEBUG_INTERACTIONS(i3x1,i3x2);
        index++;
      }
      if (!has_excluded_atoms) {
        ++excludedAtomIndex;
      }
    }
  } else if (hasForce) {
    rforce = &(*force)[0];
    for (int index1 = 0; index1 < index1_end; ++index1) {
      LOG("{} ====== top of outer loop - index1 = {}\n", __FUNCTION__, index1);
      // Skip 0 in excluded atom list that amber requires
      bool has_excluded_atoms = ((*excludedAtomIndexes)[excludedAtomIndex] >= 0);
      int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
      num_real charge11 = (*mthis->_charge_vector)[index1];
      for (int index2 = index1 + 1, index2_end(endIndex); index2 < index2_end; ++index2) {
        int maybe_excluded_atom = (*excludedAtomIndexes)[excludedAtomIndex];
        // state TOP-INNER
        if (numberOfExcludedAtomsRemaining > 0 && maybe_excluded_atom == index2) {
          ++excludedAtomIndex;
          --numberOfExcludedAtomsRemaining;
          continue;
        }
        int localindex1 = (*mthis->_iac_vec)[index1];
        int localindex2 = (*mthis->_iac_vec)[index2];
        dA = (*mthis->_cn1_vec)[(*mthis->_ico_vec)[nlocaltype * (localindex1 - 1) + localindex2 - 1] - 1];
        dC = (*mthis->_cn2_vec)[(*mthis->_ico_vec)[nlocaltype * (localindex1 - 1) + localindex2 - 1] - 1];
        num_real charge22 = (*mthis->_charge_vector)[index2];
        dQ1Q2 = calculate_dQ1Q2(1.0, dQ1Q2Scale, charge11, charge22);
        i3x1 = index1 * 3;
        i3x2 = index2 * 3;

        KERNEL_EXCLUDED_NONBOND_APPLY_ATOM_MASK(i3x1,i3x2);
        Energy = nonbond.gradient( dA, dC, dQ1Q2,
                        i3x1, i3x2,
                        position,
                        &totalEnergy,
                        rforce, // force
                        NoHessian(), // hessia
                        NULL, // dvec
                        NULL ); // hdvec
        NONBOND_DEBUG_INTERACTIONS(i3x1,i3x2);
        index++;
      }
      if (!has_excluded_atoms) {
        ++excludedAtomIndex;
      }
    }
  } else { // hasDvecAndHdvec
    rforce = &(*force)[0];
    rdvec = &(*dvec)[0];
    rhdvec = &(*hdvec)[0];
    for (int index1 = 0; index1 < index1_end; ++index1) {
      LOG("{} ====== top of outer loop - index1 = {}\n", __FUNCTION__, index1);
      // Skip 0 in excluded atom list that amber requires
      bool has_excluded_atoms = ((*excludedAtomIndexes)[excludedAtomIndex] >= 0);
      int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
      num_real charge11 = (*mthis->_charge_vector)[index1];
      for (int index2 = index1 + 1, index2_end(endIndex); index2 < index2_end; ++index2) {
        int maybe_excluded_atom = (*excludedAtomIndexes)[excludedAtomIndex];
        // state TOP-INNER
        if (numberOfExcludedAtomsRemaining > 0 && maybe_excluded_atom == index2) {
          ++excludedAtomIndex;
          --numberOfExcludedAtomsRemaining;
          continue;
        }
        int localindex1 = (*mthis->_iac_vec)[index1];
        int localindex2 = (*mthis->_iac_vec)[index2];
        dA = (*mthis->_cn1_vec)[(*mthis->_ico_vec)[nlocaltype * (localindex1 - 1) + localindex2 - 1] - 1];
        dC = (*mthis->_cn2_vec)[(*mthis->_ico_vec)[nlocaltype * (localindex1 - 1) + localindex2 - 1] - 1];
        num_real charge22 = (*mthis->_charge_vector)[index2];
        dQ1Q2 = calculate_dQ1Q2(1.0, dQ1Q2Scale, charge11, charge22);
        i3x1 = index1 * 3;
        i3x2 = index2 * 3;

        KERNEL_EXCLUDED_NONBOND_APPLY_ATOM_MASK(i3x1,i3x2);
        Energy = nonbond.hessian( dA, dC, dQ1Q2,
                        i3x1, i3x2,
                        position,
                        &totalEnergy,
                        rforce, // force
                        NoHessian(), // hessia
                        rdvec, // dvec
                        rhdvec ); // hdvec
        NONBOND_DEBUG_INTERACTIONS(i3x1,i3x2);
        index++;
      }
      if (!has_excluded_atoms) {
        ++excludedAtomIndex;
      }
    }
  }
  maybeSetEnergy(energyComponents, _sym_energyNonbondTotal, totalEnergy );
  return totalEnergy;
}

/*! The core nonbond code using pairwise terms.
    It as a template function so that template arguments can inline or elide testing code.
 */
template <class MaybeFiniteDiff>
double template_evaluateUsingTerms(EnergyNonbond_O *mthis,
                                   bool interactionsAre14,
                                   const gctools::Vec0<EnergyNonbond>& terms,
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
  mthis->maybeRebuildPairList(nvposition);
  EnergyFunction_sp energyFunction = energyFunctionNonbondParameters(score, energyScale, dielectricConstant, dQ1Q2Scale, cutoff);
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
  bool hasForce = force.notnilp();
  bool hasHessian = hessian.notnilp();
  bool hasHdAndD = (hdvec.notnilp()) && (dvec.notnilp());

  double totalNonbondEnergy = 0.0;
  DOUBLE* position = &(*nvposition)[0];
  DOUBLE* rforce = NULL;
  DOUBLE* rhessian = NULL; // &(*hessian)[0];
  DOUBLE* rdvec = NULL;
  DOUBLE* rhdvec = NULL;
  DOUBLE  Energy = 0.0;
  Nonbond_Dd_Cutoff<NoHessian> nonbond_dd_cutoff;

  // If you are going to use openmp here, you need to control access to the force and hessian
  // arrays so that only one thread updates each element at a time.
  {
    double r_switch2 = mthis->_Nonbond_r_switch*mthis->_Nonbond_r_switch;
    double r_cut2 = mthis->_Nonbond_r_cut*mthis->_Nonbond_r_cut;
    double inv_range = 1.0/(mthis->_Nonbond_r_cut-mthis->_Nonbond_r_switch); // 1.0/(r_cut - r_switch)

#define KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(I1, I2)                                                              \
  if (hasActiveAtomMask && !(bitvectorActiveAtomMask->testBit(I1 / 3) && bitvectorActiveAtomMask->testBit(I2 / 3)))  continue;

    if (!hasForce) {
      for (auto si = terms.begin(); si != terms.end(); si++ ) {
        KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(si->term.I1,si->term.I2);
        Energy = nonbond_dd_cutoff.energy(
            si->term.dA,
            si->term.dC,
            si->term.dQ1Q2,
            mthis->_Nonbond_invdd,         // dd for distance dependent dielectric
            mthis->_Nonbond_r_switch,   // r_switch
            r_switch2,                 // r_switch^2
            mthis->_Nonbond_r_cut,      // r_cut
            r_cut2,                    // r_cut^
            inv_range,                 // inv_range
            si->term.I1,
            si->term.I2,
            position,
            &totalNonbondEnergy,
            NULL,
            NoHessian(),
            NULL,
            NULL);
        const double& dA = si->term.dA;
        const double& dC = si->term.dC;
        const double& dQ1Q2 = si->term.dQ1Q2;
        NONBOND_DEBUG_INTERACTIONS(si->term.I1,si->term.I2);
      }
    } else if (hasForce) {
      rforce = &(*force)[0];
      for (auto si = terms.begin(); si != terms.end(); si++ ) {
        KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(si->term.I1,si->term.I2);
        Energy = nonbond_dd_cutoff.gradient(
            si->term.dA,
            si->term.dC,
            si->term.dQ1Q2,
            mthis->_Nonbond_invdd,         // dd for distance dependent dielectric
            mthis->_Nonbond_r_switch,   // r_switch
            r_switch2,                 // r_switch^2
            mthis->_Nonbond_r_cut,      // r_cut
            r_cut2,                    // r_cut^
            inv_range,                 // inv_range
            si->term.I1,
            si->term.I2,
            position,
            &totalNonbondEnergy,
            rforce,
            NoHessian(),
            NULL,
            NULL );
        const double& dA = si->term.dA;
        const double& dC = si->term.dC;
        const double& dQ1Q2 = si->term.dQ1Q2;
        NONBOND_DEBUG_INTERACTIONS(si->term.I1,si->term.I2);
      }
    } else { // if (hasHdAndD) {
      rforce = &(*force)[0];
      rdvec = &(*dvec)[0];
      rhdvec = &(*hdvec)[0];
      for (auto si = terms.begin(); si != terms.end(); si++ ) {
        KERNEL_TERM_NONBOND_APPLY_ATOM_MASK(si->term.I1,si->term.I2);
        Energy = nonbond_dd_cutoff.hessian(
            si->term.dA,
            si->term.dC,
            si->term.dQ1Q2,
            mthis->_Nonbond_invdd,         // dd for distance dependent dielectric
            mthis->_Nonbond_r_switch,   // r_switch
            r_switch2,                 // r_switch^2
            mthis->_Nonbond_r_cut,      // r_cut
            r_cut2,                    // r_cut^
            inv_range,                 // inv_range
            si->term.I1,
            si->term.I2,
            position,
            &totalNonbondEnergy,
            rforce,
            NoHessian(),
            rdvec,
            rhdvec);
        const double& dA = si->term.dA;
        const double& dC = si->term.dC;
        const double& dQ1Q2 = si->term.dQ1Q2;
        NONBOND_DEBUG_INTERACTIONS(si->term.I1,si->term.I2);
      }
    }
  }
  maybeSetEnergy(energyComponents, termSymbol, totalNonbondEnergy );
  return totalNonbondEnergy;
}

SYMBOL_EXPORT_SC_(ChemPkg, find_type);

double EnergyNonbond_O::evaluateAllComponent(ScoringFunction_sp score,
                                             NVector_sp pos,
                                             core::T_sp energyScale,
                                             core::T_sp energyComponents,
                                             bool calcForce,
                                             gc::Nilable<NVector_sp> force, bool calcDiagonalHessian, bool calcOffDiagonalHessian,
                                             gc::Nilable<AbstractLargeSquareMatrix_sp> hessian, gc::Nilable<NVector_sp> hdvec,
                                             gc::Nilable<NVector_sp> dvec, core::T_sp activeAtomMask,
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
  if (this->_UsesExcludedAtoms) {
    // Evaluate the nonbonds using the excluded atom list
    energy += template_evaluateUsingExcludedAtoms<NoFiniteDifference>(this, score, pos, energyScale, energyComponents, calcForce, force,
                                                                     calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec,
                                                                     dvec, activeAtomMask, debugInteractions, fails, index);
  } else {
    // Evaluate everything using terms
    energy += template_evaluateUsingTerms<NoFiniteDifference>(this, false, this->_Terms, _sym_energyNonbondTotal,
                                                              score, pos, energyScale, energyComponents, calcForce, force,
                                                             calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec,
                                                             activeAtomMask, debugInteractions, fails, index);
  }
  return energy;
}

CL_DEFMETHOD
double EnergyNonbond_O::debugAllComponent(ScoringFunction_sp score,
                                          NVector_sp pos,
                                          core::T_sp energyScale,
                                          core::T_sp energyComponents,
                                          bool calcForce,
                                          gc::Nilable<NVector_sp> force, bool calcDiagonalHessian, bool calcOffDiagonalHessian,
                                          gc::Nilable<AbstractLargeSquareMatrix_sp> hessian, gc::Nilable<NVector_sp> hdvec,
                                          gc::Nilable<NVector_sp> dvec, core::T_sp activeAtomMask, core::T_sp debugInteractions) {
  //  printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  double energy = 0.0;
  size_t fails = 0;
  size_t index = 0;
    // Evaluate the 1-4 terms
  energy += template_evaluateUsingTerms<DebugFiniteDifference>(this, true, this->_Terms14, _sym_energyNonbond14,
                                                               score, pos, energyScale, energyComponents, calcForce, force,
                                                               calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec,
                                                               activeAtomMask, debugInteractions, fails, index, true);
  if (this->_UsesExcludedAtoms) {
    // Evaluate the nonbonds using the excluded atom list
    energy += template_evaluateUsingExcludedAtoms<DebugFiniteDifference>(
        this, score, pos, energyScale, energyComponents, calcForce, force, calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec,
        activeAtomMask, debugInteractions, fails, index, true);
  } else {
    // Evaluate everything using terms
    energy += template_evaluateUsingTerms<DebugFiniteDifference>(this, false, this->_Terms, _sym_energyNonbondTotal,
                                                                 score, pos, energyScale, energyComponents, calcForce, force,
                                                                calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec,
                                                                activeAtomMask, debugInteractions, fails, index, true);
  }
  return energy;
}

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
#define NONBOND_DEBUG_INTERACTIONS(i1, i2)
#include <cando/chem/energy_functions/_Nonbond_termCode.cc>
#undef NONBOND_DEBUG_INTERACTIONS
#undef CUTOFF_SQUARED
#undef DIELECTRIC
  return Energy;
}

CL_LAMBDA((energy-nonbond chem:energy-nonbond) score pos &optional active-atom-mask);
CL_DEFMETHOD void EnergyNonbond_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(ScoringFunction_sp score,
                                                                                          NVector_sp pos,
                                                                                          core::T_sp energyScale,
                                                                                          core::T_sp activeAtomMask) {
  double energy = 0.0;
  size_t fails = 0;
  size_t index = 0;
  auto force = NVector_O::make(pos->size(), 0.0, true);
  // Evaluate the 1-4 terms
  energy += template_evaluateUsingTerms<DebugFiniteDifference>(this, true, this->_Terms14, _sym_energyNonbond14,
                                                               score, pos, energyScale, nil<core::T_O>(), true, force, false, false,
                                                               nil<core::T_O>(), nil<core::T_O>(), nil<core::T_O>(),
                                                               nil<core::T_O>(), nil<core::T_O>(), fails, index);
  if (this->_UsesExcludedAtoms) {
    // Evaluate the nonbonds using the excluded atom list
    energy += template_evaluateUsingExcludedAtoms<DebugFiniteDifference>(this, score, pos, energyScale, nil<core::T_O>(), true, force, false,
                                                                         false, nil<core::T_O>(), nil<core::T_O>(), nil<core::T_O>(),
                                                                         nil<core::T_O>(), nil<core::T_O>(), fails, index);
  } else {
    // Evaluate everything using terms
    energy += template_evaluateUsingTerms<DebugFiniteDifference>(this, false, this->_Terms, _sym_energyNonbondTotal,
                                                                 score, pos, energyScale, nil<core::T_O>(), true, force, false, false,
                                                                 nil<core::T_O>(), nil<core::T_O>(), nil<core::T_O>(),
                                                                 nil<core::T_O>(), nil<core::T_O>(), fails, index);
  }
  core::lisp_write(fmt::format("There are {} fails out of {} terms\n", fails, index));
}

bool EnergyNonbond::defineForAtomPair(core::T_sp forceField, bool is14, Atom_sp a1, Atom_sp a2, size_t a1CoordinateIndexTimes3,
                                      size_t a2CoordinateIndexTimes3, EnergyNonbond_sp energyNonbond,
                                      core::HashTable_sp atomTypes, core::T_sp keepInteraction ) {
  if (_sym_trapEnergyNonbondDefineForAtomPair->fboundp()) {
    core::eval::funcall(_sym_trapEnergyNonbondDefineForAtomPair,
                        a1, a2,
                        core::make_fixnum(a1CoordinateIndexTimes3),
                        core::make_fixnum(a2CoordinateIndexTimes3),
                        energyNonbond,
                        atomTypes,
                        keepInteraction);
  }
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
  if (is14) {
    // See AMBER manual 1-4 Non-Bonded Interaction Scaling
    // Section 15.1.6
    // vdwScale is 1/2 of vdw
    // electrostaticScale is 1.0/1.2 of electrostatic
    vdwScale = 1.0 / 2.0;
    electrostaticScale = 1.0 / 1.2;
  } else {
    vdwScale = 1.0;
    electrostaticScale = 1.0;
  }
  LOG("vdwScale = {}", (double)(vdwScale));
  LOG("electrostaticScale = {}", (double)(electrostaticScale));
  LOG(" is14={}", is14);
  FFNonbond_sp ffNonbond1 = gc::As<FFNonbond_sp>(tffNonbond1);
  FFNonbond_sp ffNonbond2 = gc::As<FFNonbond_sp>(tffNonbond2);
  {
    double rstar = ffNonbond1->getRadius_Angstroms() + ffNonbond2->getRadius_Angstroms();
    epsilonij = sqrt(ffNonbond1->getEpsilon_kcal() * ffNonbond2->getEpsilon_kcal());
    double rtothe6th = rstar*rstar*rstar*rstar*rstar*rstar;
    double rtothe12th = rtothe6th*rtothe6th;
    double parmA = epsilonij * rtothe12th; // pow(this->_RStar, 12.0);
    double parmC = 2.0 * epsilonij * rtothe6th; // pow(this->_RStar, 6.0);
    this->term.dA = parmA * vdwScale;
    this->term.dC = parmC * vdwScale;
  }
  {
    double charge1 = a1->getCharge();
    double charge2 = a2->getCharge();
    this->term.dQ1Q2 = calculate_dQ1Q2(electrostaticScale, dQ1Q2Scale, charge1, charge2 );
    // ,energyNonbond->getDielectricConstant());
    // (this->_Charge1*this->_Charge2)/energyNonbond->getDielectricConstant();
    LOG("Calc dQ1Q2 electrostaticScale= {}", (double)(electrostaticScale));
    LOG("Calc dQ1Q2 Dielectric constant = {}", (double)(energyNonbond->getDielectricConstant()));
    LOG("Calc dQ1Q2 Charge1 = {}", (double)(charge1));
    LOG("Calc dQ1Q2 Charge2 = {}", (double)(charge2));
    LOG("dQ1Q2 = {}", (double)(this->term.dQ1Q2));
  }
  this->term.I1 = a1CoordinateIndexTimes3;
  this->term.I2 = a2CoordinateIndexTimes3;
  return true;
}


void EnergyNonbond_O::dumpTerms(core::HashTable_sp atomTypes) {
  gctools::Vec0<EnergyNonbond>::iterator eni;
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
}

void EnergyNonbond_O::setupHessianPreconditioner(NVector_sp nvPosition, AbstractLargeSquareMatrix_sp m, core::T_sp activeAtomMask) {
  SIMPLE_ERROR("Nonbond term isn't used when calculating setupHessianPreconditioner but it was called!!! only the bonded "
               "components of energy are used for the precondition to keep it sparse");
}

CL_DEFMETHOD void EnergyNonbond_O::expandExcludedAtomsToTerms(ScoringFunction_sp score, core::T_sp energyScale ) {
  double dielectricConstant;
  double dQ1Q2Scale;
  double cutoff;
  EnergyFunction_sp energyFunction = energyFunctionNonbondParameters(score, energyScale, dielectricConstant, dQ1Q2Scale, cutoff);
  //  printf("%s:%d In evaluateUsingExcludedAtoms starting this->_DebugEnergy -> %d\n", __FILE__, __LINE__, this->_DebugEnergy );
  if (!this->_iac_vec) {
    SIMPLE_ERROR("The nonbonded excluded atoms parameters have not been set up");
  }
  core::SimpleVector_int32_t_sp numberOfExcludedAtoms = this->_NumberOfExcludedAtomIndexes;
  core::SimpleVector_int32_t_sp excludedAtomIndexes = this->_ExcludedAtomIndexes;
  num_real electrostaticScale = 1.0;
  LOG("Nonbond component is enabled");
  // printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  num_real dA, dC;
  int I1, I2;
  int i = 0;
  //  int endIndex = pos->length()/3;
  int endIndex = this->_AtomTable->getNumberOfAtoms();
  //  printf( "%s:%d maxindex %d\n", __FILE__, __LINE__, endIndex);
  int nlocaltype = 0;
  // Find the max local type
  for (i = 0; i < this->_iac_vec->length(); ++i) {
    if (nlocaltype < (*this->_iac_vec)[i]) {
      nlocaltype = (*this->_iac_vec)[i];
    }
  }
  //  printf("%s:%d:%s    values -> %lu\n", __FILE__, __LINE__, __FUNCTION__, this->_Terms.size());

  // Count the number of nonbond cross terms that we are going to generate
  size_t excludedAtomIndex = 0;
  size_t count = 0;
  int index1_end = endIndex - 1;
  for (int index1 = 0; index1 < index1_end; ++index1) {
    LOG("{} ====== top of outer loop - index1 = {}\n", __FUNCTION__, index1);
    // Skip 0 in excluded atom list that amber requires
    bool has_excluded_atoms = ((*excludedAtomIndexes)[excludedAtomIndex] >= 0);
    int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
    for (int index2 = index1 + 1, index2_end(endIndex); index2 < index2_end; ++index2) {
      int maybe_excluded_atom = (*excludedAtomIndexes)[excludedAtomIndex];
      if (numberOfExcludedAtomsRemaining > 0 && maybe_excluded_atom == index2) {
        LOG("    Excluding atom {}\n", index2);
        ++excludedAtomIndex;
        --numberOfExcludedAtomsRemaining;
        continue;
      }
      ++count;
    }
    if (!has_excluded_atoms) {
      // No excluded atoms for the current index2 so increment excludedAtomIndex because inner loop didn't do it.
      // see Swails document http://ambermd.org/prmtop.pdf
      ++excludedAtomIndex;
    }
  }
  this->_Terms.resize(count);
  //  printf("%s:%d:%s    values -> %lu\n", __FILE__, __LINE__, __FUNCTION__, this->_Terms.size());

  // Now fill in the terms
  excludedAtomIndex = 0;

  size_t termIndex = 0;
  for (int index1 = 0; index1 < index1_end; ++index1) {
    LOG("{} ====== top of outer loop - index1 = {}\n", __FUNCTION__, index1);
    // Skip 0 in excluded atom list that amber requires
    bool has_excluded_atoms = ((*excludedAtomIndexes)[excludedAtomIndex] >= 0);
    int numberOfExcludedAtomsRemaining = numberOfExcludedAtoms->operator[](index1);
    num_real charge11 = (*this->_charge_vector)[index1];
    //    [[maybe_unused]]num_real electrostatic_scaled_charge11 = charge11*electrostaticScale;
    for (int index2 = index1 + 1, index2_end(endIndex); index2 < index2_end; ++index2) {
      LOG("    --- top of inner loop   numberOfExcludedAtomsRemaining -> {}    index2 -> {}\n", numberOfExcludedAtomsRemaining,
          index2);
      int maybe_excluded_atom = (*excludedAtomIndexes)[excludedAtomIndex];
      if (numberOfExcludedAtomsRemaining > 0 && maybe_excluded_atom == index2) {
        LOG("    Excluding atom {}\n", index2);
        ++excludedAtomIndex;
        --numberOfExcludedAtomsRemaining;
        continue;
      }
      //   cn1i    (*this->_cn1_vec)[cn1i]
      int localindex1 = (*this->_iac_vec)[index1];
      int localindex2 = (*this->_iac_vec)[index2];
      dA = (*this->_cn1_vec)[(*this->_ico_vec)[nlocaltype * (localindex1 - 1) + localindex2 - 1] - 1];
      dC = (*this->_cn2_vec)[(*this->_ico_vec)[nlocaltype * (localindex1 - 1) + localindex2 - 1] - 1];
      //      printf("%s:%d localindex1 %d and localindex2 %d\n", __FILE__, __LINE__, localindex1, localindex2);

      //      printf("%s:%d dA     %lf and dC     %lf\n", __FILE__, __LINE__, dA, dC);
      num_real charge22 = (*this->_charge_vector)[index2];

      I1 = index1 * 3;
      I2 = index2 * 3;
      EnergyNonbond enb;
      /// Atom-vector in AtomTable is empty
      enb._Atom1_enb = this->_AtomTable->elt_atom(index1);
      enb._Atom2_enb = this->_AtomTable->elt_atom(index2);
      enb.term.dA = dA;
      enb.term.dC = dC;
      num_real dQ1Q2 = calculate_dQ1Q2(1.0, dQ1Q2Scale, charge11, charge22);
      enb.term.dQ1Q2 = dQ1Q2;
      enb.term.I1 = I1;
      enb.term.I2 = I2;
      if (termIndex < count) {
        this->_Terms[termIndex] = enb;
      } else {
        SIMPLE_ERROR("Overflowed the _Terms array with termIndex={} and count = {}\n", termIndex, count);
      }
      ++termIndex;
//      printf( "nonbond index1 name %s index2 %s\n",  this->_AtomTable->elt_atom_name(index1),
//      this->_AtomTable->elt_atom_name(index2));
    }
    if (!has_excluded_atoms) {
      // No excluded atoms for the current index2 so increment excludedAtomIndex because inner loop didn't do it.
      // see Swails document http://ambermd.org/prmtop.pdf
      ++excludedAtomIndex;
    }
  }
  //  printf("%s:%d:%s    values -> %lu\n", __FILE__, __LINE__, __FUNCTION__, this->_Terms.size());
  //  printf( "Nonbond energy vdw({}) electrostatic({})\n", (double)this->_EnergyVdw,  this->_EnergyElectrostatic );
  LOG("Nonbond energy vdw({}) electrostatic({})\n", (double)this->_EnergyVdw, this->_EnergyElectrostatic);
  LOG("Nonbond energy }\n");
}

SYMBOL_EXPORT_SC_(KeywordPkg, nonbond);
SYMBOL_EXPORT_SC_(KeywordPkg, vdwRadii);
SYMBOL_EXPORT_SC_(KeywordPkg, distance);
SYMBOL_EXPORT_SC_(KeywordPkg, force);

void EnergyNonbond_O::initialize() {
  this->Base::initialize();
}

void EnergyNonbond_O::addTerm(const EnergyNonbond &term) { this->_Terms.push_back(term); }
void EnergyNonbond_O::addTerm14(const EnergyNonbond &term) { this->_Terms14.push_back(term); }

void EnergyNonbond_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw, terms), this->_Terms);
  node->field(INTERN_(kw, terms14), this->_Terms14);
  node->field(INTERN_(kw, InteractionsKept), this->_InteractionsKept );
  node->field(INTERN_(kw, InteractionsDiscarded), this->_InteractionsDiscarded );
  node->field(INTERN_(kw, FFNonbondDb), this->_FFNonbondDb);
  node->field(INTERN_(kw, AtomTable), this->_AtomTable);
  node->field(INTERN_(kw, ntypes), this->_ntypes);
  node->field(INTERN_(kw, atom_name_vector), this->_atom_name_vector);         // atom-name-vector
  node->field(INTERN_(kw, charge_vector), this->_charge_vector);               // charge-vector
  node->field(INTERN_(kw, mass_vector), this->_mass_vector);                   // masses
  node->field(INTERN_(kw, atomic_number_vector), this->_atomic_number_vector); // vec
  node->field(INTERN_(kw, ico_vec), this->_ico_vec);                           // ico-vec
  node->field(INTERN_(kw, iac_vec), this->_iac_vec);                           // iac-vec
  node->field(INTERN_(kw, atom_type_vector), this->_atom_type_vector);         // Amber atom type names
  node->field(INTERN_(kw, local_typej_vec), this->_local_typej_vec);           // local-typej-vec
  node->field(INTERN_(kw, cn1_vec), this->_cn1_vec);
  node->field(INTERN_(kw, cn2_vec), this->_cn2_vec);
  node->field(INTERN_(kw, NumberOfExcludedAtomIndexes), this->_NumberOfExcludedAtomIndexes);
  node->field(INTERN_(kw, ExcludedAtomIndexes), this->_ExcludedAtomIndexes);
  node->field(INTERN_(kw, UsesExcludedAtoms), this->_UsesExcludedAtoms);
  node->field_if_not_nil(INTERN_(kw, Matter1), this->_Matter1);
  node->field_if_not_nil(INTERN_(kw, Matter2), this->_Matter2);
  node->field(INTERN_(kw, KeepInteractionFactory), this->_KeepInteractionFactory);
  this->Base::fields(node);
}

string EnergyNonbond_O::beyondThresholdInteractionsAsString() {
  return component_beyondThresholdInteractionsAsString<EnergyNonbond_O, EnergyNonbond>(*this);
}

void EnergyNonbond_O::construct14InteractionTerms(AtomTable_sp atomTable, Matter_sp matter, core::T_sp forceField,
                                                  core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes) {
  if (keepInteractionFactory.nilp()) return;
  core::T_sp keepInteraction = specializeKeepInteractionFactory(keepInteractionFactory,EnergyNonbond_O::staticClass());
  bool hasKeepInteractionFunction = gc::IsA<core::Function_sp>(keepInteraction);
  EnergyNonbond energyNonbond;
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
      energyNonbond.defineForAtomPair(forceField, true,
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


void EnergyNonbond_O::constructNonbondTermsFromAtomTable(AtomTable_sp atomTable, core::T_sp nbForceField,
                                                         core::HashTable_sp atomTypes, core::T_sp keepInteractionFactory,
                                                         core::T_sp tcoordinates ) {
  this->_AtomTable = atomTable;
  this->_NonbondForceField = nbForceField;
  this->_AtomTypes = atomTypes;
  this->_KeepInteractionFactory = keepInteractionFactory;
  if (tcoordinates.notnilp()) {
    this->rebuildPairList(tcoordinates);
  } else {
    // Don't calculate the pair-list if we don't have coordinates
  }
};

CL_DEFMETHOD
core::T_mv EnergyNonbond_O::maybeRebuildPairList(core::T_sp tcoordinates) {
  auto coords = gc::As<NVector_sp>(tcoordinates);
  if (this->_DisplacementBuffer.nilp()) {
    return this->rebuildPairList(tcoordinates);
  } else if (gc::IsA<NVector_sp>(this->_DisplacementBuffer)) {
    NVector_sp nvDisplacementBuffer = gc::As_unsafe<NVector_sp>(this->_DisplacementBuffer);
    if (nvDisplacementBuffer->size() != coords->size()) {
      SIMPLE_ERROR("The size of the _DispacementBuffer({}) MUST match the size of the coordinatess({})", nvDisplacementBuffer->size(), coords->size());
    }
    double skinThickness = this->_Nonbond_r_pairlist-this->_Nonbond_r_cut;
    double movedTrigger = 0.5*skinThickness;
    double movedTrigger2 = movedTrigger*movedTrigger;
    vecreal* raw_db = &(*nvDisplacementBuffer)[0];
    vecreal* raw_coords = &(*coords)[0];
    for (size_t ci=0; ci<nvDisplacementBuffer->size(); ci+=3 ) {
      const vecreal& dx = raw_db[ci];
      const vecreal& dy = raw_db[ci+1];
      const vecreal& dz = raw_db[ci+2];
      const vecreal& cx = raw_coords[ci];
      const vecreal& cy = raw_coords[ci+1];
      const vecreal& cz = raw_coords[ci+2];
      if (dx==cx && dy==cy && dz==cz ) continue;
      vecreal dist2 = (dx-cx)*(dx-cx)+ (dy-cy)*(dy-cy)+ (dz-cz)*(dz-cz);
      if (dist2>movedTrigger2) {
        return this->rebuildPairList(tcoordinates);
      }
    }
    return Values0<core::T_O>();
  }
  SIMPLE_ERROR("{}: We should never get here", __FUNCTION__);
}

CL_DEFMETHOD
core::T_mv EnergyNonbond_O::rebuildPairList(core::T_sp tcoordinates) {

  // If the nonbond pairs are between two matters then don't rebuild the pair-list
  if (this->_Matter1.notnilp()) {
    ASSERT(this->_Matter2.notnilp());
    this->rebuildPairListBetweenMatters(tcoordinates);
    return Values0<core::T_O>();
  }

  this->_UsesExcludedAtoms = false;
  this->_DisplacementBuffer = copy_nvector(gc::As<NVector_sp>(tcoordinates));
  size_t interactionsKept = 0;
  size_t interactionsDiscarded = 0;
  size_t totalInteractions = 0;
  if (this->_KeepInteractionFactory.notnilp()) {
    this->_Terms.clear();
    core::T_sp keepInteraction = specializeKeepInteractionFactory(this->_KeepInteractionFactory,EnergyNonbond_O::staticClass());
    bool hasKeepInteractionFunction = gc::IsA<core::Function_sp>(keepInteraction);
    double r_pairlist2 = this->_Nonbond_r_pairlist*this->_Nonbond_r_pairlist;

    // ------------------------------------------------------------
    //
    if (this->_AtomTable->getNumberOfAtoms() >= 2) {
      LOG("Defining NONBONDS");
      vecreal* coords = NULL;
      if (gc::IsA<NVector_sp>(tcoordinates)) {
        coords = &(*gc::As_unsafe<NVector_sp>(tcoordinates))[0];
      }
      gctools::Vec0<EnergyAtom>::iterator iea1;
      gctools::Vec0<EnergyAtom>::iterator iea2;
      size_t total_comparisons = this->_AtomTable->getNumberOfAtoms() * this->_AtomTable->getNumberOfAtoms() / 2;
      if (chem__verbose(0)) {
        core::clasp_write_string(
            fmt::format("For nonbonded interactions, about to carry out {} atom-to-atom comparisons\n", total_comparisons));
      }
      for (iea1 = this->_AtomTable->begin(); iea1 != this->_AtomTable->end() - 1; iea1++) {
        Vector3 pos1;
        if (coords) {
          pos1.set(coords[iea1->coordinateIndexTimes3()],
                   coords[iea1->coordinateIndexTimes3()+1],
                   coords[iea1->coordinateIndexTimes3()+2]);
        } else {
          pos1 = iea1->atom()->getPosition();
        }
        for (iea2 = iea1 + 1; iea2 != this->_AtomTable->end(); iea2++) {
          if (!(iea1->inBondOrAngle(iea2->atom()) || iea1->relatedBy14(iea2->atom()))) {
            totalInteractions++;
            Vector3 pos2;
            if (coords) {
              pos2.set(coords[iea2->coordinateIndexTimes3()],
                       coords[iea2->coordinateIndexTimes3()+1],
                       coords[iea2->coordinateIndexTimes3()+2]);
            } else {
              pos2 = iea2->atom()->getPosition();
            }
            Vector3 posDelta = pos1-pos2;
            double dist2 = posDelta.dotProduct(posDelta);
            if (dist2<r_pairlist2) {
              bool keep = true;
              if (hasKeepInteractionFunction) {
                core::T_sp result = core::eval::funcall(keepInteraction,
                                                        iea1->atom(),iea2->atom(),
                                                        core::make_fixnum(iea1->coordinateIndexTimes3()),
                                                        core::make_fixnum(iea2->coordinateIndexTimes3()));
                keep = result.notnilp();
              }
              if (keep) {
                EnergyNonbond energyNonbond;
                energyNonbond.defineForAtomPair(this->_NonbondForceField,
                                                false,
                                                iea1->atom(), iea2->atom(),
                                                iea1->coordinateIndexTimes3(),
                                                iea2->coordinateIndexTimes3(),
                                                this->asSmartPtr(),
                                                this->_AtomTypes,
                                                keepInteraction);
                this->addTerm(energyNonbond);
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

core::T_mv EnergyNonbond_O::rebuildPairListBetweenMatters(core::T_sp tcoordinates) {
  core::print(fmt::format("DEBUG/in rebuildPairListBetweenMatters --------------\n"));
  core::T_sp keepInteractionFactory = this->_KeepInteractionFactory;
  if (keepInteractionFactory.nilp()) return Values0<core::T_O>();
  NVector_sp coords = gc::As<NVector_sp>(tcoordinates);
  core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyNonbond_O::staticClass() );
  Matter_sp mat1 = this->_Matter1;
  Matter_sp mat2 = this->_Matter2;
  bool hasKeepInteractionFunction = gc::IsA<core::Function_sp>(keepInteraction);
  double rpairlist2 = this->_Nonbond_r_pairlist*this->_Nonbond_r_pairlist;
  auto atomTable = this->_AtomTable;
  auto nbForceField = atomTable->nonbondForceFieldForAggregate();
  auto atomTypes = this->_AtomTypes; // This was energyFunction->atomTypes() - I hope this->_AtomTypes works the same.
  {
    LOG("Defining NONBONDS");
    Loop lMat1(mat1, ATOMS);
    while (lMat1.advanceLoopAndProcess()) {
      Atom_sp a1 = lMat1.getAtom();
      size_t i3x1 = atomTable->getCoordinateIndexTimes3(a1);
      Vector3 v1(coords,i3x1,Safe());
      Loop lMat2(mat2, ATOMS);
      while (lMat2.advanceLoopAndProcess()) {
        Atom_sp a2 = lMat2.getAtom();
        size_t i3x2 = atomTable->getCoordinateIndexTimes3(a2);
        Vector3 v2(coords,i3x2,Safe());
        Vector3 vdiff;
        vdiff = v1-v2;
        bool added = false;
        double dist2 = vdiff.dotProduct(vdiff);
//        if ((i3x1==609&&i3x1==480) || (i3x1==480&&i3x1==609)) {
        core::print(fmt::format("ai{}/ai{} dist2 = {}  rpairlist2 = {}\n", i3x1, i3x2, dist2, rpairlist2));
//        }
        if (dist2 < rpairlist2) {
          EnergyNonbond energyNonbond;
          bool in14 = false; // Between residues there are never 1-4
          if (hasKeepInteractionFunction){
            core::T_sp result = core::eval::funcall(keepInteraction,a1,a2,
                                                    core::make_fixnum(i3x1),
                                                    core::make_fixnum(i3x2));
            if (result.notnilp()) {
              energyNonbond.defineForAtomPair(nbForceField, in14, a1, a2, i3x1, i3x2, this->asSmartPtr(), atomTypes, keepInteraction);
              added = true;
              this->addTerm(energyNonbond);
            }
          } else {
            energyNonbond.defineForAtomPair(nbForceField, in14, a1, a2,
                                            i3x1, i3x2,
                                            this->asSmartPtr(), atomTypes, keepInteraction );
            added = true;
            this->addTerm(energyNonbond);
          }
        }
        if (!added) {
          if (_sym_trapEnergyNonbondNotDefineForAtomPair->fboundp()) {
            core::eval::funcall(_sym_trapEnergyNonbondDefineForAtomPair,
                                a1, a2,
                                core::make_fixnum(i3x1),
                                core::make_fixnum(i3x2),
                                this->asSmartPtr(),
                                atomTypes,
                                keepInteraction);
          }
        }
      }
    }
  }
  return Values0<core::T_O>();
}






void EnergyNonbond_O::constructExcludedAtomListFromAtomTable(AtomTable_sp atomTable, core::T_sp nbForceField, core::T_sp keepInteractionFactory ) {
  // ------------------------------------------------------------
  //
  // The new code creates an excluded atom list
  //
  this->_UsesExcludedAtoms = true;
  this->_AtomTable = atomTable;
  this->_FFNonbondDb = nbForceField;
  core::T_mv values_mv = this->_AtomTable->calculate_excluded_atom_list(keepInteractionFactory);
  core::SimpleVector_int32_t_sp number_of_excluded_atoms = gc::As<core::SimpleVector_int32_t_sp>(values_mv);
  core::MultipleValues &values = core::lisp_multipleValues();
  core::SimpleVector_int32_t_sp excluded_atoms_list =
      gc::As<core::SimpleVector_int32_t_sp>(values.second(values_mv.number_of_values()));
  this->_NumberOfExcludedAtomIndexes = number_of_excluded_atoms;
  this->_ExcludedAtomIndexes = excluded_atoms_list;
}

/* Construct nonbond terms between two molecules or two residues or a residue and a molecule that are not
bonded to each other.  Generate a pair-list in _Terms and that pair-list should never change. */
CL_DEFMETHOD void EnergyNonbond_O::constructNonbondTermsBetweenMatters(Matter_sp mat1, Matter_sp mat2,
                                                                       EnergyFunction_sp energyFunction,
                                                                       core::T_sp keepInteractionFactory) {
  core::print(fmt::format("DEBUG/in constructNonbondTermsBetweenMatters -------------- {} {}\n", _rep_(mat1), _rep_(mat2)));
  this->_UsesExcludedAtoms = false;
  this->_Matter1 = mat1;
  this->_Matter2 = mat2;
  this->_KeepInteractionFactory = keepInteractionFactory;
  this->_Terms.clear();
  this->_Terms14.clear();


}

SYMBOL_EXPORT_SC_(KeywordPkg, da);
SYMBOL_EXPORT_SC_(KeywordPkg, dc);
SYMBOL_EXPORT_SC_(KeywordPkg, i1);
SYMBOL_EXPORT_SC_(KeywordPkg, i2);
SYMBOL_EXPORT_SC_(KeywordPkg, charge1);
SYMBOL_EXPORT_SC_(KeywordPkg, charge2);
SYMBOL_EXPORT_SC_(KeywordPkg, atomname1);
SYMBOL_EXPORT_SC_(KeywordPkg, atomname2);

CL_DEFMETHOD core::List_sp EnergyNonbond_O::extract_vectors_as_alist() const {
  size_t size = this->_Terms.size();
  NVector_sp da_vec = NVector_O::make(size);
  NVector_sp dc_vec = NVector_O::make(size);
  core::SimpleVector_int32_t_sp i1_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i2_vec = core::SimpleVector_int32_t_O::make(size);
  NVector_sp charge1_vec = NVector_O::make(size);
  NVector_sp charge2_vec = NVector_O::make(size);
  core::SimpleVector_sp atomname1_vec = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atomname2_vec = core::SimpleVector_O::make(size);
  for (size_t i = 0; i < size; ++i) {
    const EnergyNonbond &entry = this->_Terms[i];
    (*da_vec)[i] = entry.term.dA;
    (*dc_vec)[i] = entry.term.dC;
    (*i1_vec)[i] = entry.term.I1;
    (*i2_vec)[i] = entry.term.I2;
    (*charge1_vec)[i] = this->_AtomTable->elt_charge(entry.term.I1 / 3);
    (*charge2_vec)[i] = this->_AtomTable->elt_charge(entry.term.I2 / 3);
    (*atomname1_vec)[i] = this->_AtomTable->elt_atom_name(entry.term.I1 / 3);
    (*atomname2_vec)[i] = this->_AtomTable->elt_atom_name(entry.term.I2 / 3);
  }
  return core::Cons_O::createList(
      core::Cons_O::create(kw::_sym_da, da_vec), core::Cons_O::create(kw::_sym_dc, dc_vec),
      core::Cons_O::create(kw::_sym_i1, i1_vec), core::Cons_O::create(kw::_sym_i2, i2_vec),
      core::Cons_O::create(kw::_sym_charge1, charge1_vec), core::Cons_O::create(kw::_sym_charge2, charge2_vec),
      core::Cons_O::create(kw::_sym_atomname1, atomname1_vec), core::Cons_O::create(kw::_sym_atomname2, atomname2_vec));
}

CL_DEFMETHOD core::T_sp EnergyNonbond_O::getFFNonbondDb() {
  if (this->_FFNonbondDb.boundp())
    return this->_FFNonbondDb;
  SIMPLE_ERROR("The _FFNonbondDb of an EnergyNonbond has not been initialized");
}

SYMBOL_EXPORT_SC_(KeywordPkg, ntypes);
SYMBOL_EXPORT_SC_(KeywordPkg, atom_name_vector);
SYMBOL_EXPORT_SC_(KeywordPkg, atom_type_vector);
SYMBOL_EXPORT_SC_(KeywordPkg, charge_vector);
SYMBOL_EXPORT_SC_(KeywordPkg, mass_vector);
SYMBOL_EXPORT_SC_(KeywordPkg, atomic_number_vector);
SYMBOL_EXPORT_SC_(KeywordPkg, ico_vec);
SYMBOL_EXPORT_SC_(KeywordPkg, iac_vec);
SYMBOL_EXPORT_SC_(KeywordPkg, local_typej_vec);
SYMBOL_EXPORT_SC_(KeywordPkg, cn1_vec);
SYMBOL_EXPORT_SC_(KeywordPkg, cn2_vec);
SYMBOL_EXPORT_SC_(KeywordPkg, excluded_atoms_list);
SYMBOL_EXPORT_SC_(KeywordPkg, number_excluded_atoms);
SYMBOL_EXPORT_SC_(KeywordPkg, atom_table);

  CL_DOCSTRING(R"doc(For each nonbond term call the CALLBACK function.
The CALLBACK function must accept two arguments:
1. The first atom of the nonbond term
2. The second atom of the nonbond term
3. The Ix3 coordinate for the first atom
4. The Ix3 coordinate for the second atom)doc");
CL_DEFMETHOD void EnergyNonbond_O::callForEachTerm(core::Function_sp callback) {
    if (this->_UsesExcludedAtoms) {
      SIMPLE_ERROR("callForEachTerm does not support this nonbond component because it uses excluded atoms");
    }
  for (auto eni = this->_Terms.begin(); eni != this->_Terms.end(); eni++) {
    core::eval::funcall(callback,eni->_Atom1_enb,eni->_Atom2_enb,
                        core::make_fixnum(eni->term.I1),
                        core::make_fixnum(eni->term.I2));
  }
}

  CL_DOCSTRING(R"doc(For each nonbond term call the CALLBACK function.
The CALLBACK function must accept two arguments:
1. The first atom of the nonbond term
2. The second atom of the nonbond term
3. The Ix3 coordinate for the first atom
4. The Ix3 coordinate for the second atom)doc");
CL_DEFMETHOD void EnergyNonbond_O::callForEachTerm14(core::Function_sp callback) {
    if (this->_UsesExcludedAtoms) {
      SIMPLE_ERROR("callForEachTerm does not support this nonbond component because it uses excluded atoms");
    }
  for (auto eni = this->_Terms14.begin(); eni != this->_Terms14.end(); eni++) {
    core::eval::funcall(callback,eni->_Atom1_enb,eni->_Atom2_enb,
                        core::make_fixnum(eni->term.I1),
                        core::make_fixnum(eni->term.I2));
  }
}


CL_DEFMETHOD void EnergyNonbond_O::constructNonbondTermsFromAList(core::List_sp values) {
  this->_ntypes = translate::from_object<size_t>(safe_alist_lookup<core::T_sp>(values, kw::_sym_ntypes))._v; // ntypes
  this->_atom_name_vector = safe_alist_lookup<core::SimpleVector_sp>(values, kw::_sym_atom_name_vector);     // atom-name-vector
  this->_atom_type_vector = safe_alist_lookup<core::SimpleVector_sp>(values, kw::_sym_atom_type_vector);     // atom-type-vector
  this->_charge_vector = safe_alist_lookup<NVector_sp>(values, kw::_sym_charge_vector);                      // charge-vector
  this->_mass_vector = safe_alist_lookup<NVector_sp>(values, kw::_sym_mass_vector);                          // masses
  this->_atomic_number_vector = safe_alist_lookup<core::SimpleVector_int32_t_sp>(values, kw::_sym_atomic_number_vector); // vec
  this->_ico_vec = safe_alist_lookup<core::SimpleVector_int32_t_sp>(values, kw::_sym_ico_vec);                           // ico-vec
  this->_iac_vec = safe_alist_lookup<core::SimpleVector_int32_t_sp>(values, kw::_sym_iac_vec);                           // iac-vec
  //  this->_local_typej_vec =      safe_alist_lookup(values,kw::_sym_local_typej_vec);      // local-typej-vec
  this->_cn1_vec = safe_alist_lookup<NVector_sp>(values, kw::_sym_cn1_vec);
  this->_cn2_vec = safe_alist_lookup<NVector_sp>(values, kw::_sym_cn2_vec);
}

CL_DOCSTRING(R"dx(These are the nonbond terms as they are read from the topology file.
The AMBER charge conversion of 18.2223 has been divided out of the charge values.
The results of this call should be able to be passed to EnergyNonbond_O::constructNonbondTermsFromAList.)dx")
CL_DEFMETHOD core::List_sp EnergyNonbond_O::nonbondTermsAsAList() {
  ql::list result;
  result << core::Cons_O::create(kw::_sym_ntypes, core::make_fixnum(this->_ntypes));          // ntypes
  result << core::Cons_O::create(kw::_sym_atom_name_vector, this->_atom_name_vector);         // atom-name-vector
  result << core::Cons_O::create(kw::_sym_atom_type_vector, this->_atom_type_vector);         // atom-type-vector
  result << core::Cons_O::create(kw::_sym_charge_vector, this->_charge_vector);               // charge-vector
  result << core::Cons_O::create(kw::_sym_mass_vector, this->_mass_vector);                   // masses
  result << core::Cons_O::create(kw::_sym_atomic_number_vector, this->_atomic_number_vector); // vec
  result << core::Cons_O::create(kw::_sym_ico_vec, this->_ico_vec);                           // ico-vec
  result << core::Cons_O::create(kw::_sym_iac_vec, this->_iac_vec);                           // iac-vec
  //  this->_local_typej_vec =      safe_alist_lookup(values,kw::_sym_local_typej_vec);      // local-typej-vec
  result << core::Cons_O::create(kw::_sym_cn1_vec, this->_cn1_vec); // cn1
  result << core::Cons_O::create(kw::_sym_cn2_vec, this->_cn2_vec); // cn2
  return result.cons();
}

CL_DEFMETHOD void EnergyNonbond_O::setNonbondExcludedAtomInfo(AtomTable_sp atom_table,
                                                              core::SimpleVector_int32_t_sp excluded_atoms_list,
                                                              core::SimpleVector_int32_t_sp number_excluded_atoms) {
  this->_AtomTable = atom_table;
  this->_ExcludedAtomIndexes = excluded_atoms_list;
  this->_NumberOfExcludedAtomIndexes = number_excluded_atoms;
}


  CL_DEFMETHOD
  void EnergyNonbond_O::checkEnergyNonbond() {
    if (!this->_AtomTable) {
      SIMPLE_ERROR("EnergyNonbond _AtomTable is NULL");
    }
  }

EnergyNonbond_sp EnergyNonbond_O::copyFilter(core::T_sp keepInteractionFactory) {
  EnergyNonbond_sp copy = EnergyNonbond_O::create();
  copyEnergyComponent( copy, this->asSmartPtr() );
  core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyNonbond_O::staticClass() );

#define COPY_FIELD(xxx) { copy->xxx = this->xxx; }
  COPY_FIELD(_EnergyElectrostatic);
  COPY_FIELD(_UsesExcludedAtoms);
  COPY_FIELD(_NonbondsKept);
  COPY_FIELD(_NonbondsDiscarded);
  COPY_FIELD(_FFNonbondDb);
  COPY_FIELD(_AtomTable);
  COPY_FIELD(_ntypes);
  COPY_FIELD(_NonbondForceField);
  COPY_FIELD(_AtomTypes);
  COPY_FIELD(_atom_name_vector);
  COPY_FIELD(_charge_vector);
  COPY_FIELD(_mass_vector);
  COPY_FIELD(_atomic_number_vector);
  COPY_FIELD(_ico_vec);
  COPY_FIELD(_iac_vec);
  COPY_FIELD(_atom_type_vector);
  COPY_FIELD(_local_typej_vec);
  COPY_FIELD(_cn1_vec);
  COPY_FIELD(_cn2_vec);
  COPY_FIELD(_NumberOfExcludedAtomIndexes);
  COPY_FIELD(_ExcludedAtomIndexes);
  COPY_FIELD(_InteractionsKept);
  COPY_FIELD(_InteractionsDiscarded);
  COPY_FIELD(_Matter1);
  COPY_FIELD(_Matter2);

  // Copy the interaction factory
  copy->_KeepInteractionFactory = keepInteractionFactory;
  if (!this->_UsesExcludedAtoms) {
    if (this->_DisplacementBuffer.notnilp()) {
      copy->_DisplacementBuffer = copy_nvector(gc::As<NVector_sp>(this->_DisplacementBuffer));
    }
    // Copy the nonbond pair-list terms
    for (auto edi = this->_Terms.begin(); edi != this->_Terms.end(); edi++) {
      Atom_sp a1 = edi->_Atom1_enb;
      Atom_sp a2 = edi->_Atom2_enb;
      size_t ia1 = edi->term.I1;
      size_t ia2 = edi->term.I2;
      if (skipInteraction_EnergyNonbond(keepInteraction, a1, a2, core::make_fixnum(ia1), core::make_fixnum(ia2))) continue;
      copy->_Terms.push_back(*edi);
    }
    // Copy the 1-4 interactions
    for (auto edi = this->_Terms14.begin(); edi != this->_Terms14.end(); edi++) {
      Atom_sp a1 = edi->_Atom1_enb;
      Atom_sp a2 = edi->_Atom2_enb;
      size_t ia1 = edi->term.I1;
      size_t ia2 = edi->term.I2;
      if (skipInteraction_EnergyNonbond(keepInteraction, a1, a2, core::make_fixnum(ia1), core::make_fixnum(ia2))) continue;
      copy->_Terms14.push_back(*edi);
    }
  } else {
    if (keepInteraction==_lisp->_true()) {
      // Copy the 1-4 interactions
      for (auto edi = this->_Terms14.begin(); edi != this->_Terms14.end(); edi++) {
        copy->_Terms14.push_back(*edi);
      }
    } else {
      SIMPLE_ERROR("EnergyNonbond uses excluded atoms and we don't support keepInteraction = {}", _rep_(keepInteraction));
    }
  }
  return copy;
}

void EnergyNonbond_O::emitTestCalls(core::T_sp stream, chem::NVector_sp pos) const
{
  SIMPLE_ERROR("Update me");
#if 0
  if (!this->_UsesExcludedAtoms) {
    for ( auto si=this->_Terms.begin();
          si!=this->_Terms.end(); si++ ) {
      core::print(fmt::format("test_nonbond( errs, {}, {}, {}, {}, {}, position_size, position, force_new, force_old, hessian_new, hessian_old, dvec_new, dvec_old, hdvec_new, hdvec_old );\n",
                              si->term.dA, si->term.dC, si->term.dQ1Q2, si->term.I1, si->term.I2 ), stream);
    }
  } else {
    SIMPLE_ERROR("You cannot emit test calls for excluded atoms");
  }
#endif
}



}; // namespace chem
