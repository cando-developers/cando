/*
    File: energyDihedral.cc
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
#define	DEBUG_LEVEL_NONE


#include <cando/main/extension.h>
#include <cando/chem/energyDihedral.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <clasp/core/ql.h>
#include <cando/chem/bond.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffPtorDb.h>
#include <cando/chem/ffItorDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/wrappers.h>

#define VEC8(v) {v,v,v,v,v,v,v,v}
#define VEC4(v) {v,v,v,v}
#define VEC2(v) {v,v}

namespace chem {

CL_LISPIFY_NAME(make-energy-dihedral);
CL_DEF_CLASS_METHOD
EnergyDihedral_sp EnergyDihedral_O::make(EnergyFunction_sp energyFunction) {
  auto component = ensureComponent<EnergyDihedral_O>(energyFunction);
  return component;
}

};

// Trig-free multiple-angle recurrence (cos_n_phi/sin_n_phi) used by the
// dihedral_fast_floored kernel; shared with the dihedralKernelCompare.cc test.
#include <cando/chem/energyKernels/dihedralRecurrence.h>

namespace chem {

// ---- Dihedral kernel A/B selector ------------------------------------------
//   0 = dihedral_fast         : AD geometry derivative; uses atan2 + sin + cos
//   1 = dihedral_fast_floored : manual Blondel-Karplus gradient + trig-free
//                               multiple-angle recurrence (linear-dihedral safe)
// Set globally to A/B, e.g.  -DCANDO_DIHEDRAL_KERNEL=1  (see energyDihedral.h:
// the SAME macro selects the matching kernel header there, so dihedral_term is
// defined exactly once per translation unit).  The default (0) is also defined
// in energyDihedral.h, which is included above.
#ifndef CANDO_DIHEDRAL_KERNEL
#define CANDO_DIHEDRAL_KERNEL 0
#endif

#if CANDO_DIHEDRAL_KERNEL == 1
#include "cando/chem/energyKernels/dihedral_fast_floored.c"
template <typename T>
using Dihedral_Component = Dihedral_Fast_Floored<T>;
#else
#include "cando/chem/energyKernels/dihedral_fast.c"
template <typename T>
using Dihedral_Component = Dihedral_Fast<T>;
#endif


std::string EnergyDihedral_O::descriptionOfContents() const {
  stringstream ss;
  ss << ":enabled " << ((this->_Enabled) ? "T" : "NIL");
  ss << " number-of-terms " << this->_Terms.size();
  return ss.str();
}


std::string EnergyDihedral_O::implementation_details() const {
  Dihedral_Component<NoHessian> dihedral;

  std::stringstream ss;
  ss << dihedral.description();
  return ss.str();
}


size_t EnergyDihedral_O::runTestCalls(core::T_sp stream, chem::NVector_sp coords) const
{
  #define POS_SIZE 12
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
  size_t idx = 0;
  size_t errs = 0;
  Dihedral_Component<double*> dihedral;
  for ( auto di=this->_Terms.begin(); di!=this->_Terms.end(); ++di ) {
    position[0] = coords[di->term.i3x1];
    position[1] = coords[di->term.i3x1+1];
    position[2] = coords[di->term.i3x1+2];
    position[3] = coords[di->term.i3x2];
    position[4] = coords[di->term.i3x2+1];
    position[5] = coords[di->term.i3x2+2];
    position[6] = coords[di->term.i3x3];
    position[7] = coords[di->term.i3x3+1];
    position[8] = coords[di->term.i3x3+2];
    position[9] = coords[di->term.i3x4];
    position[10] = coords[di->term.i3x4+1];
    position[11] = coords[di->term.i3x4+2];
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    dihedral_term lt(di->term.v, di->term.n, di->term.sinphase, di->term.cosphase, 0, 3, 6, 9);
    dihedral.gradient(    lt, position, &energy_new,    force_new );
    dihedral.gradient_fd( lt, position, &energy_ground, force_ground );
    if (!test_match( stream, "dihedral_gradient", POS_SIZE,
                     force_new, force_ground,
                     0, 0,
                     0, 0 )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH dihedral_gradient #{} V = {}  DN = {}  sinPhase = {}  cosPhase = {}\n",
                              idx, di->term.v, di->term.n, di->term.sinphase, di->term.cosphase ), stream );
    }
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    dihedral_term lt2(di->term.v, di->term.n, di->term.sinphase, di->term.cosphase, 0, 3, 6, 9);
    dihedral.hessian(    lt2, position, &energy_new,    force_new,    hessian_new,    dvec_new,    hdvec_new );
    dihedral.hessian_fd( lt2, position, &energy_ground, force_ground, hessian_ground, dvec_ground, hdvec_ground );
    if (!test_match( stream, "dihedral_hessian", POS_SIZE,
                     force_new, force_ground,
                     hessian_new, hessian_ground,
                     hdvec_new, hdvec_ground )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH dihedral_hessian #{} V = {}  DN = {}  sinPhase = {}  cosPhase = {}\n",
                              idx, di->term.v, di->term.n, di->term.sinphase, di->term.cosphase ), stream );
    }
    idx++;
  }
  core::print(fmt::format("dihedral errors = {}\n", errs), stream);
  return errs;
}

}


namespace chem {

template <int N,typename Real>
struct SinCos {
  static void sinNPhiCosNPhi(Real& sinNPhi, Real& cosNPhi,
                      Real sinPhi, Real cosPhi ) {
    Real sinNm1Phi;
    Real cosNm1Phi;
    SinCos<N-1,Real>::sinNPhiCosNPhi(sinNm1Phi,cosNm1Phi,sinPhi,cosPhi);
    sinNPhi = cosPhi*sinNm1Phi+sinPhi*cosNm1Phi;
    cosNPhi = cosPhi*cosNm1Phi-sinPhi*sinNm1Phi;
  }
};

template <typename Real>
struct SinCos<1,Real> {
  static void sinNPhiCosNPhi(Real& sinNPhi, Real& cosNPhi,
                      Real sinPhi, Real cosPhi ) {
    sinNPhi = sinPhi;
    cosNPhi = cosPhi;
  }
};


template <typename Real>
void	sinNPhiCosNPhi(int n, Real& sinNPhi, Real& cosNPhi,
                       Real sinPhi, Real cosPhi )
{
  switch (n) {
  case 1:
      SinCos<1,Real>::sinNPhiCosNPhi(sinNPhi,cosNPhi,sinPhi,cosPhi);
      break;
  case 2:
      SinCos<2,Real>::sinNPhiCosNPhi(sinNPhi,cosNPhi,sinPhi,cosPhi);
      break;
  case 3:
      SinCos<3,Real>::sinNPhiCosNPhi(sinNPhi,cosNPhi,sinPhi,cosPhi);
      break;
  case 4:
      SinCos<4,Real>::sinNPhiCosNPhi(sinNPhi,cosNPhi,sinPhi,cosPhi);
      break;
  case 5:
      SinCos<5,Real>::sinNPhiCosNPhi(sinNPhi,cosNPhi,sinPhi,cosPhi);
      break;
  case 6:
      SinCos<6,Real>::sinNPhiCosNPhi(sinNPhi,cosNPhi,sinPhi,cosPhi);
      break;
  default:
      SIMPLE_ERROR("N {} should always be 1-6", n);
      break;
  }
}

#define ENERGY_FUNCTION I1, I2, I3, I4, activeAtomMask

SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralIn1);
SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralIn2);
SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralIn3);
SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralIn4);
SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralIn5);
SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralIn6);
SYMBOL_EXPORT_SC_(ChemPkg,EnergyDihedralUnknown);

core::T_sp dihedral_type(int in) {
  switch (in) {
  case 1:
      return chem::_sym_EnergyDihedralIn1;
  case 2:
      return chem::_sym_EnergyDihedralIn2;
  case 3:
      return chem::_sym_EnergyDihedralIn3;
  case 4:
      return chem::_sym_EnergyDihedralIn4;
  case 5:
      return chem::_sym_EnergyDihedralIn5;
  case 6:
      return chem::_sym_EnergyDihedralIn6;
  }
  return chem::_sym_EnergyDihedralUnknown;
}

#undef DIHEDRAL_APPLY_ATOM_MASK
#define DIHEDRAL_APPLY_ATOM_MASK(I1,I2,I3,I4) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         && bitvectorActiveAtomMask->testBit(I3/3) \
         && bitvectorActiveAtomMask->testBit(I4/3)) \
    ) goto SKIP_term;



double	_evaluateEnergyOnly_Dihedral( int I1, int I2, int I3, int I4,
                                      core::T_sp activeAtomMask,
                                     num_real x1, num_real y1, num_real z1,
                                     num_real x2, num_real y2, num_real z2,
                                     num_real x3, num_real y3, num_real z3,
                                     num_real x4, num_real y4, num_real z4,
                                     num_real V, num_real DN, int IN,
                                     num_real cosPhase, num_real sinPhase )
{
  // I1..I4, IN, and activeAtomMask are part of the public API for symmetry with
  // the gradient/hessian forms; they have no effect on a single-term energy
  // calculation, so suppress unused-parameter warnings.
  (void)I1; (void)I2; (void)I3; (void)I4; (void)activeAtomMask; (void)IN;
  double localPos[12] = { (double)x1, (double)y1, (double)z1,
                          (double)x2, (double)y2, (double)z2,
                          (double)x3, (double)y3, (double)z3,
                          (double)x4, (double)y4, (double)z4 };
  double energy_accum = 0.0;
  Dihedral_Component<NoHessian> dihedral;
  dihedral_term lt((double)V, (double)DN, (double)sinPhase, (double)cosPhase, 0, 3, 6, 9);
  return dihedral.energy(lt, localPos, &energy_accum);
}

CL_LAMBDA((energy-dihedral chem:energy-dihedral) pos &optional active-atom-mask);
CL_DEFMETHOD void	EnergyDihedral_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp 	pos,
                                                                                                 core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  constexpr size_t PS = Dihedral_Component<double*>::PositionSize;  // 12
  const double energyTol = 1.0e-7;
  const double forceTol  = 1.0e-5;
  const double hessTol   = 1.0e-3;
  double localPos[PS];
  double force_a[PS], force_fd[PS];
  double hess_a[PS*PS], hess_fd[PS*PS];
  double dvec_a[PS], dvec_fd[PS];
  double hdvec_a[PS], hdvec_fd[PS];
  Dihedral_Component<double*> dihedral;
  int fails = 0;
  int idx = 0;
  for ( auto di=this->_Terms.begin(); di!=this->_Terms.end(); di++, idx++ ) {
    int I1 = di->term.i3x1;
    int I2 = di->term.i3x2;
    int I3 = di->term.i3x3;
    int I4 = di->term.i3x4;
    if (hasActiveAtomMask &&
        !(bitvectorActiveAtomMask->testBit(I1/3) &&
          bitvectorActiveAtomMask->testBit(I2/3) &&
          bitvectorActiveAtomMask->testBit(I3/3) &&
          bitvectorActiveAtomMask->testBit(I4/3))) continue;
    localPos[0]  = pos->getElement(I1+0);
    localPos[1]  = pos->getElement(I1+1);
    localPos[2]  = pos->getElement(I1+2);
    localPos[3]  = pos->getElement(I2+0);
    localPos[4]  = pos->getElement(I2+1);
    localPos[5]  = pos->getElement(I2+2);
    localPos[6]  = pos->getElement(I3+0);
    localPos[7]  = pos->getElement(I3+1);
    localPos[8]  = pos->getElement(I3+2);
    localPos[9]  = pos->getElement(I4+0);
    localPos[10] = pos->getElement(I4+1);
    localPos[11] = pos->getElement(I4+2);
    for (size_t k = 0; k < PS; ++k) {
      force_a[k]=0.0; force_fd[k]=0.0;
      dvec_a[k]=0.0;  dvec_fd[k]=0.0;
      hdvec_a[k]=0.0; hdvec_fd[k]=0.0;
    }
    for (size_t k = 0; k < PS*PS; ++k) {
      hess_a[k]=0.0; hess_fd[k]=0.0;
    }
    double e_a = 0.0, e_fd = 0.0;
    dihedral_term lt(di->term.v, di->term.n, di->term.sinphase, di->term.cosphase, 0, 3, 6, 9);
    dihedral.hessian(   lt, localPos, &e_a,  force_a,  hess_a,  dvec_a,  hdvec_a);
    dihedral.hessian_fd(lt, localPos, &e_fd, force_fd, hess_fd, dvec_fd, hdvec_fd);
    bool termFails = false;
    if (std::fabs(e_a - e_fd) > energyTol) {
      termFails = true;
      core::lisp_write(fmt::format(
        "  energy mismatch term {}: ana={} fd={} diff={}\n",
        idx, e_a, e_fd, e_a - e_fd));
    }
    for (size_t k = 0; k < PS; ++k) {
      if (std::fabs(force_a[k] - force_fd[k]) > forceTol) {
        termFails = true;
        core::lisp_write(fmt::format(
          "  force mismatch term {} dof {}: ana={} fd={} diff={}\n",
          idx, k, force_a[k], force_fd[k], force_a[k] - force_fd[k]));
      }
    }
    for (size_t r = 0; r < PS; ++r) {
      for (size_t c = 0; c < PS; ++c) {
        double a = hess_a[r*PS+c];
        double f = hess_fd[r*PS+c];
        if (std::fabs(a - f) > hessTol) {
          termFails = true;
          core::lisp_write(fmt::format(
            "  hessian mismatch term {} ({},{}): ana={} fd={} diff={}\n",
            idx, r, c, a, f, a - f));
        }
      }
    }
    if (termFails) {
      fails++;
      core::lisp_write(fmt::format(
        "  -> term {} I1={} I2={} I3={} I4={} V={} DN={}\n",
        idx, I1, I2, I3, I4, di->term.v, di->term.n));
    }
  }
  core::lisp_write(fmt::format(
    "compareAnalyticalAndNumericalForceAndHessianTermByTerm: {} fails out of {} terms\n",
    fails, idx));
}



#define USE_DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4) \
    if (doDebugInteractions) { \
      core::eval::funcall(debugInteractions,dihedral_type(IN), \
                          mk_double_float(Energy), \
                          core::make_fixnum(I1), core::make_fixnum(I2), core::make_fixnum(I3), core::make_fixnum(I4)); \
    }

#define IGNORE_DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4) {}

core::List_sp EnergyDihedral::encode() const {
  ql::list result;
  result << core::Cons_O::create(INTERN_(kw,v), core::clasp_make_double_float(this->term.v))
         << core::Cons_O::create(INTERN_(kw,in), core::make_fixnum((int)this->term.n))
         << core::Cons_O::create(INTERN_(kw,phase), core::clasp_make_double_float(this->_PhaseRad))
         << core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.i3x1))
         << core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->term.i3x2))
         << core::Cons_O::create(INTERN_(kw,i3), core::make_fixnum(this->term.i3x3))
         << core::Cons_O::create(INTERN_(kw,i4), core::make_fixnum(this->term.i3x4))
         << core::Cons_O::create(INTERN_(kw,proper), _lisp->_boolean(this->_Proper))
         << core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1)
         << core::Cons_O::create(INTERN_(kw,atom2), this->_Atom2)
         << core::Cons_O::create(INTERN_(kw,atom3), this->_Atom3)
         << core::Cons_O::create(INTERN_(kw,atom4), this->_Atom4)
    ;
  return result.cons();
}

void EnergyDihedral::decode(core::List_sp alist) {
  SIMPLE_ERROR("Implement decode of EnergyDihedral");
}

#ifdef XML_ARCHIVE
void	EnergyDihedral::archive(core::ArchiveP node)
{
//    node->attribute("_Type1",this->_Type1);
//    node->attribute("_Type2",this->_Type2);
//    node->attribute("_Type3",this->_Type3);
//    node->attribute("_Type4",this->_Type4);
  node->attribute("_Proper",this->_Proper);
  node->attribute("_PhaseRad",this->_PhaseRad);
  node->attribute("sinPhase",this->term.sinphase);
  node->attribute("cosPhase",this->term.cosphase);
  node->attribute("V",this->term.v);
  node->attribute("DN",this->term.n);
  node->attribute("IN",(int)this->term.n);
  node->attribute("I1",this->term.i3x1);
  node->attribute("I2",this->term.i3x2);
  node->attribute("I3",this->term.i3x3);
  node->attribute("I4",this->term.i3x4);
  node->attribute("a1",this->_Atom1);
  node->attribute("a2",this->_Atom2);
  node->attribute("a3",this->_Atom3);
  node->attribute("a4",this->_Atom4);
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
  node->attributeIfDefined("calcForce",this->_calcForce,this->_calcForce);
  node->attributeIfDefined("calcDiagonalHessian",this->_calcDiagonalHessian,this->_calcDiagonalHessian);
  node->attributeIfDefined("calcOffDiagonalHessian",this->_calcOffDiagonalHessian,this->_calcOffDiagonalHessian);
#include <cando/chem/energy_functions/_Dihedral_debugEvalSerialize.cc>
#endif //]
}
#endif


void EnergyDihedral::defineFrom( int n, FFPtor_sp ffterm , EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4, double scale)
{
//    this->_Term = term;
  this->_Proper = true;
//    this->_Type1 = term->_T1;
//    this->_Type2 = term->_T2;
//    this->_Type3 = term->_T3;
//    this->_Type4 = term->_T4;
  this->_Atom1 = ea1->atom();
  this->_Atom2 = ea2->atom();
  this->_Atom3 = ea3->atom();
  this->_Atom4 = ea4->atom();
  this->_PhaseRad = ffterm->getPhaseRadians(n);
  this->term.sinphase = sin(this->_PhaseRad);
  this->term.cosphase = cos(this->_PhaseRad);
  this->term.v = ffterm->getV_kcal(n)*scale;
  this->term.n = n;
  this->term.i3x1 = ea1->coordinateIndexTimes3();
  this->term.i3x2 = ea2->coordinateIndexTimes3();
  this->term.i3x3 = ea3->coordinateIndexTimes3();
  this->term.i3x4 = ea4->coordinateIndexTimes3();
}

double	EnergyDihedral::getDihedral()
{
  Vector3	pos1, pos2, pos3, pos4;
  num_real	angle;
  pos1 = this->_Atom1->getPosition();
  pos2 = this->_Atom2->getPosition();
  pos3 = this->_Atom3->getPosition();
  pos4 = this->_Atom4->getPosition();
  angle = geom::calculateDihedral( pos1, pos2, pos3, pos4);
  return angle;
}

double	EnergyDihedral::getDihedralDeviation()
{
  Vector3	pos1, pos2, pos3, pos4;
  num_real	phi, dev;
  pos1 = this->_Atom1->getPosition();
  pos2 = this->_Atom2->getPosition();
  pos3 = this->_Atom3->getPosition();
  pos4 = this->_Atom4->getPosition();
  phi = geom::calculateDihedral( pos1, pos2, pos3, pos4);
  dev = 1.0+cos(this->term.n*phi-this->_PhaseRad);
  return dev;
}


void EnergyDihedral::defineFrom( int n, FFItor_sp term , EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4, double scale)
{
//    this->_Term = term;
  this->_Proper = false;
//    this->_Type1 = term->_T1;
//    this->_Type2 = term->_T2;
//    this->_Type3 = term->_T3;
//    this->_Type4 = term->_T4;
  this->_Atom1 = ea1->atom();
  this->_Atom2 = ea2->atom();
  this->_Atom3 = ea3->atom();
  this->_Atom4 = ea4->atom();
  this->term.n = n;
  this->term.v = term->getV_kcal(n)*scale;
  this->_PhaseRad = term->getPhaseRadians(n);
  this->term.cosphase = cos(this->_PhaseRad);
  this->term.sinphase = sin(this->_PhaseRad);
  this->term.i3x1 = ea1->coordinateIndexTimes3();
  this->term.i3x2 = ea2->coordinateIndexTimes3();
  this->term.i3x3 = ea3->coordinateIndexTimes3();
  this->term.i3x4 = ea4->coordinateIndexTimes3();
}

#if 0 //[
void	EnergyDihedral::defineMissingProper( EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, EnergyAtom *ea4 )
{
  this->_Proper = true;
//    this->_Type1 = ea1->_Atom->getType(atomTypes);
//    this->_Type2 = ea2->_Atom->getType(atomTypes);
  /    this->_Type3 = ea3->_Atom->getType(atomTypes);
  this->_Type4 = ea4->_Atom->getType(atomTypes);
}
#endif //]


#if 0
adapt::QDomNode_sp	EnergyDihedral::asXml()
{
  adapt::QDomNode_sp	node;
  Vector3	vdiff;

  node = adapt::QDomNode_O::create(env,"EnergyDihedral");
  node->addAttributeBool("proper",this->_Proper );
  node->addAttributeString("atom1Name",this->_Atom1->getName());
  node->addAttributeString("atom2Name",this->_Atom2->getName());
  node->addAttributeString("atom3Name",this->_Atom3->getName());
  node->addAttributeString("atom4Name",this->_Atom4->getName());
//    node->addAttributeString("atom1Type",this->_Type1 );
//    node->addAttributeString("atom2Type",this->_Type2 );
//    node->addAttributeString("atom3Type",this->_Type3 );
//    node->addAttributeString("atom4Type",this->_Type4 );
  node->addAttributeInt("I1",this->term.i3x1);
  node->addAttributeInt("I2",this->term.i3x2);
  node->addAttributeInt("I3",this->term.i3x3);
  node->addAttributeInt("I4",this->term.i3x4);
  node->addAttributeDoubleScientific("DN",this->term.n);
  node->addAttributeInt("IN",(int)this->term.n);
  node->addAttributeDoubleScientific("V",this->term.v);
  node->addAttributeDoubleScientific("PhaseRad",this->_PhaseRad);
#if TURN_ENERGY_FUNCTION_DEBUG_ON
  adapt::QDomNode_sp xml = adapt::QDomNode_O::create(env,"Evaluated");
  xml->addAttributeBool("calcForce",this->_calcForce );
  xml->addAttributeBool("calcDiagonalHessian",this->_calcDiagonalHessian );
  xml->addAttributeBool("calcOffDiagonalHessian",this->_calcOffDiagonalHessian );
#include <_Dihedral_debugEvalXml.cc>
  node->addChild(xml);
#endif
  return node;
}

void	EnergyDihedral::parseFromXmlUsingAtomTable(adapt::QDomNode_sp	xml,
                                                   AtomTable_sp at )
{
  this->term.i3x1 = xml->getAttributeInt("I1");
  this->term.i3x2 = xml->getAttributeInt("I2");
  this->term.i3x3 = xml->getAttributeInt("I3");
  this->term.i3x4 = xml->getAttributeInt("I4");
  this->_Atom1 = at->findEnergyAtomWithCoordinateIndex(this->term.i3x1)->atom();
  this->_Atom2 = at->findEnergyAtomWithCoordinateIndex(this->term.i3x2)->atom();
  this->_Atom3 = at->findEnergyAtomWithCoordinateIndex(this->term.i3x3)->atom();
  this->_Atom4 = at->findEnergyAtomWithCoordinateIndex(this->term.i3x4)->atom();
  this->term.n = xml->getAttributeDouble("DN");
  this->term.v = xml->getAttributeDouble("V");
  this->_PhaseRad = xml->getAttributeDouble("PhaseRad");
  this->_Proper = xml->getAttributeBool("proper");
}
#endif


void EnergyDihedral_O::addTerm(const EnergyDihedral& term)
{
  if (this->_Terms.size() == this->_Terms.capacity()-1) {
    this->_Terms.reserve(this->_Terms.size()*2);
  }
  this->_Terms.push_back(term);
}


void	EnergyDihedral_O::dumpTerms(core::HashTable_sp atomTypes)
{
  gctools::Vec0<EnergyDihedral>::iterator	edi;
  string				as1,as2,as3,as4;
  string				str1, str2, str3, str4, type;
  for (edi=this->_Terms.begin();edi!=this->_Terms.end();edi++)
  {
    if ( edi->term.v == 0.0 ) continue;
    if ( edi->_Proper ) {
      type = "1PROPER  ";
    } else {
      type = "9IMPROPER";
    }
    as1 = atomLabel(edi->_Atom1);
    as2 = atomLabel(edi->_Atom2);
    as3 = atomLabel(edi->_Atom3);
    as4 = atomLabel(edi->_Atom4);
    if ( edi->_Proper ) {
      if ( as1 < as4 ) {
        str1 = as1;
        str2 = as2;
        str3 = as3;
        str4 = as4;
      } else {
        str4 = as1;
        str3 = as2;
        str2 = as3;
        str1 = as4;
      }
    } else {
      stringstream ss;
      ss.str("");
      ss << as3 << "#Cr";
      str1 = ss.str();
      ss.str("");
      ss << as1 << "#A1";
      str2 = ss.str();
      ss.str("");
      ss << as2 << "#A2";
      str3 = ss.str();
      ss.str("");
      ss << as4 << "#A4";
      str4 = ss.str();
    }
    core::clasp_write_string(fmt::format("TERM 3DIH {} {:<9} - {:<9} - {:<9} - {:<9} {:8.2f} {:8.2f} {:2.0f}\n" , type , str1 , str2 , str3 , str4 , edi->term.v , edi->_PhaseRad , edi->term.n));
  }
}


CL_DEFMETHOD
core::List_sp	EnergyDihedral_O::lookupDihedralTerms(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, core::HashTable_sp atomTypes )
{
  gctools::Vec0<EnergyDihedral>::iterator	edi;
  string				as1,as2,as3,as4;
  string				str1, str2, str3, str4, type;
  ql::list  result;
  core::T_sp tia1 = atomTable->_AtomTableIndexes->gethash(a1);
  core::T_sp tia2 = atomTable->_AtomTableIndexes->gethash(a2);
  core::T_sp tia3 = atomTable->_AtomTableIndexes->gethash(a3);
  core::T_sp tia4 = atomTable->_AtomTableIndexes->gethash(a4);
  if (!tia1.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a1));
  if (!tia2.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a2));
  if (!tia3.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a3));
  if (!tia4.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a4));
  for (edi=this->_Terms.begin();edi!=this->_Terms.end();edi++) {
    if ((edi->_Atom1==a1 &&
         edi->_Atom2==a2 &&
         edi->_Atom3==a3 &&
         edi->_Atom4==a4)
        || (edi->_Atom4==a1 &&
            edi->_Atom3==a2 &&
            edi->_Atom2==a3 &&
            edi->_Atom1==a4)) {
      ql::list oneResult;
      oneResult << INTERN_(kw,type1) << edi->_Atom1->getType(atomTypes)
                << INTERN_(kw,type2) << edi->_Atom2->getType(atomTypes)
                << INTERN_(kw,type3) << edi->_Atom3->getType(atomTypes)
                << INTERN_(kw,type4) << edi->_Atom4->getType(atomTypes)
                << INTERN_(kw,proper) << _lisp->_boolean(edi->_Proper)
                << INTERN_(kw,multiplicity) << core::clasp_make_fixnum((int)edi->term.n)
                << INTERN_(kw,phase_degrees) << core::clasp_make_double_float(edi->_PhaseRad)
                << INTERN_(kw,v) << core::clasp_make_double_float(edi->term.v);
      result << oneResult.result();
    }
  }
  return result.result();
}



void	EnergyDihedral_O::setupHessianPreconditioner(
                                                     NVector_sp nvPosition,
                                                     AbstractLargeSquareMatrix_sp m,
                                                     core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  Dihedral_Component<double*> dihedral;
  constexpr size_t PS = Dihedral_Component<double*>::PositionSize;  // 12
  double localPos[PS];
  double localForce[PS];
  double localHess[PS*PS];
  double localDvec[PS];
  double localHdvec[PS];
  for ( gctools::Vec0<EnergyDihedral>::iterator di=this->_Terms.begin();
        di!=this->_Terms.end(); di++ ) {
    int I1 = di->term.i3x1;
    int I2 = di->term.i3x2;
    int I3 = di->term.i3x3;
    int I4 = di->term.i3x4;
    if (hasActiveAtomMask &&
        !(bitvectorActiveAtomMask->testBit(I1/3) &&
          bitvectorActiveAtomMask->testBit(I2/3) &&
          bitvectorActiveAtomMask->testBit(I3/3) &&
          bitvectorActiveAtomMask->testBit(I4/3))) continue;
    localPos[0]  = (*nvPosition)[I1+0];
    localPos[1]  = (*nvPosition)[I1+1];
    localPos[2]  = (*nvPosition)[I1+2];
    localPos[3]  = (*nvPosition)[I2+0];
    localPos[4]  = (*nvPosition)[I2+1];
    localPos[5]  = (*nvPosition)[I2+2];
    localPos[6]  = (*nvPosition)[I3+0];
    localPos[7]  = (*nvPosition)[I3+1];
    localPos[8]  = (*nvPosition)[I3+2];
    localPos[9]  = (*nvPosition)[I4+0];
    localPos[10] = (*nvPosition)[I4+1];
    localPos[11] = (*nvPosition)[I4+2];
    // Pack the four atoms' 3 coords into a contiguous 12-vector at offsets 0, 3, 6 and 9.
    // Zero local accumulator buffers (kernel uses += into them).
    for (size_t k = 0; k < PS; ++k) {
      localForce[k] = 0.0; localDvec[k] = 0.0; localHdvec[k] = 0.0;
    }
    for (size_t k = 0; k < PS*PS; ++k) localHess[k] = 0.0;
    double energyAccum = 0.0;
    dihedral_term lt(di->term.v, di->term.n, di->term.sinphase, di->term.cosphase, 0, 3, 6, 9);
    dihedral.hessian(lt, localPos, &energyAccum, localForce, localHess,
                     localDvec, localHdvec);
    // Scatter the 12x12 local Hessian into the global sparse matrix.
    // Local offsets 0..2 map to global I1+0..I1+2, 3..5 to I2, 6..8 to I3 and 9..11 to I4.
    // The kernel writes BOTH (r,c) and (c,r) for off-diagonal entries
    // (energyComponent.h:269-270), so accumulate only the lower triangle here.
    int globalIdx[PS] = { I1+0, I1+1, I1+2, I2+0, I2+1, I2+2,
                          I3+0, I3+1, I3+2, I4+0, I4+1, I4+2 };
    for (size_t r = 0; r < PS; ++r) {
      for (size_t c = 0; c <= r; ++c) {
        double v = localHess[r*PS + c];
        if (v != 0.0) {
          m->addToElement(globalIdx[r], globalIdx[c], v);
        }
      }
    }
  }
}
double	EnergyDihedral_O::evaluateAllComponentSingle(
    gctools::Vec0<EnergyDihedral>::iterator di_start,
    gctools::Vec0<EnergyDihedral>::iterator di_end,
    ScoringFunction_sp score,
    NVector_sp 	pos,
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
#undef DIHEDRAL_APPLY_ATOM_MASK
#define DIHEDRAL_APPLY_ATOM_MASK(I1,I2,I3,I4) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         && bitvectorActiveAtomMask->testBit(I3/3) \
         && bitvectorActiveAtomMask->testBit(I4/3)) \
    ) goto SKIP_term_and_angle_test;
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  double totalEnergy = 0.0;

  auto evalType = determineEnergyComponentEvalType(force,hdvec,dvec);

  double termEnergy = 0.0;

  int	I1, I2, I3, I4;
  num_real sinPhase, cosPhase, SinNPhi, CosNPhi;
  gctools::Vec0<EnergyDihedral>::iterator di;

  int i = 0;

#define DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4) USE_DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4)

  DOUBLE* position = &(*pos)[0];
  DOUBLE* rforce = NULL;
  DOUBLE* rhessian = NULL; // &(*hessian)[0];
  DOUBLE* rdvec = NULL;
  DOUBLE* rhdvec = NULL;
  double  Energy;
  Dihedral_Component<NoHessian> dihedral;
#define KERNEL_DIHEDRAL_APPLY_ATOM_MASK(I1,I2,I3,I4) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         && bitvectorActiveAtomMask->testBit(I3/3) \
         && bitvectorActiveAtomMask->testBit(I4/3)) \
    ) continue;
  if (evalType==energyEval) {
    for ( di=di_start; di!= di_end; di++ ) {
      KERNEL_DIHEDRAL_APPLY_ATOM_MASK(di->term.i3x1, di->term.i3x2, di->term.i3x3, di->term.i3x4);
      Energy = dihedral.energy(di->term, position, &termEnergy);
      const INT IN = (int)di->term.n;
      DIHEDRAL_DEBUG_INTERACTIONS(di->term.i3x1, di->term.i3x2, di->term.i3x3, di->term.i3x4);
    }
  } else if (evalType==gradientEval) {
    rforce = &(*force)[0];
    for ( di=di_start; di!= di_end; di++ ) {
      KERNEL_DIHEDRAL_APPLY_ATOM_MASK(di->term.i3x1, di->term.i3x2, di->term.i3x3, di->term.i3x4);
      Energy = dihedral.gradient(di->term, position, &termEnergy, rforce);
      const INT IN = (int)di->term.n;
      DIHEDRAL_DEBUG_INTERACTIONS(di->term.i3x1, di->term.i3x2, di->term.i3x3, di->term.i3x4);
    }
  } else {
      rforce = &(*force)[0];
      rdvec = &(*dvec)[0];
      rhdvec = &(*hdvec)[0];
    for ( di=di_start; di!= di_end; di++ ) {
      KERNEL_DIHEDRAL_APPLY_ATOM_MASK(di->term.i3x1, di->term.i3x2, di->term.i3x3, di->term.i3x4);
      Energy = dihedral.hessian(di->term, position, &termEnergy, rforce, NoHessian(), rdvec, rhdvec);
      const INT IN = (int)di->term.n;
      DIHEDRAL_DEBUG_INTERACTIONS(di->term.i3x1, di->term.i3x2, di->term.i3x3, di->term.i3x4);
    }
  }

#undef DIHEDRAL_DEBUG_INTERACTIONS

  // Near-linear (degenerate) torsions are handled inside the kernel by the
  // linear_dihedral_clamp (see dihedral_fast.c / energy-kernels.lisp), which
  // floors t1^2+t2^2 so the derivatives stay finite.  No separate error path.

  SKIP_term_and_angle_test: (void)0;

  return termEnergy;
}







void EnergyDihedral_O::initialize()
{
  this->Base::initialize();
}

void EnergyDihedral_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}

SYMBOL_EXPORT_SC_(KeywordPkg,v);
SYMBOL_EXPORT_SC_(KeywordPkg,in);
SYMBOL_EXPORT_SC_(KeywordPkg,dn);
SYMBOL_EXPORT_SC_(KeywordPkg,phase);
SYMBOL_EXPORT_SC_(KeywordPkg,i1);
SYMBOL_EXPORT_SC_(KeywordPkg,i2);
SYMBOL_EXPORT_SC_(KeywordPkg,i3);
SYMBOL_EXPORT_SC_(KeywordPkg,i4);
SYMBOL_EXPORT_SC_(KeywordPkg,proper);
SYMBOL_EXPORT_SC_(KeywordPkg,atom1);
SYMBOL_EXPORT_SC_(KeywordPkg,atom2);
SYMBOL_EXPORT_SC_(KeywordPkg,atom3);
SYMBOL_EXPORT_SC_(KeywordPkg,atom4);
SYMBOL_EXPORT_SC_(KeywordPkg,coordinates);
SYMBOL_EXPORT_SC_(KeywordPkg,indexes);
SYMBOL_EXPORT_SC_(KeywordPkg,simd_index);

CL_DEFMETHOD core::List_sp EnergyDihedral_O::extract_vectors_as_alist() const{
  size_t size = this->_Terms.size();
  NVector_sp v_vec = NVector_O::make(size);
  core::SimpleVector_int32_t_sp in_vec = core::SimpleVector_int32_t_O::make(size);
  NVector_sp phase_vec = NVector_O::make(size);
  core::SimpleVector_int32_t_sp i1_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i2_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i3_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i4_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_sp proper_vec = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom1_vec = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom2_vec = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom3_vec = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom4_vec = core::SimpleVector_O::make(size);
  for (size_t i=0;i<size;++i){
    const EnergyDihedral& entry = this->_Terms[i];
    (*v_vec)[i] = entry.term.v;
    (*in_vec)[i] = (int)entry.term.n;
    (*phase_vec)[i] = entry._PhaseRad;
    (*i1_vec)[i] = entry.term.i3x1;
    (*i2_vec)[i] = entry.term.i3x2;
    (*i3_vec)[i] = entry.term.i3x3;
    (*i4_vec)[i] = entry.term.i3x4;
    (*proper_vec)[i] = _lisp->_boolean(entry._Proper);
    (*atom1_vec)[i] = entry._Atom1;
    (*atom2_vec)[i] = entry._Atom2;
    (*atom3_vec)[i] = entry._Atom3;
    (*atom4_vec)[i] = entry._Atom4;
  }
  ql::list result;
  result
    << core::Cons_O::create(kw::_sym_v, v_vec)
    << core::Cons_O::create(kw::_sym_in, in_vec)
    << core::Cons_O::create(kw::_sym_phase, phase_vec)
    << core::Cons_O::create(kw::_sym_i1, i1_vec)
    << core::Cons_O::create(kw::_sym_i2, i2_vec)
    << core::Cons_O::create(kw::_sym_i3, i3_vec)
    << core::Cons_O::create(kw::_sym_i4, i4_vec)
    << core::Cons_O::create(kw::_sym_proper, proper_vec)
    << core::Cons_O::create(kw::_sym_atom1, atom1_vec)
    << core::Cons_O::create(kw::_sym_atom2, atom2_vec)
    << core::Cons_O::create(kw::_sym_atom3, atom3_vec)
    << core::Cons_O::create(kw::_sym_atom4, atom4_vec);
  return result.cons();
}

CL_DEFMETHOD void EnergyDihedral_O::fill_from_vectors_in_alist(core::List_sp vectors)
{
  NVector_sp v_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_v));
  NVector_sp dn_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_dn));
  NVector_sp phase_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_phase));
  core::SimpleVector_int32_t_sp i1_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i1));
  core::SimpleVector_int32_t_sp i2_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i2));
  core::SimpleVector_int32_t_sp i3_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i3));
  core::SimpleVector_int32_t_sp i4_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i4));
  core::SimpleVector_sp proper_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_proper));
  core::SimpleVector_sp atom1_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom1));
  core::SimpleVector_sp atom2_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom2));
  core::SimpleVector_sp atom3_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom3));
  core::SimpleVector_sp atom4_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom4));
  this->_Terms.resize(v_vec->length());
  for (size_t i=0, iEnd(v_vec->length());i<iEnd;++i){
    EnergyDihedral& entry = this->_Terms[i];
    entry.term.v = (*v_vec)[i];
    entry.term.n = (*dn_vec)[i];
    entry._PhaseRad = (*phase_vec)[i];
    entry.term.cosphase = cos((*phase_vec)[i]);
    entry.term.sinphase = sin((*phase_vec)[i]);
    entry.term.i3x1 = (*i1_vec)[i];
    entry.term.i3x2 = (*i2_vec)[i];
    entry.term.i3x3 = (*i3_vec)[i];
    entry.term.i3x4 = (*i4_vec)[i];
    entry._Proper =  ((*proper_vec)[i]).isTrue();
    entry._Atom1 = gc::As_unsafe<Atom_sp>((*atom1_vec)[i]);
    entry._Atom2 = gc::As_unsafe<Atom_sp>((*atom2_vec)[i]);
    entry._Atom3 = gc::As_unsafe<Atom_sp>((*atom3_vec)[i]);
    entry._Atom4 = gc::As_unsafe<Atom_sp>((*atom4_vec)[i]);
  }
}

CL_DEFMETHOD void EnergyDihedral_O::addDihedralTerm(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4,
                                                    double phase, bool proper, double v, int multiplicity)
{
  EnergyAtom* ea1 = atomTable->getEnergyAtomPointer(a1);
  EnergyAtom* ea2 = atomTable->getEnergyAtomPointer(a2);
  EnergyAtom* ea3 = atomTable->getEnergyAtomPointer(a3);
  EnergyAtom* ea4 = atomTable->getEnergyAtomPointer(a4);
  EnergyDihedral energyDihedral(a1,a2,a3,a4,phase, proper, 
                                ea1->coordinateIndexTimes3(),
                                ea2->coordinateIndexTimes3(),
                                ea3->coordinateIndexTimes3(),
                                ea4->coordinateIndexTimes3(),
                                std::sin(phase),
                                std::cos(phase),
                                v,
                                1.0*multiplicity,
                                multiplicity);
  this->addTerm(energyDihedral);
}

EnergyComponent_sp EnergyDihedral_O::copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc) {
  core::T_sp keepInteraction = specializeKeepInteractionFactory(keepInteractionFactory,EnergyDihedral_O::staticClass());
  EnergyDihedral_sp copy = EnergyDihedral_O::create();
  copyEnergyComponent( copy, this->asSmartPtr() );
  for ( auto edi=this->_Terms.begin(); edi!=this->_Terms.end(); edi++ ) {
    Atom_sp a1 = edi->_Atom1;
    Atom_sp a2 = edi->_Atom2;
    Atom_sp a3 = edi->_Atom3;
    Atom_sp a4 = edi->_Atom4;
    if ( skipInteraction( keepInteraction, a1, a2, a3, a4 ) ) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}


#ifndef _TARGET_OS_DARWIN
#define VREAL8_WIDTH 8
typedef double real;
typedef real real8 __attribute__((vector_size(8*VREAL8_WIDTH))) __attribute__((aligned(4)));
typedef int64_t int8 __attribute__((vector_size(8*VREAL8_WIDTH))) __attribute__((aligned(4)));

double EnergyDihedral_O::evaluateAllComponentSimd8(
    gctools::Vec0<EnergyDihedral>::iterator di_start8,
    gctools::Vec0<EnergyDihedral>::iterator di_end8,
    ScoringFunction_sp          score,
    NVector_sp 	                pos,
    bool 		        calcForce,
    gc::Nilable<NVector_sp> 	force,
    bool		        calcDiagonalHessian,
    bool		        calcOffDiagonalHessian,
    gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
    gc::Nilable<NVector_sp>	hdvec,
    gc::Nilable<NVector_sp>     dvec,
    core::T_sp activeAtomMask )
{
  IMPLEMENT_ME();
}

//
//
//  -------------------------------------------------------------
//
//  Simd4 function
//
//

#define VREAL4_WIDTH 4
typedef double real;
typedef real real4 __attribute__((vector_size(8*VREAL4_WIDTH)));
typedef int64_t int4 __attribute__((vector_size(8*VREAL4_WIDTH)));

double EnergyDihedral_O::evaluateAllComponentSimd4(
    gctools::Vec0<EnergyDihedral>::iterator di_start4,
    gctools::Vec0<EnergyDihedral>::iterator di_end4,
    ScoringFunction_sp          score,
    NVector_sp 	                pos,
    bool 		        calcForce,
    gc::Nilable<NVector_sp> 	force,
    bool		        calcDiagonalHessian,
    bool		        calcOffDiagonalHessian,
    gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
    gc::Nilable<NVector_sp>	hdvec,
    gc::Nilable<NVector_sp>     dvec,
    core::T_sp activeAtomMask )
{
  IMPLEMENT_ME();
}

//
//
//  -------------------------------------------------------------
//
//  Simd2 function
//
//

#define VREAL2_WIDTH 2
typedef double real;
typedef real real2 __attribute__((vector_size(8*VREAL2_WIDTH)));
typedef int64_t int2 __attribute__((vector_size(8*VREAL2_WIDTH)));

double EnergyDihedral_O::evaluateAllComponentSimd2(
    gctools::Vec0<EnergyDihedral>::iterator di_start2,
    gctools::Vec0<EnergyDihedral>::iterator di_end2,
    ScoringFunction_sp          score,
    NVector_sp 	                pos,
    bool 		        calcForce,
    gc::Nilable<NVector_sp> 	force,
    bool		        calcDiagonalHessian,
    bool		        calcOffDiagonalHessian,
    gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
    gc::Nilable<NVector_sp>	hdvec,
    gc::Nilable<NVector_sp>     dvec,
    core::T_sp activeAtomMask )
{
  IMPLEMENT_ME();
}
#endif // !_TARGET_OS_DARWIN

double EnergyDihedral_O::evaluateAllComponent(ScoringFunction_sp          score,
                                              NVector_sp 	                pos,
                                              core::T_sp energyScale,
                                              core::T_sp               energyComponents,
                                              bool 		        calcForce,
                                              gc::Nilable<NVector_sp> 	force,
                                              bool		        calcDiagonalHessian,
                                              bool		        calcOffDiagonalHessian,
                                              gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                              gc::Nilable<NVector_sp>	hdvec,
                                              gc::Nilable<NVector_sp>     dvec,
                                              core::T_sp activeAtomMask,
                                              core::T_sp debugInteractions )
{
  num_real  energy = 0.0;
  this->_Evaluations++;
#ifdef _TARGET_OS_DARWIN
    energy += this->evaluateAllComponentSingle(this->_Terms.begin(),
                                               this->_Terms.end(),
                                               score,
                                               pos,
                                               calcForce,
                                               force,
                                               calcDiagonalHessian,
                                               calcOffDiagonalHessian,
                                               hessian,
                                               hdvec,
                                               dvec,
                                               activeAtomMask,
                                               debugInteractions );
#else
  if (cando::global_simd_width ==1 || debugInteractions.notnilp()) {
    energy += this->evaluateAllComponentSingle(this->_Terms.begin(),
                                               this->_Terms.end(),
                                               score,
                                               pos,
                                               calcForce,
                                               force,
                                               calcDiagonalHessian,
                                               calcOffDiagonalHessian,
                                               hessian,
                                               hdvec,
                                               dvec,
                                               activeAtomMask,
                                               debugInteractions );
  } else if (cando::global_simd_width == 8 ) {
    gctools::Vec0<EnergyDihedral>::iterator di_vector_end8 = this->_Terms.begin()+((int)(std::distance(this->_Terms.begin(),this->_Terms.end())/VREAL8_WIDTH))*VREAL8_WIDTH;
    energy += this->evaluateAllComponentSimd8(this->_Terms.begin(),
                                              di_vector_end8,
                                              score,
                                              pos,
                                              calcForce,
                                              force,
                                              calcDiagonalHessian,
                                              calcOffDiagonalHessian,
                                              hessian,
                                              hdvec,
                                              dvec,
                                              activeAtomMask );
    energy += this->evaluateAllComponentSingle(di_vector_end8,
                                               this->_Terms.end(),
                                               score,
                                               pos,
                                               calcForce,
                                               force,
                                               calcDiagonalHessian,
                                               calcOffDiagonalHessian,
                                               hessian,
                                               hdvec,
                                               dvec,
                                               activeAtomMask,
                                               nil<core::T_O>());
  } else if (cando::global_simd_width == 4 ) {
    gctools::Vec0<EnergyDihedral>::iterator di_vector_end4 = this->_Terms.begin()+((int)(std::distance(this->_Terms.begin(),this->_Terms.end())/VREAL4_WIDTH))*VREAL4_WIDTH;
    energy += this->evaluateAllComponentSimd4(this->_Terms.begin(),
                                              di_vector_end4,
                                              score,
                                              pos,
                                              calcForce,
                                              force,
                                              calcDiagonalHessian,
                                              calcOffDiagonalHessian,
                                              hessian,
                                              hdvec,
                                              dvec,
                                              activeAtomMask );
    energy += this->evaluateAllComponentSingle(di_vector_end4,
                                               this->_Terms.end(),
                                               score,
                                               pos,
                                               calcForce,
                                               force,
                                               calcDiagonalHessian,
                                               calcOffDiagonalHessian,
                                               hessian,
                                               hdvec,
                                               dvec,
                                               activeAtomMask,
                                               nil<core::T_O>());
  } else if (cando::global_simd_width == 2 ) {
    gctools::Vec0<EnergyDihedral>::iterator di_vector_end2 = this->_Terms.begin()+((int)(std::distance(this->_Terms.begin(),this->_Terms.end())/VREAL2_WIDTH))*VREAL2_WIDTH;
    energy += this->evaluateAllComponentSimd2(this->_Terms.begin(),
                                              di_vector_end2,
                                              score,
                                              pos,
                                              calcForce,
                                              force,
                                              calcDiagonalHessian,
                                              calcOffDiagonalHessian,
                                              hessian,
                                              hdvec,
                                              dvec,
                                              activeAtomMask );
    energy += this->evaluateAllComponentSingle(di_vector_end2,
                                               this->_Terms.end(),
                                               score,
                                               pos,
                                               calcForce,
                                               force,
                                               calcDiagonalHessian,
                                               calcOffDiagonalHessian,
                                               hessian,
                                               hdvec,
                                               dvec,
                                               activeAtomMask,
                                               nil<core::T_O>() );
  }
#endif
  maybeSetEnergy( energyComponents, EnergyDihedral_O::static_classSymbol(), energy );
  return energy;
};


void EnergyDihedral_O::emitTestCalls(core::T_sp stream, chem::NVector_sp pos) const
{
  SIMPLE_ERROR("Update EnergyDihedral_O::emitTestCalls");
  #if 0
  for ( auto si=this->_Terms.begin();
        si!=this->_Terms.end(); si++ ) {
    core::print(fmt::format("test_dihedral( errs, {}, {}, {}, {}, {}, {}, {}, {}, position_size, position, force_new, force_old, hessian_new, hessian_old, dvec_new, dvec_old, hdvec_new, hdvec_old );\n",
                            si->term.v, si->term.n, si->term.sinphase, si->term.cosphase, si->term.i3x1, si->term.i3x2, si->term.i3x3, si->term.i3x4 ),stream);
  }
  #endif
}



};
