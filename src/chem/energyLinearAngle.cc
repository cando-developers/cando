/*
    File: energyLinearAngle.cc
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

#include <clasp/core/foundation.h>
#include <clasp/core/numerics.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/ql.h>
#include <cando/chem/energyLinearAngle.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/bond.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>

namespace chem {

CL_LISPIFY_NAME(make-energy-linear-angle);
CL_DEF_CLASS_METHOD
EnergyLinearAngle_sp EnergyLinearAngle_O::make(EnergyFunction_sp energyFunction) {
  auto component = ensureComponent<EnergyLinearAngle_O>(energyFunction);
  return component;
}

};


namespace chem {
// ----------------------------------------------------------------------
// include the generated kernel: defines Linear_Angle<HESSIAN>.
// The linear-bend energy E = kt*(1+cos theta) is written purely in cos theta,
// so the kernel is smooth at linearity -- no acos, no clamp, and no
// LinearAngleError machinery (unlike energyAngle.cc).
// ----------------------------------------------------------------------

#include "cando/chem/energyKernels/linear_angle.c"

};


namespace chem {

double	_evaluateEnergyOnly_LinearAngle(
    int I1, int I2, int I3, core::T_sp activeAtomMask,
		num_real x1, num_real y1, num_real z1,
		num_real x2, num_real y2, num_real z2,
		num_real x3, num_real y3, num_real z3,
		num_real kt )
{
  // I1, I2, I3, and activeAtomMask are part of the public API for symmetry with
  // the gradient/hessian forms; they have no effect on a single-term energy
  // calculation, so suppress unused-parameter warnings.
  (void)I1; (void)I2; (void)I3; (void)activeAtomMask;
  double localPos[9] = { (double)x1, (double)y1, (double)z1,
                         (double)x2, (double)y2, (double)z2,
                         (double)x3, (double)y3, (double)z3 };
  double energy_accum = 0.0;
  Linear_Angle<NoHessian> linear_angle;
  linear_angle_term lt((double)kt, 0, 3, 6);
  return linear_angle.energy(lt, localPos, &energy_accum);
}

std::string EnergyLinearAngle_O::descriptionOfContents() const {
  stringstream ss;
  ss << ":enabled " << ((this->_Enabled) ? "T" : "NIL");
  ss << " number-of-terms " << this->_Terms.size();
  return ss.str();
}

CL_LAMBDA((energy-linear-angle chem:energy-linear-angle) pos &optional activeAtomMask);
CL_DEFMETHOD void	EnergyLinearAngle_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(chem::NVector_sp 	pos, core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  constexpr size_t PS = Linear_Angle<double*>::PositionSize;  // 9
  const double energyTol = 1.0e-7;
  const double forceTol  = 1.0e-5;
  const double hessTol   = 1.0e-3;
  double localPos[PS];
  double force_a[PS], force_fd[PS];
  double hess_a[PS*PS], hess_fd[PS*PS];
  double dvec_a[PS], dvec_fd[PS];
  double hdvec_a[PS], hdvec_fd[PS];
  Linear_Angle<double*> linear_angle;
  int fails = 0;
  int idx = 0;
  for ( auto ai=this->_Terms.begin(); ai!=this->_Terms.end(); ai++, idx++ ) {
    int I1 = ai->term.i3x1;
    int I2 = ai->term.i3x2;
    int I3 = ai->term.i3x3;
    if (hasActiveAtomMask &&
        !(bitvectorActiveAtomMask->testBit(I1/3) &&
          bitvectorActiveAtomMask->testBit(I2/3) &&
          bitvectorActiveAtomMask->testBit(I3/3))) continue;
    localPos[0] = pos->getElement(I1+0);
    localPos[1] = pos->getElement(I1+1);
    localPos[2] = pos->getElement(I1+2);
    localPos[3] = pos->getElement(I2+0);
    localPos[4] = pos->getElement(I2+1);
    localPos[5] = pos->getElement(I2+2);
    localPos[6] = pos->getElement(I3+0);
    localPos[7] = pos->getElement(I3+1);
    localPos[8] = pos->getElement(I3+2);
    for (size_t k = 0; k < PS; ++k) {
      force_a[k]=0.0; force_fd[k]=0.0;
      dvec_a[k]=0.0;  dvec_fd[k]=0.0;
      hdvec_a[k]=0.0; hdvec_fd[k]=0.0;
    }
    for (size_t k = 0; k < PS*PS; ++k) {
      hess_a[k]=0.0; hess_fd[k]=0.0;
    }
    double e_a = 0.0, e_fd = 0.0;
    linear_angle_term lt(ai->term.kt, 0, 3, 6);
    linear_angle.hessian(   lt, localPos, &e_a,  force_a,  hess_a,  dvec_a,  hdvec_a);
    linear_angle.hessian_fd(lt, localPos, &e_fd, force_fd, hess_fd, dvec_fd, hdvec_fd);
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
        "  -> term {} I1={} I2={} I3={} kt={}\n",
        idx, I1, I2, I3, ai->term.kt));
    }
  }
  core::lisp_write(fmt::format(
    "compareAnalyticalAndNumericalForceAndHessianTermByTerm: {} fails out of {} terms\n",
    fails, idx));
}




#define LINEAR_ANGLE_DEBUG_INTERACTIONS(I1,I2,I3) \
    if (doDebugInteractions) { \
      core::eval::funcall(debugInteractions,EnergyLinearAngle_O::static_classSymbol(), \
                          mk_double_float(Energy), \
                          core::make_fixnum(I1), core::make_fixnum(I2), core::make_fixnum(I3)); \
    }

core::List_sp EnergyLinearAngle::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,kt),core::clasp_make_double_float(this->term.kt)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.i3x1)),
                                  core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->term.i3x2)),
                                  core::Cons_O::create(INTERN_(kw,i3), core::make_fixnum(this->term.i3x3)),
                                  core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1),
                                  core::Cons_O::create(INTERN_(kw,atom2), this->_Atom2),
                                  core::Cons_O::create(INTERN_(kw,atom3), this->_Atom3));
}

void EnergyLinearAngle::decode(core::List_sp alist) {
  SIMPLE_ERROR("Implement decode of EnergyLinearAngle");
}

double EnergyLinearAngle::getT()
{
    Vector3	pos1, pos2, pos3;
num_real	angle;
    pos1 = this->_Atom1->getPosition();
    pos2 = this->_Atom2->getPosition();
    pos3 = this->_Atom3->getPosition();
    angle = geom::calculateAngle( pos1, pos2, pos3 );
    return angle;
}


void EnergyLinearAngle_O::addTerm(const EnergyLinearAngle& term)
{
    this->_Terms.push_back(term);
}



void	EnergyLinearAngle_O::dumpTerms(core::HashTable_sp atomTypes)
{
gctools::Vec0<EnergyLinearAngle>::iterator	eai;
string				as1,as2,as3;
string				str1, str2, str3;
    for (eai=this->_Terms.begin();eai!=this->_Terms.end();eai++)
    {
	as1 = atomLabel(eai->_Atom1);
	as2 = atomLabel(eai->_Atom2);
	as3 = atomLabel(eai->_Atom3);
	if ( as1 < as3 ) {
	    str1 = as1;
	    str2 = as2;
	    str3 = as3;
	} else {
	    str3 = as1;
	    str2 = as2;
	    str1 = as3;
	}
        core::clasp_write_string(fmt::format("TERM 2LINANG {:<9} - {:<9} - {:<9} {:8.2f}\n" , str1.c_str() , str2.c_str() , str3.c_str() , eai->term.kt ));
    }
}




void	EnergyLinearAngle_O::setupHessianPreconditioner(
					chem::NVector_sp nvPosition,
					chem::AbstractLargeSquareMatrix_sp m,
                                        core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  Linear_Angle<double*> linear_angle;
  constexpr size_t PS = Linear_Angle<double*>::PositionSize;  // 9
  double localPos[PS];
  double localForce[PS];
  double localHess[PS*PS];
  double localDvec[PS];
  double localHdvec[PS];
  for ( gctools::Vec0<EnergyLinearAngle>::iterator ai=this->_Terms.begin();
        ai!=this->_Terms.end(); ai++ ) {
    int I1 = ai->term.i3x1;
    int I2 = ai->term.i3x2;
    int I3 = ai->term.i3x3;
    if (hasActiveAtomMask &&
        !(bitvectorActiveAtomMask->testBit(I1/3) &&
          bitvectorActiveAtomMask->testBit(I2/3) &&
          bitvectorActiveAtomMask->testBit(I3/3))) continue;
    localPos[0] = (*nvPosition)[I1+0];
    localPos[1] = (*nvPosition)[I1+1];
    localPos[2] = (*nvPosition)[I1+2];
    localPos[3] = (*nvPosition)[I2+0];
    localPos[4] = (*nvPosition)[I2+1];
    localPos[5] = (*nvPosition)[I2+2];
    localPos[6] = (*nvPosition)[I3+0];
    localPos[7] = (*nvPosition)[I3+1];
    localPos[8] = (*nvPosition)[I3+2];
    // Pack the three atoms' 3 coords into a contiguous 9-vector at offsets 0, 3 and 6.
    // Zero local accumulator buffers (kernel uses += into them).
    for (size_t k = 0; k < PS; ++k) {
      localForce[k] = 0.0; localDvec[k] = 0.0; localHdvec[k] = 0.0;
    }
    for (size_t k = 0; k < PS*PS; ++k) localHess[k] = 0.0;
    double energyAccum = 0.0;
    linear_angle_term lt(ai->term.kt, 0, 3, 6);
    linear_angle.hessian(lt, localPos, &energyAccum, localForce, localHess,
                         localDvec, localHdvec);
    // Scatter the 9x9 local Hessian into the global sparse matrix.
    // Local offsets 0..2 map to global I1+0..I1+2, 3..5 to I2+0..I2+2 and 6..8 to I3+0..I3+2.
    // The kernel writes BOTH (r,c) and (c,r) for off-diagonal entries
    // (energyComponent.h:269-270), so accumulate only the lower triangle here.
    int globalIdx[PS] = { I1+0, I1+1, I1+2, I2+0, I2+1, I2+2, I3+0, I3+1, I3+2 };
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


double EnergyLinearAngle_O::evaluateAllComponent( ScoringFunction_sp score,
                                                  chem::NVector_sp 	pos,
                                                  core::T_sp energyScale,
                                                  core::T_sp energyComponents,
                                                  bool 		calcForce,
                                                  gc::Nilable<chem::NVector_sp> 	force,
                                                  bool		calcDiagonalHessian,
                                                  bool		calcOffDiagonalHessian,
                                                  gc::Nilable<chem::AbstractLargeSquareMatrix_sp>	hessian,
                                                  gc::Nilable<chem::NVector_sp> hdvec,
                                                  gc::Nilable<chem::NVector_sp> dvec,
                                                  core::T_sp activeAtomMask,
                                                  core::T_sp debugInteractions )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  num_real termEnergy = 0.0;
  this->_Evaluations++;

  auto evalType = determineEnergyComponentEvalType(force,hdvec,dvec);

  double  Energy = 0.0;  // referenced by LINEAR_ANGLE_DEBUG_INTERACTIONS macro
  double* position = &(*pos)[0];
  double* rforce = NULL;
  double* rdvec = NULL;
  double* rhdvec = NULL;
  Linear_Angle<NoHessian> linear_angle;
  gctools::Vec0<EnergyLinearAngle>::iterator ai;

#define KERNEL_LINEAR_ANGLE_APPLY_ATOM_MASK(I1,I2,I3) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         && bitvectorActiveAtomMask->testBit(I3/3) \
         ) \
    ) continue;

  if (evalType==energyEval) {
    for ( ai=this->_Terms.begin(); ai!=this->_Terms.end(); ai++ ) {
      KERNEL_LINEAR_ANGLE_APPLY_ATOM_MASK(ai->term.i3x1,ai->term.i3x2,ai->term.i3x3);
      Energy = linear_angle.energy(ai->term, position, &termEnergy);
      LINEAR_ANGLE_DEBUG_INTERACTIONS(ai->term.i3x1, ai->term.i3x2, ai->term.i3x3);
    }
  } else if (evalType==gradientEval) {
    rforce = &(*force)[0];
    for ( ai=this->_Terms.begin(); ai!=this->_Terms.end(); ai++ ) {
      KERNEL_LINEAR_ANGLE_APPLY_ATOM_MASK(ai->term.i3x1,ai->term.i3x2,ai->term.i3x3);
      Energy = linear_angle.gradient(ai->term, position, &termEnergy, rforce);
      LINEAR_ANGLE_DEBUG_INTERACTIONS(ai->term.i3x1, ai->term.i3x2, ai->term.i3x3);
    }
  } else { // hessianEval
    rforce = &(*force)[0];
    rdvec = &(*dvec)[0];
    rhdvec = &(*hdvec)[0];
    for ( ai=this->_Terms.begin(); ai!=this->_Terms.end(); ai++ ) {
      KERNEL_LINEAR_ANGLE_APPLY_ATOM_MASK(ai->term.i3x1,ai->term.i3x2,ai->term.i3x3);
      Energy = linear_angle.hessian(ai->term, position, &termEnergy, rforce, NoHessian(), rdvec, rhdvec);
      LINEAR_ANGLE_DEBUG_INTERACTIONS(ai->term.i3x1, ai->term.i3x2, ai->term.i3x3);
    }
  }
  maybeSetEnergy( energyComponents, EnergyLinearAngle_O::static_classSymbol(), termEnergy );
  return termEnergy;
}



SYMBOL_EXPORT_SC_(KeywordPkg,linear_angle);


void EnergyLinearAngle_O::initialize()
{
    this->Base::initialize();
}

void EnergyLinearAngle_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}


SYMBOL_EXPORT_SC_(KeywordPkg,kt);
SYMBOL_EXPORT_SC_(KeywordPkg,i1);
SYMBOL_EXPORT_SC_(KeywordPkg,i2);
SYMBOL_EXPORT_SC_(KeywordPkg,i3);
SYMBOL_EXPORT_SC_(KeywordPkg,atom1);
SYMBOL_EXPORT_SC_(KeywordPkg,atom2);
SYMBOL_EXPORT_SC_(KeywordPkg,atom3);

CL_DEFMETHOD core::List_sp EnergyLinearAngle_O::extract_vectors_as_alist() const {
  size_t size = this->_Terms.size();
  NVector_sp kt_vec = NVector_O::make(size);
  core::SimpleVector_int32_t_sp i1_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i2_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i3_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_sp atom1_vec    = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom2_vec    = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom3_vec    = core::SimpleVector_O::make(size);

  for (size_t i=0;i<size;++i) {
    const EnergyLinearAngle& entry = this->_Terms[i];
    (*kt_vec)[i] = entry.term.kt;
    (*i1_vec)[i] = entry.term.i3x1;
    (*i2_vec)[i] = entry.term.i3x2;
    (*i3_vec)[i] = entry.term.i3x3;
    (*atom1_vec)[i] = entry._Atom1;
    (*atom2_vec)[i] = entry._Atom2;
    (*atom3_vec)[i] = entry._Atom3;
  }
  return core::Cons_O::createList(core::Cons_O::create(kw::_sym_kt,kt_vec),
                                  core::Cons_O::create(kw::_sym_i1,i1_vec),
                                  core::Cons_O::create(kw::_sym_i2,i2_vec),
                                  core::Cons_O::create(kw::_sym_i3,i3_vec),
                                  core::Cons_O::create(kw::_sym_atom1,atom1_vec),
                                  core::Cons_O::create(kw::_sym_atom2,atom2_vec),
                                  core::Cons_O::create(kw::_sym_atom3,atom3_vec));

}

CL_DEFMETHOD void EnergyLinearAngle_O::fill_from_vectors_in_alist(core::List_sp vectors)
{
  NVector_sp kt_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_kt));
  core::SimpleVector_int32_t_sp i1_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i1));
  core::SimpleVector_int32_t_sp i2_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i2));
  core::SimpleVector_int32_t_sp i3_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i3));
  core::SimpleVector_sp atom1_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom1));
  core::SimpleVector_sp atom2_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom2));
  core::SimpleVector_sp atom3_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom3));
  this->_Terms.resize(kt_vec->length());
  for (size_t i=0, iEnd(kt_vec->length());i<iEnd;++i){
       EnergyLinearAngle& entry = this->_Terms[i];
       entry.term.kt = (*kt_vec)[i];
       entry.term.i3x1 = (*i1_vec)[i];
       entry.term.i3x2 = (*i2_vec)[i];
       entry.term.i3x3 = (*i3_vec)[i];
       entry._Atom1 = gc::As_unsafe<Atom_sp>((*atom1_vec)[i]);
       entry._Atom2 = gc::As_unsafe<Atom_sp>((*atom2_vec)[i]);
       entry._Atom3 = gc::As_unsafe<Atom_sp>((*atom3_vec)[i]);
  }
}

CL_DEFMETHOD void EnergyLinearAngle_O::addLinearAngleTerm(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, Atom_sp a3, double kt) {
  EnergyAtom* ea1 = atomTable->getEnergyAtomPointer(a1);
  EnergyAtom* ea2 = atomTable->getEnergyAtomPointer(a2);
  EnergyAtom* ea3 = atomTable->getEnergyAtomPointer(a3);
  EnergyLinearAngle energyLinearAngle(a1,a2,a3,
                                      ea1->coordinateIndexTimes3(),
                                      ea2->coordinateIndexTimes3(),
                                      ea3->coordinateIndexTimes3(),
                                      kt);
  this->addTerm(energyLinearAngle);
}




CL_DEFMETHOD
core::List_sp	EnergyLinearAngle_O::lookupLinearAngleTerms(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, Atom_sp a3, core::HashTable_sp atomTypes )
{
  ql::list  result;
  core::T_sp tia1 = atomTable->_AtomTableIndexes->gethash(a1);
  core::T_sp tia2 = atomTable->_AtomTableIndexes->gethash(a2);
  core::T_sp tia3 = atomTable->_AtomTableIndexes->gethash(a3);
  if (!tia1.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a1));
  if (!tia2.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a2));
  if (!tia3.fixnump()) SIMPLE_ERROR("Could not find {} in energy function" , _rep_(a3));
  for (auto edi=this->_Terms.begin();edi!=this->_Terms.end();edi++) {
    if ((edi->_Atom1==a1 &&
         edi->_Atom2==a2 &&
         edi->_Atom3==a3)
        || (edi->_Atom3==a1 &&
            edi->_Atom2==a2 &&
            edi->_Atom1==a3)) {
      ql::list oneResult;
      oneResult << INTERN_(kw,type1) << edi->_Atom1->getType(atomTypes)
                << INTERN_(kw,type2) << edi->_Atom2->getType(atomTypes)
                << INTERN_(kw,type3) << edi->_Atom3->getType(atomTypes)
                << INTERN_(kw,kt) << core::clasp_make_double_float(edi->term.kt);
      result << oneResult.result();
    }
  }
  return result.result();
}

EnergyComponent_sp EnergyLinearAngle_O::copyFilter(core::T_sp keepInteractionFactory , SetupAccumulator& setupAcc) {
  EnergyLinearAngle_sp copy = EnergyLinearAngle_O::create();
  copyEnergyComponent( copy, this->asSmartPtr() );
  core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyLinearAngle_O::staticClass() );
  for ( auto edi=this->_Terms.begin(); edi!=this->_Terms.end(); edi++ ) {
    Atom_sp a1 = edi->_Atom1;
    Atom_sp a2 = edi->_Atom2;
    Atom_sp a3 = edi->_Atom3;
    if ( skipInteraction( keepInteraction, a1, a2, a3 ) ) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}

void EnergyLinearAngle_O::emitTestCalls(core::T_sp stream, chem::NVector_sp ) const
{
  SIMPLE_ERROR("Update me");
}

size_t EnergyLinearAngle_O::runTestCalls(core::T_sp stream, chem::NVector_sp coords) const
{
#define POS_SIZE 9
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
  Linear_Angle<double*> linear_angle;
  for ( auto si=this->_Terms.begin();
        si!=this->_Terms.end(); si++ ) {
    position[0] = coords[si->term.i3x1];
    position[1] = coords[si->term.i3x1+1];
    position[2] = coords[si->term.i3x1+2];
    position[3] = coords[si->term.i3x2];
    position[4] = coords[si->term.i3x2+1];
    position[5] = coords[si->term.i3x2+2];
    position[6] = coords[si->term.i3x3];
    position[7] = coords[si->term.i3x3+1];
    position[8] = coords[si->term.i3x3+2];
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    linear_angle_term lt(si->term.kt, 0, 3, 6);
    linear_angle.gradient(    lt, position, &energy_new,    force_new );
    linear_angle.gradient_fd( lt, position, &energy_ground, force_ground );
    if (!test_match(stream, "linear_angle_gradient", POS_SIZE,
                    force_new, force_ground,
                    0, 0,
                    0, 0 )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH linear_angle_gradient #{} kt = {}\n",
                              idx, si->term.kt ), stream );
    }
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    linear_angle_term lt2(si->term.kt, 0, 3, 6);
    linear_angle.hessian(    lt2, position, &energy_new,    force_new,    hessian_new,    dvec_new,    hdvec_new );
    linear_angle.hessian_fd( lt2, position, &energy_ground, force_ground, hessian_ground, dvec_ground, hdvec_ground );
    if (!test_match( stream, "linear_angle_hessian", POS_SIZE,
                     force_new, force_ground,
                     hessian_new, hessian_ground,
                     hdvec_new, hdvec_ground )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH linear_angle_hessian #{} kt = {}\n",
                              idx, si->term.kt ), stream );
    }
    idx++;
  }
  core::print(fmt::format("linear_angle errors = {}\n", errs), stream);
  return errs;
}

};
