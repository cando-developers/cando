/*
    File: energyAngle.cc
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
#include <cando/chem/energyAngle.h>
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
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>

namespace chem {

CL_LISPIFY_NAME(make-energy-angle);
CL_DEF_CLASS_METHOD
EnergyAngle_sp EnergyAngle_O::make(EnergyFunction_sp energyFunction) {
  auto component = ensureComponent<EnergyAngle_O>(energyFunction);
  return component;
}

};


namespace chem {
// ----------------------------------------------------------------------
// include energy functions
// ----------------------------------------------------------------------

#define VERYSMALL 1e-8

class LinearAngleError {
public:
  LinearAngleError() {};
};

void linear_angle_error() {
  throw LinearAngleError();
}

inline double linear_angle_clamp(double cos_theta) {
  return copysign(1.0-VERYSMALL,cos_theta);
}

#include "cando/chem/energyKernels/angle.c"

};


namespace chem {

#define ENERGY_FUNCTION I1, I2, I3, activeAtomMask

#define ANGLE_APPLY_ATOM_MASK(I1,I2,I3) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         && bitvectorActiveAtomMask->testBit(I3/3) \
         ) \
    ) goto SKIP_term;


double	_evaluateEnergyOnly_Angle(
    int I1, int I2, int I3, core::T_sp activeAtomMask,
		num_real x1, num_real y1, num_real z1,
		num_real x2, num_real y2, num_real z2,
		num_real x3, num_real y3, num_real z3,
		num_real t0, num_real kt )
{
  // I1, I2, I3, and activeAtomMask are part of the public API for symmetry with
  // the gradient/hessian forms; they have no effect on a single-term energy
  // calculation, so suppress unused-parameter warnings.
  (void)I1; (void)I2; (void)I3; (void)activeAtomMask;
  double localPos[9] = { (double)x1, (double)y1, (double)z1,
                         (double)x2, (double)y2, (double)z2,
                         (double)x3, (double)y3, (double)z3 };
  double energy_accum = 0.0;
  Angle<NoHessian> angle;
  angle_term lt((double)kt, (double)t0, 0, 3, 6);
  return angle.energy(lt, localPos, &energy_accum);
}

std::string EnergyAngle_O::descriptionOfContents() const {
  stringstream ss;
  ss << ":enabled " << ((this->_Enabled) ? "T" : "NIL");
  ss << " number-of-terms " << this->_Terms.size();
  return ss.str();
}

CL_LAMBDA((energy-angle chem:energy-angle) pos &optional activeAtomMask);
CL_DEFMETHOD void	EnergyAngle_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(chem::NVector_sp 	pos, core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  constexpr size_t PS = Angle<double*>::PositionSize;  // 9
  const double energyTol = 1.0e-7;
  const double forceTol  = 1.0e-5;
  const double hessTol   = 1.0e-3;
  double localPos[PS];
  double force_a[PS], force_fd[PS];
  double hess_a[PS*PS], hess_fd[PS*PS];
  double dvec_a[PS], dvec_fd[PS];
  double hdvec_a[PS], hdvec_fd[PS];
  Angle<double*> angle;
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
    angle_term lt(ai->term.kt, ai->term.t0, 0, 3, 6);
    angle.hessian(   lt, localPos, &e_a,  force_a,  hess_a,  dvec_a,  hdvec_a);
    angle.hessian_fd(lt, localPos, &e_fd, force_fd, hess_fd, dvec_fd, hdvec_fd);
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
        "  -> term {} I1={} I2={} I3={} kt={} t0={}\n",
        idx, I1, I2, I3, ai->term.kt, ai->term.t0));
    }
  }
  core::lisp_write(fmt::format(
    "compareAnalyticalAndNumericalForceAndHessianTermByTerm: {} fails out of {} terms\n",
    fails, idx));
}





#define ANGLE_DEBUG_INTERACTIONS(I1,I2,I3) \
    if (doDebugInteractions) { \
      core::eval::funcall(debugInteractions,EnergyAngle_O::static_classSymbol(), \
                          mk_double_float(Energy), \
                          core::make_fixnum(I1), core::make_fixnum(I2), core::make_fixnum(I3)); \
    }

core::List_sp EnergyAngle::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,kt),core::clasp_make_double_float(this->term.kt)),
                                  core::Cons_O::create(INTERN_(kw,t0_deg),core::clasp_make_double_float(this->term.t0/0.0174533)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.i3x1)),
                                  core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->term.i3x2)),
                                  core::Cons_O::create(INTERN_(kw,i3), core::make_fixnum(this->term.i3x3)),
                                  core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1),
                                  core::Cons_O::create(INTERN_(kw,atom2), this->_Atom2),
                                  core::Cons_O::create(INTERN_(kw,atom3), this->_Atom3));
}

void EnergyAngle::decode(core::List_sp alist) {
  SIMPLE_ERROR("Implement decode of EnergyAngle");
}



void EnergyAngle::defineFrom( FFAngle_sp term , EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3, double scale)
{

//    this->_Term = term;
//    this->_Type1 = term->_Type1;
//    this->_Type2 = term->_Type2;
//    this->_Type3 = term->_Type3;
    this->term.t0 = term->_AngRad;
    this->term.kt = term->getK2_kcalPerRadianSquared()*scale;
//    this->_K3 = term->_K3;
//    this->_K4 = term->_K4;
//    this->_Ub_k = term->_Ub_k;
//    this->_Ub_len = term->_Ub_len;
    this->_Atom1 = ea1->atom();
    this->_Atom2 = ea2->atom();
    this->_Atom3 = ea3->atom();
    this->term.i3x1 = ea1->coordinateIndexTimes3();
    this->term.i3x2 = ea2->coordinateIndexTimes3();
    this->term.i3x3 = ea3->coordinateIndexTimes3();
}



#if 0 //[
void	EnergyAngle::defineMissing( EnergyAtom *ea1, EnergyAtom *ea2, EnergyAtom *ea3 )
{
    this->_Type1 = ea1->_Atom->getType(atomTypes);
    this->_Type2 = ea2->_Atom->getType(atomTypes);
    this->_Type3 = ea3->_Atom->getType(atomTypes);
}
#endif //]

double EnergyAngle::getT()
{
    Vector3	pos1, pos2, pos3;
num_real	angle;
    pos1 = this->_Atom1->getPosition();
    pos2 = this->_Atom2->getPosition();
    pos3 = this->_Atom3->getPosition();
    angle = geom::calculateAngle( pos1, pos2, pos3 );
    return angle;
}

#if 0
    adapt::QDomNode_sp	EnergyAngle::asXml()
{
    adapt::QDomNode_sp	node;
    Vector3	vdiff;

    node = adapt::QDomNode_O::create(env,"EnergyAngle");
    node->addAttributeString("atom1Name",this->_Atom1->getName());
    node->addAttributeString("atom2Name",this->_Atom2->getName());
    node->addAttributeString("atom3Name",this->_Atom3->getName());
//    node->addAttributeString("atom1Type",this->_Type1 );
//    node->addAttributeString("atom2Type",this->_Type2 );
//    node->addAttributeString("atom3Type",this->_Type3 );
    node->addAttributeInt("I1",this->term.i3x1);
    node->addAttributeInt("I2",this->term.i3x2);
    node->addAttributeInt("I3",this->term.i3x3);
    node->addAttributeDoubleScientific("kt",this->term.kt);
    node->addAttributeDoubleScientific("t0",this->term.t0);
#if TURN_ENERGY_FUNCTION_DEBUG_ON
    adapt::QDomNode_sp xml = adapt::QDomNode_O::create(env,"Evaluated");
    xml->addAttributeBool("calcForce",this->_calcForce );
    xml->addAttributeBool("calcDiagonalHessian",this->_calcDiagonalHessian );
    xml->addAttributeBool("calcOffDiagonalHessian",this->_calcOffDiagonalHessian );
#include <_Angle_debugEvalXml.cc>
    node->addChild(xml);
#endif
    return node;
}

    void	EnergyAngle::parseFromXmlUsingAtomTable(adapt::QDomNode_sp	xml,
					AtomTable_sp atomTable )
{
    int i1 = xml->getAttributeInt("idx1");
    int i2 = xml->getAttributeInt("idx2");
    int i3 = xml->getAttributeInt("idx3");
    this->_Atom1 = atomTable->findEnergyAtomWithCoordinateIndex(i1)->atom();
    this->_Atom2 = atomTable->findEnergyAtomWithCoordinateIndex(i2)->atom();
    this->_Atom3 = atomTable->findEnergyAtomWithCoordinateIndex(i3)->atom();
    this->term.kt = xml->getAttributeDouble("kt");
    this->term.t0 = xml->getAttributeDouble("t0");
    this->term.i3x1 = xml->getAttributeInt("I1");
    this->term.i3x2 = xml->getAttributeInt("I2");
    this->term.i3x3 = xml->getAttributeInt("I3");
}
#endif


void EnergyAngle_O::addTerm(const EnergyAngle& term)
{
    this->_Terms.push_back(term);
}



void	EnergyAngle_O::dumpTerms(core::HashTable_sp atomTypes)
{
gctools::Vec0<EnergyAngle>::iterator	eai;
string				as1,as2,as3,as4;
string				str1, str2, str3, str4;
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
        core::clasp_write_string(fmt::format("TERM 2ANG {:<9} - {:<9} - {:<9} {:8.2f} {:8.2f}\n" , str1.c_str() , str2.c_str() , str3.c_str() , eai->term.kt , core::degrees(eai->term.t0) ));
    }
}





void	EnergyAngle_O::setupHessianPreconditioner(
					chem::NVector_sp nvPosition,
					chem::AbstractLargeSquareMatrix_sp m,
                                        core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  Angle<double*> angle;
  constexpr size_t PS = Angle<double*>::PositionSize;  // 9
  double localPos[PS];
  double localForce[PS];
  double localHess[PS*PS];
  double localDvec[PS];
  double localHdvec[PS];
  for ( gctools::Vec0<EnergyAngle>::iterator ai=this->_Terms.begin();
        ai!=this->_Terms.end(); ai++ ) {
    {
      int I1 = ai->term.i3x1;
      int I2 = ai->term.i3x2;
      int I3 = ai->term.i3x3;
      ANGLE_APPLY_ATOM_MASK(I1, I2, I3);
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
      angle_term lt(ai->term.kt, ai->term.t0, 0, 3, 6);
      angle.hessian(lt, localPos, &energyAccum, localForce, localHess,
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
    SKIP_term:;
  }
}


double EnergyAngle_O::evaluateAllComponent( ScoringFunction_sp score,
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

  EnergyAngle* badAngle = NULL;

  //
  // Copy from implementAmberFunction::evaluateAll
  //
#define ANGLE_CALC_FORCE
#define ANGLE_CALC_DIAGONAL_HESSIAN
#define ANGLE_CALC_OFF_DIAGONAL_HESSIAN
#undef ANGLE_SET_PARAMETER
#define ANGLE_SET_PARAMETER(x)	{x=ai->term.x;}
#undef ANGLE_SET_POSITION
#define ANGLE_SET_POSITION(x,ii,of)	{x=pos->element(ii+of);}
#undef ANGLE_ENERGY_ACCUMULATE
#define ANGLE_ENERGY_ACCUMULATE(e) { termEnergy += (e); }
#undef	ANGLE_FORCE_ACCUMULATE
#undef	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE
#undef	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANGLE_FORCE_ACCUMULATE 		ForceAcc
#define	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc



  {
    num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,kt,t0; //,angleScale;
    //	double DotAbCb;
    int I1, I2, I3,i;
    DOUBLE* position = &(*pos)[0];
    DOUBLE* rforce = NULL;
    DOUBLE* rhessian = NULL; // &(*hessian)[0];
    DOUBLE* rdvec = NULL;
    DOUBLE* rhdvec = NULL;
    double  Energy;
    Angle<NoHessian> angle;
    gctools::Vec0<EnergyAngle>::iterator ai;

#define KERNEL_ANGLE_APPLY_ATOM_MASK(I1,I2,I3) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         && bitvectorActiveAtomMask->testBit(I3/3) \
         ) \
    ) continue;

    if (evalType==energyEval) {
      for ( i=0,ai=this->_Terms.begin(); ai!=this->_Terms.end(); ai++,i++ ) {
        KERNEL_ANGLE_APPLY_ATOM_MASK(ai->term.i3x1,ai->term.i3x2,ai->term.i3x3);
        try {
          Energy = angle.energy(ai->term, position, &termEnergy);
          ANGLE_DEBUG_INTERACTIONS(ai->term.i3x1, ai->term.i3x2, ai->term.i3x3);
        } catch (LinearAngleError(err)) {
          badAngle = ai;
          goto ERROR_LINEAR_ANGLE;
        }
      }
    } else if (evalType==gradientEval) {
      rforce = &(*force)[0];
      for ( i=0,ai=this->_Terms.begin(); ai!=this->_Terms.end(); ai++,i++ ) {
        KERNEL_ANGLE_APPLY_ATOM_MASK(ai->term.i3x1,ai->term.i3x2,ai->term.i3x3);
        try {
          Energy = angle.gradient(ai->term, position, &termEnergy, rforce);
          ANGLE_DEBUG_INTERACTIONS(ai->term.i3x1, ai->term.i3x2, ai->term.i3x3);
        } catch (LinearAngleError(err)) {
          badAngle = ai;
          goto ERROR_LINEAR_ANGLE;
        }
      }
    } else { // if (hasHdAndD) {
      rforce = &(*force)[0];
      rdvec = &(*dvec)[0];
      rhdvec = &(*hdvec)[0];
      for ( i=0,ai=this->_Terms.begin(); ai!=this->_Terms.end(); ai++,i++ ) {
        KERNEL_ANGLE_APPLY_ATOM_MASK(ai->term.i3x1,ai->term.i3x2,ai->term.i3x3);
        try {
          Energy = angle.hessian(ai->term, position, &termEnergy, rforce, NoHessian(), rdvec, rhdvec);
          ANGLE_DEBUG_INTERACTIONS(ai->term.i3x1, ai->term.i3x2, ai->term.i3x3);
        } catch (LinearAngleError(err)) {
          badAngle = ai;
          goto ERROR_LINEAR_ANGLE;
        }
      }
    }
  }
  maybeSetEnergy( energyComponents, EnergyAngle_O::static_classSymbol(), termEnergy );
  return termEnergy;
 ERROR_LINEAR_ANGLE:
  ERROR(chem::_sym_LinearAngleError,core::Cons_O::createList(kw::_sym_atoms,core::Cons_O::createList(badAngle->_Atom1,badAngle->_Atom2,badAngle->_Atom3),
                                                             kw::_sym_coordinates,pos,
                                                             kw::_sym_indexes,core::Cons_O::createList(core::make_fixnum(badAngle->term.i3x1), core::make_fixnum(badAngle->term.i3x2), core::make_fixnum(badAngle->term.i3x3))));
}



SYMBOL_EXPORT_SC_(KeywordPkg,angle);
SYMBOL_EXPORT_SC_(KeywordPkg,atoms);
SYMBOL_EXPORT_SC_(KeywordPkg,t0);
SYMBOL_EXPORT_SC_(KeywordPkg,angle_deviation);


void EnergyAngle_O::initialize()
{
    this->Base::initialize();
}

void EnergyAngle_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}


SYMBOL_EXPORT_SC_(KeywordPkg,kt);
SYMBOL_EXPORT_SC_(KeywordPkg,t0);
SYMBOL_EXPORT_SC_(KeywordPkg,i1);
SYMBOL_EXPORT_SC_(KeywordPkg,i2);
SYMBOL_EXPORT_SC_(KeywordPkg,i3);
SYMBOL_EXPORT_SC_(KeywordPkg,atom1);
SYMBOL_EXPORT_SC_(KeywordPkg,atom2);
SYMBOL_EXPORT_SC_(KeywordPkg,atom3);

CL_DEFMETHOD core::List_sp EnergyAngle_O::extract_vectors_as_alist() const {
  size_t size = this->_Terms.size();
  NVector_sp kt_vec = NVector_O::make(size);
  NVector_sp t0_vec = NVector_O::make(size);
  core::SimpleVector_int32_t_sp i1_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i2_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_int32_t_sp i3_vec = core::SimpleVector_int32_t_O::make(size);
  core::SimpleVector_sp atom1_vec    = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom2_vec    = core::SimpleVector_O::make(size);
  core::SimpleVector_sp atom3_vec    = core::SimpleVector_O::make(size);

  for (size_t i=0;i<size;++i) {
    const EnergyAngle& entry = this->_Terms[i];
    (*kt_vec)[i] = entry.term.kt;
    (*t0_vec)[i] = entry.term.t0;
    (*i1_vec)[i] = entry.term.i3x1;
    (*i2_vec)[i] = entry.term.i3x2;
    (*i3_vec)[i] = entry.term.i3x3;
    (*atom1_vec)[i] = entry._Atom1;
    (*atom2_vec)[i] = entry._Atom2;
    (*atom3_vec)[i] = entry._Atom3;
  }
  return core::Cons_O::createList(core::Cons_O::create(kw::_sym_kt,kt_vec),
                                  core::Cons_O::create(kw::_sym_t0,t0_vec),
                                  core::Cons_O::create(kw::_sym_i1,i1_vec),
                                  core::Cons_O::create(kw::_sym_i2,i2_vec),
                                  core::Cons_O::create(kw::_sym_i3,i3_vec),
                                  core::Cons_O::create(kw::_sym_atom1,atom1_vec),
                                  core::Cons_O::create(kw::_sym_atom2,atom2_vec),
                                  core::Cons_O::create(kw::_sym_atom3,atom3_vec));

}

CL_DEFMETHOD void EnergyAngle_O::fill_from_vectors_in_alist(core::List_sp vectors)
{
  NVector_sp kt_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_kt));
  NVector_sp t0_vec = (safe_alist_lookup<NVector_sp>(vectors,kw::_sym_t0));
  core::SimpleVector_int32_t_sp i1_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i1));
  core::SimpleVector_int32_t_sp i2_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i2));
  core::SimpleVector_int32_t_sp i3_vec = (safe_alist_lookup<core::SimpleVector_int32_t_sp>(vectors,kw::_sym_i3));
  core::SimpleVector_sp atom1_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom1));
  core::SimpleVector_sp atom2_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom2));
  core::SimpleVector_sp atom3_vec = (safe_alist_lookup<core::SimpleVector_sp>(vectors,kw::_sym_atom3));
  this->_Terms.resize(kt_vec->length());
  for (size_t i=0, iEnd(kt_vec->length());i<iEnd;++i){
       EnergyAngle& entry = this->_Terms[i];
       entry.term.kt = (*kt_vec)[i];
       entry.term.t0 = (*t0_vec)[i];
       entry.term.i3x1 = (*i1_vec)[i];
       entry.term.i3x2 = (*i2_vec)[i];
       entry.term.i3x3 = (*i3_vec)[i];
       entry._Atom1 = gc::As_unsafe<Atom_sp>((*atom1_vec)[i]);
       entry._Atom2 = gc::As_unsafe<Atom_sp>((*atom2_vec)[i]);
       entry._Atom3 = gc::As_unsafe<Atom_sp>((*atom3_vec)[i]);
  }
}

CL_DEFMETHOD void EnergyAngle_O::addAngleTerm(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, Atom_sp a3, double kt, double t0) {
  EnergyAtom* ea1 = atomTable->getEnergyAtomPointer(a1);
  EnergyAtom* ea2 = atomTable->getEnergyAtomPointer(a2);
  EnergyAtom* ea3 = atomTable->getEnergyAtomPointer(a3);
  EnergyAngle energyAngle(a1,a2,a3,
                          ea1->coordinateIndexTimes3(),
                          ea2->coordinateIndexTimes3(),
                          ea3->coordinateIndexTimes3(),
                          kt,t0);
  this->addTerm(energyAngle);
}




CL_DEFMETHOD
core::List_sp	EnergyAngle_O::lookupAngleTerms(AtomTable_sp atomTable, Atom_sp a1, Atom_sp a2, Atom_sp a3, core::HashTable_sp atomTypes )
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
                << INTERN_(kw,kt) << core::clasp_make_double_float(edi->term.kt)
                << INTERN_(kw,t0) << core::clasp_make_double_float(edi->term.t0/0.0174533);
      result << oneResult.result();
    }
  }
  return result.result();
}

EnergyComponent_sp EnergyAngle_O::copyFilter(core::T_sp keepInteractionFactory , SetupAccumulator& setupAcc) {
  EnergyAngle_sp copy = EnergyAngle_O::create();
  copyEnergyComponent( copy, this->asSmartPtr() );
  core::T_sp keepInteraction = specializeKeepInteractionFactory( keepInteractionFactory, EnergyAngle_O::staticClass() );
  for ( auto edi=this->_Terms.begin(); edi!=this->_Terms.end(); edi++ ) {
    Atom_sp a1 = edi->_Atom1;
    Atom_sp a2 = edi->_Atom2;
    Atom_sp a3 = edi->_Atom3;
    if ( skipInteraction( keepInteraction, a1, a2, a3 ) ) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}

void EnergyAngle_O::emitTestCalls(core::T_sp stream, chem::NVector_sp ) const
{
  SIMPLE_ERROR("Update me");
#if 0
  for ( auto si=this->_Terms.begin();
        si!=this->_Terms.end(); si++ ) {
    core::print(fmt::format("test_angle( errs, {}, {}, {}, {}, {}, position_size, position, force_new, force_old, hessian_new, hessian_old, dvec_new, dvec_old, hdvec_new, hdvec_old );\n",
                            si->term.kt, si->term.t0, si->term.i3x1, si->term.i3x2, si->term.i3x3 ),stream);
  }
#endif
}

size_t EnergyAngle_O::runTestCalls(core::T_sp stream, chem::NVector_sp coords) const
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
  Angle<double*> angle;
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
    try {
      angle_term lt(si->term.kt, si->term.t0, 0, 3, 6);
      angle.gradient(    lt, position, &energy_new,    force_new );
      angle.gradient_fd( lt, position, &energy_ground, force_ground );
    } catch (LinearAngleError err) {
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("In angle_gradient hit linear-angle-error index = {}\n", idx), stream);
      errs++;
    }
    if (!test_match(stream, "force", POS_SIZE,
                    force_new, force_ground,
                    0, 0,
                    0, 0 )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH angle_gradient #{} kt = {}  t0 = {}\n",
                              idx,
                              si->term.kt, si->term.t0 ), stream );
    }
    energy_new = 0.0;
    energy_ground = 0.0;
    test_zero( POS_SIZE,
               force_new, force_ground,
               hessian_new, hessian_ground,
               dvec_new, dvec_ground,
               hdvec_new, hdvec_ground );
    try {
      angle_term lt2(si->term.kt, si->term.t0, 0, 3, 6);
      angle.hessian(    lt2, position, &energy_new,    force_new,    hessian_new,    dvec_new,    hdvec_new );
      angle.hessian_fd( lt2, position, &energy_ground, force_ground, hessian_ground, dvec_ground, hdvec_ground );
    } catch (LinearAngleError err) {
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("In angle_hessian hit linear-angle-error index = {}\n", idx), stream);
      errs++;
    }
    if (!test_match( stream, "hessian", POS_SIZE,
                     force_new, force_ground,
                     hessian_new, hessian_ground,
                     hdvec_new, hdvec_ground )) {
      errs++;
      test_position( stream, POS_SIZE, position );
      core::print(fmt::format("MISMATCH angle_hessian #{} kt = {}  t0 = {}\n",
                              idx,
                              si->term.kt, si->term.t0 ), stream );
    }
    idx++;
  }
  core::print(fmt::format("angle errors = {}\n", errs), stream);
  return errs;
}

};
