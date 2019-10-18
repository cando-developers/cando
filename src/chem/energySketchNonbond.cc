/*
    File: energySketchNonbond.cc


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
#define	DEBUG_LEVEL_FULL

#include <clasp/core/foundation.h>
#include <clasp/core/bformat.h>
#include <cando/chem/energySketchNonbond.h>
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


#define DEBUG_NONBOND_TERM 1
#define LOG_ENERGY(x)
//#define LOG_ENERGY core::write_bf_stream

namespace chem
{

core::List_sp EnergySketchNonbond::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,constant),core::clasp_make_double_float(this->_Constant)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->I1)),
                                  core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->I2)),
                                  core::Cons_O::create(INTERN_(kw,flags), core::make_fixnum(this->_FreezeFlags)));
}

void EnergySketchNonbond::decode(core::List_sp alist) {
  SIMPLE_ERROR(BF("Implement decode of EnergySketchNonbond"));
}

CL_DEFMETHOD void EnergySketchNonbond_O::setScaleSketchNonbond(double d) { this->_ScaleSketchNonbond = d; };
CL_DEFMETHOD double	EnergySketchNonbond_O::getScaleSketchNonbond()	{return this->_ScaleSketchNonbond; };

CL_DEFMETHOD void EnergySketchNonbond_O::setFreezeFlags(size_t flags) { this->_FreezeFlags = flags; };

CL_DEFMETHOD
void EnergySketchNonbond_O::addSketchNonbondTerm(size_t coordinate1IndexTimes3,
                                                 size_t coordinate2IndexTimes3,
                                                 size_t freezeFlags,
                                                 double constant)
{
  EnergySketchNonbond term(freezeFlags,constant,coordinate1IndexTimes3,coordinate2IndexTimes3);
  this->_Terms.push_back(term);
}

void	EnergySketchNonbond_O::zeroEnergy()
{
  this->Base::zeroEnergy();
  this->_Energy = 0.0;
}

double	EnergySketchNonbond_O::getEnergy()
{
  double	e;
  e = this->_Energy;
  return e;
}


void	EnergySketchNonbond_O::setupHessianPreconditioner(
                                                          NVector_sp nvPosition,
                                                          AbstractLargeSquareMatrix_sp m )
{
  SIMPLE_ERROR(BF("Nonbond term isn't used when calculating setupHessianPreconditioner but it was called!!!"));
}

double	EnergySketchNonbond_O::evaluateAllComponent( ScoringFunction_sp score,
                                            NVector_sp 	pos,
                                            bool 		calcForce,
                                            gc::Nilable<NVector_sp> 	force,
                                            bool		calcDiagonalHessian,
                                            bool		calcOffDiagonalHessian,
                                            gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                            gc::Nilable<NVector_sp>	hdvec, 
                                            gc::Nilable<NVector_sp> 	dvec )
{
    // Evaluate everything using terms
  this->evaluateTerms(pos,calcForce,force,calcDiagonalHessian,
                      calcOffDiagonalHessian,hessian,hdvec,dvec);
  return this->_TotalEnergy;
}
    
    

void	EnergySketchNonbond_O::evaluateTerms(NVector_sp 	pos,
                                             bool 		calcForce,
                                             gc::Nilable<NVector_sp> 	force,
                                             bool		calcDiagonalHessian,
                                             bool		calcOffDiagonalHessian,
                                             gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                             gc::Nilable<NVector_sp>	hdvec, 
                                             gc::Nilable<NVector_sp> 	dvec )
{
//  printf("%s:%d:%s Entering\n", __FILE__, __LINE__, __FUNCTION__ );
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
#define Log(x) log(x)
#define EREP_CALC_FORCE
#define EREP_CALC_DIAGONAL_HESSIAN
#define EREP_CALC_OFF_DIAGONAL_HESSIAN
#undef	EREP_SET_PARAMETER
#define	EREP_SET_PARAMETER(x)
#undef	EREP_SET_POSITION
#define	EREP_SET_POSITION(x,ii,of)	{x=coordinates_ptr[ii+of];};
#undef	EREP_ENERGY_ACCUMULATE
#define	EREP_ENERGY_ACCUMULATE(e) {this->_Energy+=(e);}
#undef	EREP_FORCE_ACCUMULATE
#undef	EREP_DIAGONAL_HESSIAN_ACCUMULATE
#undef	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	EREP_FORCE_ACCUMULATE(II,IO,val) 		{force_ptr[II+IO] += val;};
#define	EREP_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc
#define MAYBE_BAIL(val) {} if (val <=0.1 ) goto TOO_CLOSE; if ( val >= this->_LongDistanceCutoff) goto TOO_FAR;
  
#include <cando/chem/energy_functions/_Erep_termDeclares.cc>

  if ( this->isEnabled() ) 
  {
    double* coordinates_ptr = (double*)(pos->rowMajorAddressOfElement_(0));
    double* force_ptr = (double*)(force->rowMajorAddressOfElement_(0));
    double x1,y1,z1,x2,y2,z2,crep;
    double dx, dy, dz;
    double dsq, ERepDistance;
    double crep_over_dsq;
    double cutoff_sq = this->_LongDistanceCutoff*this->_LongDistanceCutoff;
    for ( size_t index = 0; index<this->_Terms.size(); ++index ) {
      EnergySketchNonbond& ea = this->_Terms[index];
      if (this->_FreezeFlags&ea._FreezeFlags) continue;
      int I1 = ea.I1;
      int I2 = ea.I2;
      crep = ea._Constant*this->_ScaleSketchNonbond;
      EREP_SET_POSITION(x1,I1,0);
      EREP_SET_POSITION(y1,I1,1);
      EREP_SET_POSITION(z1,I1,2);
      EREP_SET_POSITION(x2,I2,0);
      EREP_SET_POSITION(y2,I2,1);
      EREP_SET_POSITION(z2,I2,2);
      dx = (x2-x1);
      dy = (y2-y1);
      dz = (z2-z1);
      if (fabs(dx)>this->_LongDistanceCutoff) goto TOO_FAR;
      if (fabs(dy)>this->_LongDistanceCutoff) goto TOO_FAR;
      if (fabs(dz)>this->_LongDistanceCutoff) goto TOO_FAR;
      dsq = dx*dx+dy*dy+dz*dz;
      ErepDistance = sqrt(dsq);
      if (ErepDistance <= (0.1)) goto TOO_CLOSE;
      if (ErepDistance > this->_LongDistanceCutoff) goto TOO_FAR;
      crep_over_dsq = crep/dsq;
      fx1 = -dx*crep_over_dsq;
      fy1 = -dy*crep_over_dsq;
      fz1 = -dz*crep_over_dsq;
      fx2 = -fx1;
      fy2 = -fy1;
      fz2 = -fz1;
      EREP_FORCE_ACCUMULATE(I1,0,fx1);
      EREP_FORCE_ACCUMULATE(I1,1,fy1);
      EREP_FORCE_ACCUMULATE(I1,2,fz1);
      EREP_FORCE_ACCUMULATE(I2,0,fx2);
      EREP_FORCE_ACCUMULATE(I2,1,fy2);
      EREP_FORCE_ACCUMULATE(I2,2,fz2);
//#include <cando/chem/energy_functions/_Erep_termCode.cc>
#if 0
        if (chem__verbose(2)) {
          core::write_bf_stream(BF("(defparameter v1 (geom:vec %f %f %f))\n") % x1 % y1 % z1 );
          core::write_bf_stream(BF("(defparameter v2 (geom:vec %f %f %f))\n") % x2 % y2 % z2 );
          core::write_bf_stream(BF("crep = %f\n") % crep );
          core::write_bf_stream(BF("ErepDistance = %f\n") % ErepDistance );
          core::write_bf_stream(BF("Erep = %f\n") % Erep );
          core::write_bf_stream(BF("fx1 fy1 fz1 = %f %f %f\n") % fx1 % fy1 % fz1 );
          core::write_bf_stream(BF("fx2 fy2 fz2 = %f %f %f\n") % fx2 % fy2 % fz2 );
        }
#endif
        goto CONTINUE;
      TOO_FAR:
#if 0
        if (chem__verbose(2)) {
          core::write_bf_stream(BF("x1 y1 z1 = %f %f %f\n") % x1 % y1 % z1 );
          core::write_bf_stream(BF("x2 y2 z2 = %f %f %f\n") % x2 % y2 % z2 );
          core::write_bf_stream(BF("Too far ErepDistance = %f\n") % ErepDistance );
        }
#endif
        // Do nothing
        goto CONTINUE;
      TOO_CLOSE:
#if 0
        if (chem__verbose(2)) {
          core::write_bf_stream(BF("x1 y1 z1 = %f %f %f\n") % x1 % y1 % z1 );
          core::write_bf_stream(BF("x2 y2 z2 = %f %f %f\n") % x2 % y2 % z2 );
          core::write_bf_stream(BF("Too close ErepDistance = %f\n") % ErepDistance );
        }
#endif
        EREP_FORCE_ACCUMULATE(I1,0,0.5+core::randomNumber01());
        EREP_FORCE_ACCUMULATE(I1,1,0.5+core::randomNumber01());
        EREP_FORCE_ACCUMULATE(I1,2,0.5+core::randomNumber01());
        EREP_FORCE_ACCUMULATE(I2,0,0.5+core::randomNumber01());
        EREP_FORCE_ACCUMULATE(I2,1,0.5+core::randomNumber01());
        EREP_FORCE_ACCUMULATE(I2,2,0.5+core::randomNumber01());
      CONTINUE:
        (void)0;
    }
  }
}


void EnergySketchNonbond_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,LongDistanceCutoff), this->_LongDistanceCutoff );
  node->field( INTERN_(kw,ScaleSketchNonbond), this->_ScaleSketchNonbond );
  node->field( INTERN_(kw,terms), this->_Terms);
  this->Base::fields(node);
}

void EnergySketchNonbond_O::initialize()
{
  this->Base::initialize();
}


CL_DEFMETHOD
void EnergySketchNonbond_O::setLongDistanceCutoff(float cutoff) {
  this->_LongDistanceCutoff = cutoff;
}

CL_DOCSTRING(R"doc(Invoke a callback for every term in the energy-sketch-nonbond.
The callback takes the lambda-list (index atom1 atom2 atom1-index*3 atom2-index*3 constant).
Constant is a single-float to avoid consing.
Use the modify-sketch-nonbond-constant to modify the parameter.)doc");
CL_DEFMETHOD void EnergySketchNonbond_O::walkSketchNonbondTerms(core::T_sp callback)
{
  for (size_t i=0;i<this->_Terms.size();++i) {
    const EnergySketchNonbond& entry = this->_Terms[i];
    core::eval::funcall(callback,core::make_fixnum(i),
                        core::make_fixnum(entry._FreezeFlags),
                        core::make_fixnum(entry.I1),
                        core::make_fixnum(entry.I2),
                        core::make_single_float(entry._Constant));
  };
};

CL_DEFMETHOD void EnergySketchNonbond_O::modifySketchNonbondTermConstant(size_t index, float constant) {
  if (index<this->_Terms.size()) {
    EnergySketchNonbond& entry = this->_Terms[index];
    entry._Constant = constant;
    return;
  }
  SIMPLE_ERROR(BF("index %d is out of bounds as a energy-sketch-nonbond term index (#entries %d)") % index % this->_Terms.size() );
}

void EnergySketchNonbond_O::reset()
{
  this->_Terms.clear();
}


};

