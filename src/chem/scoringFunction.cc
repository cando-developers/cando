/*
    File: energyFunction.cc
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
//
// (C) 2004 Christian E. Schafmeister
//

/*
__BEGIN_DOC(classes.EnergyFunction,section,EnergyFunction)
EnergyFunction objects are used to calculate Amber or Gaff molecular 
mechanics energies and first and second derivatives.
__END_DOC
*/

/*
 *	energyFunction.cc
 *
 */

#include <clasp/core/common.h>
#include <clasp/core/bformat.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/loop.h>
#include <cando/adapt/indexedObjectBag.h>
#include <clasp/core/environment.h>
#include <cando/chem/minimizerLog.h>
#include <cando/chem/restraint.h>
#include <cando/chem/iterateRestraints.h>
#include <cando/adapt/iterateCons.h>
#include <cando/chem/ringFinder.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/atom.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyStretch.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/chemInfo.h>
#include <clasp/core/wrappers.h>


#ifdef	darwin
#include <xmmintrin.h>
#endif


namespace chem
{

ScoringFunction_O::ScoringFunction_O() : _VelocityScale(Vector3(1.0,1.0,1.0)) {};

CL_DEFMETHOD void ScoringFunction_O::setf_velocity_scale(double xscale, double yscale, double zscale) {
  this->_VelocityScale.set(xscale,yscale,zscale);
}

void ScoringFunction_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,Debug),this->_Debug);
  node->field_if_not_default(INTERN_(kw,VelocityScale),this->_VelocityScale,Vector3());
}

CL_LISPIFY_NAME("evaluateEnergy");
CL_DEFMETHOD double	ScoringFunction_O::evaluateEnergy( NVector_sp pos )
{
  double		energy;
  energy = this->evaluateAll(pos,
                             false,
                             nil<core::T_O>(),
                             false, false,
                             nil<core::T_O>(),
                             nil<core::T_O>(),
                             nil<core::T_O>() );
  return energy;
}



CL_LISPIFY_NAME("evaluateEnergyForce");
CL_DEFMETHOD double	ScoringFunction_O::evaluateEnergyForce( NVector_sp pos, bool calcForce, NVector_sp force )
{
  double	energy;
  gc::Nilable<NVector_sp>	rawGrad;
  if ( calcForce ) {
    rawGrad = force;
  } else {
    rawGrad = nil<core::T_O>();
  }
  energy = this->evaluateAll(pos,calcForce,
                             rawGrad,
                             false,
                             false,
                             nil<core::T_O>(),
                             nil<core::T_O>(),
                             nil<core::T_O>() );
  return energy;
}


CL_LISPIFY_NAME("evaluateEnergyForceFullHessian");
CL_DEFMETHOD double	ScoringFunction_O::evaluateEnergyForceFullHessian(
                                                                          NVector_sp pos,
                                                                          bool calcForce, NVector_sp force,
                                                                          bool calcDiagonalHessian,
                                                                          bool calcOffDiagonalHessian,
                                                                          AbstractLargeSquareMatrix_sp hessian )
{
  double	energy;
  gc::Nilable<NVector_sp> rawGrad;
  if ( calcForce ) rawGrad = force;
  else rawGrad = nil<core::T_O>();
  energy = this->evaluateAll( pos,
                              calcForce, rawGrad,
                              calcDiagonalHessian,
                              calcOffDiagonalHessian,
                              hessian,
                              nil<core::T_O>(),
                              nil<core::T_O>() );
  return energy;
}


CL_LISPIFY_NAME("evaluateEnergyForceFullHessianForDebugging");
CL_DEFMETHOD double	ScoringFunction_O::evaluateEnergyForceFullHessianForDebugging()
{
  NVector_sp	pos, force;
  AbstractLargeSquareMatrix_sp	hessian;
  double		energy;
  pos = NVector_O::create(this->getNVectorSize());
  force = NVector_O::create(this->getNVectorSize());
  hessian = FullLargeSquareMatrix_O::create(this->getNVectorSize(),SymmetricDiagonalLower);
  this->loadCoordinatesIntoVector(pos);
  energy = this->evaluateAll(pos,
                             true,
                             force,
                             true,
                             true,
                             hessian,
                             nil<core::T_O>(),
                             nil<core::T_O>() );
  return energy;
}




CL_LISPIFY_NAME("calculateEnergy");
CL_DEFMETHOD double	ScoringFunction_O::calculateEnergy( )
{
  NVector_sp	pos;
  pos = NVector_O::create(this->getNVectorSize());
  this->loadCoordinatesIntoVector(pos);
  return this->evaluateEnergy(pos);
}


CL_LISPIFY_NAME("calculateEnergyAndForce");
CL_DEFMETHOD core::T_mv ScoringFunction_O::calculateEnergyAndForce( )
{
  NVector_sp	pos;
  NVector_sp	force;
  double		energy;
  pos = NVector_O::create(this->getNVectorSize());
  force = NVector_O::create(this->getNVectorSize());
  this->loadCoordinatesIntoVector(pos);
  energy = this->evaluateEnergyForce(pos,true,force);
    	// To calculate the force magnitude use force->magnitude();
    	// To calculate the force rmsMagnitude use force->rmsMagnitude();
//  this->writeForceToAtoms(force);
  return Values(core::clasp_make_double_float(energy),force,pos);
}


CL_LISPIFY_NAME("checkForBeyondThresholdInteractions");
CL_DEFMETHOD core::List_sp	ScoringFunction_O::checkForBeyondThresholdInteractions(double threshold)
{
  SUBIMP();
}



CL_DOCSTRING(R"dx(Velocity-verlet-step moves the atoms one velocity-verlet-step.  
If tfrozen is a simple-bit-vector then it is used to determine frozen atoms.)dx")
CL_LISPIFY_NAME("velocity-verlet-step");
CL_DEFUN void chem__velocity_verlet_step(ScoringFunction_sp scoringFunc,
                                         NVector_sp position,
                                         NVector_sp velocity,
                                         NVector_sp force,
                                         NVector_sp force_dt,
                                         NVector_sp delta_t_over_mass,
                                         double delta_t,
                                         core::T_sp tfrozen)
{
  core::SimpleBitVector_sp frozen;
  if (gc::IsA<core::SimpleBitVector_sp>(tfrozen)) {
    frozen = gc::As_unsafe<core::SimpleBitVector_sp>(tfrozen);
    if (frozen->length() != (position->length())) {
      SIMPLE_ERROR(BF("frozen must be a simple-bit-vector of length %d or NIL - got %s") % (position->length()) % _rep_(tfrozen));
    }
  } else if (tfrozen.notnilp()) {
    SIMPLE_ERROR(BF("frozen must be a simple-bit-vector or NIL"));
  }
  double delta_tsquared = delta_t*delta_t;
  double delta_tsquared_div2 = delta_tsquared/2.0;
  NVector_sp position_dt = NVector_O::create(position->size());
  size_t atom_idx = 0;
  if (frozen) {
    for ( size_t idx = 0; idx<position->size(); idx += 3) {
      if (frozen->testBit(idx)==0) {
        double offsetx = delta_t*(*velocity)[idx+0] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+0];
        (*position)[idx+0] = (*position)[idx+0] + offsetx;
      }
      if (frozen->testBit(idx+1)==0) {
        double offsety = delta_t*(*velocity)[idx+1] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+1];
        (*position)[idx+1] = (*position)[idx+1] + offsety;
      }
      if (frozen->testBit(idx+2)==0) {
        double offsetz = delta_t*(*velocity)[idx+2] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+2];
        (*position)[idx+2] = (*position)[idx+2] + offsetz;
      }
      atom_idx++;
    }
    scoringFunc->evaluateEnergyForce(position,true,force_dt);
    atom_idx = 0;
    for ( size_t idx = 0; idx<position->size(); idx+=3 ) {
      if (frozen->testBit(idx)==0) {
        (*velocity)[idx+0] = ((*velocity)[idx+0] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+0]+(*force_dt)[idx+0]))*scoringFunc->_VelocityScale.getX();
        (*force)[idx+0] = (*force_dt)[idx+0];
      }
      if (frozen->testBit(idx+1)==0) {
        (*velocity)[idx+1] = ((*velocity)[idx+1] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+1]+(*force_dt)[idx+1]))*scoringFunc->_VelocityScale.getY();
        (*force)[idx+1] = (*force_dt)[idx+1];
      }
      if (frozen->testBit(idx+2)==0) {
        (*velocity)[idx+2] = ((*velocity)[idx+2] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+2]+(*force_dt)[idx+2]))*scoringFunc->_VelocityScale.getZ();
        (*force)[idx+2] = (*force_dt)[idx+2];
      }
      atom_idx++;
    }
  } else {
    for ( size_t idx = 0; idx<position->size(); idx += 3) {
      double offsetx = delta_t*(*velocity)[idx+0] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+0];
      double offsety = delta_t*(*velocity)[idx+1] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+1];
      double offsetz = delta_t*(*velocity)[idx+2] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+2];
      (*position)[idx+0] = (*position)[idx+0] + offsetx;
      (*position)[idx+1] = (*position)[idx+1] + offsety;
      (*position)[idx+2] = (*position)[idx+2] + offsetz;
      atom_idx++;
    }
    scoringFunc->evaluateEnergyForce(position,true,force_dt);
    atom_idx = 0;
    for ( size_t idx = 0; idx<position->size(); idx+=3 ) {
      (*velocity)[idx+0] = ((*velocity)[idx+0] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+0]+(*force_dt)[idx+0]))*scoringFunc->_VelocityScale.getX();
      (*velocity)[idx+1] = ((*velocity)[idx+1] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+1]+(*force_dt)[idx+1]))*scoringFunc->_VelocityScale.getY();
      (*velocity)[idx+2] = ((*velocity)[idx+2] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+2]+(*force_dt)[idx+2]))*scoringFunc->_VelocityScale.getZ();
      (*force)[idx+0] = (*force_dt)[idx+0];
      (*force)[idx+1] = (*force_dt)[idx+1];
      (*force)[idx+2] = (*force_dt)[idx+2];
      atom_idx++;
    }
  }
}
 


CL_DOCSTRING(R"dx(Like velocity-verlet-step but limits displacement of atoms in the x,y,z directions using the limit_displacement vector)dx")
CL_DOCSTRING_LONG(R"dx(Return the number of atoms whose displacement was limited.)dx")
CL_LISPIFY_NAME("velocity-verlet-step-limit-displacement");
CL_DEFUN size_t chem__velocity_verlet_step_limit_displacement(ScoringFunction_sp scoringFunc,
                                                              NVector_sp position,
                                                              NVector_sp velocity,
                                                              NVector_sp force,
                                                              NVector_sp force_dt,
                                                              NVector_sp delta_t_over_mass,
                                                              double delta_t,
                                                              core::T_sp tfrozen,
                                                              const Vector3& limit_displacement)
{
  core::SimpleBitVector_sp frozen;
  if (gc::IsA<core::SimpleBitVector_sp>(tfrozen)) {
    frozen = gc::As_unsafe<core::SimpleBitVector_sp>(tfrozen);
    if (frozen->length() != (position->length())) {
      SIMPLE_ERROR(BF("frozen must be a simple-bit-vector of length %d or NIL - got %s") % (position->length()) % _rep_(tfrozen));
    }
  } else if (tfrozen.notnilp()) {
    SIMPLE_ERROR(BF("frozen must be a simple-bit-vector or NIL"));
  }
  double delta_tsquared = delta_t*delta_t;
  double delta_tsquared_div2 = delta_tsquared/2.0;
  NVector_sp position_dt = NVector_O::create(position->size());
  size_t atom_idx = 0;
  size_t atoms_limited = 0;
  for ( size_t idx = 0; idx<position->size(); idx += 3) {
    double offsetx = delta_t*(*velocity)[idx+0] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+0];
    double offsety = delta_t*(*velocity)[idx+1] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+1];
    double offsetz = delta_t*(*velocity)[idx+2] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+2];
    if ((offsetx>limit_displacement.getX()) ||
        (offsety>limit_displacement.getY()) ||
        (offsetz>limit_displacement.getZ()) ||
        (offsetx<-limit_displacement.getX()) ||
        (offsety<-limit_displacement.getY()) ||
        (offsetz<-limit_displacement.getZ())) {
      atoms_limited++;
      if (offsetx>limit_displacement.getX()) offsetx = limit_displacement.getX();
      if (offsety>limit_displacement.getY()) offsety = limit_displacement.getY();
      if (offsetz>limit_displacement.getZ()) offsetz = limit_displacement.getZ();
      if (offsetx<-limit_displacement.getX()) offsetx = -limit_displacement.getX();
      if (offsety<-limit_displacement.getY()) offsety = -limit_displacement.getY();
      if (offsetz<-limit_displacement.getZ()) offsetz = -limit_displacement.getZ();
    }
    if (!frozen || frozen->testBit(idx  )==0) (*position)[idx+0] = (*position)[idx+0] + offsetx;
    if (!frozen || frozen->testBit(idx+1)==0) (*position)[idx+1] = (*position)[idx+1] + offsety;
    if (!frozen || frozen->testBit(idx+2)==0) (*position)[idx+2] = (*position)[idx+2] + offsetz;
    atom_idx++;
  }
  scoringFunc->evaluateEnergyForce(position,true,force_dt);
  atom_idx = 0;
  for ( size_t idx = 0; idx<position->size(); idx+=3 ) {
    if (!frozen || frozen->testBit(idx)==0) {
      (*velocity)[idx+0] = ((*velocity)[idx+0] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+0]+(*force_dt)[idx+0]))*scoringFunc->_VelocityScale.getX();
      (*force)[idx+0] = (*force_dt)[idx+0];
    }
    if (!frozen || frozen->testBit(idx+1)==0) {
      (*velocity)[idx+1] = ((*velocity)[idx+1] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+1]+(*force_dt)[idx+1]))*scoringFunc->_VelocityScale.getY();
      (*force)[idx+1] = (*force_dt)[idx+1];
    }
    if (!frozen || frozen->testBit(idx+2)==0) {
      (*velocity)[idx+2] = ((*velocity)[idx+2] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+2]+(*force_dt)[idx+2]))*scoringFunc->_VelocityScale.getZ();
      (*force)[idx+2] = (*force_dt)[idx+2];
    }
    atom_idx++;
  }
  return atoms_limited;
}


CL_DEFUN core::List_sp chem__all_components(ScoringFunction_sp scoring_function) {
  return scoring_function->allComponents();
}
};
