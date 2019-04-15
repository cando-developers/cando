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

CL_LISPIFY_NAME("evaluateEnergy");
CL_DEFMETHOD double	ScoringFunction_O::evaluateEnergy( NVector_sp pos )
{
  double		energy;
  energy = this->evaluateAll(pos,
                             false,
                             _Nil<core::T_O>(),
                             false, false,
                             _Nil<core::T_O>(),
                             _Nil<core::T_O>(),
                             _Nil<core::T_O>() );
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
    rawGrad = _Nil<core::T_O>();
  }
  energy = this->evaluateAll(pos,calcForce,
                             rawGrad,
                             false,
                             false,
                             _Nil<core::T_O>(),
                             _Nil<core::T_O>(),
                             _Nil<core::T_O>() );
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
  else rawGrad = _Nil<core::T_O>();
  energy = this->evaluateAll( pos,
                              calcForce, rawGrad,
                              calcDiagonalHessian,
                              calcOffDiagonalHessian,
                              hessian,
                              _Nil<core::T_O>(),
                              _Nil<core::T_O>() );
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
                             _Nil<core::T_O>(),
                             _Nil<core::T_O>() );
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
CL_DEFMETHOD double ScoringFunction_O::calculateEnergyAndForce( )
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
  return energy;
}


CL_LISPIFY_NAME("checkForBeyondThresholdInteractions");
CL_DEFMETHOD core::List_sp	ScoringFunction_O::checkForBeyondThresholdInteractions(double threshold)
{
  SUBIMP();
}



};
