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
#include <clasp/core/numbers.h>
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


#include <cando/chem/rigidBodyEnergyFunction.h>

#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/chemInfo.h>
#include <clasp/core/wrappers.h>


/* Nothing */

namespace chem {

CL_LISPIFY_NAME("CHEM:MAKE-RIGID-BODY-ENERGY-FUNCTION");
CL_DEFUN RigidBodyEnergyFunction_sp RigidBodyEnergyFunction_O::make(size_t number_of_rigid_bodies) {
  GC_ALLOCATE_VARIADIC(RigidBodyEnergyFunction_O, energy, number_of_rigid_bodies);
  energy->_SavedCoordinates = NVector_O::create(number_of_rigid_bodies*7);
  return energy;
}

string RigidBodyEnergyFunction_O::energyTermsEnabled() {
  stringstream ss;
  size_t termCount = 0;
  for ( auto cur : this->_Terms ) {
    EnergyRigidBodyComponent_sp term = gc::As<EnergyRigidBodyComponent_sp>(CONS_CAR(cur));
    if ( termCount ) ss << " ";
    ss << core::lisp_instance_class(term)->classNameAsString();
  }
  return ss.str();
}

void RigidBodyEnergyFunction_O::enableDebug() {
  for ( auto cur : this->_Terms ) {
    EnergyRigidBodyComponent_sp term = gc::As<EnergyRigidBodyComponent_sp>(CONS_CAR(cur));
    term->enableDebug();
  }
}

void RigidBodyEnergyFunction_O::disableDebug() {
  for ( auto cur : this->_Terms ) {
    EnergyRigidBodyComponent_sp term = gc::As<EnergyRigidBodyComponent_sp>(CONS_CAR(cur));
    term->disableDebug();
  }
}

void RigidBodyEnergyFunction_O::setupHessianPreconditioner( NVector_sp pos, AbstractLargeSquareMatrix_sp hessian) {
  IMPLEMENT_ME();
}

void	RigidBodyEnergyFunction_O::loadCoordinates(NVector_sp pos)
{
  IMPLEMENT_ME();
}

void RigidBodyEnergyFunction_O::saveCoordinates(NVector_sp pos)
{
  IMPLEMENT_ME();
}

void RigidBodyEnergyFunction_O::saveCoordinatesAndForces(NVector_sp pos, NVector_sp force)
{
  IMPLEMENT_ME();
}


double	RigidBodyEnergyFunction_O::evaluateRaw( NVector_sp pos, NVector_sp force ) {

  IMPLEMENT_ME();
}

//    virtual double	evaluate( NVector_sp pos, NVector_sp force, bool calculateForce ) ;
  adapt::QDomNode_sp	identifyTermsBeyondThreshold();
//    uint	countBadVdwInteractions(double scaleSumOfVdwRadii, geom::DisplayList_sp displayIn);

ForceMatchReport_sp RigidBodyEnergyFunction_O::checkIfAnalyticalForceMatchesNumericalForce( NVector_sp pos, NVector_sp force ) {
  IMPLEMENT_ME();
}

void	RigidBodyEnergyFunction_O::useDefaultSettings()
{
  IMPLEMENT_ME();
}

    /*! Set a single options */
void	RigidBodyEnergyFunction_O::setOption( core::Symbol_sp option, core::T_sp val)
{
  IMPLEMENT_ME();
}



    /*! Set the energy function options. List the options as a flat list of keyword/value pairs */
void	RigidBodyEnergyFunction_O::setOptions( core::List_sp options )
{
  IMPLEMENT_ME();
}


double	RigidBodyEnergyFunction_O::evaluateAll( 	NVector_sp pos,
                                                        bool calcForce,
                                                        gc::Nilable<NVector_sp> force,
                                                        bool calcDiagonalHessian,
                                                        bool calcOffDiagonalHessian,
                                                        gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                        gc::Nilable<NVector_sp> hdvec,
                                                        gc::Nilable<NVector_sp> dvec	)
{_G()
  bool	hasForce = force.notnilp();
  bool   hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());

#ifdef	DEBUG_ON //[
	// Summarize entry state for debugging
  LOG(BF("calcForce = %d") % calcForce  );
  LOG(BF("calcDiagonalHessian = %d") % calcDiagonalHessian  );
  LOG(BF("calcOffDiagonalHessian = %d") % calcOffDiagonalHessian  );
  LOG(BF("hasForce = %d") % hasForce  );
  LOG(BF("hasHdAndD = %d") % hasHdAndD  );
  if ( hasForce && force->size() < pos->size() ) {
    SIMPLE_ERROR(BF("Force does not have the necessary dimensions"));
  }
#endif //]

  if ( !calcForce && ( calcDiagonalHessian || calcOffDiagonalHessian ) ) {
    SIMPLE_ERROR(BF("Inconsistant arguments: if you want to calcDiagonalHessian or calcOffDiagonalHessian you must calcForce"));
  }
  if ( !calcDiagonalHessian & calcOffDiagonalHessian ) {
    SIMPLE_ERROR(BF("Inconsistant arguments: if you want to calcOffDiagonalHessian you must calcDiagonalHessian"));
  }

////    _lisp->profiler().pushTimerStates();
  if ( hasForce ) force->zero();
  if ( hasHessian ) hessian->zero();
  if ( hasHdAndD ) {
    LOG(BF("Zeroing hdvec") );
    hdvec->zero();	// Zero the result
  }

    LOG(BF("Starting evaluation of energy") );
    double totalEnergy = 0.0;
    for ( auto cur : this->_Terms ) {
      EnergyRigidBodyComponent_sp term = gc::As<EnergyRigidBodyComponent_sp>(CONS_CAR(cur));
      totalEnergy += term->evaluateAll(pos,calcForce,force,
                                       calcDiagonalHessian,
                                       calcOffDiagonalHessian,
                                       hessian, hdvec, dvec);
    }
    return totalEnergy;
}


double	RigidBodyEnergyFunction_O::evaluateEnergy( NVector_sp pos )
{
  IMPLEMENT_ME();
};
double	RigidBodyEnergyFunction_O::evaluateEnergyForce( NVector_sp pos, bool calcForce, NVector_sp force )
{
  IMPLEMENT_ME();
};

void	RigidBodyEnergyFunction_O::dealWithProblem(core::Symbol_sp error_symbol, core::T_sp arguments)
{
  IMPLEMENT_ME();
};

void RigidBodyEnergyFunction_O::setPosition(size_t index, double a, double b, double c, double d, double x, double y, double z) {
  if (index>= this->_RigidBodies) {
    SIMPLE_ERROR(BF("set-position at index %d out of range <= %d for coordinate") % index % this->_RigidBodies);
  }
  double*  p = &(*this->_SavedCoordinates)[0];
  size_t base = index*7;
  p[base+0] = a;
  p[base+1] = b;
  p[base+2] = c;
  p[base+3] = d;
  p[base+4] = x;
  p[base+5] = y;
  p[base+6] = z;
}
core::T_mv RigidBodyEnergyFunction_O::getPosition(size_t index) {
  if (index>= this->_RigidBodies) {
    SIMPLE_ERROR(BF("set-position at index %d out of range <= %d for coordinate") % index % this->_RigidBodies);
  }
  double*  p = &(*this->_SavedCoordinates)[0];
  size_t base = index*7;
  return Values(core::clasp_make_double_float(p[base+0]),
                core::clasp_make_double_float(p[base+1]),
                core::clasp_make_double_float(p[base+2]),
                core::clasp_make_double_float(p[base+3]),
                core::clasp_make_double_float(p[base+4]),
                core::clasp_make_double_float(p[base+5]),
                core::clasp_make_double_float(p[base+6]) );
}

  


void RigidBodyEnergyFunction_O::normalizePosition()
{
  // Normalize the quaternions
  double* dpos = &(*pos)[0];
  for ( size_t i(0), iEnd(pos->length()); i<iEnd; i += 7 ) {
    double qlen = sqrt(dpos[i+0]*dpos[i+0] +
                       dpos[i+1]*dpos[i+1] +
                       dpos[i+2]*dpos[i+2] +
                       dpos[i+3]*dpos[i+3]);
    ASSERT(qlen!=0.0);
    double qinvlen = 1.0/qlen;
    dpos[i+0] *= qinvlen;
    dpos[i+1] *= qinvlen;
    dpos[i+2] *= qinvlen;
    dpos[i+3] *= qinvlen;
  }
}



};
