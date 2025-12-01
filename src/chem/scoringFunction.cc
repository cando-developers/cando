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
#include <clasp/core/evaluator.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/loop.h>
#include <cando/adapt/indexedObjectBag.h>
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


namespace chem {

double energyScaleDielectricConstant(core::T_sp energyScale) {
  if (gc::IsA<EnergyScale_sp>(energyScale)) {
    return gc::As_unsafe<EnergyScale_sp>(energyScale)->getDielectricConstant();
  }
    return 1.0;
}

double energyScaleVdwScale(core::T_sp energyScale) {
  if (gc::IsA<EnergyScale_sp>(energyScale)) {
    return gc::As_unsafe<EnergyScale_sp>(energyScale)->getVdwScale();
  }
    return 1.0;
}

double energyScaleElectrostaticScale(core::T_sp energyScale) {
  if (gc::IsA<EnergyScale_sp>(energyScale)) {
    return gc::As_unsafe<EnergyScale_sp>(energyScale)->getElectrostaticScale();
  }
    return 1.0;
}

double energyScaleNonbondCutoff(core::T_sp energyScale) {
  if (gc::IsA<EnergyScale_sp>(energyScale)) {
    return gc::As_unsafe<EnergyScale_sp>(energyScale)->getNonbondCutoff();
  }
    return 1.0;
}

};

namespace chem
{

// set_energy is a generic function for setting energies
SYMBOL_EXPORT_SC_(ChemPkg,set_energy);


ScoringFunction_O::ScoringFunction_O() : _VelocityScale(Vector3(1.0,1.0,1.0)), _AtomTypes(core::HashTable_O::createEq()) {};

CL_DEFMETHOD void ScoringFunction_O::setf_velocity_scale(double xscale, double yscale, double zscale) {
  this->_VelocityScale.set(xscale,yscale,zscale);
}

void ScoringFunction_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,Debug),this->_Debug);
  node->field_if_not_default(INTERN_(kw,VelocityScale),this->_VelocityScale,Vector3());
  node->field_if_not_unbound(INTERN_(kw,AtomTypes),this->_AtomTypes);
}

CL_DEFMETHOD 
core::HashTable_sp ScoringFunction_O::atomTypes() const {
  if (this->_AtomTypes.boundp()) return this->_AtomTypes;
  SIMPLE_ERROR("The atom-types are unbound");
}

CL_DEFMETHOD void ScoringFunction_O::setAtomTypes(core::HashTable_sp atomTypes) {
  this->_AtomTypes = atomTypes;
}

CL_LISPIFY_NAME("evaluateEnergy");
CL_LAMBDA((scoring-function chem:scoring-function) positions &key energy-scale component-energy active-atom-mask debug-interactions disable-restraints);
CL_DEFMETHOD double	ScoringFunction_O::evaluateEnergy( NVector_sp pos,
                                                           core::T_sp energyScale,
                                                           core::T_sp componentEnergy,
                                                           core::T_sp activeAtomMask,
                                                           core::T_sp debugInteractions,
                                                           bool disableRestraints)
{
  double		energy;
  energy = this->evaluateAll(pos,
                             energyScale,
                             componentEnergy,
                             false,
                             nil<core::T_O>(),
                             false, false,
                             nil<core::T_O>(),
                             nil<core::T_O>(),
                             nil<core::T_O>(),
                             activeAtomMask,
                             debugInteractions,
                             disableRestraints);
  return energy;
}



CL_LISPIFY_NAME("evaluateEnergyForce");
CL_LAMBDA((scoring-function chem:scoring-function) positions &key energy-scale calc-force force active-atom-mask debug-interactions disable-restraints);
CL_DEFMETHOD double	ScoringFunction_O::evaluateEnergyForce( NVector_sp pos,
                                                                core::T_sp energyScale,
                                                                bool calcForce,
                                                                NVector_sp force,
                                                                core::T_sp activeAtomMask,
                                                                core::T_sp debugInteractions,
                                                                bool disableRestraints )
{
  double	energy;
  gc::Nilable<NVector_sp>	rawGrad;
  if ( calcForce ) {
    rawGrad = force;
  } else {
    rawGrad = nil<core::T_O>();
  }
  energy = this->evaluateAll(pos,
                             energyScale,
                             nil<core::T_O>(),
                             calcForce,
                             rawGrad,
                             false,
                             false,
                             nil<core::T_O>(),
                             nil<core::T_O>(),
                             nil<core::T_O>(),
                             activeAtomMask,
                             debugInteractions,
                             disableRestraints );
  return energy;
}

CL_DEFMETHOD NVector_sp ScoringFunction_O::makeCoordinates() const
{
  NVector_sp pos = NVector_O::make(this->getNVectorSize(),0.0,true);
  return pos;
};


CL_LISPIFY_NAME("evaluateEnergyForceFullHessian");
CL_LAMBDA((scoring-function chem:scoring-function) positions &key energy-scale calc-force force calc-diagonal-hessian calc-off-diagonal-hessian hessian active-atom-mask debug-interactions disable-restraints);
CL_DEFMETHOD double	ScoringFunction_O::evaluateEnergyForceFullHessian(
                                                                          NVector_sp pos,
                                                                          core::T_sp energyScale,
                                                                          bool calcForce, NVector_sp force,
                                                                          bool calcDiagonalHessian,
                                                                          bool calcOffDiagonalHessian,
                                                                          AbstractLargeSquareMatrix_sp hessian,
                                                                          core::T_sp activeAtomMask,
                                                                          core::T_sp debugInteractions,
                                                                          bool disableRestraints )
{
  double	energy;
  gc::Nilable<NVector_sp> rawGrad;
  if ( calcForce ) rawGrad = force;
  else rawGrad = nil<core::T_O>();
  energy = this->evaluateAll( pos,
                              energyScale,
                              nil<core::T_O>(),
                              calcForce, rawGrad,
                              calcDiagonalHessian,
                              calcOffDiagonalHessian,
                              hessian,
                              nil<core::T_O>(),
                              nil<core::T_O>(),
                              activeAtomMask,
                              debugInteractions,
                              disableRestraints );
  return energy;
}


#if 0
CL_LISPIFY_NAME("evaluateEnergyForceFullHessianForDebugging");
CL_LAMBDA((scoring-function chem:scoring-function) &key energy-scale active-atom-mask debug-interactions disable-restraints);
CL_DEFMETHOD double	ScoringFunction_O::evaluateEnergyForceFullHessianForDebugging(core::T_sp energyScale,
                                                                                      core::T_sp activeAtomMask,
                                                                                      core::T_sp debugInteractions,
                                                                                      bool disableRestraints )
{
  NVector_sp	pos, force;
  AbstractLargeSquareMatrix_sp	hessian;
  double		energy;
  pos = NVector_O::create(this->getNVectorSize());
  force = NVector_O::create(this->getNVectorSize());
  hessian = FullLargeSquareMatrix_O::create(this->getNVectorSize(),SymmetricDiagonalLower);
  this->loadCoordinatesIntoVector(pos);
  energy = this->evaluateAll(pos,
                             energyScale,
                             nil<core::T_O>(),
                             true,
                             force,
                             true,
                             true,
                             hessian,
                             nil<core::T_O>(),
                             nil<core::T_O>(),
                             activeAtomMask,
                             debugInteractions,
                             disableRestraints );
  return energy;
}
#endif



#if 0
CL_LISPIFY_NAME("calculateEnergy");
CL_LAMBDA((scoring-function chem:scoring-function) &optional active-atom-mask debug-interactions disable-restraints);
CL_DEFMETHOD double	ScoringFunction_O::calculateEnergy(core::T_sp activeAtomMask, core::T_sp debugInteractions, bool disableRestraints )
{
  NVector_sp	pos;
  pos = NVector_O::create(this->getNVectorSize());
  this->loadCoordinatesIntoVector(pos);
  return this->evaluateEnergy(pos,activeAtomMask,debugInteractions,disableRestraints);
}


CL_LISPIFY_NAME("calculateEnergyAndForce");
CL_LAMBDA((scoring-function chem:scoring-function) &key energy-scale active-atom-mask debug-interactions disable-restraints);
CL_DEFMETHOD core::T_mv ScoringFunction_O::calculateEnergyAndForce( core::T_sp energyScale, core::T_sp activeAtomMask, core::T_sp debugInteractions, bool disableRestraints )
{
  NVector_sp	pos;
  NVector_sp	force;
  double		energy;
  pos = NVector_O::create(this->getNVectorSize());
  force = NVector_O::create(this->getNVectorSize());
  this->loadCoordinatesIntoVector(pos);
  energy = this->evaluateEnergyForce(pos,energyScale,true,force,activeAtomMask,debugInteractions,disableRestraints);
    	// To calculate the force magnitude use force->magnitude();
    	// To calculate the force rmsMagnitude use force->rmsMagnitude();
//  this->writeForceToAtoms(force);
  return Values(core::clasp_make_double_float(energy),force,pos);
}
#endif

CL_LAMBDA((scoring-function chem:scoring-function) coords &key energy-scale force (delta 0.00001) active-atom-mask);
CL_DEFMETHOD void ScoringFunction_O::evaluateFiniteDifferenceForce(NVector_sp coords, core::T_sp energyScale, NVector_sp force, double delta, core::T_sp activeAtomMask ) {
{
  double x, ylow, yhigh, grad;
  double	deltaDiv2 = delta/2.0;
//  NVector_sp tpos = NVector_O::make(coords->size(),0.0,false,coords->size(),(Vector_real*)&(*coords)[0]);
  for ( size_t ii=0; ii<force->size(); ii++ ) {
    NVector_sp tpos = NVector_O::make(coords->size(),0.0,false,coords->size(),(Vector_real*)&(*coords)[0]);
    x = tpos->element(ii);
    tpos->setElement(ii,x-deltaDiv2);
    ylow = this->evaluateEnergy(tpos,energyScale,nil<core::T_O>(),activeAtomMask);
    tpos->setElement(ii,x+deltaDiv2);
    yhigh = this->evaluateEnergy(tpos,energyScale,nil<core::T_O>(),activeAtomMask);
    //tpos->setElement(ii,x); // restore
    grad = (yhigh-ylow)/delta;
    force->setElement(ii, -grad);
  }
}


  
}

CL_LISPIFY_NAME("checkForBeyondThresholdInteractions");
CL_DEFMETHOD core::List_sp	ScoringFunction_O::checkForBeyondThresholdInteractions(double threshold)
{
  SUBIMP();
}



CL_DOCSTRING(R"dx(Velocity-verlet-step moves the atoms one velocity-verlet-step.
If tunfrozen is a simple-bit-vector then it is used to determine unfrozen atoms.)dx")
CL_LISPIFY_NAME("velocity-verlet-step");
DOCGROUP(cando);
CL_DEFUN void chem__velocity_verlet_step(ScoringFunction_sp scoringFunc,
                                         NVector_sp position,
                                         core::T_sp energyScale,
                                         NVector_sp velocity,
                                         NVector_sp force,
                                         NVector_sp force_dt,
                                         NVector_sp delta_t_over_mass,
                                         double delta_t,
                                         core::T_sp tunfrozen)
{
//  SIMPLE_WARN("FIXactiveAtomMask is tunfrozen and activeAtomMask the same?");
  core::SimpleBitVector_sp unfrozen;
  if (gc::IsA<core::SimpleBitVector_sp>(tunfrozen)) {
    unfrozen = gc::As_unsafe<core::SimpleBitVector_sp>(tunfrozen);
    if (unfrozen->length() != (position->length()/3)) {
      SIMPLE_ERROR("unfrozen must be a simple-bit-vector of length {} or NIL - got {}" , (position->length()/3) , _rep_(tunfrozen));
    }
  } else if (tunfrozen.notnilp()) {
    SIMPLE_ERROR("unfrozen must be a simple-bit-vector or NIL");
  }
  size_t atom_idx = 0;
  if (unfrozen) {
    {
      for ( size_t idx = 0; idx<position->size(); idx += 3 ) {
        if (unfrozen->testBit(atom_idx)==1) {
          double offsetx = delta_t*(*velocity)[idx+0] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+0];
          (*position)[idx+0] = (*position)[idx+0] + offsetx;
          double offsety = delta_t*(*velocity)[idx+1] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+1];
          (*position)[idx+1] = (*position)[idx+1] + offsety;
          double offsetz = delta_t*(*velocity)[idx+2] + delta_t*(*delta_t_over_mass)[atom_idx]*(*force)[idx+2];
          (*position)[idx+2] = (*position)[idx+2] + offsetz;
        }
        atom_idx++;
      }
    }
    scoringFunc->evaluateEnergyForce(position,energyScale,true,force_dt,tunfrozen);
    atom_idx = 0;
    for ( size_t idx = 0; idx<position->size(); idx+=3 ) {
      if (unfrozen->testBit(atom_idx)==1) {
        (*velocity)[idx+0] = ((*velocity)[idx+0] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+0]+(*force_dt)[idx+0]))*scoringFunc->_VelocityScale.getX();
        (*force)[idx+0] = (*force_dt)[idx+0];
        (*velocity)[idx+1] = ((*velocity)[idx+1] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+1]+(*force_dt)[idx+1]))*scoringFunc->_VelocityScale.getY();
        (*force)[idx+1] = (*force_dt)[idx+1];
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
    scoringFunc->evaluateEnergyForce(position,energyScale,true,force_dt,tunfrozen);
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
 


CL_DOCSTRING(R"dx(Like velocity-verlet-step but limits displacement of atoms in the x,y,z directions using the limit_displacement vector)dx");
CL_DOCSTRING_LONG(R"dx(Return the number of atoms whose displacement was limited.)dx");
CL_LISPIFY_NAME("velocity-verlet-step-limit-displacement");
DOCGROUP(cando);
CL_DEFUN size_t chem__velocity_verlet_step_limit_displacement(ScoringFunction_sp scoringFunc,
                                                              NVector_sp position,
                                                              core::T_sp energyScale,
                                                              NVector_sp velocity,
                                                              NVector_sp force,
                                                              NVector_sp force_dt,
                                                              NVector_sp delta_t_over_mass,
                                                              double delta_t,
                                                              core::T_sp tunfrozen,
                                                              const Vector3& limit_displacement )
{
//  SIMPLE_WARN("FIXactiveAtomMask Is tunfrozen and activeAtomMask the same");
  core::SimpleBitVector_sp unfrozen;
  if (gc::IsA<core::SimpleBitVector_sp>(tunfrozen)) {
    unfrozen = gc::As_unsafe<core::SimpleBitVector_sp>(tunfrozen);
    if (unfrozen->length() != (position->length()/3)) {
      SIMPLE_ERROR("unfrozen must be a simple-bit-vector of length {} or NIL - got {}" , (position->length()/3) , _rep_(tunfrozen));
    }
  } else if (tunfrozen.notnilp()) {
    SIMPLE_ERROR("unfrozen must be a simple-bit-vector or NIL");
  }
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
    if (!unfrozen || unfrozen->testBit(atom_idx  )==1) {
      (*position)[idx+0] = (*position)[idx+0] + offsetx;
      (*position)[idx+1] = (*position)[idx+1] + offsety;
      (*position)[idx+2] = (*position)[idx+2] + offsetz;
    }
    atom_idx++;
  }
  scoringFunc->evaluateEnergyForce(position,energyScale,true,force_dt,tunfrozen);
  atom_idx = 0;
  for ( size_t idx = 0; idx<position->size(); idx+=3 ) {
    if (!unfrozen || unfrozen->testBit(atom_idx)==1) {
      (*velocity)[idx+0] = ((*velocity)[idx+0] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+0]+(*force_dt)[idx+0]))*scoringFunc->_VelocityScale.getX();
      (*force)[idx+0] = (*force_dt)[idx+0];
      (*velocity)[idx+1] = ((*velocity)[idx+1] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+1]+(*force_dt)[idx+1]))*scoringFunc->_VelocityScale.getY();
      (*force)[idx+1] = (*force_dt)[idx+1];
      (*velocity)[idx+2] = ((*velocity)[idx+2] + (*delta_t_over_mass)[atom_idx]*0.5*((*force)[idx+2]+(*force_dt)[idx+2]))*scoringFunc->_VelocityScale.getZ();
      (*force)[idx+2] = (*force_dt)[idx+2];
    }
    atom_idx++;
  }
  return atoms_limited;
}


DOCGROUP(cando);
CL_DEFUN core::List_sp chem__all_components(ScoringFunction_sp scoring_function) {
  return scoring_function->allComponents();
}

CL_DEFUN EnergyComponents_sp chem__make_energy_components() {
  auto obj = gctools::GC<EnergyComponents_O>::allocate_with_default_constructor();
  return obj;
}


void EnergyComponents_O::setEnergy( core::T_sp name, double energy ) {
  core::Cons_sp pair = core::Cons_O::create(name, mk_double_float(energy));
  this->_Components = core::Cons_O::create(pair,this->_Components);
}


void maybeSetEnergy( core::T_sp componentEnergy, core::T_sp name, double energy ) {
  if (componentEnergy.notnilp()) {
    gc::As<EnergyComponents_sp>(componentEnergy)->setEnergy( name, energy );
  }
}

string EnergyComponents_O::__repr__() const {
  stringstream ss;
  ss << "#<EnergyComponents ";
  for (auto cur : (core::List_sp)this->_Components ) {
    core::Cons_sp entry = gc::As<core::Cons_sp>(CONS_CAR(cur));
    ss << _rep_(CONS_CAR(entry)) << " " << _rep_(CONS_CDR(entry)) << "\n";
  }
  ss << ">";
  return ss.str();
}
};
