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
#include <clasp/core/numbers.h>
#include <cando/adapt/indexedObjectBag.h>
#include <cando/chem/minimizerLog.h>
#include <cando/chem/restraint.h>
#include <cando/chem/iterateRestraints.h>
#include <cando/adapt/iterateCons.h>
#include <cando/chem/ringFinder.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/energyRigidBodyNonbond.h>
#include <cando/chem/energyRigidBodyStaple.h>
#include <cando/chem/atom.h>
#include <cando/chem/virtualAtom.h>


#include <cando/chem/rigidBodyEnergyFunction.h>

#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/chemInfo.h>
#include <clasp/core/array.h>
#include <clasp/core/wrappers.h>


/* Nothing */

namespace chem {

CL_LISPIFY_NAME("CHEM:MAKE-RIGID-BODY-ENERGY-FUNCTION");
DOCGROUP(cando);
CL_DEFUN RigidBodyEnergyFunction_sp RigidBodyEnergyFunction_O::make(size_t number_of_rigid_bodies, BoundingBox_sp boundingBox) {
  auto  energy = gctools::GC<RigidBodyEnergyFunction_O>::allocate( number_of_rigid_bodies, boundingBox);
  energy->_SavedCoordinates = NVector_O::create(number_of_rigid_bodies*7);
  return energy;
}

string RigidBodyEnergyFunction_O::energyTermsEnabled() {
  stringstream ss;
  size_t termCount = 0;
  for ( auto cur : this->_Terms ) {
    EnergyRigidBodyComponent_sp term = gc::As<EnergyRigidBodyComponent_sp>(CONS_CAR(cur));
    if ( termCount ) ss << " ";
    ss << core::lisp_instance_class(term)->_classNameAsString();
  }
  return ss.str();
}

void RigidBodyEnergyFunction_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,RigidBodies),this->_RigidBodies);
  node->field_if_not_unbound(INTERN_(kw,SavedCoordinates),this->_SavedCoordinates);
  node->field(INTERN_(kw,Terms),this->_Terms);
  node->field_if_not_unbound(INTERN_(kw,BoundingBox),this->_BoundingBox);
  this->Base::fields(node);
}

CL_DOCSTRING(R"dx(Return the bounding-box for the atom-table.)dx");
CL_LISPIFY_NAME(rigidBodyEnergyFunction-bounding-box);
CL_DEFMETHOD BoundingBox_sp RigidBodyEnergyFunction_O::boundingBox() const
{
  if (this->_BoundingBox.boundp()) {
    return this->_BoundingBox;
  }
  SIMPLE_ERROR("The bounding-box slot is unbound");
}

CL_DOCSTRING(R"dx(Return T if the bounding-box is bound)dx");
CL_LISPIFY_NAME(rigidBodyEnergyFunction-bounding-box-bound-p);
CL_DEFMETHOD bool RigidBodyEnergyFunction_O::boundingBoxBoundP() const
{
  return this->_BoundingBox.boundp();
}

CL_DOCSTRING(R"dx(Set the bounding-box)dx");
CL_LISPIFY_NAME(rigidBodyEnergyFunction-set-bounding-box);
CL_DEFMETHOD void RigidBodyEnergyFunction_O::setBoundingBox(BoundingBox_sp boundingBox) {
  this->_BoundingBox = boundingBox;
}

CL_DOCSTRING(R"dx(Make the bounding-box unbound)dx");
CL_LISPIFY_NAME(rigidBodyEnergyFunction-mak-unbound-bounding-box);
CL_DEFMETHOD void RigidBodyEnergyFunction_O::makUnboundBoundingBox() {
  this->_BoundingBox = unbound<BoundingBox_O>();
}

CL_DEFMETHOD void RigidBodyEnergyFunction_O::set_coordinates(NVector_sp pos) {
  if (pos->length() != this->_RigidBodies*7) {
    SIMPLE_ERROR("The coordinates you pass do not have enough components - there are {} rigid bodies with 7 (a,b,c,d,x,y,z) coordinates each - total {} are needed" , this->_RigidBodies , (this->_RigidBodies*7));
  }
  this->_SavedCoordinates = pos;
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

void RigidBodyEnergyFunction_O::setupHessianPreconditioner( NVector_sp pos, AbstractLargeSquareMatrix_sp hessian,
                                                            core::T_sp activeAtomMask ) {
  IMPLEMENT_ME();
}

void	RigidBodyEnergyFunction_O::loadCoordinatesIntoVector(NVector_sp pos)
{
  ASSERT(this->_SavedCoordinates && this->_SavedCoordinates->length() == pos->length());
  core::core__copy_subarray(pos,core::make_fixnum(0),this->_SavedCoordinates,core::make_fixnum(0),core::make_fixnum(pos->length()));
}

void RigidBodyEnergyFunction_O::saveCoordinatesFromVector(NVector_sp pos)
{
  ASSERT(this->_SavedCoordinates && this->_SavedCoordinates->length() == pos->length());
  core::core__copy_subarray(this->_SavedCoordinates,core::make_fixnum(0),pos,core::make_fixnum(0),core::make_fixnum(pos->length()));
}

void RigidBodyEnergyFunction_O::saveCoordinatesAndForcesFromVectors(NVector_sp pos, NVector_sp force)
{
  ASSERT(this->_SavedCoordinates && this->_SavedCoordinates->length() == pos->length());
  this->saveCoordinatesFromVector(pos);
  // Do the forces
}


double	RigidBodyEnergyFunction_O::evaluateRaw( NVector_sp pos, NVector_sp force ) {

  IMPLEMENT_ME();
}

//    virtual double	evaluate( NVector_sp pos, NVector_sp force, bool calculateForce ) ;
  adapt::QDomNode_sp	identifyTermsBeyondThreshold();
//    uint	countBadVdwInteractions(double scaleSumOfVdwRadii, geom::DisplayList_sp displayIn);

ForceMatchReport_sp RigidBodyEnergyFunction_O::checkIfAnalyticalForceMatchesNumericalForce( NVector_sp pos, NVector_sp force, core::T_sp activeAtomMask ) {
  IMPLEMENT_ME();
}

void	RigidBodyEnergyFunction_O::useDefaultSettings()
{
  IMPLEMENT_ME();
}

    /*! Set a single options */
void	RigidBodyEnergyFunction_O::setOption( core::Symbol_sp option, core::T_sp val)
{
  if ( option == _sym_nonbondTerm)
  {
    return;
  }
  SIMPLE_ERROR("Unknown option");
}



    /*! Set the energy function options. List the options as a flat list of keyword/value pairs */
void	RigidBodyEnergyFunction_O::setOptions( core::List_sp options )
{
  IMPLEMENT_ME();
}


double	RigidBodyEnergyFunction_O::evaluateAll( NVector_sp pos,
                                                core::T_sp componentEnergy,
                                                bool calcForce,
                                                gc::Nilable<NVector_sp> force,
                                                bool calcDiagonalHessian,
                                                bool calcOffDiagonalHessian,
                                                gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                gc::Nilable<NVector_sp> hdvec,
                                                gc::Nilable<NVector_sp> dvec,
                                                core::T_sp activeAtomMask,
                                                core::T_sp debugInteractions )
{
  bool	hasForce = force.notnilp();
  bool  hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
#ifdef	DEBUG_ON //[
	// Summarize entry state for debugging
  LOG("calcForce = {}" , calcForce  );
  LOG("calcDiagonalHessian = {}" , calcDiagonalHessian  );
  LOG("calcOffDiagonalHessian = {}" , calcOffDiagonalHessian  );
  LOG("hasForce = {}" , hasForce  );
  LOG("hasHdAndD = {}" , hasHdAndD  );
  if ( hasForce && force->size() < pos->size() ) {
    SIMPLE_ERROR("Force does not have the necessary dimensions");
  }
#endif //]

  if ( !calcForce && ( calcDiagonalHessian || calcOffDiagonalHessian ) ) {
    SIMPLE_ERROR("Inconsistant arguments: if you want to calcDiagonalHessian or calcOffDiagonalHessian you must calcForce");
  }
  if ( !calcDiagonalHessian & calcOffDiagonalHessian ) {
    SIMPLE_ERROR("Inconsistant arguments: if you want to calcOffDiagonalHessian you must calcDiagonalHessian");
  }

////    _lisp->profiler().pushTimerStates();
  if ( hasForce ) {
    force->zero();
  }
  if ( hasHessian ) hessian->zero();
  if ( hasHdAndD ) {
    LOG("Zeroing hdvec" );
    hdvec->zero();	// Zero the result
  }

  LOG("Starting evaluation of energy" );
  double totalEnergy = 0.0;
  for ( auto cur : this->_Terms ) {
    EnergyRigidBodyComponent_sp term = gc::As<EnergyRigidBodyComponent_sp>(CONS_CAR(cur));
    if (term->isEnabled()) {
      totalEnergy += term->evaluateAllComponent(this->asSmartPtr(),
                                                pos,
                                                componentEnergy,
                                                calcForce,
                                                force,
                                                calcDiagonalHessian,
                                                calcOffDiagonalHessian,
                                                hessian,
                                                hdvec,
                                                dvec,
                                                activeAtomMask,
                                                nil<core::T_O>());
    }
  }
  return totalEnergy;
}



void	RigidBodyEnergyFunction_O::dealWithProblem(core::Symbol_sp error_symbol, core::T_sp arguments)
{
  IMPLEMENT_ME();
};

CL_DEFMETHOD
size_t RigidBodyEnergyFunction_O::numberOfRigidBodies() const {
  return this->_RigidBodies;
}

void RigidBodyEnergyFunction_O::setPosition(size_t index, double a, double b, double c, double d, double x, double y, double z) {
  if (index>= this->_RigidBodies) {
    SIMPLE_ERROR("set-position at index {} out of range <= {} for coordinate" , index , this->_RigidBodies);
  }
  Vector_real*  p = &(*this->_SavedCoordinates)[0];
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
    SIMPLE_ERROR("set-position at index {} out of range <= {} for coordinate" , index , this->_RigidBodies);
  }
  Vector_real*  p = &(*this->_SavedCoordinates)[0];
  size_t base = index*7;
  return Values(core::clasp_make_double_float(p[base+0]),
                core::clasp_make_double_float(p[base+1]),
                core::clasp_make_double_float(p[base+2]),
                core::clasp_make_double_float(p[base+3]),
                core::clasp_make_double_float(p[base+4]),
                core::clasp_make_double_float(p[base+5]),
                core::clasp_make_double_float(p[base+6]) );
}

  


void RigidBodyEnergyFunction_O::normalizePosition(NVector_sp pos)
{
  // Normalize the quaternions
  Vector_real* dpos = &(*pos)[0];
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

string RigidBodyEnergyFunction_O::energyComponentsAsString() {
  stringstream ss;
  for ( auto cur : this->_Terms ) {
    EnergyRigidBodyComponent_sp term = gc::As<EnergyRigidBodyComponent_sp>(CONS_CAR(cur));
    ss << _rep_(core::lisp_instance_class(term)) << " ";
  }
  return ss.str();
}
 

void RigidBodyEnergyFunction_O::dumpTerms()
{
  for ( auto cur : this->_Terms ) {
    EnergyRigidBodyComponent_sp term = gc::As<EnergyRigidBodyComponent_sp>(CONS_CAR(cur));
    term->dumpTerms(this->_AtomTypes);
  }
}


CL_DOCSTRING(R"dx(Like velocity-verlet-step)dx");
CL_DOCSTRING_LONG(R"dx(but limits displacement of rigid body quaternions and centers using the limit_displacement vector.
Return the number of rigid bodies whose displacement was limited.  The limit_displacement nvector contains seven double precision values, 
the a,b,c,d quaternion and x,y,z displacements.  The quaternion displacement idea is experimental.  I'm not sure how it will
effect the dynamics.)dx")
CL_LISPIFY_NAME("rigid-body-velocity-verlet-step-limit-displacement");
DOCGROUP(cando);
CL_DEFUN size_t chem__rigid_body_velocity_verlet_step_limit_displacement(ScoringFunction_sp scoringFunc,
                                                                         NVector_sp position,
                                                                         NVector_sp velocity,
                                                                         NVector_sp force,
                                                                         NVector_sp force_dt,
                                                                         NVector_sp delta_t_over_mass,
                                                                         double delta_t,
                                                                         core::T_sp tunfrozen,
                                                                         NVector_sp velocity_scale,
                                                                         NVector_sp limit_displacement )
{
//  SIMPLE_WARN("FIXactiveAtomMask Should tunfrozen and activeAtomMask be the same");
  core::SimpleBitVector_sp unfrozen;
  if (gc::IsA<core::SimpleBitVector_sp>(tunfrozen)) {
    unfrozen = gc::As_unsafe<core::SimpleBitVector_sp>(tunfrozen);
    if (unfrozen->length() != (position->length()/7)) {
      SIMPLE_ERROR("frozen must be a simple-bit-vector of length {} or NIL - it is {}" , (position->length()/7) , _rep_(position));
    }
  } else if (tunfrozen.notnilp()) {
    SIMPLE_ERROR("unfrozen must be a simple-bit-vector or NIL");
  }
  size_t body_idx = 0;
  size_t body_limited = 0;
  for ( size_t idx = 0; idx<position->size(); idx += 7) {
    if (!unfrozen || unfrozen->testBit(body_idx)==1) {
      double offseta = delta_t*(*velocity)[idx+0] + delta_t*(*delta_t_over_mass)[body_idx]*(*force)[idx+0];
      double offsetb = delta_t*(*velocity)[idx+1] + delta_t*(*delta_t_over_mass)[body_idx]*(*force)[idx+1];
      double offsetc = delta_t*(*velocity)[idx+2] + delta_t*(*delta_t_over_mass)[body_idx]*(*force)[idx+2];
      double offsetd = delta_t*(*velocity)[idx+3] + delta_t*(*delta_t_over_mass)[body_idx]*(*force)[idx+3];
      double offsetx = delta_t*(*velocity)[idx+4] + delta_t*(*delta_t_over_mass)[body_idx]*(*force)[idx+4];
      double offsety = delta_t*(*velocity)[idx+5] + delta_t*(*delta_t_over_mass)[body_idx]*(*force)[idx+5];
      double offsetz = delta_t*(*velocity)[idx+6] + delta_t*(*delta_t_over_mass)[body_idx]*(*force)[idx+6];
      if ((offseta>(*limit_displacement)[0]) ||
          (offsetb>(*limit_displacement)[1]) ||
          (offsetc>(*limit_displacement)[2]) ||
          (offsetd>(*limit_displacement)[3]) ||
          (offsetx>(*limit_displacement)[4]) ||
          (offsety>(*limit_displacement)[5]) ||
          (offsetz>(*limit_displacement)[6])) {
        body_limited++;
        if (offseta>(*limit_displacement)[0]) offseta = (*limit_displacement)[0];
        if (offsetb>(*limit_displacement)[1]) offsetb = (*limit_displacement)[1];
        if (offsetc>(*limit_displacement)[2]) offsetc = (*limit_displacement)[2];
        if (offsetd>(*limit_displacement)[3]) offsetd = (*limit_displacement)[3];
        if (offsetx>(*limit_displacement)[4]) offsetx = (*limit_displacement)[4];
        if (offsety>(*limit_displacement)[5]) offsety = (*limit_displacement)[5];
        if (offsetz>(*limit_displacement)[6]) offsetz = (*limit_displacement)[6];
        if (offseta<(*limit_displacement)[0]) offseta = -(*limit_displacement)[0];
        if (offsetb<(*limit_displacement)[1]) offsetb = -(*limit_displacement)[1];
        if (offsetc<(*limit_displacement)[2]) offsetc = -(*limit_displacement)[2];
        if (offsetd<(*limit_displacement)[3]) offsetd = -(*limit_displacement)[3];
        if (offsetx<(*limit_displacement)[4]) offsetx = -(*limit_displacement)[4];
        if (offsety<(*limit_displacement)[5]) offsety = -(*limit_displacement)[5];
        if (offsetz<(*limit_displacement)[6]) offsetz = -(*limit_displacement)[6];
      }
      (*position)[idx+0] = (*position)[idx+0] + offseta;
      (*position)[idx+1] = (*position)[idx+1] + offsetb;
      (*position)[idx+2] = (*position)[idx+2] + offsetc;
      (*position)[idx+3] = (*position)[idx+3] + offsetd;
      (*position)[idx+4] = (*position)[idx+4] + offsetx;
      (*position)[idx+5] = (*position)[idx+5] + offsety;
      (*position)[idx+6] = (*position)[idx+6] + offsetz;
    }
    body_idx++;
  }
  scoringFunc->evaluateEnergyForce(position,true,force_dt,tunfrozen);
  body_idx = 0;
  for ( size_t idx = 0; idx<position->size(); idx+=7 ) {
    if (!unfrozen || unfrozen->testBit(body_idx)==1) {
      (*velocity)[idx+0] = ((*velocity)[idx+0] + (*delta_t_over_mass)[body_idx]*0.5*((*force)[idx+0]+(*force_dt)[idx+0]))*(*velocity_scale)[0];
      (*velocity)[idx+1] = ((*velocity)[idx+1] + (*delta_t_over_mass)[body_idx]*0.5*((*force)[idx+1]+(*force_dt)[idx+1]))*(*velocity_scale)[1];
      (*velocity)[idx+2] = ((*velocity)[idx+2] + (*delta_t_over_mass)[body_idx]*0.5*((*force)[idx+2]+(*force_dt)[idx+2]))*(*velocity_scale)[2];
      (*velocity)[idx+3] = ((*velocity)[idx+3] + (*delta_t_over_mass)[body_idx]*0.5*((*force)[idx+3]+(*force_dt)[idx+3]))*(*velocity_scale)[3];
      (*velocity)[idx+4] = ((*velocity)[idx+4] + (*delta_t_over_mass)[body_idx]*0.5*((*force)[idx+4]+(*force_dt)[idx+4]))*(*velocity_scale)[4];
      (*velocity)[idx+5] = ((*velocity)[idx+5] + (*delta_t_over_mass)[body_idx]*0.5*((*force)[idx+5]+(*force_dt)[idx+5]))*(*velocity_scale)[5];
      (*velocity)[idx+6] = ((*velocity)[idx+6] + (*delta_t_over_mass)[body_idx]*0.5*((*force)[idx+6]+(*force_dt)[idx+6]))*(*velocity_scale)[6];
      (*force)[idx+0] = (*force_dt)[idx+0];
      (*force)[idx+1] = (*force_dt)[idx+1];
      (*force)[idx+2] = (*force_dt)[idx+2];
      (*force)[idx+3] = (*force_dt)[idx+3];
      (*force)[idx+4] = (*force_dt)[idx+4];
      (*force)[idx+5] = (*force_dt)[idx+5];
      (*force)[idx+6] = (*force_dt)[idx+6];
    }
    body_idx++;
  }
  return body_limited;
}


CL_DEFMETHOD Aggregate_sp RigidBodyEnergyFunction_O::buildPseudoAggregate(NVector_sp pos) const {
  core::SimpleVector_float_sp coords = this->buildPseudoAggregateCoordinates(pos);
  size_t idx = 0;
  EnergyRigidBodyNonbond_sp nonbondTerm = unbound<EnergyRigidBodyNonbond_O>();
  EnergyRigidBodyStaple_sp stapleTerm = unbound<EnergyRigidBodyStaple_O>();
  for ( auto cur : this->_Terms ) {
    core::T_sp term = CONS_CAR(cur);
    if (gc::IsA<EnergyRigidBodyNonbond_sp>(term)) {
      nonbondTerm = gc::As<EnergyRigidBodyNonbond_sp>(term);
    } else if (gc::IsA<EnergyRigidBodyStaple_sp>(term)) {
      stapleTerm = gc::As<EnergyRigidBodyStaple_sp>(term);
    }
  }
  Aggregate_sp aggregate = Aggregate_O::make(nil<core::T_O>());
  if (nonbondTerm.boundp()) {
    size_t prevAtomInfoEnd = 0;
    for ( size_t nbodyi = 0; nbodyi < nonbondTerm->_RigidBodyEndAtom->length(); nbodyi++ ) {
      Molecule_sp molecule = Molecule_O::make(nil<core::T_O>());
      Residue_sp residue = Residue_O::make(nil<core::T_O>());
      Atom_sp prevAtom = unbound<Atom_O>();
      aggregate->addMatter(molecule);
      molecule->addMatter(residue);
      size_t atomInfoStart = prevAtomInfoEnd;
      size_t atomInfoEnd = (*nonbondTerm->_RigidBodyEndAtom)[nbodyi];
      for ( size_t ai = atomInfoStart; ai<atomInfoEnd; ai++ ) {
        RigidBodyAtomInfo& atomInfo = nonbondTerm->_AtomInfoTable[ai];
        Atom_sp part_atom;
        if (gc::IsA<Atom_sp>(atomInfo._Object)) {
          Atom_sp atm = gc::As_unsafe<Atom_sp>(atomInfo._Object);
          part_atom = Atom_O::make(atm->getName(),atm->getElement());
        } else {
          part_atom = Atom_O::make(nil<core::T_O>(),element_C);
        }
        float x = (*coords)[idx++];
        float y = (*coords)[idx++];
        float z = (*coords)[idx++];
        Vector3 pos( x, y, z );
        part_atom->setPosition(pos);
        residue->addMatter(part_atom);
        if (prevAtom.boundp()) {
          prevAtom->bondTo(part_atom,singleBond);
        }
        prevAtom = part_atom;
      }
      prevAtomInfoEnd = atomInfoEnd;
    }
  }
  if (stapleTerm.boundp()) {
    for ( size_t ii = 0; ii < stapleTerm->_Terms.size(); ii++ ) {
      Molecule_sp molecule = Molecule_O::make(nil<core::T_O>());
      Residue_sp residue = Residue_O::make(nil<core::T_O>());
      aggregate->addMatter(molecule);
      molecule->addMatter(residue);
      Atom_sp atomK = Atom_O::make(nil<core::T_O>(),element_N);
      float xK = (*coords)[idx++];
      float yK = (*coords)[idx++];
      float zK = (*coords)[idx++];
      Vector3 posK( xK, yK, zK );
      atomK->setPosition(posK);
      Atom_sp atomL = Atom_O::make(nil<core::T_O>(),element_N);
      float xL = (*coords)[idx++];
      float yL = (*coords)[idx++];
      float zL = (*coords)[idx++];
      Vector3 posL( xL, yL, zL );
      atomL->setPosition(posL);
      atomK->bondTo(atomL,singleBond);
      residue->addMatter(atomK);
      residue->addMatter(atomL);
    }
  }
  return aggregate;
}

CL_DEFMETHOD core::SimpleVector_float_sp RigidBodyEnergyFunction_O::buildPseudoAggregateCoordinates(NVector_sp pos) const {
  core::SimpleVector_float_sp coords = core::SimpleVector_float_O::make(this->numberOfPoints()*3);
  EnergyRigidBodyNonbond_sp nonbondTerm = unbound<EnergyRigidBodyNonbond_O>();
  EnergyRigidBodyStaple_sp stapleTerm = unbound<EnergyRigidBodyStaple_O>();
  for ( auto cur : this->_Terms ) {
    core::T_sp term = CONS_CAR(cur);
    if (gc::IsA<EnergyRigidBodyNonbond_sp>(term)) {
      nonbondTerm = gc::As<EnergyRigidBodyNonbond_sp>(term);
    } else if (gc::IsA<EnergyRigidBodyStaple_sp>(term)) {
      stapleTerm = gc::As<EnergyRigidBodyStaple_sp>(term);
    }
  }
  size_t idx = 0;
  if (nonbondTerm.boundp()) {
    idx = nonbondTerm->partsCoordinates(pos,idx,coords);
  }
  if (stapleTerm.boundp()) {
    idx = stapleTerm->partsCoordinates(pos,idx,coords);
  }
  return coords;
}


CL_DEFMETHOD size_t RigidBodyEnergyFunction_O::numberOfPoints() const {
  EnergyRigidBodyNonbond_sp nonbondTerm = unbound<EnergyRigidBodyNonbond_O>();
  EnergyRigidBodyStaple_sp stapleTerm = unbound<EnergyRigidBodyStaple_O>();
  for ( auto cur : this->_Terms ) {
    core::T_sp term = CONS_CAR(cur);
    if (gc::IsA<EnergyRigidBodyNonbond_sp>(term)) {
      nonbondTerm = gc::As<EnergyRigidBodyNonbond_sp>(term);
    } else if (gc::IsA<EnergyRigidBodyStaple_sp>(term)) {
      stapleTerm = gc::As<EnergyRigidBodyStaple_sp>(term);
    }
  }
  size_t numberOfPoints = 0;
  if (nonbondTerm.boundp()) {
    numberOfPoints += nonbondTerm->_AtomInfoTable.size();
  }
  if (stapleTerm.boundp()) {
    numberOfPoints += 2* stapleTerm->_Terms.size();
  }
  return numberOfPoints;
}
  


#if 0
uint RigidBodyEnergyFunction_O::checkForBeyondThresholdInteractions()
{
  return 0;
}
#endif

};
