/*
    File: bondedJoint.cc
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
#define DEBUG_LEVEL_FULL

#define DEBUG_UPDATEXYZCOORDS

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/lisp.h>
#include <cando/geom/matrix.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/atomId.h>
#include <clasp/core/numerics.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/jumpJoint.h>
#include <cando/kinematics/bondedJoint.h>
#include <cando/kinematics/xyzJoint.h>

namespace kinematics
{
void BondedJoint_O::fields(core::Record_sp record) {
  static_assert(BondedJoint_O::MaxChildren==5,"BondedJoint_O::MaxChildren has changed from 5 - update the code below");
  record->field_if_not_unbound(INTERN_(kw,child4),this->_Children[4]);
  record->field_if_not_unbound(INTERN_(kw,child3),this->_Children[3]);
  record->field_if_not_unbound(INTERN_(kw,child2),this->_Children[2]);
  record->field_if_not_unbound(INTERN_(kw,child1),this->_Children[1]);
  record->field_if_not_unbound(INTERN_(kw,child0),this->_Children[0]);
  record->field(INTERN_(kw,num_children),this->_NumberOfChildren);
  switch (record->stage()) {
  case core::Record_O::initializing:
  case core::Record_O::loading:
    {
      double thetaDegrees;
      double phiDegrees;
      record->field(INTERN_(kw,phi_deg),phiDegrees);
      record->field(INTERN_(kw,theta_deg),thetaDegrees);
      this->_Phi = phiDegrees*0.0174533;
      this->_Theta = thetaDegrees*0.0174533;
    }
    break;
  case core::Record_O::saving: {
    double thetaDegrees = this->_Theta/0.0174533;
    double phiDegrees = this->_Phi/0.0174533;
    record->field(INTERN_(kw,phi_deg),phiDegrees);
    record->field(INTERN_(kw,theta_deg),thetaDegrees);
  }
      break;
  default:
      break;
  }
  record->field(INTERN_(kw,distance),this->_Distance);
  record->field(INTERN_(kw,dof_propagates),this->_DofChangePropagatesToYoungerSiblings);
  this->Base::fields(record);
}


void BondedJoint_O::initialize() {
  for (int i=0; i<MaxChildren; ++i ) {
    this->_Children[i] = unbound<Joint_O>();
  }
}

	/*! Return the stubJoint1 = point c */
CL_DEFMETHOD Joint_sp BondedJoint_O::inputStubJoint0() const { return this->parent();}

	/*! Return the stubJoint2 = point b */
CL_DEFMETHOD Joint_sp BondedJoint_O::inputStubJoint1() const { return this->parent()->parent();};

	/*! Return the stubJoint3 = point a */
CL_DEFMETHOD Joint_sp BondedJoint_O::inputStubJoint2() const { return this->parent()->parent()->parent(); };


CL_LAMBDA(atom-id name atom-table);
CL_LISPIFY_NAME("make_BondedJoint");
CL_DEF_CLASS_METHOD
BondedJoint_sp BondedJoint_O::make(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable) {
  return gctools::GC<BondedJoint_O>::allocate(atomId,name,atomTable);
}

bool BondedJoint_O::definedp() const {
  return (!std::isnan(this->_Phi));
}

void BondedJoint_O::_appendChild(Joint_sp c)
{
  size_t index = this->_NumberOfChildren++;
  LOG(" Appending to node {} child {} at index {}\n" , _rep_(this->asSmartPtr()) , _rep_(c) , index);
  this->_Children[index] = c;
}
  
void BondedJoint_O::_insertChild(int before, Joint_sp child)
{
  if ( this->_numberOfChildren()>this->_maxNumberOfChildren() )
  {
    THROW_HARD_ERROR("You exceeded the maximum[%d] number of children allowed for a BondedJoint" , this->_maxNumberOfChildren());
  }
  for ( int i=this->_numberOfChildren(); i>before; i-- )
    this->_Children[i] = this->_Children[i-1];
  this->_Children[before] = child;
  this->_NumberOfChildren++;
  child->_Parent = this->asSmartPtr();
}

void BondedJoint_O::_releaseChild(int idx)
{
  if ( this->_numberOfChildren() == 0 )
  {
    THROW_HARD_ERROR(("There are no children to delete"));
  }
  int num = this->_numberOfChildren() - 1;
  for ( int i=idx; i < num; i++ )
  {
    this->_Children[i] = this->_Children[i+1];
  }
  this->_NumberOfChildren--;
  this->_Children[this->_NumberOfChildren] = unbound<Joint_O>();
}

void BondedJoint_O::_releaseAllChildren()
{
  if ( this->_numberOfChildren() == 0 ) return;
  int num = this->_numberOfChildren();
  for ( int idx=0; idx < num; idx++ )
  {
    this->_Children[idx] = unbound<Joint_O>();
  }
  this->_NumberOfChildren = 0;
}

void BondedJoint_O::_updateInternalCoord(chem::NVector_sp coords)
{
  if (!this->position(coords).isDefined()) return;
  KIN_LOG(" <<< {}\n" , _rep_(this->asSmartPtr()));
//	using numeric::x_rotation_matrix_radians;
//	using numerioc::z_rotation_matrix_radians;
//	using numeric::constants::d::pi;
  Joint_sp jC = this->parent();
  Vector3 C = jC->position(coords);
  // If C is not defined then we simply return
  if (!C.isDefined()) return;
  KIN_LOG("Calculated _Distance = {}\n" , this->_Distance );
  if (gc::IsA<BondedJoint_sp>(jC) && gc::IsA<BondedJoint_sp>(jC->parent())) { // don't move past JumpJoint_O nodes
    auto distance = geom::calculateDistance(this->position(coords),C);
    KIN_LOG("!gc::IsA<JumpJoint_sp>(jC)   jC = {}\n" , _rep_(jC));
    Joint_sp jB = jC->parent();
    Vector3 B = jB->position(coords);
    if (!B.isDefined()) {
      SIMPLE_ERROR("The B joint {} has undefined position", _rep_(jB));
    }
    auto theta = geom::calculateAngle(this->position(coords),C,B); // Must be from incoming direction
    KIN_LOG("_Theta = {}\n" , (this->_Theta/0.0174533));
    Joint_sp jA = jB->parent();
    Vector3 A = jA->position(coords);
    if (!A.isDefined()) {
      SIMPLE_ERROR("The A joint {} has undefined position", _rep_(jB));
    }
    double phi = geom::calculateDihedral(this->position(coords),C,B,A);
    this->_Distance = distance;
    this->_Theta = theta;
    this->setPhi(phi);
    KIN_LOG("_Phi = {}\n", (this->_Phi/0.0174533));
    return;
  }
  Vector3 vec = this->position(coords);
  if (!vec.isDefined()) return;
  geom::internalCoordinatesFromPointAndStub(vec,
                                            this->getInputStub(coords)._Transform,
                                            this->_Distance, this->_Theta, this->_Phi );
}

bool BondedJoint_O::keepDofFixed(DofType dof) const
{
  IMPLEMENT_ME();
#if 0
  if ( dof == DofType::distance ) {
    return false;
  } else if ( dof == DofType::theta ) {
    return ( gc::IsA<JumpJoint_sp>(this->parent())
             && this->id() == this->parent()->stubJoint2Id() );
  } else if ( dof == DofType::phi ) {
    Joint_sp parent = this->parent();
    if ( gc::IsA<JumpJoint_sp>(parent) &&
         ( this->id() == parent->stubJoint2Id() ||
           this->id() == parent->stubJoint3Id(at) ) ) return true;
    if ( parent->parent().unboundp() ) return false;
    Joint_sp grandParent = parent->parent();
    return ( gc::IsA<JumpJoint_sp>(grandParent)
             && this->id() == grandParent->stubJoint3Id(at));
  } else {
    SIMPLE_ERROR("BondedJoint_O::keepDofFixed: BAD_DOF: {}" , dof.asString() );
  }
  return false;
#endif
}


string BondedJoint_O::asString() const
{
  stringstream ss;
  ss << this->Joint_O::asString();
  ss << fmt::format("  _Distance[{}]  _Theta[{}]/deg  _Phi[{}]/deg"
                     , this->_Distance
                     , (this->_Theta/0.0174533)
                     , (this->_Phi/0.0174533) ) << std::endl;
  return ss.str();
}

CL_DEFMETHOD void BondedJoint_O::setPhi(double dihedral) {
  if (isnan(dihedral)) {
    SIMPLE_ERROR("Hit NAN dihedral");
  }
  this->_Phi = dihedral;
}

/*! There are three possible situations
    C is always a BondedJoint_sp
    (1)    <BondedJoint_sp>(B) <BondedJoint_sp>(A)
    (2)    <BondedJoint_sp>(B) <JumpJoint_sp>(A)
    (3)    <JumpJoint_sp>(A)
*/
Stub BondedJoint_O::getInputStub(chem::NVector_sp coords) const
{
  Stub stub;
  auto parent = this->parent();
  if (gc::IsA<StubJoint_sp>(parent)) {
    auto stubParent = gc::As_unsafe<StubJoint_sp>(parent);
    geom::stubFromThreePoints( stub._Transform,
                               stubParent->position(coords),
                               stubParent->transformedParentPos(),
                               stubParent->transformedGrandParentPos() );
    return stub;
  }
  auto grandParent = parent->parent();
  if (gc::IsA<StubJoint_sp>(grandParent)) {
    auto stubGrandParent = gc::As_unsafe<StubJoint_sp>(grandParent);
    geom::stubFromThreePoints( stub._Transform,
                               parent->position(coords),
                               stubGrandParent->position(coords),
                               stubGrandParent->transformedParentPos() );
    return stub;
  }
  auto greatGrandParent = grandParent->parent();
  if (gc::IsA<StubJoint_sp>(greatGrandParent)) {
    auto stubGreatGrandParent = gc::As_unsafe<StubJoint_sp>(greatGrandParent);
    geom::stubFromThreePoints( stub._Transform,
                               parent->position(coords),
                               grandParent->position(coords),
                               stubGreatGrandParent->position(coords) );
    return stub;
  }
  geom::stubFromThreePoints(stub._Transform,
                      parent->position(coords),
                      grandParent->position(coords),
                      greatGrandParent->position(coords));
  KIN_LOG("for {} stub = {}\n", _rep_(this->_Name), stub._Transform.asString());
  return stub;
}


void BondedJoint_O::_updateChildrenXyzCoords(chem::NVector_sp coords) {
  if (this->_numberOfChildren()>0) {
    int firstNonJumpIndex = this->firstNonJumpChildIndex();
    for ( int ii=0; ii < firstNonJumpIndex; ii++) {
      Stub jstub = this->_child(ii)->getInputStub(coords);
    // I should ratchet the newStub around the X axis and use relative dihedral
      this->_child(ii)->_updateXyzCoord(coords,jstub);
    // ratchet newStub
//    this->_DofChangePropagatesToYoungerSiblings = false;
      this->noteXyzUpToDate();
    }
    Stub stub = this->_child(firstNonJumpIndex)->getInputStub(coords);
    for ( int ii=firstNonJumpIndex; ii < this->_numberOfChildren(); ii++) {
    // I should ratchet the stub around the X axis and use relative dihedral
      this->_child(ii)->_updateXyzCoord(coords,stub);
//    this->_DofChangePropagatesToYoungerSiblings = false;
      this->noteXyzUpToDate();
    }
    for ( int ii=0; ii < this->_numberOfChildren(); ii++) {
      this->_child(ii)->_updateChildrenXyzCoords(coords);
    }
  }
}

void BondedJoint_O::_updateXyzCoord(chem::NVector_sp coords, Stub& stub)
{
      // https://math.stackexchange.com/questions/133177/finding-a-unit-vector-perpendicular-to-another-vector
  ASSERT(this->definedp());
  KIN_LOG("name = {} stub = \n{}\n", _rep_(this->_Name), stub._Transform.asString());
  KIN_LOG("_Distance = {}  _Theta = {} deg   _Phi = {} deg\n", this->_Distance , (this->_Theta/0.0174533) , (this->_Phi/0.0174533) );
  double bcTheta = this->_Theta;
  Vector3 d2;
//  printf("%s:%d:%s Calculating position for joint %s\n", __FILE__, __LINE__, __FUNCTION__, _rep_(this->_Name).c_str());
//  printf("%s:%d:%s distance = %lf  angle_deg = %lf   dihedral_deg = %lf\n", __FILE__, __LINE__, __FUNCTION__, this->_Distance, this->_Theta/0.0174533, this->_Phi/0.0174533 );
  Vector3 newpos = geom::pointFromStubAndInternalCoordinates(stub._Transform,this->_Distance, bcTheta, this->_Phi, d2 );
  if (!newpos.isDefined()) SIMPLE_ERROR("newpos could not be determined for {}", _rep_(this->asSmartPtr()));
  this->setPosition(coords,newpos);
#if 0
  // The following code prints everything for the calculation of newpos for a particular joint
  if (this->_PositionIndexX3 == 465) {
    core::clasp_write_string(fmt::format("{}:{}:{} distance = {}  angle_deg = {}   dihedral_deg = {}\n", __FILE__, __LINE__, __FUNCTION__, this->_Distance, this->_Theta/0.0174533, this->_Phi/0.0174533 ));
    core::clasp_write_string(fmt::format("{}:{}:{}  ggp = {}\n", __FILE__, __LINE__, __FUNCTION__, this->parent()->parent()->parent()->position(coords).asString() ));
    core::clasp_write_string(fmt::format("{}:{}:{}   gp = {}\n", __FILE__, __LINE__, __FUNCTION__, this->parent()->parent()->position(coords).asString() ));
    core::clasp_write_string(fmt::format("{}:{}:{}    p = {}\n", __FILE__, __LINE__, __FUNCTION__, this->parent()->position(coords).asString() ));
    core::clasp_write_string(fmt::format("{}:{}:{}  _PositionIndex3 {} stub = {}\n", __FILE__, __LINE__, __FUNCTION__, this->_PositionIndexX3, stub._Transform.asString() ));
    core::clasp_write_string(fmt::format("{}:{}:{} newp = {}\n", __FILE__, __LINE__, __FUNCTION__, newpos.asString() ));
    core::clasp_write_string(fmt::format("{}:{}:{}   d2 = {}, {}, {}\n", __FILE__, __LINE__, __FUNCTION__, d2.getX(), d2.getY(), d2.getZ() ));
    double cdih = geom::calculateDihedral(
        this->parent()->parent()->parent()->position(coords),
        this->parent()->parent()->position(coords),
        this->parent()->position(coords),
        newpos);
    core::clasp_write_string(fmt::format("{}:{}:{} calculated dihedral = {}\n", __FILE__, __LINE__, __FUNCTION__, cdih/0.0174533 ));
  }
#endif
  
}

CL_DEFMETHOD void BondedJoint_O::updateXyzCoord(chem::NVector_sp coords) {
  Stub stub = this->getInputStub(coords);
  this->_updateXyzCoord(coords,stub);
}


double BondedJoint_O::dof(DofType const& dof) const
{
  if ( dof == DofType::phi )
  {
    return this->_Phi;
  } else if ( dof == DofType::theta )
  {
    return this->_Theta;
  } else if ( dof == DofType::distance )
  {
    return this->_Distance;
  }
  SIMPLE_ERROR("Illegal dof request for BondedJoint - I can only handle internal dofs not rigid body");
}


core::Symbol_sp BondedJoint_O::typeSymbol() const { return _sym_bonded;};




};
