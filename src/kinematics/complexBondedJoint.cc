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
#include <clasp/core/lisp.h>
#include <cando/geom/matrix.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/atomId.h>
#include <clasp/core/numerics.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/jumpJoint.h>
#include <cando/kinematics/complexBondedJoint.h>

namespace kinematics
{
void ComplexBondedJoint_O::fields(core::Record_sp record) {
  static_assert(ComplexBondedJoint_O::MaxInputStubJoints==2,"ComplexBondedJoint_O::MaxInputStubJoints has changed from 2 - update the code below");
  record->field_if_not_unbound(INTERN_(kw,inputStub1),this->_InputStubJoints[0]);
  record->field_if_not_unbound(INTERN_(kw,inputStub2),this->_InputStubJoints[1]);
  this->Base::fields(record);
}


void ComplexBondedJoint_O::initialize() {
  for (int i=0; i<MaxInputStubJoints; ++i ) {
    this->_InputStubJoints[i] = unbound<Joint_O>();
  }
}

CL_LAMBDA(atom-id name atom-table);
CL_LISPIFY_NAME("make_ComplexBondedJoint");
CL_DEF_CLASS_METHOD
ComplexBondedJoint_sp ComplexBondedJoint_O::make(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable ) {
  return gctools::GC<ComplexBondedJoint_O>::allocate(atomId, name, atomTable);
}



Stub ComplexBondedJoint_O::getInputStub(chem::NVector_sp coords) const
{
  Stub stub;
  if (this->inputStubJoint1().unboundp()) {
    if (gc::IsA<JumpJoint_sp>(this->parent())) {
      JumpJoint_sp jumpJoint = gc::As_unsafe<JumpJoint_sp>(this->parent());
//      Matrix flipped = jumpJoint->transform().flipXY();
      stub = jumpJoint->getInputStub(coords); // .flipXY();
      KIN_LOG("for {} stub = {}\n", _rep_(this->_Name), stub._Transform.asString());
      return stub;
    }
    SIMPLE_ERROR("Illegal to getInputStub with only a parent and that parent is not a JumpJoint_sp - it's a {}" , _rep_(this->parent()));
  }
  if (this->inputStubJoint2().unboundp()) {
    if (gc::IsA<JumpJoint_sp>(this->parent())) {
      Vector3 origin = this->parent()->position(coords);
      Vector3 axisX = Vector3(1.0,0.0,0.0);
      Vector3 axisY = Vector3(0.0,1.0,0.0);
      Vector3 axisZ = Vector3(0.0,0.0,1.0);
      Vector3 originPX = origin.add(axisX);
      Vector3 originPY = origin.add(axisY);
//      stub.fromFourPoints(origin, originPX, origin, originPY);
      stub.fromThreePoints(origin, originPX, originPY);
      KIN_LOG("for {} stub = {}\n", _rep_(this->_Name), stub._Transform.asString());
      return stub;
    }
    stub.fromCenterAndRotation( this->inputStubJoint0()->position(coords), gc::As<JumpJoint_sp>(this->inputStubJoint1())->transform().flipXY());
//    stub.fromCenterAndRotation( this->inputStubJoint0()->getPosition(), gc::As<JumpJoint_sp>(this->inputStubJoint1())->transform());
    KIN_LOG("for {} stub = {}\n", _rep_(this->_Name), stub._Transform.asString());
    return stub;
  }
  stub.fromThreePoints(this->inputStubJoint0()->position(coords),
                       this->inputStubJoint1()->position(coords),
                       this->inputStubJoint2()->position(coords));
  KIN_LOG("for {} stub = {}\n", _rep_(this->_Name), stub._Transform.asString());
  return stub;
}

CL_DEFMETHOD bool ComplexBondedJoint_O::inputStubJoint1BoundP() const {
  return !this->_InputStubJoints[0].unboundp();
}

CL_DEFMETHOD void ComplexBondedJoint_O::setInputStubJoint1(Joint_sp joint) {
  this->_InputStubJoints[0] = joint;
}

CL_DEFMETHOD void ComplexBondedJoint_O::makeUnboundInputStubJoint1() {
  this->_InputStubJoints[0] = unbound<Joint_O>();
}

CL_DEFMETHOD bool ComplexBondedJoint_O::inputStubJoint2BoundP() const {
  return !this->_InputStubJoints[1].unboundp();
}

CL_DEFMETHOD void ComplexBondedJoint_O::setInputStubJoint2(Joint_sp joint) {
  this->_InputStubJoints[1] = joint;
}

CL_DEFMETHOD void ComplexBondedJoint_O::makeUnboundInputStubJoint2() {
  this->_InputStubJoints[1] = unbound<Joint_O>();
}

void ComplexBondedJoint_O::_updateChildrenXyzCoords(chem::NVector_sp coords) {
  this->Joint_O::_updateChildrenXyzCoords(coords);
}

void ComplexBondedJoint_O::_updateInternalCoord(chem::NVector_sp coords)
{
  KIN_LOG(" <<< {}\n", _rep_(this->asSmartPtr()));
//	using numeric::x_rotation_matrix_radians;
//	using numerioc::z_rotation_matrix_radians;
//	using numeric::constants::d::pi;
  Joint_sp jC = this->parent();
  Vector3 C = jC->position(coords);
  this->_Distance = geom::calculateDistance(this->position(coords),C);
  KIN_LOG("Calculated _Distance = {}\n", this->_Distance );
  if (!this->inputStubJoint1().unboundp()) {
    Joint_sp jB = this->inputStubJoint1();
    Vector3 B = jB->position(coords);
//    this->_Theta = PREPARE_ANGLE(geom::calculateAngle(this->position(coords),C,B)); // Must be from incoming direction
    this->_Theta = geom::calculateAngle(this->position(coords),C,B); // Must be from incoming direction
    KIN_LOG("_Theta = {}\n", (this->_Theta/0.0174533));
    if (!this->inputStubJoint2().unboundp()) {
      Joint_sp jA = this->inputStubJoint2();
      Vector3 A = jA->position(coords);
      this->_Phi = geom::calculateDihedral(this->position(coords),C,B,A);
      KIN_LOG("_Phi = {}\n", (this->_Phi/0.0174533));
      return;
    }
  } else if (!this->inputStubJoint2().unboundp()) {
    Joint_sp jB = this->inputStubJoint2();
    Vector3 B = jB->position(coords);
//    this->_Theta = PREPARE_ANGLE(geom::calculateAngle(this->position(coords),C,B)); // Must be from incoming direction
    this->_Theta = geom::calculateAngle(this->position(coords),C,B); // Must be from incoming direction
    KIN_LOG("_Theta = {}\n", (this->_Theta/0.0174533));
  }
}

};
