/*
    File: coordinateCalculators.cc
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
#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/geom/matrix.h>
#include <clasp/core/symbolTable.h>
#include <cando/kinematics/joint.h>
#include <cando/kinematics/bondedJoint.h>
#include <cando/kinematics/jumpJoint.h>
#include <clasp/core/translators.h>
#include <cando/kinematics/coordinateCalculators.h>


namespace kinematics {

void noop_to_internal(Joint_sp joint) {
  printf("%s:%d:%s\n",  __FILE__, __LINE__, __FUNCTION__ );
};

void noop_to_external(Joint_sp joint) {
  printf("%s:%d:%s\n",  __FILE__, __LINE__, __FUNCTION__ );
};

void origin_to_internal(Joint_sp joint) {
  printf("%s:%d Here do something about modifying the origin joint transform\n", __FILE__, __LINE__ );
};


void origin_to_external(chem::NVector_sp coords, Joint_sp joint) {
#if 0
  // Apply the transform from the joint atom to the tempPos.
  if (!joint->isJump()) {
    SIMPLE_ERROR("joint must be a JumpJoint_sp {}" , _rep_(joint));
  }
  JumpJoint_sp jj = gc::As<JumpJoint_sp>(joint);
  Jump& jump = jj->_Jump;
  Vector3 tempPos = Vector3(0.0,0.0,0.0);
  joint->setPosition(coords, jump._Transform.multiplyByVector3(tempPos));
#endif
};

void bond_to_internal(chem::NVector_sp coords, Joint_sp joint ) {

  BondedJoint_sp bjoint = gc::As<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  bjoint->_Distance = geom::calculateDistance(joint->position(coords),bondJoint->position(coords));
  printf("%s:%d Here do something about modifying the bondJoint transform\n", __FILE__, __LINE__ );
};


void bond_to_external(chem::NVector_sp coords, Joint_sp joint) {
#if 0
  BondedJoint_sp bjoint = gc::As<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  Vector3 tempPos = geom::geom__build_using_bond(bjoint->_Distance,bondJoint->position(coords));
  // Apply the transform from the bondJoint atom to the tempPos.
  if (!bondJoint->isJump()) {
    SIMPLE_ERROR("bondJoint must be a JumpJoint_sp {}" , _rep_(bondJoint));
  }
  JumpJoint_sp jj = gc::As<JumpJoint_sp>(bondJoint);
  Jump& jump = jj->_Jump;
  bjoint->setPosition( coords, jump._Transform.multiplyByVector3(tempPos));
#endif
};


void bond_angle_to_internal(chem::NVector_sp coords, Joint_sp joint) {
  BondedJoint_sp bjoint = gc::As<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  Joint_sp angleJoint = bondJoint->_Parent;
  Vector3 joint_position = joint->position(coords);
  Vector3 bondJoint_position = bondJoint->position(coords);
  Vector3 angleJoint_position = angleJoint->position(coords);
  bjoint->_Distance = geom::calculateDistance(joint_position,bondJoint_position);
  bjoint->_Theta = geom::calculateAngle(joint_position,bondJoint_position,angleJoint_position);
  printf("%s:%d Here do something about modifying the angleJoint transform\n", __FILE__, __LINE__ );
};


void bond_angle_to_external(Joint_sp joint) {
  printf("%s:%d:%s Calculating\n", __FILE__, __LINE__, __FUNCTION__ );
#if 0
  BondedJoint_sp bjoint = gc::As<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  Joint_sp angleJoint = bondJoint->_Parent;
  Vector3 tempPos = geom::geom__build_using_bond_angle(bjoint->_Distance,bondJoint->position(coords),
                                                       bjoint->_Theta,angleJoint->position(coords));
  // Apply the transform from the angleJoint atom to the tempPos.
  if (!angleJoint->isJump()) {
    SIMPLE_ERROR("angleJoint must be a JumpJoint_sp {}" , _rep_(angleJoint));
  }
  JumpJoint_sp jj = gc::As<JumpJoint_sp>(angleJoint);
  Jump& jump = jj->_Jump;
  joint->setPosition( coords, jump._Transform.multiplyByVector3(tempPos));
#endif
}

void bond_angle_dihedral_to_internal(chem::NVector_sp coords, Joint_sp joint ) {
  BondedJoint_sp bjoint = gc::As<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  Joint_sp angleJoint = bondJoint->_Parent;
  Joint_sp dihedralJoint = angleJoint->_Parent;
  Vector3 joint_Position = joint->position(coords);
  Vector3 bondJoint_Position = bondJoint->position(coords);
  Vector3 angleJoint_Position = angleJoint->position(coords);
  Vector3 dihedralJoint_Position = dihedralJoint->position(coords);
  bjoint->_Distance = geom::calculateDistance(joint_Position,bondJoint_Position);
  bjoint->_Theta = geom::calculateAngle(joint_Position,bondJoint_Position,angleJoint_Position);
  bjoint->_Phi = geom::calculateDihedral(joint_Position,bondJoint_Position,angleJoint_Position,dihedralJoint_Position);
};


void bond_angle_dihedral_to_external(chem::NVector_sp coords, Joint_sp joint) {
  BondedJoint_sp bjoint = gc::As<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  Joint_sp angleJoint = bondJoint->_Parent;
  Joint_sp dihedralJoint = angleJoint->_Parent;
  joint->setPosition( coords, geom::geom__build_using_bond_angle_dihedral(bjoint->_Distance,bondJoint->position(coords),
                                                                          bjoint->_Theta,angleJoint->position(coords),
                                                                          bjoint->_Phi,dihedralJoint->position(coords)));
};


void general_to_internal(chem::NVector_sp coords, Joint_sp joint) {
  if (joint->_Parent.unboundp()) {
    origin_to_internal(joint);
    return;
  }
  BondedJoint_sp bjoint = gc::As<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  if (bondJoint->_Parent.unboundp()) {
    bond_to_internal(coords,joint);
    return;
  }
  Joint_sp angleJoint = bondJoint->_Parent;
  if (angleJoint->_Parent.unboundp()) {
    bond_angle_to_internal(coords,joint);
    return;
  }
  Joint_sp dihedralJoint = angleJoint->_Parent;
  bjoint->_Distance = geom::calculateDistance(joint->position(coords),bondJoint->position(coords));
  bjoint->_Theta = geom::calculateAngle(joint->position(coords),bondJoint->position(coords),angleJoint->position(coords));
  bjoint->_Phi = geom::calculateDihedral(joint->position(coords),bondJoint->position(coords),angleJoint->position(coords),dihedralJoint->position(coords));
}

void general_to_external(chem::NVector_sp coords,Joint_sp joint) {
  if (joint->_Parent.unboundp()) {
    origin_to_internal(joint);
    return;
  }
  BondedJoint_sp bjoint = gc::As<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  if (bondJoint->_Parent.unboundp()) {
    bond_to_internal(coords,joint);
    return;
  }
  Joint_sp angleJoint = bondJoint->_Parent;
  if (angleJoint->_Parent.unboundp()) {
    bond_angle_to_internal(coords,joint);
    return;
  }
  Joint_sp dihedralJoint = angleJoint->_Parent;
  joint->setPosition( coords, geom::geom__build_using_bond_angle_dihedral(bjoint->_Distance,bondJoint->position(coords),
                                                                          bjoint->_Theta,angleJoint->position(coords),
                                                                          bjoint->_Phi,dihedralJoint->position(coords)));
}



};
