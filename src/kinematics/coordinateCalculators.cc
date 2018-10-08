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
#include <cando/kinematics/atom.h>
#include <cando/kinematics/bondedAtom.h>
#include <cando/kinematics/jumpAtom.h>
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


void origin_to_external(Joint_sp joint) {
#if 0
  // Apply the transform from the joint atom to the tempPos.
  if (!joint->isJump()) {
    SIMPLE_ERROR(BF("joint must be a JumpJoint_sp %s") % _rep_(joint));
  }
  JumpJoint_sp jj = gc::As_unsafe<JumpJoint_sp>(joint);
  Jump& jump = jj->_Jump;
  Vector3 tempPos = Vector3(0.0,0.0,0.0);
  joint->_Position = jump._Transform.multiplyByVector3(tempPos);
#endif
};

void bond_to_internal(Joint_sp joint) {

  BondedJoint_sp bjoint = gc::As_unsafe<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  bjoint->_Distance = geom::calculateDistance(joint->_Position,bondJoint->_Position);
  printf("%s:%d Here do something about modifying the bondJoint transform\n", __FILE__, __LINE__ );
};


void bond_to_external(Joint_sp joint) {
#if 0
  BondedJoint_sp bjoint = gc::As_unsafe<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  Vector3 tempPos = geom::buildUsingBond(bjoint->_Distance,bondJoint->_Position);
  // Apply the transform from the bondJoint atom to the tempPos.
  if (!bondJoint->isJump()) {
    SIMPLE_ERROR(BF("bondJoint must be a JumpJoint_sp %s") % _rep_(bondJoint));
  }
  JumpJoint_sp jj = gc::As_unsafe<JumpJoint_sp>(bondJoint);
  Jump& jump = jj->_Jump;
  bjoint->_Position = jump._Transform.multiplyByVector3(tempPos);
#endif
};


void bond_angle_to_internal(Joint_sp joint) {
  BondedJoint_sp bjoint = gc::As_unsafe<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  Joint_sp angleJoint = bondJoint->_Parent;
  bjoint->_Distance = geom::calculateDistance(joint->_Position,bondJoint->_Position);
  bjoint->_Theta = geom::calculateAngle(joint->_Position,bondJoint->_Position,angleJoint->_Position);
  printf("%s:%d Here do something about modifying the angleJoint transform\n", __FILE__, __LINE__ );
};


void bond_angle_to_external(Joint_sp joint) {
#if 0
  BondedJoint_sp bjoint = gc::As_unsafe<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  Joint_sp angleJoint = bondJoint->_Parent;
  Vector3 tempPos = geom::buildUsingBondAngle(bjoint->_Distance,bondJoint->_Position,
                                              bjoint->_Theta,angleJoint->_Position);
  // Apply the transform from the angleJoint atom to the tempPos.
  if (!angleJoint->isJump()) {
    SIMPLE_ERROR(BF("angleJoint must be a JumpJoint_sp %s") % _rep_(angleJoint));
  }
  JumpJoint_sp jj = gc::As_unsafe<JumpJoint_sp>(angleJoint);
  Jump& jump = jj->_Jump;
  joint->_Position = jump._Transform.multiplyByVector3(tempPos);
#endif
}

void bond_angle_dihedral_to_internal(Joint_sp joint) {
  BondedJoint_sp bjoint = gc::As_unsafe<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  Joint_sp angleJoint = bondJoint->_Parent;
  Joint_sp dihedralJoint = angleJoint->_Parent;
  bjoint->_Distance = geom::calculateDistance(joint->_Position,bondJoint->_Position);
  bjoint->_Theta = geom::calculateAngle(joint->_Position,bondJoint->_Position,angleJoint->_Position);
  bjoint->_Phi = geom::calculateDihedral(joint->_Position,bondJoint->_Position,angleJoint->_Position,dihedralJoint->_Position);
};


void bond_angle_dihedral_to_external(Joint_sp joint) {
  BondedJoint_sp bjoint = gc::As_unsafe<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  Joint_sp angleJoint = bondJoint->_Parent;
  Joint_sp dihedralJoint = angleJoint->_Parent;
  joint->_Position = geom::buildUsingBondAngleDihedral(bjoint->_Distance,bondJoint->_Position,
                                                       bjoint->_Theta,angleJoint->_Position,
                                                       bjoint->_Phi,dihedralJoint->_Position);
};


void general_to_internal(Joint_sp joint) {
  if (joint->_Parent.unboundp()) {
    origin_to_internal(joint);
    return;
  }
  BondedJoint_sp bjoint = gc::As_unsafe<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  if (bondJoint->_Parent.unboundp()) {
    bond_to_internal(joint);
    return;
  }
  Joint_sp angleJoint = bondJoint->_Parent;
  if (angleJoint->_Parent.unboundp()) {
    bond_angle_to_internal(joint);
    return;
  }
  Joint_sp dihedralJoint = angleJoint->_Parent;
  bjoint->_Distance = geom::calculateDistance(joint->_Position,bondJoint->_Position);
  bjoint->_Theta = geom::calculateAngle(joint->_Position,bondJoint->_Position,angleJoint->_Position);
  bjoint->_Phi = geom::calculateDihedral(joint->_Position,bondJoint->_Position,angleJoint->_Position,dihedralJoint->_Position);
}

void general_to_external(Joint_sp joint) {
  if (joint->_Parent.unboundp()) {
    origin_to_internal(joint);
    return;
  }
  BondedJoint_sp bjoint = gc::As_unsafe<BondedJoint_sp>(joint);
  Joint_sp bondJoint = joint->_Parent;
  if (bondJoint->_Parent.unboundp()) {
    bond_to_internal(joint);
    return;
  }
  Joint_sp angleJoint = bondJoint->_Parent;
  if (angleJoint->_Parent.unboundp()) {
    bond_angle_to_internal(joint);
    return;
  }
  Joint_sp dihedralJoint = angleJoint->_Parent;
  joint->_Position = geom::buildUsingBondAngleDihedral(bjoint->_Distance,bondJoint->_Position,
                                                       bjoint->_Theta,angleJoint->_Position,
                                                       bjoint->_Phi,dihedralJoint->_Position);
}



};
