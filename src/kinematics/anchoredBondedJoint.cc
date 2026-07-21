/*
    File: anchoredBondedJoint.cc
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
#include <clasp/core/evaluator.h>
#include <cando/geom/vector3.h>
#include <cando/geom/ovector3.h>
#include <cando/geom/matrix.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/atomId.h>
#include <clasp/core/numerics.h>
#include <cando/kinematics/joint.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/jumpJoint.h>
#include <cando/kinematics/anchoredBondedJoint.h>

namespace kinematics
{

Stub AnchoredBondedJoint_O::getInputStub(chem::NVector_sp coords) const {
    Stub stub;
    geom::stubFromThreePoints(stub._Transform,
                              this->transformedParentPos(),
                              this->transformedGrandParentPos(),
                              this->transformedGreatGrandParentPos());
    return stub;               
  }
    
  void AnchoredBondedJoint_O::updateAnchorPositionsFromCoords(chem::NVector_sp coords) {
    geom::geom__vec_extract(this->_ParentAtomPosition,           coords, this->_ParentAtomPositionIndex3);
    geom::geom__vec_extract(this->_GrandParentAtomPosition,      coords, this->_GrandParentAtomPositionIndex3);
    geom::geom__vec_extract(this->_GreatGrandParentAtomPosition, coords, this->_GreatGrandParentAtomPositionIndex3);
  }   

CL_DEFMETHOD
Vector3 AnchoredBondedJoint_O::transformedParentPos() const {
  const Vector3& pos = this->_ParentAtomPosition;
  Matrix m = gc::As<geom::OMatrix_sp>(core::eval::funcall(_sym_orientation_transform, nil<core::T_O>()))->ref();
  Vector3 tpos = m.multiplyByVector3(pos);
  return tpos;
}

CL_DEFMETHOD
Vector3 AnchoredBondedJoint_O::transformedGrandParentPos() const {
  const Vector3& pos = this->_GrandParentAtomPosition;
  Matrix m = gc::As<geom::OMatrix_sp>(core::eval::funcall(_sym_orientation_transform, nil<core::T_O>()))->ref();
  Vector3 tpos = m.multiplyByVector3(pos);
  return tpos;
}

CL_DEFMETHOD
Vector3 AnchoredBondedJoint_O::transformedGreatGrandParentPos() const {
  const Vector3& pos = this->_GreatGrandParentAtomPosition;
  Matrix m = gc::As<geom::OMatrix_sp>(core::eval::funcall(_sym_orientation_transform, nil<core::T_O>()))->ref();
  Vector3 tpos = m.multiplyByVector3(pos);
  return tpos;
}

///CL_LAMBDA(atom-id name atom-table parent-pos grand-parent-pos great-grand-parent-pos);
CL_LISPIFY_NAME("make_AnchoredBondedJoint");
CL_DEF_CLASS_METHOD
AnchoredBondedJoint_sp AnchoredBondedJoint_O::make(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable
                               ,int parentAtomPositionIndex3
                               ,int grandParentAtomPositionIndex3
                               ,int greatGrandParentAtomPositionIndex3 ) {
  auto obj = gctools::GC<AnchoredBondedJoint_O>::allocate( atomId, name, atomTable, parentAtomPositionIndex3, grandParentAtomPositionIndex3, greatGrandParentAtomPositionIndex3 );
  return obj;
}

void AnchoredBondedJoint_O::initialize() {

}


void AnchoredBondedJoint_O::fields(core::Record_sp node) {
  node->field_if_not_default(INTERN_(kw,ggpatomPos),this->_GreatGrandParentAtomPosition,Vector3());
  node->field_if_not_default(INTERN_(kw,gpatomPos),this->_GrandParentAtomPosition,Vector3());
  node->field_if_not_default(INTERN_(kw,patomPos),this->_ParentAtomPosition,Vector3());
  node->field(INTERN_(kw,ggpatomPosIndex3),this->_GreatGrandParentAtomPositionIndex3);
  node->field(INTERN_(kw,gpatomPosIndex3),this->_GrandParentAtomPositionIndex3);
  node->field(INTERN_(kw,patomPosIndex3),this->_ParentAtomPositionIndex3);
  this->Base::fields(node);
}

};


