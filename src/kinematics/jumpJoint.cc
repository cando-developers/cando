/*
    File: jumpJoint.cc
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



#include <clasp/core/foundation.h>
#include <cando/kinematics/jointTree.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/symbolTable.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/jumpJoint.h>

namespace kinematics
{

void JumpJoint_O::fields(core::Record_sp node) {
  node->field_if_not_empty(INTERN_(kw,children),this->_Children);
  this->Base::fields(node);
}

void JumpJoint_O::_insertChild(int idx, Joint_sp c)
{
  Joint_sp empty;
  this->_Children.ensure_initialized();
  this->_Children.insert(this->_Children.begin()+idx,empty);
  this->_Children[idx] = c;
}


void JumpJoint_O::_appendChild(Joint_sp c)
{_OF();
  Joint_sp empty;
  size_t index = this->_Children.size();
  this->_Children.push_back(empty);
  this->_Children[index] = c;
  LOG(BF(" Appending to node %s child %s at index %lu\n") % _rep_(this->asSmartPtr()) % _rep_(c) % index);
}


void JumpJoint_O::_releaseChild(int idx)
{
  this->_Children[idx] = unbound<Joint_O>();
  this->_Children.erase(this->_Children.begin()+idx);
}


void JumpJoint_O::_releaseAllChildren()
{
  this->_Children.clear();
}


void JumpJoint_O::_updateInternalCoord()
{_OF();
  KIN_LOG(BF(" <<< %s\n") % _rep_(this->asSmartPtr()));
  Vector3 O = this->_Position;
  if (this->_numberOfChildren()>=2) {
    Vector3 A = this->_child(0)->_Position;
    Vector3 B = this->_child(1)->_Position;
    Vector3 OA = A - O;
    double lengthOA = OA.length();
    if (lengthOA<1e-6) SIMPLE_ERROR(BF("You are about to divide by zero"));
    Vector3 a = OA*(1.0/lengthOA);
    Vector3 OB = B - O;
    Vector3 OBxa = OB.crossProduct(a);
    double lengthOBxa = OBxa.length();
    if (lengthOBxa<1e-6) SIMPLE_ERROR(BF("You are about to divide by zero"));
    Vector3 c = OBxa*(1.0/lengthOBxa);
    Vector3 b = c.crossProduct(a);
    Matrix labFrame;
    labFrame.colX(a);
    labFrame.colY(b);
    labFrame.colZ(c);
    labFrame.setTranslate(O);
    this->_LabFrame = labFrame;
    JumpJoint_sp jjParent = gc::As<JumpJoint_sp>(this->parent());
    Matrix parentInverted = jjParent->_LabFrame.invertTransform();
    this->_ParentRelativeTransform = labFrame*parentInverted;
    KIN_LOG(BF("relativeTransform = \n%s\n") % this->_ParentRelativeTransform.asString());
  } else if (this->_numberOfChildren()==1) {
    Vector3 A = this->_child(0)->_Position;
    Vector3 OA = A - O;
    double lengthOA = OA.length();
    if (lengthOA<1e-6) SIMPLE_ERROR(BF("You are about to divide by zero"));
    Vector3 a = OA*(1.0/lengthOA);
      // We need a vector orthogonal to a - used the following reference
      // https://math.stackexchange.com/questions/133177/finding-a-unit-vector-perpendicular-to-another-vector
    Vector3& x = a;  // alias a with x to match algorithm
    Vector3 y(0.0,0.0,0.0);
    int m = 0;
    int n = 0;
    if (std::fabs(x[0]) > 1e-3) { m = 0; n = 1; }
    else if (std::fabs(x[1]) > 1e-3 ) { m = 1; n = 2; }
    else { m = 2; n = 0; };
    y[n] = x[m];
    y[m] = -x[n];
    double lengthy = y.length();
    if (lengthy<1e-6) SIMPLE_ERROR(BF("You are about to divide by zero"));
    y = y * (1.0/lengthy);
    // Now y should be orthogonal to a
    Vector3& b = y;
    Vector3 c = a.crossProduct(b);
    Matrix labFrame;
    labFrame.colX(a);
    labFrame.colY(b);
    labFrame.colZ(c);
    labFrame.setTranslate(O);
    this->_LabFrame = labFrame;
    JumpJoint_sp jjParent = gc::As<JumpJoint_sp>(this->parent());
    Matrix parentInverted = jjParent->_LabFrame.invertTransform();
    this->_ParentRelativeTransform = labFrame*parentInverted;
    KIN_LOG(BF("relativeTransform = \n%s\n") % this->_ParentRelativeTransform.asString());
  }
}

void JumpJoint_O::updateInternalCoords(bool const recursive,
                                       JointTree_sp at) {
  this->_updateInternalCoord();
  if ( recursive )
  {
    for (int childIdx=0; childIdx<this->_numberOfChildren(); childIdx++ ) {
      this->_child(childIdx).get()->updateInternalCoords(true,at);
    }
  }
}



Joint_sp JumpJoint_O::stubJoint1() const
{
  if (this->stubDefined()) {
    return this->asSmartPtr();
  }
  return this->parent();
}

Joint_sp JumpJoint_O::stubJoint2() const
{
  if ( this->stubDefined() )
  {
    return this->getNonJumpJoint(0);
  } else
  {
    return this->parent()->stubJoint2();
  }
}

Joint_sp JumpJoint_O::stubJoint3(JointTree_sp at) const
{_OF();
  if ( this->stubDefined() )
  {
    Joint_sp first(this->getNonJumpJoint(0));
    ASSERT(first.boundp());
    Joint_sp second(first->getNonJumpJoint(0));
    if ( second.boundp() )
    {
      return second;
    } else
    {
      return this->getNonJumpJoint(1);
    }
  } else
  {
    return this->parent()->stubJoint3(at);
  }
}

bool JumpJoint_O::keepDofFixed(DofType dof) const
{
  return (this->parent().boundp());
}


void JumpJoint_O::_updateXyzCoord(Stub& stub)
{_OF();
  ASSERTF(stub.isOrthogonal(1e-3),BF("Stub is not orthogonal - stub:\n%s") % stub.asString());
  this->_LabFrame = stub._Transform.multiplyByMatrix(this->_ParentRelativeTransform);
  this->position(this->_LabFrame.getTranslation());
  KIN_LOG(BF("LabFrame.getTranslation() = %s\n") % this->_LabFrame.getTranslation().asString());
}


    /*! Update the external coordinates using the input stub */
void JumpJoint_O::_updateXyzCoords(Stub& stub)
{_OF();
  this->_updateXyzCoord(stub);
  Stub newStub;
  newStub._Transform = this->_LabFrame;
  for ( int ii=0; ii<this->_numberOfChildren(); ii++ ) {
    this->_child(ii).get()->_updateXyzCoords(newStub);
  }
  this->noteXyzUpToDate();
}

Stub JumpJoint_O::getStub() const {
  Stub jj;
  jj._Transform = this->_LabFrame;
  return jj;
}



void JumpJoint_O::updateXyzCoords()
{
  // parent must be an OriginJumpJoint
  Stub stub = this->parent()->getStub();
  this->JumpJoint_O::_updateXyzCoords(stub);
}



double JumpJoint_O::dof(DofType const& dof) const
{_OF();
  if ( dof.isRigidBodyDof() )
  {
    SIMPLE_ERROR(BF("Do something for dof"));
//    return _Jump.getRigidBodyDelta(dof,RigidBodyDirection::n2c);
  }
  SIMPLE_ERROR(BF("Illegal dof for JumpJoint - I can only return rigid body dofs"));
}


core::Symbol_sp JumpJoint_O::typeSymbol() const {_OF(); return _sym_jump;};



};
