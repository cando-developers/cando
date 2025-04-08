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
#include <clasp/core/lispStream.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/evaluator.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/jumpJoint.h>

namespace kinematics
{

SYMBOL_EXPORT_SC_(KinPkg,orientation_selfFrame);
SYMBOL_EXPORT_SC_(KinPkg,orientation_labFrame);
SYMBOL_EXPORT_SC_(KinPkg,orientation_transform);

JumpJoint_O::JumpJoint_O( const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable ) : Joint_O(atomId,name,atomTable) {
};

bool JumpJoint_O::definedp(chem::NVector_sp internals) const
{
  return true;
}
CL_LAMBDA(atom-id name atom-table);
CL_LISPIFY_NAME("make_JumpJoint");
CL_DEF_CLASS_METHOD
    JumpJoint_sp JumpJoint_O::make(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable ) {
  return gctools::GC<JumpJoint_O>::allocate( atomId, name, atomTable );
}

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
{
  Joint_sp empty;
  size_t index = this->_Children.size();
  this->_Children.push_back(empty);
  this->_Children[index] = c;
  LOG(" Appending to node {} child {} at index {}\n" , _rep_(this->asSmartPtr()) , _rep_(c) , index);
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


void JumpJoint_O::_updateInternalCoord(chem::NVector_sp internals, chem::NVector_sp coords)
{
  KIN_LOG(" <<< {}\n", _rep_(this->asSmartPtr()));
}


bool JumpJoint_O::keepDofFixed(DofType dof) const
{
  return (this->parent().boundp());
}


void JumpJoint_O::_updateXyzCoord(chem::NVector_sp internals, chem::NVector_sp coords, Stub& stub)
{
//  ASSERTF(stub.isOrthogonal(1e-3),("Stub is not orthogonal - stub:\n{}") , stub.asString());
  this->setPosition(coords,stub._Transform.getTranslation());
  KIN_LOG("orientation {}\n", stub.asString() );
}

void JumpJoint_O::updateXyzCoord(chem::NVector_sp internals, chem::NVector_sp coords)
{
  Stub newStub = this->getInputStub(coords);
  this->_updateXyzCoord(internals, coords,newStub);
}

    /*! Update the external coordinates using the input stub */
void JumpJoint_O::_updateChildrenXyzCoords(chem::NVector_sp internals, chem::NVector_sp coords)
{
  Stub newStub;
  newStub._Transform = this->transform();
  for ( int ii=0; ii < this->_numberOfChildren(); ii++) {
    this->_child(ii)->_updateXyzCoord(internals, coords,newStub);
    // ratchet newStub
//    this->_DofChangePropagatesToYoungerSiblings = false;
    this->noteXyzUpToDate();
  }
  for ( int ii=0; ii < this->_numberOfChildren(); ii++) {
    this->_child(ii)->_updateChildrenXyzCoords(internals, coords);
  }
}

 Matrix JumpJoint_O::transform() const {
   // Invoke kin:orientation-transform to return the transform
   // for this jump-joint
   Matrix m = gc::As<geom::OMatrix_sp>(core::eval::funcall(_sym_orientation_transform, this->asSmartPtr()))->ref();
   return m;
 }

Stub JumpJoint_O::getInputStub(chem::NVector_sp coords) const {
  Stub stub;
  stub._Transform = this->transform();
  KIN_LOG("for {} stub = {}\n", _rep_(this->_Name), stub._Transform.asString());
  return stub;
}


double JumpJoint_O::dof(DofType const& dof) const
{
  if ( dof.isRigidBodyDof() )
  {
    SIMPLE_ERROR("Do something for dof");
//    return _Jump.getRigidBodyDelta(dof,RigidBodyDirection::n2c);
  }
  SIMPLE_ERROR("Illegal dof for JumpJoint - I can only return rigid body dofs");
}


core::Symbol_sp JumpJoint_O::typeSymbol() const { return _sym_jump;};



};
