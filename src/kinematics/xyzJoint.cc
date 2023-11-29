/*
    File: xyzJoint.cc
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
#include <cando/geom/matrix.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/atomId.h>
#include <clasp/core/numerics.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/jumpJoint.h>
#include <cando/kinematics/xyzJoint.h>

namespace kinematics
{
void XyzJoint_O::fields(core::Record_sp record) {
  static_assert(XyzJoint_O::MaxChildren==5,"XyzJoint_O::MaxChildren has changed from 5 - update the code below");
  record->field_if_not_unbound(INTERN_(kw,child4),this->_Children[4]);
  record->field_if_not_unbound(INTERN_(kw,child3),this->_Children[3]);
  record->field_if_not_unbound(INTERN_(kw,child2),this->_Children[2]);
  record->field_if_not_unbound(INTERN_(kw,child1),this->_Children[1]);
  record->field_if_not_unbound(INTERN_(kw,child0),this->_Children[0]);
  record->field(INTERN_(kw,num_children),this->_NumberOfChildren);
  record->field(INTERN_(kw,x),this->_X);
  record->field(INTERN_(kw,y),this->_Y);
  record->field(INTERN_(kw,z),this->_Z);
  this->Base::fields(record);
}


void XyzJoint_O::initialize() {
  for (int i=0; i<MaxChildren; ++i ) {
    this->_Children[i] = unbound<Joint_O>();
  }
}

	/*! Return the stubJoint1 */
CL_DEFMETHOD Joint_sp XyzJoint_O::inputStubJoint0() const { return this->parent();}

	/*! Return the stubJoint2 */
CL_DEFMETHOD Joint_sp XyzJoint_O::inputStubJoint1() const { return this->parent()->parent();};

	/*! Return the stubJoint3 */
CL_DEFMETHOD Joint_sp XyzJoint_O::inputStubJoint2() const { return this->parent()->parent()->parent(); };


CL_LAMBDA(atom-id name atom-table);
CL_LISPIFY_NAME("make_XyzJoint");
CL_DEF_CLASS_METHOD
XyzJoint_sp XyzJoint_O::make(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable, core::T_sp orientation) {
  auto obj = gctools::GC<XyzJoint_O>::allocate(atomId,name,atomTable,orientation);
  return obj;
}


void XyzJoint_O::_appendChild(Joint_sp c)
{
  size_t index = this->_NumberOfChildren++;
  LOG(" Appending to node {} child {} at index {}\n" , _rep_(this->asSmartPtr()) , _rep_(c) , index);
  this->_Children[index] = c;
}
  
void XyzJoint_O::_insertChild(int before, Joint_sp child)
{
  if ( this->_numberOfChildren()>this->_maxNumberOfChildren() )
  {
    THROW_HARD_ERROR("You exceeded the maximum[%d] number of children allowed for a XyzJoint" , this->_maxNumberOfChildren());
  }
  for ( int i=this->_numberOfChildren(); i>before; i-- )
    this->_Children[i] = this->_Children[i-1];
  this->_Children[before] = child;
  this->_NumberOfChildren++;
  child->_Parent = this->asSmartPtr();
}

void XyzJoint_O::_releaseChild(int idx)
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

void XyzJoint_O::_releaseAllChildren()
{
  if ( this->_numberOfChildren() == 0 ) return;
  int num = this->_numberOfChildren();
  for ( int idx=0; idx < num; idx++ )
  {
    this->_Children[idx] = unbound<Joint_O>();
  }
  this->_NumberOfChildren = 0;
}

void XyzJoint_O::_updateInternalCoord(chem::NVector_sp coords)
{
  KIN_LOG(" <<< {}\n" , _rep_(this->asSmartPtr()));
//	using numeric::x_rotation_matrix_radians;
//	using numerioc::z_rotation_matrix_radians;
//	using numeric::constants::d::pi;
  if (this->_Orientation.nilp()) {
    auto x = (*coords)[this->_PositionIndexX3];
    auto y = (*coords)[this->_PositionIndexX3+1];
    auto z = (*coords)[this->_PositionIndexX3+2];
    this->_X = x;
    this->_Y = y;
    this->_Z = z;
  } else {
    auto x = (*coords)[this->_PositionIndexX3];
    auto y = (*coords)[this->_PositionIndexX3+1];
    auto z = (*coords)[this->_PositionIndexX3+2];
    this->_X = x;
    this->_Y = y;
    this->_Z = z;
  }
}

bool XyzJoint_O::keepDofFixed(DofType dof) const
{
  IMPLEMENT_ME();
}


string XyzJoint_O::asString() const
{
  stringstream ss;
  ss << this->Joint_O::asString();
  ss << fmt::format("  _X[{}]  _Y[{}]  _Z[{}]"
                    , this->_X, this->_Y, this->_Z );
  return ss.str();
}

Stub XyzJoint_O::getInputStub(chem::NVector_sp coords) const
{
  Stub stub;
  if (this->_Orientation.nilp()) {
    stub._Transform.setToIdentity();
  } else {
    Matrix m = gc::As<geom::OMatrix_sp>(core::eval::funcall(_sym_orientation_transform, this->_Orientation))->ref();
    stub._Transform = m;
  }
  return stub;
}


void XyzJoint_O::_updateChildrenXyzCoords(chem::NVector_sp coords) {
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

void XyzJoint_O::_updateXyzCoord(chem::NVector_sp coords, Stub& stub)
{
      // https://math.stackexchange.com/questions/133177/finding-a-unit-vector-perpendicular-to-another-vector
  Vector3 d2;
  d2.set(this->_X,this->_Y,this->_Z);
  if (this->_Orientation.nilp()) {
    this->setPosition(coords,d2);
  } else {
   Matrix m = gc::As<geom::OMatrix_sp>(core::eval::funcall(_sym_orientation_transform, this->_Orientation))->ref();
   this->setPosition(coords,m*d2);
  }
}

CL_DEFMETHOD void XyzJoint_O::updateXyzCoord(chem::NVector_sp coords) {
  Stub stub = this->getInputStub(coords);
  this->_updateXyzCoord(coords,stub);
}


double XyzJoint_O::dof(DofType const& dof) const
{
  IMPLEMENT_ME();
}


SYMBOL_EXPORT_SC_(KinPkg,xyz);

core::Symbol_sp XyzJoint_O::typeSymbol() const { return _sym_xyz;};




};
