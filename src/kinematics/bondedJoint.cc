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
#include <cando/kinematics/jointTree.h>
#include <cando/kinematics/jumpJoint.h>
#include <cando/kinematics/bondedJoint.h>

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
  record->field(INTERN_(kw,phi),this->_Phi);
  record->field(INTERN_(kw,theta),this->_Theta);
  record->field(INTERN_(kw,distance),this->_Distance);
  record->field(INTERN_(kw,dof_propagates),this->_DofChangePropagatesToYoungerSiblings);
  this->Base::fields(record);
}


void BondedJoint_O::initialize() {
  for (int i=0; i<MaxChildren; ++i ) {
    this->_Children[i] = unbound<Joint_O>();
  }
}

CL_LAMBDA(atom-id &optional name);
CL_LISPIFY_NAME("make_BondedJoint");
CL_DEF_CLASS_METHOD
BondedJoint_sp BondedJoint_O::make(const chem::AtomId& atomId, core::T_sp name) {
  return gctools::GC<BondedJoint_O>::allocate(atomId,name);
}


void BondedJoint_O::_appendChild(Joint_sp c)
{_OF();
  size_t index = this->_NumberOfChildren++;
  LOG(BF(" Appending to node %s child %s at index %lu\n") % _rep_(this->asSmartPtr()) % _rep_(c) % index);
  this->_Children[index] = c;
}
  
    void BondedJoint_O::_insertChild(int before, Joint_sp child)
    {
	if ( this->_numberOfChildren()>this->_maxNumberOfChildren() )
	{
	    THROW_HARD_ERROR(BF("You exceede the maximum[%d] number of children allowed for a BondedJoint") % this->_maxNumberOfChildren());
	}
	for ( int i=this->_numberOfChildren(); i>before; i-- )
	    this->_Children[i] = this->_Children[i-1];
	this->_Children[before] = child;
	this->_NumberOfChildren++;
    }

    void BondedJoint_O::_releaseChild(int idx)
    {
	if ( this->_numberOfChildren() == 0 )
	{
	    THROW_HARD_ERROR(BF("There are no children to delete"));
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



void BondedJoint_O::_updateInternalCoord()
{_OF();
  KIN_LOG(BF(" <<< %s\n") % _rep_(this->asSmartPtr()));
//	using numeric::x_rotation_matrix_radians;
//	using numerioc::z_rotation_matrix_radians;
//	using numeric::constants::d::pi;
  Joint_sp jC = this->parent();
  Vector3 C = jC->position();
  this->_Distance = geom::calculateDistance(this->_Position,C);
  KIN_LOG(BF("Calculated _Distance = %lf\n") % this->_Distance );
  if (gc::IsA<BondedJoint_sp>(jC) && gc::IsA<BondedJoint_sp>(jC->parent())) { // don't move past JumpJoint_O nodes
    KIN_LOG(BF("!gc::IsA<JumpJoint_sp>(jC)   jC = %s\n") % _rep_(jC));
    Joint_sp jB = jC->parent();
    Vector3 B = jB->position();
    this->_Theta = PREPARE_ANGLE(geom::calculateAngle(this->_Position,C,B)); // Must be from incoming direction
    KIN_LOG(BF("_Theta = %lf\n") % (this->_Theta/0.0174533));
    Joint_sp jA = jB->parent();
    Vector3 A = jA->position();
    this->_Phi = geom::calculateDihedral(this->_Position,C,B,A);
    KIN_LOG(BF("_Phi = %lf\n") % (this->_Phi/0.0174533));
    return;
  }
#if 1
  internalCoordinatesFromPointAndCoordinateSystem(this->getPosition(),this->getInputStub()._Transform,
                                                   this->_Distance, this->_Theta, this->_Phi );
#else
  KIN_LOG(BF("gc::IsA<JumpJoint_sp>(jC)   jC = %s\n") % _rep_(jC));
  Stub stub = jC->getStub();
  KIN_LOG(BF("stub = \n%s\n") % stub._Transform.asString());
  Vector3 x = stub._Transform.colX();
  Vector3 y = stub._Transform.colY();
  Vector3 z = stub._Transform.colZ();
  KIN_LOG(BF("x = %s\n") % x.asString());
  KIN_LOG(BF("y = %s\n") % y.asString());
  KIN_LOG(BF("z = %s\n") % z.asString());
  Vector3 D = this->getPosition();
  Vector3 CD = D - C;
  double lengthCD = CD.length();
  if (lengthCD<SMALL_NUMBER) SIMPLE_ERROR(BF("About to divide by zero"));
  Vector3 d = CD*(1.0/lengthCD);
  KIN_LOG(BF("d = %s\n") % d.asString());
  double dx = d.dotProduct(x);
  double dy = d.dotProduct(y);
  double dz = d.dotProduct(z);
  KIN_LOG(BF("dx = %lf  dy = %lf  dz = %lf\n") % dx % dy % dz );
  this->_Phi = geom::geom__planeVectorAngle(dy,dz);
  KIN_LOG(BF("  dy = %lf   dz = %lf\n") % dy % dz );
  KIN_LOG(BF("_Phi = %lf deg\n") % (this->_Phi/0.0174533));
  Vector3 dox(1.0,0.0,0.0);
  Vector3 dop(dx,dy,dz);
  KIN_LOG(BF("dop = %s\n") % dop.asString());
  KIN_LOG(BF("dop.dotProduct(dox) = %lf\n") % dop.dotProduct(dox));
  if (dop.dotProduct(dox) > (1.0-SMALL_NUMBER)) {
    this->_Theta = 0.0;
    return;
  }
  Vector3 doz = dox.crossProduct(dop);
  doz = doz.normalized();
  KIN_LOG(BF("doz = %s\n") % doz.asString());
  Vector3 doy = doz.crossProduct(dox);
  KIN_LOG(BF("doy = %s\n") % doy.asString());
  double eox = dop.dotProduct(dox);
  double eoy = dop.dotProduct(doy);
  KIN_LOG(BF("eox = %lf  eoy = %lf\n") % eox % eoy );
//  double eoz = dop.dotProduct(doz); // Must be 0.0
  this->_Theta = geom::geom__planeVectorAngle(eox,eoy);
  KIN_LOG(BF("    this->_Theta = %lf deg\n") % (this->_Theta/0.0174533));
#endif
}

void BondedJoint_O::updateInternalCoords( bool recursive,
                                          JointTree_sp at)
{_OF();
  this->_updateInternalCoord();
  for ( int it=0; it<this->_numberOfChildren(); it++ ) {
    this->_child(it)->updateInternalCoords( true, at );
  }
  return;
}

bool BondedJoint_O::keepDofFixed(DofType dof,JointTree_sp at) const
{_OF();
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
    SIMPLE_ERROR(BF("BondedJoint_O::keepDofFixed: BAD_DOF: %s") % dof.asString() );
  }
  return false;
#endif
}


    string BondedJoint_O::asString() const
    {
	stringstream ss;
	ss << this->Joint_O::asString();
	ss << (BF("  _Distance[%lf]  _Theta[%lf]/deg  _Phi[%lf]/deg")
	       % this->_Distance
	       % (this->_Theta/0.0174533)
	       % (this->_Phi/0.0174533) ).str() << std::endl;
	return ss.str();
    }

/*! There are three possible situations
    C is always a BondedJoint_sp
    (1)    <BondedJoint_sp>(B) <BondedJoint_sp>(A)
    (2)    <BondedJoint_sp>(B) <JumpJoint_sp>(A)
    (3)    <JumpJoint_sp>(A)
*/
Stub BondedJoint_O::getInputStub() const
{_OF();
  Stub stub;
  stub.fromFourPoints(this->inputStubJoint0()->position(),
                      this->inputStubJoint1()->position(),
                      this->inputStubJoint0()->position(),
                      this->inputStubJoint2()->position());
  return stub;
#if 0
  KIN_LOG(BF("getStub() -> %s\n") % _rep_(this->asSmartPtr()));
  Stub stub;
  BondedJoint_sp jbC = this->asSmartPtr();
  Vector3 C = jbC->position();
  KIN_LOG(BF("C position = %s\n") % C.asString());
  if (gc::IsA<BondedJoint_sp>(jbC->parent())) {
    // jB is a BondedJoint - we have situation (1) or (2) - we only need the position of jC
    KIN_LOG(BF("Situation (1) or (2) calculating stub\n"));
    KIN_LOG(BF("Got jbC->%s\n") % _rep_(jbC));
    BondedJoint_sp jbB = gc::As_unsafe<BondedJoint_sp>(this->parent());
    KIN_LOG(BF("Got jbB->%s\n") % _rep_(jbB));
    Vector3 B = jbB->position();
    KIN_LOG(BF("B position = %s\n") % B.asString());
    Joint_sp jA = jbB->parent();
    KIN_LOG(BF("Got jA->%s\n") % _rep_(jA));
    Vector3 A = jA->position();
    KIN_LOG(BF("A position = %s\n") % A.asString());
    Vector3 BC = (C - B);
    double lengthBC = BC.length();
    if (lengthBC<SMALL_NUMBER) SIMPLE_ERROR(BF("About to divide by zero lengthBC"));
    Vector3 bc = BC.multiplyByScalar(1.0/lengthBC);
    KIN_LOG(BF("bc = %s\n") % bc.asString());
    Vector3 AB = (B - A);
    KIN_LOG(BF("AB = %s\n") % AB.asString());
    Vector3 ABxbc = AB.crossProduct(bc);
    KIN_LOG(BF("ABxbc = %s\n") % ABxbc.asString());
    double lengthABxbc = ABxbc.length();
    KIN_LOG(BF("lengthABxbc = %lf\n") % lengthABxbc );
    if (lengthABxbc<SMALL_NUMBER) {
      SIMPLE_ERROR(BF("About to divide by zero when calculating stub for %s")% _rep_(jbC));
    }
    Vector3 n = ABxbc.multiplyByScalar(1.0/lengthABxbc);
    KIN_LOG(BF("n = %s\n") % n.asString());
    Vector3 my = n.crossProduct(bc);
    KIN_LOG(BF("Got my->%s\n") % my.asString());
        // Define the columns of the stub
        // stub._Transform = [ bc, my, n, C ]
    stub._Transform.colX(bc);
    stub._Transform.colY(my);
    stub._Transform.colZ(n);
    stub._Transform.setTranslate(C);
    KIN_LOG(BF("Returning stub:\n%s\n") % stub._Transform.asString());
    return stub;
  } else if (gc::IsA<JumpJoint_sp>(jbC->parent())) {
    // Situation (3) - we need to handle specially
    KIN_LOG(BF("Situation (3) calculating stub\n"));
    JumpJoint_sp jjB = gc::As_unsafe<JumpJoint_sp>(jbC->parent());
    KIN_LOG(BF("B is jump-joint %s\n") % _rep_(jbC));
    Matrix Btransform = jjB->_LabFrame;
    Vector3 Bx = Btransform.colX();
    Vector3 By = Btransform.colY();
    Vector3 Bz = Btransform.colZ();
    Vector3 X = C-jjB->getPosition();
    Vector3 x = X.normalized();
    stub._Transform.colX(x);
    if (Bz.dotProduct(x)<0.1) {
      Vector3 y = Bz.crossProduct(x);
      y = y.normalized();
      Vector3 z = x.crossProduct(y);
      stub._Transform.colY(y);
      stub._Transform.colZ(z);
    } else {
      Vector3 y = By.crossProduct(x);
      y = y.normalized();
      Vector3 z = x.crossProduct(y);
      stub._Transform.colY(y);
      stub._Transform.colZ(z);
    }
    stub._Transform.setTranslate(C);
    KIN_LOG(BF("Returning stub:\n%s\n") % stub._Transform.asString());
    return stub;
  }
  SIMPLE_ERROR(BF("Parent jB %s must be a jump-joint or a bonded-joint ") % _rep_(jbC->parent()));
#endif
}

void BondedJoint_O::_updateXyzCoord(Stub& stub)
{
      // https://math.stackexchange.com/questions/133177/finding-a-unit-vector-perpendicular-to-another-vector
  KIN_LOG(BF("stub = \n%s\n") % stub._Transform.asString());
  KIN_LOG(BF("_Distance = %lf  _Theta = %lf deg   _Phi = %lf deg\n")
          % this->_Distance
          % (this->_Theta/0.0174533)
          % (this->_Phi/0.0174533) );
  double bcTheta = FINAL_ANGLE(this->_Theta);
#if 1
  Vector3 d2;
  this->_Position = pointFromMatrixAndInternalCoordinates(stub._Transform,this->_Distance, bcTheta, this->_Phi, d2 );
#else
  KIN_LOG(BF(" rtTheta = %lf deg\n") % (bcTheta/0.0174533));
  double cosTheta = std::cos(bcTheta);
  double sinTheta = std::sin(bcTheta);
  double cosPhi = std::cos(this->_Phi);
  double sinPhi = std::sin(this->_Phi);
  Vector3 d2(this->_Distance*cosTheta,this->_Distance*cosPhi*sinTheta,this->_Distance*sinPhi*sinTheta);
  KIN_LOG(BF("d2 = %s\n") % d2.asString());
  this->_Position = stub._Transform * d2;
  KIN_LOG(BF("this->_Position = %s\n") % this->_Position.asString());
#endif
}

void BondedJoint_O::_updateXyzCoords(Stub& stub)
{
  IMPLEMENT_ME();
#if 0
  this->_updateXyzCoord(stub);
  Stub newStub = this->getStub();
  for ( int ii=0; ii < this->_numberOfChildren(); ii++) {
    this->_child(ii)->_updateXyzCoords(newStub);
    this->_DofChangePropagatesToYoungerSiblings = false;
    this->noteXyzUpToDate();
  }
#endif
}


void BondedJoint_O::updateXyzCoords()
{_OF();
	/// dof_change_propagates_to_younger_siblings_ will be set to false inside
	/// update_xyz_coords -- keep a local copy.
  bool local_dof_change_propagates_to_younger_siblings( this->_DofChangePropagatesToYoungerSiblings );
 
	/// Ancestral coordinates are up-to-date since this node
	/// is the root of a subtree that needs refolding.
	/// The stub is passed to update_xyz_coords, and this atom will modify it;
	/// after which the stub is ready to be passed to the younger siblings.
  Stub stub = this->getInputStub();
  this->_updateXyzCoords(stub);
#if 0  
  this->BondedJoint_O::_updateXyzCoords(stub);
  if ( local_dof_change_propagates_to_younger_siblings )
  {
    ASSERTF(this->_Parent.boundp(),BF("Parent is not defined"));
    Joint_sp parent = this->_Parent;
    int ii = 0;
	    /// you had better find yourself in your parent's atom list.
    while ( parent->_child(ii) != this->asSmartPtr() )
    {
      ++ii;
      ASSERTF(ii != parent.get()->_numberOfChildren(),
              BF("While iterating over all younger siblings I hit the end"
                 " - this should never happen, I should hit myself"));
    }
    while ( ii != parent.get()->_numberOfChildren() )
    {
      parent->_child(ii)->_updateXyzCoords(stub,at);
      ++ii;
    }
  }
#endif
}


    double BondedJoint_O::dof(DofType const& dof) const
    {_OF();
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
	SIMPLE_ERROR(BF("Illegal dof request for BondedJoint - I can only handle internal dofs not rigid body"));
    }


    core::Symbol_sp BondedJoint_O::typeSymbol() const {_OF(); return _sym_bonded;};




};
