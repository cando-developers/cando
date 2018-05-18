/*
    File: bondedAtom.cc
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
#define DEBUG_LEVEL_NONE



#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/atomId.h>
#include <clasp/core/numerics.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/atomTree.h>
#include <cando/kinematics/bondedAtom.h>

namespace kinematics
{
void BondedJoint_O::initialize() {
  for (int i=0; i<MaxChildren; ++i ) {
    this->_Children[i].reset_();
  }
}

    void BondedJoint_O::_insertChild(int before, Joint_sp child)
    {
	if ( this->_numberOfChildren()>this->_maxNumberOfChildren() )
	{
	    THROW_HARD_ERROR(BF("You exceede the maximum[%d] number of children allowed for a BondedAtom") % this->_maxNumberOfChildren());
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
        this->_Children[this->_NumberOfChildren].reset_();
    }


    void BondedJoint_O::_releaseAllChildren()
    {
	if ( this->_numberOfChildren() == 0 ) return;
	int num = this->_numberOfChildren();
	for ( int idx=0; idx < num; idx++ )
	{
          this->_Children[idx].reset_();
	}
	this->_NumberOfChildren = 0;
    }



    Joint_sp BondedJoint_O::stubAtom3(AtomTree_sp at) const
    {_OF();
      if ( gc::IsA<JumpJoint_sp>(this->parent()))
	{
	    ASSERT(this->parent().get()->stubDefined());
	    Joint_sp p_stub2( this->parent().get()->stubAtom2() );
	    chem::AtomId const & p_stub2_id( p_stub2.get()->id() );
	    if ( this->id() == p_stub2_id )
	    {
		// very special case!!
		return this->parent().get()->stubAtom3(at);
	    } else {
		return p_stub2;
	    }
	} else {
	    return this->parent().get()->stubAtom2();
	}
    }

    void BondedJoint_O::updateInternalCoords( Stub & stub,
					   bool recursive,
					   AtomTree_sp at)
    {_OF();
	LOG(BF("Incoming stub:\n%s") % stub.asString() );
//	using numeric::x_rotation_matrix_radians;
//	using numeric::z_rotation_matrix_radians;
//	using numeric::constants::d::pi;
 
	ASSERTF( stub.isOrthogonal( 1e-3 ), BF("The stub is not orthogonal: %s") % stub.asString() );
	LOG(BF("My position = %s") % this->position().asString() );
	Vector3 w(this->position()-stub.translation());
	this->_Distance  = w.length();
	LOG(BF("Calculated _Distance=%lf") % this->_Distance );
	bool flip_stub( false );
	if ( this->_Distance < 1e-2 )
	{
	    LOG(BF("Distance is very small, setting _Phi and _Theta to zero"));
	    // phi, theta dont make much sense
	    //  std::cerr << "WARNING:: very small d= " << d_ << ' ' << id() << std::endl;
	    this->_Phi = 0.0;
	    this->_Theta = 0.0;
	} else
	{
	    LOG(BF("Distance is normal, calculating _Theta and _Phi"));
	    w = w.normalized();
	    Real const x( w.dotProduct(stub.colX() ) );
	    Real const y( w.dotProduct(stub.colY() ) );
	    Real const z( w.dotProduct(stub.colZ() ) );
 
	    Real const tol( 1e-6 );
	    if ( x < (-1.0 + tol) )
	    {
		
		// very special case:
		// confirm that we are the stub_atom2 of a jump:
		if ( this->keepDofFixed(DofType::theta,at ) )
		{
		    LOG(BF("Keeping DOF(theta) fixed and flipping stub"));
		    this->_Theta = core::numerics::pi;
		    this->_Phi = 0.0;
		    flip_stub = true; // very special case
		} else {
		    LOG(BF("theta is not fixed"));
		    this->_Theta = core::numerics::pi;
		    this->_Phi = 0.0;
		}
	    } else if ( x > (1.0 - tol) )
	    {
		LOG(BF("x > 1.0-tol"));
		//std::cout << "WARNING:: update_internal_coords: exactly parallel? " << id() << std::endl;
		this->_Theta = 0.0;
		this->_Phi = 0.0;
	    } else {
		this->_Theta = std::acos( x ); // DANGER
		LOG(BF("Calculated _Theta = %lf degrees") % (this->_Theta/0.0174533) );
		if ( this->_Theta < 1e-2 || core::numerics::pi - this->_Theta < 1e-2 )
		{
		    // less than 0.57 degrees
		    //std::cout << "WARNING:: small theta but we are calculating phi: " <<
		    //  theta_ << std::endl;
		}
		this->_Phi  = std::atan2( z, y );
		LOG(BF("Calculated _Phi=%lf degrees") % (this->_Phi/0.0174533) );
	    } // small theta
	} // small d
	stub.multiplyRotationPart(XRotationMatrixRadians(this->_Phi));
	if ( recursive )
	{
	    Matrix zRotationTheta(false);
	    zRotationTheta.rotationZ(this->_Theta);
	    Matrix rotationTrans(stub.rotation());
	    rotationTrans.setTranslate(this->position());
	    Stub newStub(rotationTrans);
	    if ( flip_stub ) {
		// special case if I'm stub_atom2 of my parent (who is a jump)
		newStub.multiplyRotationPart(XRotationMatrixRadians(core::numerics::pi));
	    }
 
	    for ( int it=0; it<this->_numberOfChildren(); it++ )
	    {
		this->_child(it)->updateInternalCoords( newStub, true, at );
	    }
	}
    }


    bool BondedJoint_O::keepDofFixed(DofType dof,AtomTree_sp at) const
    {_OF();
	if ( dof == DofType::distance )
	{
	    return false;
	} else if ( dof == DofType::theta )
	{
          return ( gc::IsA<JumpJoint_sp>(this->parent())
		     && this->id() == this->parent()->stubAtom2Id() );
	} else if ( dof == DofType::phi )
	{
          Joint_sp parent = this->parent();
          if ( gc::IsA<JumpJoint_sp>(parent) &&
               ( this->id() == parent->stubAtom2Id() ||
                 this->id() == parent->stubAtom3Id(at) ) ) return true;
          if ( parent->parent().unboundp() ) return false;
          Joint_sp grandParent = parent->parent();
          return ( gc::IsA<JumpJoint_sp>(grandParent)
                   && this->id() == grandParent->stubAtom3Id(at));
	} else
	{
          SIMPLE_ERROR(BF("BondedJoint_O::keepDofFixed: BAD_DOF: %s") % dof.asString() );
	}
	return false;
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


    void BondedJoint_O::updateXyzCoords(Stub& stub,AtomTree_sp at)
    {_OF();
	ASSERTF(stub.isOrthogonal(1e-3),BF("The Stub is not orthogonal - stub:\n%s") % stub.asString());
	stub.multiplyRotationPart(XRotationMatrixRadians(this->_Phi));
	Stub newStub(stub);
	newStub.multiplyRotationPart(ZRotationMatrixRadians(this->_Theta));
	if ( std::abs(this->_Theta-core::numerics::pi) < 1e-6 )
	{
	    // very special case 
	    if ( this->keepDofFixed(DofType::theta,at))
	    {
		newStub.multiplyRotationPart(XRotationMatrixRadians(core::numerics::pi));
	    }
	}
	newStub.addToTranslation(newStub.rotationColX()*this->_Distance);
	this->position(newStub.translation());
	for ( int ii=0; ii < this->_numberOfChildren(); ii++)
	{
	    this->_child(ii)->updateXyzCoords(newStub,at);
	}
	this->_DofChangePropagatesToYoungerSiblings = false;
	this->noteXyzUpToDate();
    }




    void BondedJoint_O::updateXyzCoords(AtomTree_sp at)
    {_OF();
	/// dof_change_propagates_to_younger_siblings_ will be set to false inside
	/// update_xyz_coords -- keep a local copy.
	bool local_dof_change_propagates_to_younger_siblings( this->_DofChangePropagatesToYoungerSiblings );
 
	/// Ancestral coordinates are up-to-date since this node
	/// is the root of a subtree that needs refolding.
	/// The stub is passed to update_xyz_coords, and this atom will modify it;
	/// after which the stub is ready to be passed to the younger siblings.
	Stub stub( this->getInputStub(at) );
	this->BondedJoint_O::updateXyzCoords(stub,at);
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
		parent->_child(ii)->updateXyzCoords(stub,at);
		++ii;
	    }
	}
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
	SIMPLE_ERROR(BF("Illegal dof request for BondedAtom - I can only handle internal dofs not rigid body"));
    }


    core::Symbol_sp BondedJoint_O::typeSymbol() const {_OF(); return _sym_bonded;};




};
