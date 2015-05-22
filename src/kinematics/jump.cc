#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "jump.h"
#include "stub.h"
#include "core/wrappers.h"
namespace kinematics
{



    Real RigidBodyDelta::dof(DofType const& dof) const
    {
	switch (dof.dof())
	{
	case DofType::rbx: return _X;
	case DofType::rby: return _Y;
	case DofType::rbz: return _Z;
	case DofType::rbalpha: return _Alpha;
	case DofType::rbbeta: return _Beta;
	case DofType::rbgamma: return _Gamma;
	default:
	    THROW_HARD_ERROR(BF("Illegal dof for RigidBodyDelta[%s]") % dof.asString() );
	}
    };

    void Jump::fromStubs(Stub const& stub1, Stub const& stub2, core::Lisp_sp lisp)
    {
	Matrix	toCanonical = stub1.getTransform().invertTransform();
	Matrix fromCanonical = stub2.getTransform();
	this->_Transform = fromCanonical*toCanonical;
    }


    void Jump::makeJump(Stub const& stub1, Stub& stub2 ) const
    {
	// make temporary local copy of our rotation-translation
	Matrix tmp_rt(this->_Transform);
 
#if 0
	// n2c
	tmp_rt.fold_in_rb_deltas( rb_delta[1], rb_center[1] );
	tmp_rt.reverse();
 
	// c2n
	tmp_rt.fold_in_rb_deltas( rb_delta[2], rb_center[2] );
 
	// back to original direction
	tmp_rt.reverse();
#endif
//	tmp_rt.make_jump( stub1, stub2 );
	stub2 = stub1.getTransform() * tmp_rt;
    }





    Real Jump::getRigidBodyDelta(DofType const& dof, RigidBodyDirection const& dir) const
    {
	if ( dir == RigidBodyDirection::n2c )
	{
	    return _N2cDelta.dof(dof);
	}
	return _C2nDelta.dof(dof);
    }


    






}; /* kinematics */
