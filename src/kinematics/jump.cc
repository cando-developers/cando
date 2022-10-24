/*
    File: jump.cc
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
#define	DEBUG_LEVEL_NONE

#include <clasp/core/common.h>
#include <cando/kinematics/jump.h>
#include <cando/kinematics/stub.h>
#include <clasp/core/wrappers.h>
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
	    THROW_HARD_ERROR(("Illegal dof for RigidBodyDelta[%s]") , dof.asString() );
	}
    };

    void Jump::fromStubs(Stub const& stub1, Stub const& stub2)
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
