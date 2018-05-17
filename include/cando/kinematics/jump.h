/*
    File: jump.h
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
#ifndef	kinematics_jump_H
#define kinematics_jump_H

#include <clasp/core/foundation.h>
#include <cando/geom/matrix.h>
#include <cando/kinematics/stub.fwd.h>
#include <cando/kinematics/dofType.h>
#include <cando/kinematics/kinFoundation.h>

namespace kinematics
{

    class RigidBodyDirection
    {
    public:
	enum Direction { n2c, c2n };
    public:
    RigidBodyDirection(Direction const& d) : _Dir(d) {};
    RigidBodyDirection() : _Dir(n2c) {};

	bool operator==(Direction const& dir) const { return this->_Dir == dir;};
	inline string asString() const
	{
	    if ( _Dir == n2c ) return "n2c";
	    return "c2n";
	}
    private:
	Direction _Dir;
    };


    class RigidBodyDelta
    {
    private:
	Real	_X;
	Real	_Y;
	Real	_Z;
	Real	_Alpha;
	Real	_Beta;
	Real	_Gamma;

    public:
    RigidBodyDelta() : _X(0.0), _Y(0.0), _Z(0.0), _Alpha(0.0), _Beta(0.0), _Gamma(0.0) {};

	Real dof(DofType const& dof) const;
    };


    class Jump
    {
    protected:
	Matrix		_Transform;
	RigidBodyDelta	_N2cDelta;
	RigidBodyDelta	_C2nDelta;
	Vector3		_N2cCenter;
	Vector3		_C2nCenter;
    public:
	/*! Calculate the Jump from the stubs */
	void fromStubs(Stub const& stub1, Stub const& stub2);

	/*! Dump my transform */
	string asString() const { return this->_Transform.asStringFormatted();};


	/* Make a jump from stub1 to stub2 */
	void makeJump(Stub const& stub1, Stub& stub2 ) const;


	/*! Return the dof */
	Real getRigidBodyDelta(DofType const& dof, RigidBodyDirection const& dir) const;

    };



};
#endif
