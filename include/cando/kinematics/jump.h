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
	void fromStubs(Stub const& stub1, Stub const& stub2, core::Lisp_sp lisp);

	/*! Dump my transform */
	string asString() const { return this->_Transform.asStringFormatted();};


	/* Make a jump from stub1 to stub2 */
	void makeJump(Stub const& stub1, Stub& stub2 ) const;


	/*! Return the dof */
	Real getRigidBodyDelta(DofType const& dof, RigidBodyDirection const& dir) const;

    };



};
#endif
