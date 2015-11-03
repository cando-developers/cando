#ifndef kinematics_dofType_H
#define kinematics_dofType_H

namespace kinematics
{

    /*! Degree of freedom dofType include internal coordinates (distance,theta/angle,phi/dihedral)
      and rigid body translations (rbx,rby,rbz) and rigid body euler rotations (rbalpha,rbbeta,rbgamma) */

    class DofType
    {
    public:
	enum Dof { distance, theta, phi, rbx, rby, rbz, rbalpha, rbbeta, rbgamma, undefined };
    public:

    DofType() : _Dof(undefined) {};
    DofType(Dof const& d) : _Dof(d) {};


	bool operator==(DofType const& dof) const
	{
	    return this->_Dof == dof._Dof;
	}

	bool operator!=(DofType const& dof) const
	{
	    return this->_Dof != dof._Dof;
	}

	bool operator<(DofType const& dof) const
	{
	    return this->_Dof < dof._Dof;
	}

	/*! Return true if we are a rigid body type */
	bool isRigidBodyDof() const;

	string asString() const;

	Dof dof() const { return this->_Dof;};
	
    private:
	Dof	_Dof;
    };
};



#endif // kinematics_dofType_H
