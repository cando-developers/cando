#ifndef	_kinematics_Stub_H
#define _kinematics_Stub_H

#include <clasp/core/common.h>
#include <cando/geom/matrix.h>
#include <cando/kinematics/stub.fwd.h>
#include <cando/kinematics/kinematicsPackage.h>

namespace kinematics
{


class Stub
{
private: // instance variables here
	Matrix	_Transform;

public: // Functions here

	Stub() { this->_Transform.setToIdentity();};

	Stub(const Matrix& mat) { this->_Transform = mat;};

	Stub(Vector3 const& v1,
	     Vector3 const& v2,
	     Vector3 const& v3,
	     Vector3 const& v4);

	/*! Return true if the stub is orthogonal to the tolerance (tol) */
	bool isOrthogonal(double tol) const;

	Matrix const& getTransform() const { return this->_Transform;};

	/*! Define a stub from four points
	  A center and three points a,b,c
	  Orthogonal coordinate frame contains three unit vectors by column,
	  the first one is the unit vector from b pointing to a
	  the second is the unit vector in the plane defined by vector b->a and perpendicular to b->a
	  and the third is the cross product of the first two. */
	void fromFourPoints(Vector3 const& center,
			    Vector3 const& a,
			    Vector3 const& b,
			    Vector3 const& c );

	/*! Return the translation part of the Stub */
	Vector3 translation() const { return this->_Transform.getTranslation();};

	/*! Return the rotation part of the stub */
	Matrix rotation() const { return this->_Transform.extract3x3();};

	/*! Multiply the rotation part by the given matrix */
	void multiplyRotationPart(const Matrix& multiplier);

	/*! Return the colX of the Stub */
	Vector3 colX() const { return this->_Transform.colX();};
	/*! Return the colY of the Stub */
	Vector3 colY() const { return this->_Transform.colY();};
	/*! Return the colZ of the Stub */
	Vector3 colZ() const { return this->_Transform.colZ();};


	/*! Add to the translation part */
	void addToTranslation(const Vector3& offset);


	/*! Return the rotation part colX */
	Vector3 rotationColX() const;




	string asString() const { return this->_Transform.asStringFormatted(); };
};

}; /* kinematics */


#endif /* _kinematics_Stub_H */


