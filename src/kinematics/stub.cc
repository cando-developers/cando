#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "stub.h"
#include "core/wrappers.h"
namespace kinematics
{

// ----------------------------------------------------------------------
//

    Stub::Stub(Vector3 const& center,
	       Vector3 const& a,
	       Vector3 const& b,
	       Vector3 const& c)
    {_G();
	this->fromFourPoints(center,a,b,c);
    }



    bool Stub::isOrthogonal(double tol) const
    {
	return this->_Transform.is3x3Orthogonal(tol);
    }


    void Stub::fromFourPoints(Vector3 const& center,
			      Vector3 const& a,
			      Vector3 const& b,
			      Vector3 const& c)
    {_G();
	this->_Transform.setToIdentity();
	LOG(BF("center = %s") % center.asString());
	LOG(BF("a = %s") % a.asString());
	LOG(BF("b = %s") % b.asString());
	LOG(BF("c = %s") % c.asString());
	Vector3 e1(a-b);
	e1 = e1.normalized();
	LOG(BF("e1 = (a-b).normalized : %s") % e1.asString() );
	Vector3 e3(e1.crossProduct(c-b));
	e3 = e3.normalized();
	LOG(BF("e3 = (e1.crossProduct(c-b)).normalized : %s") % e3.asString() );
	Vector3 e2(e3.crossProduct(e1));
	LOG(BF("e2 = (e3.crossProduct(e1)): %s") % e2.asString() );
	LOG(BF("Stub before being set: %s") % this->_Transform.asStringFormatted());
	this->_Transform.colX(e1).colY(e2).colZ(e3).setTranslate(center);
	LOG(BF("Stub = \n%s") % this->_Transform.asStringFormatted());
    }


    void Stub::multiplyRotationPart(const Matrix& multiplier)
    {
	this->_Transform.multiply3x3SaveTranslation(multiplier);
    }
	    

    void Stub::addToTranslation(const Vector3& offset)
    {
	this->_Transform.addToTranslation(offset);
    }

    Vector3 Stub::rotationColX() const
    {
	return this->_Transform.colX();
    }


}; /* kinematics */
