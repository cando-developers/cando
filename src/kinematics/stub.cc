/*
    File: stub.cc
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
#include <clasp/core/environment.h>
#include <cando/kinematics/stub.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//

    Stub::Stub(Vector3 const& center,
	       Vector3 const& a,
	       Vector3 const& b,
	       Vector3 const& c)
    {
	this->fromFourJoints(center,a,b,c);
    }



    bool Stub::isOrthogonal(double tol) const
    {
	return this->_Transform.is3x3Orthogonal(tol);
    }


    void Stub::fromFourJoints(Vector3 const& center,
			      Vector3 const& a,
			      Vector3 const& b,
			      Vector3 const& c)
    {
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
