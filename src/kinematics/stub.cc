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
#include <cando/kinematics/stub.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//

#if 0
Stub::Stub(Vector3 const& center,
           Vector3 const& a,
           Vector3 const& b,
           Vector3 const& c)
{
  this->fromFourPoints(center,a,b,c);
}
#endif


bool Stub::isOrthogonal(double tol) const
{
  return this->_Transform.is3x3Orthogonal(tol);
}


/** Generate a Stub using four points centered on center.
    The Stub axes are e1_ = normalized(b-c).
                      e3_ = normalized(crossProduct(e1_,(a-c)))
                      e2_ = crossProduct(e3_,e1_)
    Stub is colX(e1).colY(e2).colZ(e3).setTranslate(c)
    c = parent
    b = grand-parent
    a = great-grand-parent
    SEE DIHEDRAL_DEFINITION
*/
void Stub::fromThreePoints(Vector3 const& parent,
                           Vector3 const& grand_parent,
                           Vector3 const& great_grand_parent)
{
  geom::stubFromThreePoints(this->_Transform, parent, grand_parent, great_grand_parent);
}

void Stub::fromCenterAndRotation(const Vector3& center,
                                 const Matrix& transform) {
  this->_Transform.colX(transform.colX())
      .colY(transform.colY())
      .colZ(transform.colZ())
      .setTranslate(center);
}

void Stub::fromMatrix(const Matrix& transform) {
  this->_Transform = transform;
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
