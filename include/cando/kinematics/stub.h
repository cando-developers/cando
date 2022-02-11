
/*
    File: stub.h
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
#ifndef	_kinematics_Stub_H
#define _kinematics_Stub_H

#include <clasp/core/common.h>
#include <cando/geom/matrix.h>
#include <cando/geom/omatrix.h>
#include <cando/kinematics/stub.fwd.h>
#include <cando/kinematics/kinematicsPackage.h>

namespace kinematics
{


class Stub
{
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

	/*! Define a stub from four Joints
	  A center and three Joints a,b,c
	  Orthogonal coordinate frame contains three unit vectors by column,
	  the first one is the unit vector from b Jointing to a
	  the second is the unit vector in the plane defined by vector b->a and perpendicular to b->a
	  and the third is the cross product of the first two. */
  void fromFourPoints(Vector3 const& center,
                      Vector3 const& a,
                      Vector3 const& b,
                      Vector3 const& c );

  /*! Define a stub from a center and a rotation matrix from the transform.
*/
  void fromCenterAndRotation(const Vector3& center,
                             const Matrix& transform);
  
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

public: // instance variables here
  Matrix	_Transform;
};

}; /* kinematics */

namespace translate
{

template <>
struct	from_object<kinematics::Stub>
{
  typedef	kinematics::Stub		ExpectedType;
  typedef	kinematics::Stub		DeclareType;
  DeclareType _v;
  from_object(core::T_sp o) : _v(o.as<geom::OMatrix_O>()->ref()) {};
};

template <>
struct	from_object<const kinematics::Stub&>
{
  typedef	kinematics::Stub		ExpectedType;
  typedef	kinematics::Stub		DeclareType;
  DeclareType _v;
  from_object(core::T_sp o) : _v(o.as<geom::OMatrix_O>()->ref()) {};
};


template <>
struct	to_object<kinematics::Stub>
{
  typedef	geom::OMatrix_sp ExpectedType;
  typedef	geom::OMatrix_sp DeclareType;
  static core::T_sp convert(kinematics::Stub pos)
  {_G();
    geom::OMatrix_sp ov = geom::OMatrix_O::create();
    ov->setAll(pos.getTransform());
    return (ov);
  }
};




template <>
struct	to_object<const kinematics::Stub& >
{
  typedef	geom::OMatrix_sp ExpectedType;
  typedef	geom::OMatrix_sp DeclareType;
  static core::T_sp convert(const kinematics::Stub& pos)
  {_G();
    geom::OMatrix_sp ov = geom::OMatrix_O::create();
    ov->setAll(pos.getTransform());
    return (ov);
  }
};



template <>
struct	to_object<kinematics::Stub& >
{
  typedef	geom::OMatrix_sp ExpectedType;
  typedef	geom::OMatrix_sp DeclareType;
  static core::T_sp convert(kinematics::Stub& pos)
  {_G();
    geom::OMatrix_sp ov = geom::OMatrix_O::create();
    ov->setAll(pos.getTransform());
    return (ov);
  }
};


};


#endif /* _kinematics_Stub_H */
