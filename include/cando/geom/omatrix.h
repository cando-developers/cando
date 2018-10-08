/*
    File: omatrix.h
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
#ifndef	OMATRIX_H
#define	OMATRIX_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <cando/geom/geomPackage.fwd.h>
#include <cando/geom/matrix.h>

namespace geom {

  SMART(OMatrix );
  class OMatrix_O : public core::CxxObject_O
  {
    LISP_CLASS(geom,GeomPkg,OMatrix_O,"M4",core::CxxObject_O);
  public:
    Matrix	_Value;
  public:
    virtual core::List_sp encode();
    virtual void decode(core::List_sp);
  public:
        //! Create a 4x4 matrix uninitialized (identity=false) or identity if true
    static OMatrix_sp make(bool identity);
  public:
    string __repr__() const;

	//! Create a clone of this matrix
    OMatrix_sp clone() const;

    void setFromVector(core::Array_sp vector);
    core::Array_sp asVector() const;
    
    void setToIdentity();
    void setToIndexAsValue();
    void translate(Vector3 pos);

    void setFromString(const string& s);
    void setFromStringFast(const string& s);
    string asString();
    string asStringFast();

    CL_DEFMETHOD void set_from_quaternion(double qw, double qx, double qy, double qz, double tx, double ty, double tz) {
      quaternion_to_matrix(this->_Value,qw,qx,qy,qz,tx,ty,tz);
    }

    CL_DEFMETHOD void set_from_normalized_quaternion(double qw, double qx, double qy, double qz, double tx, double ty, double tz) {
      normalized_quaternion_to_matrix(this->_Value,qw,qx,qy,qz,tx,ty,tz);
    }

    CL_DEFMETHOD core::T_mv rotation_to_quaternion() const {
      double qw, qx, qy, qz;
      rotation_matrix_to_quaternion(qw,qx,qy,qz,this->_Value);
      return Values(core::DoubleFloat_O::create(qw),
                    core::DoubleFloat_O::create(qx),
                    core::DoubleFloat_O::create(qy),
                    core::DoubleFloat_O::create(qz));
    }

    void rotationX(double radians);
    void rotationY(double radians);
    void rotationZ(double radians);

    Matrix m_TIMES_m(const Matrix& other) const;
    Vector3 m_TIMES_v(const Vector3& vec) const;

CL_LISPIFY_NAME("extract3x3");
CL_DEFMETHOD     Matrix extract3x3() const { return this->_Value.extract3x3();};
CL_LISPIFY_NAME("transposed3x3");
CL_DEFMETHOD     Matrix transposed3x3() const { return this->_Value.transposed3x3();};

CL_LISPIFY_NAME("atRowColGet");
CL_DEFMETHOD     double atRowColGet(int row, int col) { return this->_Value.atRowCol(row,col);};
CL_LISPIFY_NAME("atRowColPut");
CL_DEFMETHOD     void atRowColPut(int row, int col, double v) { return this->_Value.atRowColPut(row,col,v);};
    Matrix getMatrix() { return this->_Value; };
    void setAll(const Matrix& m);
    const Matrix& ref() const { return this->_Value;};

	/*! Return true if the upper 3x3 rotation matrix is orthogonal */
CL_LISPIFY_NAME("is3x3Orthogonal");
CL_DEFMETHOD     bool is3x3Orthogonal(double tol) const { return this->_Value.is3x3Orthogonal(tol);};

    OMatrix_O() {};
  OMatrix_O(bool identity) : _Value(identity) {};
    virtual ~OMatrix_O() {};
  };



}; // namespace geom



namespace translate
{

    template <>
    struct	from_object<Matrix>
    {
	typedef	Matrix		ExpectedType;
	typedef	Matrix		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o) : _v(o.as<geom::OMatrix_O>()->ref()) {};
    };

    template <>
    struct	from_object<const Matrix&>
    {
	typedef	Matrix		ExpectedType;
	typedef	const Matrix&		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o) : _v(o.as<geom::OMatrix_O>()->ref()) {};
#if 0
	{
	    geom::OMatrix_sp v = o.as<geom::OMatrix_O>();
	    _v.setAll(v->getMatrix());
	}
#endif
    };


    template <>
    struct	to_object<Matrix>
    {
	typedef	geom::OMatrix_sp ExpectedType;
	typedef	geom::OMatrix_sp DeclareType;
	static core::T_sp convert(Matrix pos)
	{_G();
	    geom::OMatrix_sp ov = geom::OMatrix_O::create();
	    ov->setAll(pos);
	    return (ov);
	}
    };




    template <>
    struct	to_object<const Matrix& >
    {
	typedef	geom::OMatrix_sp ExpectedType;
	typedef	geom::OMatrix_sp DeclareType;
	static core::T_sp convert(const Matrix& pos)
	{_G();
	    geom::OMatrix_sp ov = geom::OMatrix_O::create();
	    ov->setAll(pos);
	    return (ov);
	}
    };



    template <>
    struct	to_object<Matrix& >
    {
	typedef	geom::OMatrix_sp ExpectedType;
	typedef	geom::OMatrix_sp DeclareType;
	static core::T_sp convert(Matrix& pos)
	{_G();
	    geom::OMatrix_sp ov = geom::OMatrix_O::create();
	    ov->setAll(pos);
	    return (ov);
	}
    };


};





#endif
