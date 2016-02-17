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
    string __str__() ;

	//! Create a clone of this matrix
    OMatrix_sp clone() const;

    void setToIdentity();
    void setToIndexAsValue();
    void translate(Vector3 pos);

    void setFromString(const string& s);
    void setFromStringFast(const string& s);
    string asString();
    string asStringFast();

    void rotationX(double radians);
    void rotationY(double radians);
    void rotationZ(double radians);

    Matrix mul_matrix(const Matrix& other) const;
    Vector3 mul_v3(const Vector3& vec) const;

CL_NAME("extract3x3");
CL_DEFMETHOD     Matrix extract3x3() const { return this->_Value.extract3x3();};
CL_NAME("transposed3x3");
CL_DEFMETHOD     Matrix transposed3x3() const { return this->_Value.transposed3x3();};

CL_NAME("atRowColGet");
CL_DEFMETHOD     double atRowColGet(int row, int col) { return this->_Value.atRowCol(row,col);};
CL_NAME("atRowColPut");
CL_DEFMETHOD     void atRowColPut(int row, int col, double v) { return this->_Value.atRowColPut(row,col,v);};
    Matrix getMatrix() { return this->_Value; };
    void setAll(const Matrix& m);
    const Matrix& ref() const { return this->_Value;};

	/*! Return true if the upper 3x3 rotation matrix is orthogonal */
CL_NAME("is3x3Orthogonal");
CL_DEFMETHOD     bool is3x3Orthogonal(double tol) const { return this->_Value.is3x3Orthogonal(tol);};

    OMatrix_O() {};
  OMatrix_O(bool identity) : _Value(identity) {};
    virtual ~OMatrix_O() {};
  };




}; // namespace geom
TRANSLATE(geom::OMatrix_O);



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


