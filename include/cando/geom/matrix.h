/*
    File: matrix.h
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
       
       
//
// (C) 2004 Christian E. Schafmeister
//

//
// Matrix.h
//


#ifndef	MATRIX_H
#define	MATRIX_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <iostream>
#include <string>
#include <vector>
#include <boost/array.hpp>
#include <cando/chem/nVector.fwd.h>
#include <cando/geom/vector3.h>
#include <cando/geom/vector4.h>

class Matrix {

public:

    /*! ctor that initializes the Matrix with identity */
	Matrix();
	/*! ctor that initializes Matrix with identity if passed true and leaves it uninitialized if false */
	Matrix(bool identity);

	virtual ~Matrix();
	Matrix operator*( const Matrix& m ) const;
//	Matrix operator*( Matrix m );
	Vector3 operator*( const Vector3& v ) const;
//	Matrix operator*( double v );
	Matrix operator*( double v ) const;
	Matrix operator+( const Matrix& m );
	Matrix operator-( const Matrix& m );
	Matrix multiplyByMatrix( const Matrix& m ) { return (*this)*m; };
	Vector3 multiplyByVector3( const Vector3& m ) const{ return (*this)*m; };
//	Vector3 multiply3x3ByVector3( Vector3 m );

	Matrix multiply3x3(Matrix const& other) const;
	/*! Multiply the 3x3 matrix in ourselves by other and copy our
	  translation into the result */
	Matrix multiply3x3SaveTranslation(Matrix const& other) const;

	/*! In place add to the translation part */
	void addToTranslation(const Vector3& offset);

  void transform_nvector_point(double& xdestination, double& ydestination, double& zdestiniation, chem::NVector_sp source, size_t source_index) const;

	Matrix	multiplyByScalar( double m ) { return (*this)*m; };

	/*! Return a transposed matrix */
	Matrix transpose() const;
	Matrix transposed() const { return this->transpose();};

	/*! Transpose the upper-left 3x3 and ignore the rest */
	Matrix transposed3x3() const;


	/*! Return the upper 3x3 part */
	Matrix extract3x3() const;
	core::VectorDoubles asOpenGLMatrix();
	core::VectorDoubles asVectorDoubles();

	/*! Set this matrix to the identity matrix */
	void setToIdentity();
	/*! Set the value of each element to its index - for testing */
	void setToIndexAsValue();

	/*! Invert the transform */
	Matrix invertTransform() const;

	/*! Return a matrix that has had the identity matrix subtracted from it */
	Matrix subtractIdentity() const;

	/*! Set the translation part of this matrix */
	void translate(const Vector3& v);
	void scale(double s);
	void rotationAxis( double radians, Vector3* vunit);
	void rotationX(double radians);
	void rotationY(double radians);
	void rotationZ(double radians);
	void rightHandedRotationX(double radians);
	void rightHandedRotationY(double radians);
	void rightHandedRotationZ(double radians);
	double trace();

        inline double& operator[](size_t i) { return this->elements[i]; };
        inline double const & operator[](size_t i) const { return this->elements[i]; };
        
	/*! Return true if the upper 3x3 rotation matrix is orthogonal within tol */
	bool is3x3Orthogonal(double tol) const;

		/*!
		 * Return the index of the first element of the row
		 */
	inline uint rowIndex(int r) const { return (r<<2);};
		/*!
		 * Return the index of the element at row,col
		 */
	inline uint index(int r, int c) const { return (r<<2)+c;};
	inline double atRowCol( int r, int c ) const {return this->elements[this->index(r,c)];};
	inline void atRowColPut( int r, int c, double d ) {this->elements[this->index(r,c)]=d;};
	inline void atColRowPut( int c, int r, double d ) {this->elements[this->index(r,c)]=d;};
	inline void atRowColAdd( int r, int c, double d ) {this->elements[this->index(r,c)]+=d;};
	Vector4	getCol(int c) { return Vector4( this->elements[c],
					this->elements[this->index(1,c)],
					this->elements[this->index(2,c)],
					this->elements[this->index(3,c)]
			    );};
	Vector4	getRow(int r) { uint rs = this->rowIndex(r);
					return Vector4(
					this->elements[rs+0],
					this->elements[rs+1],
					this->elements[rs+2],
					this->elements[rs+3]); };
	void dump();
    void writeToStream( string prefix, std::ostream& out);
    void writeToXMLStream( string prefix, std::ostream& out);
	void setFromDoubleVector( vector<double> vd );
	void setFromDoubleArray( int iElements, double vals[] );
	void setFromQuaternion( Vector4& quat );
	void setValue( double val );
	Vector3 getTranslation() const;

	/*! Set the X column */
	Matrix& colX(Vector3 const& v);
	/*! Set the Y column */
	Matrix& colY(Vector3 const& v);
	/*! Set the Z column */
	Matrix& colZ(Vector3 const& v);

	/*! Set the translation column, return ourselves */
	Matrix& setTranslate(Vector3 const& v);
	
	/*! Set the translation column from the other matrix */
	void setTranslate(const Matrix& m);

	/*! Return the X column */
	Vector3 colX() const;
	/*! Return the Y column */
	Vector3 colY() const;
	/*! Return the Z column */
	Vector3 colZ() const;

//	double	*getElements() { return &this->elements.data(); };

		//
		// Eigen vectors are returned in the rows of evecs
		//
	void	eigenSystem(Vector4& evals, Matrix& evecs );

	void setFromString(const string& s);
//	void setFromStringFast(const string& s,const core::Lisp_sp& lisp);

	string	asString() const;
	string	asStringFormatted() const;
//	string	asStringFast() const;

	void	setAll(const Matrix & p);

#ifdef	USEBOOSTPYTHON
	void		boost_set( boost::python::list names );
//	void		boost_eigenSystem( boost::python::list evals, Matrix& evecs );
#endif

public:
	typedef boost::array<double,16>::iterator	iterator;
// Data
public:
	boost::array<double,16> elements;


};



inline Matrix XRotationMatrixRadians(const double& rad)
{
    Matrix m(false);
    m.rotationX(rad);
    return m;
}

inline Matrix YRotationMatrixRadians(const double& rad)
{
    Matrix m(false);
    m.rotationY(rad);
    return m;
}

inline Matrix ZRotationMatrixRadians(const double& rad)
{
    Matrix m(false);
    m.rotationZ(rad);
    return m;
}








class ArrayOfMatrices
{
private:
	vector<Matrix>	_Matrices;
private:
	void push_back(Matrix& m) { this->_Matrices.push_back(m); };
	void clear() { this->_Matrices.clear(); };
public:
	Matrix& get(uint i);
	void set(uint i, const Matrix& m) { HARD_ASSERT(i<this->_Matrices.size()); this->_Matrices[i]=m;}
	void setSize(uint sz) { this->_Matrices.resize(sz); };
	ArrayOfMatrices();
	virtual ~ArrayOfMatrices();
};


void quaternion_to_matrix(Matrix& matrix, double w, double x, double y, double z, double tx, double ty, double tz );
void normalized_quaternion_to_matrix(Matrix& matrix, double w, double x, double y, double z, double tx, double ty, double tz );
void rotation_matrix_to_quaternion(double& w, double& x, double& y, double& z, const Matrix& m);


Vector3 pointFromMatrixAndInternalCoordinates(const Matrix& coordinateSystem,
                                              double distance,
                                              double angle,
                                              double dihedral,
                                              Vector3& d2);

void internalCoordinatesFromPointAndCoordinateSystem(const Vector3& D, const Matrix& coordinateSystem, double& distance, double& theta, double& phi);


#endif
