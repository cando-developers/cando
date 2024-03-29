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

This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University
at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */

//
// (C) 2004 Christian E. Schafmeister
//

//
// Matrix.h
//

#ifndef MATRIX_H
#define MATRIX_H

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
  Matrix& operator=(const Matrix& other);
  Matrix operator*(const Matrix &m) const;
  //	Matrix operator*( Matrix m );
  Vector3 operator*(const Vector3 &v) const;
  //	Matrix operator*( double v );
  Matrix operator*(double v) const;
  Matrix operator+(const Matrix &m) const;
  const Matrix &operator+=(const Matrix &m);
  Matrix operator-(const Matrix &m) const;
  Matrix multiplyByMatrix(const Matrix &m) { return (*this) * m; };
  Vector3 multiplyByVector3(const Vector3 &m) const { return (*this) * m; };
  //	Vector3 multiply3x3ByVector3( Vector3 m );

  Matrix multiply3x3(Matrix const &other) const;
  /*! Multiply the 3x3 matrix in ourselves by other and copy our
    translation into the result */
  Matrix multiply3x3SaveTranslation(Matrix const &other) const;

  /*! In place add to the translation part */
  void addToTranslation(const Vector3 &offset);

  void transform_nvector_point(chem::Vector_real &xdestination,
                               chem::Vector_real &ydestination, 
                               chem::Vector_real &zdestiniation,
                               chem::NVector_sp source,
                               size_t source_index) const;

  Matrix multiplyByScalar(double m) { return (*this) * m; };

  /*! Return a transposed matrix */
  Matrix transpose() const;
  Matrix transposed() const { return this->transpose(); };

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
  void translate(const Vector3 &v);
  void scale(double s);
  void rotationAxis(double radians, Vector3 *vunit);
  void rotationX(double radians);
  void rotationY(double radians);
  void rotationZ(double radians);
  void rightHandedRotationX(double radians);
  void rightHandedRotationY(double radians);
  void rightHandedRotationZ(double radians);
  double trace();

  Matrix flipXY() const;
  inline double &operator[](size_t i) { return this->elements[i]; };
  inline double const &operator[](size_t i) const { return this->elements[i]; };

  /*! Return true if the upper 3x3 rotation matrix is orthogonal within tol */
  bool is3x3Orthogonal(double tol) const;

  /*!
   * Return the index of the first element of the row
   */
  inline uint rowIndex(int r) const { return (r << 2); };
  /*!
   * Return the index of the element at row,col
   */
  inline uint index(int r, int c) const { return (r << 2) + c; };
  inline double at(int r, int c) const { return this->elements[this->index(r, c)]; };
  inline double &at(int r, int c) { return this->elements[this->index(r, c)]; };
  Vector4 getCol(int c) {
    return Vector4(this->elements[c], this->elements[this->index(1, c)], this->elements[this->index(2, c)],
                   this->elements[this->index(3, c)]);
  };
  Vector4 getRow(int r) {
    uint rs = this->rowIndex(r);
    return Vector4(this->elements[rs + 0], this->elements[rs + 1], this->elements[rs + 2], this->elements[rs + 3]);
  };
  void setRow(int r, double v0, double v1, double v2, double v3 ) {
    uint rs = this->rowIndex(r);
    this->elements[rs+0] = v0;
    this->elements[rs+1] = v1;
    this->elements[rs+2] = v2;
    this->elements[rs+3] = v3;
  };
  void setCol(int c, double v0, double v1, double v2, double v3 ) {
    this->elements[c] = v0;
    this->elements[this->index(1,c)] = v1;
    this->elements[this->index(2,c)] = v2;
    this->elements[this->index(3,c)] = v3;
  };
  void dump();
  void writeToStream(string prefix, std::ostream &out);
  void writeToXMLStream(string prefix, std::ostream &out);
  void setFromDoubleVector(vector<double> vd);
  void setFromDoubleArray(int iElements, double vals[]);
  void setFromQuaternion(Vector4 &quat);
  void setValue(double val);
  Vector3 getTranslation() const;

  /*! Set the X column */
  Matrix &colX(Vector3 const &v);
  /*! Set the Y column */
  Matrix &colY(Vector3 const &v);
  /*! Set the Z column */
  Matrix &colZ(Vector3 const &v);

  /*! Set the translation column, return ourselves */
  Matrix &setTranslate(Vector3 const &v);

  /*! Set the translation column from the other matrix */
  void setTranslate(const Matrix &m);

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
  void eigenSystem(Vector4 &evals, Matrix &evecs);

  void setFromString(const string &s);
  //	void setFromStringFast(const string& s,const core::LispPtr& lisp);

  string asString() const;
  string asStringFormatted() const;
  //	string	asStringFast() const;

  void setAll(const Matrix &p);

  core::List_sp encode() const;
  void decode(core::List_sp alist);

#ifdef USEBOOSTPYTHON
  void boost_set(boost::python::list names);
//	void		boost_eigenSystem( boost::python::list evals, Matrix& evecs );
#endif

public:
  typedef boost::array<double, 16>::iterator iterator;
  // Data
public:
  boost::array<double, 16> elements;
};

inline Matrix XRotationMatrixRadians(const double &rad) {
  Matrix m(false);
  m.rotationX(rad);
  return m;
}

inline Matrix YRotationMatrixRadians(const double &rad) {
  Matrix m(false);
  m.rotationY(rad);
  return m;
}

inline Matrix ZRotationMatrixRadians(const double &rad) {
  Matrix m(false);
  m.rotationZ(rad);
  return m;
}

class ArrayOfMatrices {
private:
  vector<Matrix> _Matrices;

private:
  void push_back(Matrix &m) { this->_Matrices.push_back(m); };
  void clear() { this->_Matrices.clear(); };

public:
  Matrix &get(uint i);
  void set(uint i, const Matrix &m) {
    HARD_ASSERT(i < this->_Matrices.size());
    this->_Matrices[i] = m;
  }
  void setSize(uint sz) { this->_Matrices.resize(sz); };
  ArrayOfMatrices();
  virtual ~ArrayOfMatrices();
};

void quaternion_to_matrix(Matrix &matrix, double w, double x, double y, double z, double tx, double ty, double tz);
void normalized_quaternion_to_matrix(Matrix &matrix, double w, double x, double y, double z, double tx, double ty, double tz);
void rotation_matrix_to_quaternion(double &w, double &x, double &y, double &z, const Matrix &m);

void rotor3_to_matrix(Matrix &matrix, double scalar, double xy, double yz, double zx, double tx, double ty, double tz);
void normalized_rotor3_to_matrix(Matrix &matrix, double scalar, double xy, double yz, double zx, double tx, double ty, double tz);
void rotation_matrix_to_rotor3(double &scalar, double &xy, double &yz, double &zx, const Matrix &m);

namespace geom {
Vector3 pointFromStubAndInternalCoordinates(const Matrix &stub, double distance, double angle, double dihedral,
                                              Vector3 &d2);

void internalCoordinatesFromPointAndStub(const Vector3 &D, const Matrix &stub, double &distance,
                                                     double &theta, double &phi);

void stubFromThreePoints(
    Matrix& transform,
    Vector3 const& parent,
    Vector3 const& grand_parent,
    Vector3 const& great_grand_parent );


Vector3	geom__build_origin();
Vector3	geom__build_using_bond( double distance, const Vector3& vb );
//! Build a vector at distance from vb and angle from v
Vector3 geom__build_using_bond_angle( double distance, const Vector3& vb,
                                      double angleRad, const Vector3& va );
Vector3 geom__build_using_bond_angle_dihedral( double distance, const Vector3& vc,
                                               double angleRad, const Vector3& vb,
                                               double dihedralRad, const Vector3& va );

double	calculateDihedral( const Vector3& va,
                           const Vector3& vb,
                           const Vector3& vc,
                           const Vector3& vd );
double	calculateAngle( const Vector3& va,
                        const Vector3& vb,
                        const Vector3& vc );

double  geom__planeVectorAngle(double dx, double dy);

};

#endif
