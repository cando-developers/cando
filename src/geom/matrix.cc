/*
    File: matrix.cc
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
#define DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/karma_real.hpp>
#include <boost/array.hpp>

#include <string.h>
#include <cando/geom/matrix.h>
#include <iostream>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/lisp.h>
#include <cando/chem/nVector.h>
#include <clasp/core/ql.h>
//#include "core/archiveNode.h"

Matrix benchmarkMatrixMultiplications(int num);

#define nrerror(s) printf("%s\n", s);

void jacobi(Matrix &a, int n, double d[], Matrix &v, int *nrot)
// Computes all eigenvalues and eigenectors of a real
// symmetric matrix a[1..n][1..n]. On
// output, elements of a above the diagonal are destroyed.
// d[1..n] returns the eigenvalues of a.
// v[1..n][1..n] is a matrix whose columns contain,
// on output, the normalized eigenvectors of
// a. nrot returns the number of Jacobi rotations that were required.
{
#define ROTATE(a, i, j, k, l)                                                                                                      \
  {                                                                                                                                \
    g = a[a.index(i, j)];                                                                                                          \
    h = a[a.index(k, l)];                                                                                                          \
    a[a.index(i, j)] = g - s * (h + g * tau);                                                                                      \
    a[a.index(k, l)] = h + s * (g - h * tau);                                                                                      \
  }
  int j, iq, ip, i;
#define MAXMATRIX 10
  double tresh, theta, tau, t, sm, s, h, g, c, b[MAXMATRIX], z[MAXMATRIX];
  if (n > MAXMATRIX) {
    nrerror("Matrix size is too large");
  }
  //    b=vector(1,n);
  //    z=vector(1,n);
  for (ip = 0; ip < n; ip++) { // Initialize to the identity matrix.
    for (iq = 0; iq < n; iq++)
      v[v.index(ip, iq)] = 0.0;
    v[v.index(ip, ip)] = 1.0;
  }
  for (ip = 0; ip < n; ip++) { // Initialize b and d to the diagonal of a.
    b[ip] = d[ip] = a[a.index(ip, ip)];
    z[ip] = 0.0; // This vector will accumulate terms of
                 // the form tapq as in equation (11.1.14).
  }
  *nrot = 0;
  for (i = 1; i <= 50; i++) {
    sm = 0.0;
    for (ip = 0; ip < n - 1; ip++) { // Sum o-diagonal elements.
      for (iq = ip + 1; iq < n; iq++)
        sm += fabs(a[a.index(ip, iq)]);
    }
    if (sm == 0.0) { // The normal return, which relies
                     // on quadratic convergence to machine underflow.
                     //	    free_vector(z,1,n);
                     //	    free_vector(b,1,n);
      return;
    }
    if (i < 4)
      tresh = 0.2 * sm / (n * n); // ...on the first three sweeps.
    else
      tresh = 0.0; // ...thereafter.
    for (ip = 0; ip < n - 1; ip++) {
      for (iq = ip + 1; iq < n; iq++) {
        g = 100.0 * fabs(a[a.index(ip, iq)]);
        // After four sweeps, skip the rotation
        // if the o-diagonal element is small.
        if (i > 4 && (double)(fabs(d[ip]) + g) == (double)fabs(d[ip]) && (double)(fabs(d[iq]) + g) == (double)fabs(d[iq]))
          a[a.index(ip, iq)] = 0.0;

        else if (fabs(a[a.index(ip, iq)]) > tresh) {
          h = d[iq] - d[ip];
          if ((double)(fabs(h) + g) == (double)fabs(h))
            t = (a[a.index(ip, iq)]) / h; // t = 1=(2)
          else {
            theta = 0.5 * h / (a[a.index(ip, iq)]); // Equation (11.1.10).
            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0)
              t = -t;
          }
          c = 1.0 / sqrt(1 + t * t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * a[a.index(ip, iq)];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[a.index(ip, iq)] = 0.0;
          for (j = 0; j <= ip - 1; j++) { // Case of rotations 1 <= j < p.
            ROTATE(a, j, ip, j, iq)
          }
          for (j = ip + 1; j <= iq - 1; j++) { // Case of rotations p < j < q.
            ROTATE(a, ip, j, j, iq)
          }
          for (j = iq + 1; j < n; j++) { // Case of rotations q < j n.
            ROTATE(a, ip, j, iq, j)
          }
          for (j = 0; j < n; j++) {
            ROTATE(v, j, ip, j, iq)
          }
          ++(*nrot);
        }
      }
    }
    for (ip = 0; ip < n; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip]; // Update d with the sum of tapq,
      z[ip] = 0.0;   // and reinitialize z.
    }
  }
  THROW_HARD_ERROR(("Too many iterations in routine jacobi"));
}

Matrix::Matrix() { this->setToIdentity(); }

Matrix::Matrix(bool identity) {
  if (identity)
    setToIdentity();
  else {
    for (Matrix::iterator mi = elements.begin(); mi != elements.end(); mi++) {
      *mi = 0.0;
    }
  }
}

Matrix::~Matrix() {}




core::List_sp Matrix::encode() const
{
  ql::list ll;
  for ( size_t ii=0; ii<16; ii++ ) {
    ll << core::make_single_float(this->elements[ii]);
  }
  return core::Cons_O::create(INTERN_(kw,terms),ll.cons());
}

void Matrix::decode(core::List_sp alist)
{
  size_t zi = 0;
  for ( auto cur : (core::List_sp)core::oCdr(alist) ) {
    this->elements[zi] = CONS_CAR(cur).unsafe_single_float();
    zi++;
  }
}


void Matrix::setToIdentity() {
  Matrix::iterator mi = this->elements.begin();
  // Row 0
  *mi = 1.0;
  mi++;
  *mi = 0.0;
  mi++;
  *mi = 0.0;
  mi++;
  *mi = 0.0;
  mi++;
  // Row 1
  *mi = 0.0;
  mi++;
  *mi = 1.0;
  mi++;
  *mi = 0.0;
  mi++;
  *mi = 0.0;
  mi++;
  // Row 2
  *mi = 0.0;
  mi++;
  *mi = 0.0;
  mi++;
  *mi = 1.0;
  mi++;
  *mi = 0.0;
  mi++;
  // Row 3
  *mi = 0.0;
  mi++;
  *mi = 0.0;
  mi++;
  *mi = 0.0;
  mi++;
  *mi = 1.0;
  mi++;
}

void Matrix::scale(double x) {
  Matrix::iterator mi = this->elements.begin();
  // Row 0
  *mi = x;
  mi++;
  *mi = 0.0;
  mi++;
  *mi = 0.0;
  mi++;
  *mi = 0.0;
  mi++;
  // Row 1
  *mi = 0.0;
  mi++;
  *mi = x;
  mi++;
  *mi = 0.0;
  mi++;
  *mi = 0.0;
  mi++;
  // Row 2
  *mi = 0.0;
  mi++;
  *mi = 0.0;
  mi++;
  *mi = x;
  mi++;
  *mi = 0.0;
  mi++;
  // Row 3
  *mi = 0.0;
  mi++;
  *mi = 0.0;
  mi++;
  *mi = 0.0;
  mi++;
  *mi = x;
  mi++;
}

Matrix Matrix::subtractIdentity() const {
  Matrix m(*this);
  m[0 + 0] -= 1.0;
  m[4 + 1] -= 1.0;
  m[8 + 2] -= 1.0;
  m[12 + 3] -= 1.0;
  return m;
}

core::VectorDoubles Matrix::asVectorDoubles() {
  vector<double> vd;
  vd.reserve(this->elements.size());
  std::back_insert_iterator<vector<double>> it(vd);
  Matrix::iterator el;
  el = this->elements.begin();
  for (uint zz = 0; zz < this->elements.size(); zz++) {
    *it++ = *el;
    el++;
  }
  return vd;
}

//
//	asOpenGLMatrix
//
//	Return a vector of doubles transposed
core::VectorDoubles Matrix::asOpenGLMatrix() {
  vector<double> vd;
  vd.reserve(this->elements.size());
  std::back_insert_iterator<vector<double>> it(vd);
  uint r, c;
  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      *it++ = this->elements[this->index(c, r)];
    }
  }
  return vd;
}

void Matrix::setAll(const Matrix &m) { this->elements = m.elements; }

void Matrix::transform_nvector_point(chem::Vector_real &xdest,
                                     chem::Vector_real &ydest,
                                     chem::Vector_real &zdest,
                                     chem::NVector_sp source,
                                     size_t source_index) const {
  double acc;
  acc = at(0, 0) * (*source)[source_index];
  acc += at(0, 1) * (*source)[source_index + 1];
  acc += at(0, 2) * (*source)[source_index + 2];
  acc += at(0, 3);
  xdest = acc;

  acc = at(1, 0) * (*source)[source_index];
  acc += at(1, 1) * (*source)[source_index + 1];
  acc += at(1, 2) * (*source)[source_index + 2];
  acc += at(1, 3);
  ydest = acc;

  acc = at(2, 0) * (*source)[source_index];
  acc += at(2, 1) * (*source)[source_index + 1];
  acc += at(2, 2) * (*source)[source_index + 2];
  acc += at(2, 3);
  zdest = acc;
}

namespace geom {
CL_DOCSTRING(R"dx(Transform a point in an nvector by the matrix and put the result in the destination)dx");
DOCGROUP(cando);
CL_DEFUN
void geom__inplace_transform_nvector_point(Vector3 &destination, const Matrix &transform, chem::NVector_sp source,
                                           size_t source_index) {
  transform.transform_nvector_point(destination.getX(), destination.getY(), destination.getZ(), source, source_index);
}
}; // namespace geom

Matrix& Matrix::operator=(const Matrix& other) {
  if (this==&other) return *this;
  this->elements = other.elements;
  return *this;
}

Matrix Matrix::operator*(const Matrix &m) const {
  Matrix res;
  double acc;
  int r, c, x;
  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      acc = 0.0;
      for (x = 0; x < 4; x++) {
        acc += at(r, x) * m.at(x, c);
      }
      res.at(r, c) = acc;
    }
  }
  return res;
}

Matrix Matrix::operator+(const Matrix &m) const {
  Matrix res;

  for (int r = 0; r < 4; r++) {
    for (int c = 0; c < 4; c++) {
      res.at(r, c) = at(r, c) + m.at(r, c);
    }
  }

  return res;
}

const Matrix& Matrix::operator+=(const Matrix &m) {
  Matrix res;
  int r, c;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      at(r, c) += m.at(r, c);
    }
  }
  return *this;
}

Matrix Matrix::operator-(const Matrix &m) const {
  Matrix res;
  int r, c;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      res.at(r, c) = at(r, c) - m.at(r, c);
    }
  }
  return res;
}

Vector3 Matrix::operator*(const Vector3 &v) const {
  int r, c;
  double res[4], vv[4];
  Vector3 vRes;

  res[0] = 0.0;
  res[1] = 0.0;
  res[2] = 0.0;
  res[3] = 0.0;
  vv[0] = v.getX();
  vv[1] = v.getY();
  vv[2] = v.getZ();
  vv[3] = 1.0;
  for (r = 0; r < 3; r++) {
    for (c = 0; c < 4; c++) {
      res[r] += at(r, c) * vv[c];
    }
  }
  vRes.set(res[0], res[1], res[2]);
  return vRes;
}

Matrix Matrix::operator*(double m) const {
  Matrix res;
  int r, c;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      res.at(r, c) = at(r, c) * m;
    }
  }
  return res;
}

Matrix Matrix::invertTransform() const {
  if (!this->is3x3Orthogonal(1e-3)) {
    THROW_HARD_ERROR(("The rotation part must be orthogonal to invert the transform"));
  }
  Matrix invertRotation = this->transposed3x3();
  Vector3 trans = this->getTranslation();
  Matrix invertTranslation;
  Vector3 invTrans = trans.multiplyByScalar(-1.0);
  invertTranslation.translate(invTrans);
  Matrix invertTransform = invertRotation * invertTranslation;
  return invertTransform;
}

namespace geom {
CL_DOCSTRING(R"dx(Return a matrix that will invert the transform of the argument)dx");
DOCGROUP(cando);
CL_DEFUN
Matrix geom__invertTransform(const Matrix &matrix) { return matrix.invertTransform(); }

}; // namespace geom

Matrix Matrix::transpose() const {
  Matrix res;
  int r, c;

  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      res.at(r, c) = at(c, r);
    }
  }
  return res;
}

//
//  0  1  2  3
//  4  5  6  7
//  8  9 10 11
// 12 13 14 15
//
Matrix &Matrix::colX(Vector3 const &v) {
  this->elements[0] = v.getX();
  this->elements[4] = v.getY();
  this->elements[8] = v.getZ();
  this->elements[12] = 0.0;
  return *this;
}

Matrix &Matrix::colY(Vector3 const &v) {
  this->elements[0 + 1] = v.getX();
  this->elements[4 + 1] = v.getY();
  this->elements[8 + 1] = v.getZ();
  this->elements[12 + 1] = 0.0;
  return *this;
}
Matrix &Matrix::colZ(Vector3 const &v) {
  this->elements[0 + 2] = v.getX();
  this->elements[4 + 2] = v.getY();
  this->elements[8 + 2] = v.getZ();
  this->elements[12 + 2] = 0.0;
  return *this;
}
Matrix &Matrix::setTranslate(Vector3 const &v) {
  this->elements[0 + 3] = v.getX();
  this->elements[4 + 3] = v.getY();
  this->elements[8 + 3] = v.getZ();
  this->elements[12 + 3] = 1.0;
  return *this;
}

void Matrix::setTranslate(const Matrix &m) {
  this->elements[0 + 3] = m.elements[0 + 3];
  this->elements[4 + 3] = m.elements[4 + 3];
  this->elements[8 + 3] = m.elements[8 + 3];
  this->elements[12 + 3] = m.elements[12 + 3];
}

Vector3 Matrix::colX() const { return Vector3(this->elements[0], this->elements[4], this->elements[8]); }
Vector3 Matrix::colY() const { return Vector3(this->elements[0 + 1], this->elements[4 + 1], this->elements[8 + 1]); }
Vector3 Matrix::colZ() const { return Vector3(this->elements[0 + 2], this->elements[4 + 2], this->elements[8 + 2]); }

Matrix Matrix::transposed3x3() const {
  Matrix res;
  int r, c;
  for (r = 0; r < 3; r++) {
    for (c = 0; c < 3; c++) {
      res.at(r, c) = at(c, r);
    }
  }
  return res;
}

void Matrix::setToIndexAsValue() {
  for (int i = 0; i < 16; i++) {
    this->elements[i] = i;
  }
}

Matrix Matrix::extract3x3() const {
  Matrix res(false);
  res.elements[0] = this->elements[0];
  res.elements[1] = this->elements[1];
  res.elements[2] = this->elements[2];
  res.elements[3] = 0.0;
  res.elements[0 + 4] = this->elements[0 + 4];
  res.elements[1 + 4] = this->elements[1 + 4];
  res.elements[2 + 4] = this->elements[2 + 4];
  res.elements[3 + 4] = 0.0;
  res.elements[0 + 8] = this->elements[0 + 8];
  res.elements[1 + 8] = this->elements[1 + 8];
  res.elements[2 + 8] = this->elements[2 + 8];
  res.elements[3 + 8] = 0.0;
  res.elements[0 + 12] = 0.0;
  res.elements[1 + 12] = 0.0;
  res.elements[2 + 12] = 0.0;
  res.elements[3 + 12] = 1.0;
  return res;
}

Matrix Matrix::multiply3x3(Matrix const &m) const {
  Matrix res;
  double acc;
  int r, c, x;
  for (r = 0; r < 3; r++) {
    for (c = 0; c < 3; c++) {
      acc = 0.0;
      for (x = 0; x < 3; x++) {
        acc += this->elements[this->index(r, x)] * m.elements[this->index(x, c)];
      }
      res.elements[this->index(r, c)] = acc;
    }
  }
  return res;
}

void Matrix::addToTranslation(const Vector3 &offset) {
  this->elements[0 + 3] += offset.getX();
  this->elements[4 + 3] += offset.getY();
  this->elements[8 + 3] += offset.getZ();
}

Matrix Matrix::multiply3x3SaveTranslation(Matrix const &m) const {
  Matrix res;
  double acc;
  int r, c, x;
  for (r = 0; r < 3; r++) {
    for (c = 0; c < 3; c++) {
      acc = 0.0;
      for (x = 0; x < 3; x++) {
        acc += this->elements[this->index(r, x)] * m.elements[this->index(x, c)];
      }
      res.elements[this->index(r, c)] = acc;
    }
  }
  res.setTranslate(*this);
  return res;
}

//
//	translate
//
// Create a translation matrix
//
void Matrix::translate(const Vector3 &v) {
  this->setToIdentity();
  at(0, 3) = v.getX();
  at(1, 3) = v.getY();
  at(2, 3) = v.getZ();
  at(3, 3) = 1.0;
}

Vector3 Matrix::getTranslation() const {
  Vector3 v(this->elements[3], this->elements[7], this->elements[11]);
  return v;
}

void Matrix::dump() {
  int r, c;
  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      printf("%8.4lf ", at(r, c));
    }
    std::cout << std::endl;
  }
}

void Matrix::writeToStream(string prefix, std::ostream &out) {
  int r, c;
  for (c = 0; c < 4; c++) {
    for (r = 0; r < 4; r++) {
      out << prefix << at(r, c) << std::endl;
    }
  }
}

void Matrix::writeToXMLStream(string prefix, std::ostream &out) {
  int r, c;
  out << prefix << "<matrix4x4> ";
  for (c = 0; c < 4; c++) {
    for (r = 0; r < 4; r++) {
      out << at(r, c) << " ";
    }
  }
  out << "</matrix4x4>" << std::endl;
}

void Matrix::setFromDoubleVector(vector<double> vd) {
  int r, c, i;
  int sz;
  i = 0;
  sz = vd.size();
  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      at(r, c) = vd[i];
      i++;
      if (i >= sz)
        return; // Only use as many elements of the list
                // as provided by the caller
    }
  }
}

void Matrix::setFromDoubleArray(int iElements, double vals[]) {
  int r, c, i;
  int sz;
  i = 0;
  sz = iElements;
  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      at(r, c) = vals[i];
      i++;
      if (i >= sz)
        return; // Only use as many elements of the list
                // as provided by the caller
    }
  }
}

double Matrix::trace() {
  return this->elements[0] + this->elements[this->index(1, 1)] + this->elements[this->index(2, 2)] +
         this->elements[this->index(3, 3)];
}

bool Matrix::is3x3Orthogonal(double tol) const {
  Matrix m(this->extract3x3());
  Matrix delta((m.multiply3x3(m.transposed3x3())).subtractIdentity());
  return ((delta.colX().length() + delta.colY().length() + delta.colZ().length()) < tol);
}

string Matrix::asString() const {
  char cstr[1024];
  string out;
  int r, c;
  double val;
  stringstream ss;
  for (r = 0; r < 4; r++) {
    if (r == 0) {
      ss << "#(";
    } else {
      ss << std::endl << "  ";
    }
    for (c = 0; c < 4; c++) {
      val = at(r, c);
      if (val == 0.0) {
        strcpy(cstr, "    0.0");
      } else if (val == 1.0) {
        strcpy(cstr, "    1.0");
      } else if (val == -1.0) {
        strcpy(cstr, "   -1.0");
      } else {
        sprintf(cstr, " %6.3lf", val);
      }
      ss << cstr << " ";
    }
  }
  ss << ")";
  return ss.str();
}

string Matrix::asStringFormatted() const {
  char cstr[1024];
  string out;
  int r, c;
  double val;
  stringstream ss;
  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      val = at(r, c);
      strcpy(cstr, fmt::format("{:8.2f}", val).c_str());
      ss << cstr << " ";
    }
    ss << std::endl;
  }
  return ss.str();
}

void Matrix::setFromString(const string &str) {
  double vals[16];
  vector<string> tokens;
  vector<string>::iterator token;
  int i;
  double val;
  core::tokenize(str, tokens, " \n\t");
  if (tokens.size() < 16) {
    stringstream serr;
    serr << "Trying to convert string(" << str << ") to matrix and there aren't enough tokens" << std::endl;
    for (uint ei = 0; ei < tokens.size(); ei++) {
      serr << "Token #" << ei << " = " << tokens[ei] << std::endl;
    }
    serr << "---done dump";
    SIMPLE_ERROR("{}", serr.str());
  }
  {
    for (i = 0, token = tokens.begin(); token != tokens.end(); token++, i++) {
      val = atof((*token).c_str());
      if (i > 15) {
        SIMPLE_ERROR("Too many numbers: {}", str);
      }
      vals[i] = val;
    }
  }
  this->setFromDoubleArray(16, vals);
}

void Matrix::setValue(double val) {
  for (uint zz = 0; zz < this->elements.size(); zz++) {
    this->elements[zz] = val;
  }
}

void Matrix::setFromQuaternion(Vector4 &q) {
  this->elements[this->index(0, 0)] = q.getZ();
  this->elements[this->index(0, 1)] = -q.getY();
  this->elements[this->index(0, 2)] = q.getX();
  this->elements[this->index(0, 3)] = q.getW();
  this->elements[this->index(1, 0)] = q.getY();
  this->elements[this->index(1, 1)] = q.getZ();
  this->elements[this->index(1, 2)] = -q.getW();
  this->elements[this->index(1, 3)] = q.getX();
  this->elements[this->index(2, 0)] = -q.getX();
  this->elements[this->index(2, 1)] = q.getW();
  this->elements[this->index(2, 2)] = q.getZ();
  this->elements[this->index(2, 3)] = q.getY();
  this->elements[this->index(3, 0)] = -q.getW();
  this->elements[this->index(3, 1)] = -q.getX();
  this->elements[this->index(3, 2)] = -q.getY();
  this->elements[this->index(3, 3)] = q.getZ();
}

std::ostream &operator<<(std::ostream &out, Matrix m) {
  int r, c;

  out << "[" << std::endl;
  for (r = 0; r < 4; r++) {
    out << " [ ";
    for (c = 0; c < 4; c++) {
      out << m.at(r, c) << " ";
    }
    out << "]" << std::endl;
  }
  out << "]" << std::endl;
  return out;
}

//
//	eigenSystem
//
//	Return the eigenvectors in the rows of evecs.
//
void Matrix::eigenSystem(Vector4 &evals, Matrix &evecs) {
  Matrix evecsTemp;
  int nrot;
  double ev[4];

  jacobi(*this, 4, ev, evecs, &nrot);
  evals.set(ev[0], ev[1], ev[2], ev[3]);
}

//
//	rotationAxis
//
//	Initialize the matrix as a rotation around the unit vector
//	specified.
//	The rotation is clockwise.
//
void Matrix::rotationAxis(double a, Vector3 *v) {
  double x = v->getX();
  double y = v->getY();
  double z = v->getZ();

  // Row 0
  at(0, 0) = 1.0 + (1.0 - cos(a)) * (x * x - 1.0);
  at(0, 1) = -z * sin(a) + (1.0 - cos(a)) * x * y;
  at(0, 2) = y * sin(a) + (1.0 - cos(a)) * x * z;
  at(0, 3) = 0.0;
  // Row 1
  at(1, 0) = z * sin(a) + (1.0 - cos(a)) * x * y;
  at(1, 1) = 1.0 + (1.0 - cos(a)) * (y * y - 1.0);
  at(1, 2) = -x * sin(a) + (1.0 - cos(a)) * y * z;
  at(1, 3) = 0.0;
  // Row 2
  at(2, 0) = -y * sin(a) + (1.0 - cos(a)) * x * z;
  at(2, 1) = x * sin(a) + (1.0 - cos(a)) * y * z;
  at(2, 2) = 1.0 + (1.0 - cos(a)) * (z * z - 1.0);
  at(2, 3) = 0.0;
  // Row 3
  at(3, 0) = 0.0;
  at(3, 1) = 0.0;
  at(3, 2) = 0.0;
  at(3, 3) = 1.0;
}

//
//	rotationX
//
//	Initialize the matrix as a rotation around the X axis matrix.
//	The rotation is clockwise looking down the positive X axis.
//
void Matrix::rotationX(double a) {
  // Row 0
  at(0, 0) = 1.0;
  at(0, 1) = 0.0;
  at(0, 2) = 0.0;
  at(0, 3) = 0.0;
  // Row 1
  at(1, 0) = 0.0;
  at(1, 1) = cos(a);
  at(1, 2) = -sin(a);
  at(1, 3) = 0.0;
  // Row 2
  at(2, 0) = 0.0;
  at(2, 1) = sin(a);
  at(2, 2) = cos(a);
  at(2, 3) = 0.0;
  // Row 3
  at(3, 0) = 0.0;
  at(3, 1) = 0.0;
  at(3, 2) = 0.0;
  at(3, 3) = 1.0;
}

//
//	rotationY
//
//	Initialize the matrix as a rotation around the Y axis matrix.
//	The rotation is clockwise looking down the positive Y axis.
//
void Matrix::rotationY(double a) {
  // Row 0
  at(0, 0) = cos(a);
  at(0, 1) = 0.0;
  at(0, 2) = sin(a);
  at(0, 3) = 0.0;
  // Row 1
  at(1, 0) = 0.0;
  at(1, 1) = 1.0;
  at(1, 2) = 0.0;
  at(1, 3) = 0.0;
  // Row 2
  at(2, 0) = -sin(a);
  at(2, 1) = 0.0;
  at(2, 2) = cos(a);
  at(2, 3) = 0.0;
  // Row 3
  at(3, 0) = 0.0;
  at(3, 1) = 0.0;
  at(3, 2) = 0.0;
  at(3, 3) = 1.0;
}

//
//	rotationZ
//
//	Initialize the matrix as a rotation around the Y axis matrix.
//	The rotation is clockwise looking down the positive Y axis.
//
void Matrix::rotationZ(double a) {
  // Row 0
  at(0, 0) = cos(a);
  at(0, 1) = -sin(a);
  at(0, 2) = 0.0;
  at(0, 3) = 0.0;
  // Row 1
  at(1, 0) = sin(a);
  at(1, 1) = cos(a);
  at(1, 2) = 0.0;
  at(1, 3) = 0.0;
  // Row 2
  at(2, 0) = 0.0;
  at(2, 1) = 0.0;
  at(2, 2) = 1.0;
  at(2, 3) = 0.0;
  // Row 3
  at(3, 0) = 0.0;
  at(3, 1) = 0.0;
  at(3, 2) = 0.0;
  at(3, 3) = 1.0;
}

//
//	rightHandedRotationX
//
//	Initialize the matrix as a rightHandedRotation around the X axis matrix.
//	The rightHandedRotation is clockwise looking down the positive X axis.
//
void Matrix::rightHandedRotationX(double a) {
  // Row 0
  at(0, 0) = 1.0;
  at(0, 1) = 0.0;
  at(0, 2) = 0.0;
  at(0, 3) = 0.0;
  // Row 1
  at(1, 0) = 0.0;
  at(1, 1) = cos(-a);
  at(1, 2) = sin(-a);
  at(1, 3) = 0.0;
  // Row 2
  at(2, 0) = 0.0;
  at(2, 1) = -sin(-a);
  at(2, 2) = cos(-a);
  at(2, 3) = 0.0;
  // Row 3
  at(3, 0) = 0.0;
  at(3, 1) = 0.0;
  at(3, 2) = 0.0;
  at(3, 3) = 1.0;
}

//
//	rightHandedRotationY
//
//	Initialize the matrix as a rightHandedRotation around the Y axis matrix.
//	The rightHandedRotation is clockwise looking down the positive Y axis.
//
void Matrix::rightHandedRotationY(double a) {
  // Row 0
  at(0, 0) = cos(-a);
  at(0, 1) = 0.0;
  at(0, 2) = -sin(-a);
  at(0, 3) = 0.0;
  // Row 1
  at(1, 0) = 0.0;
  at(1, 1) = 1.0;
  at(1, 2) = 0.0;
  at(1, 3) = 0.0;
  // Row 2
  at(2, 0) = sin(-a);
  at(2, 1) = 0.0;
  at(2, 2) = cos(-a);
  at(2, 3) = 0.0;
  // Row 3
  at(3, 0) = 0.0;
  at(3, 1) = 0.0;
  at(3, 2) = 0.0;
  at(3, 3) = 1.0;
}

Matrix Matrix::flipXY() const {
  Matrix rotz;
  rotz.rotationZ(180.0 * 0.0174533);
  Matrix flipped = this->multiply3x3(rotz);
  return flipped;
}

#ifdef XML_ARCHIVE
void Matrix::archive(core::ArchiveP node) {
  if (node->loading()) {
    const string &values = node->characters();
    LOG("Loading matrix with values: [{}]", (values.c_str()));
    this->setFromString(values, node->lisp());
    LOG("Set myself using that string and the result is: {}", (this->asString().c_str()));
  } else {
    string values = this->asString();
    node->setCharacters(values);
  }
};
#endif

//
//	rightHandedRotationZ
//
//	Initialize the matrix as a rightHandedRotation around the Y axis matrix.
//	The rightHandedRotation is clockwise looking down the positive Y axis.
//
void Matrix::rightHandedRotationZ(double a) {
  // Row 0
  at(0, 0) = cos(-a);
  at(0, 1) = sin(-a);
  at(0, 2) = 0.0;
  at(0, 3) = 0.0;
  // Row 1
  at(1, 0) = -sin(-a);
  at(1, 1) = cos(-a);
  at(1, 2) = 0.0;
  at(1, 3) = 0.0;
  // Row 2
  at(2, 0) = 0.0;
  at(2, 1) = 0.0;
  at(2, 2) = 1.0;
  at(2, 3) = 0.0;
  // Row 3
  at(3, 0) = 0.0;
  at(3, 1) = 0.0;
  at(3, 2) = 0.0;
  at(3, 3) = 1.0;
}

Matrix benchmarkMatrixMultiplications(int num) {
  Matrix m1, m2, mr;
  int i;
  m1.setToIdentity();
  m2.setToIdentity();
  for (i = 0; i < num; i++) {
    mr = m1 * m2;
    m1 = m2;
  }
  return mr;
}

ArrayOfMatrices::ArrayOfMatrices() {}

ArrayOfMatrices::~ArrayOfMatrices() {}

#ifdef XML_ARCHIVE
void ArrayOfMatrices::archive(core::ArchiveP node) {
  core::VectorStrings splitData;
  Matrix m;
  uint i;
  if (node->loading()) {
    this->_Matrices.clear();
    const string &rawData = node->characters();
    core::tokenize(rawData, splitData, "|");
    for (i = 0; i < splitData.size(); i++) {
      m.setFromString(splitData[i], _lisp);
      this->_Matrices.push_back(m);
    }
  } else {
    stringstream rawData;
    for (i = 0; i < this->_Matrices.size(); i++) {
      if (i != 0) {
        rawData << "|" << std::endl;
      }
      rawData << this->_Matrices[i].asString();
    }
    node->setCharacters(rawData.str());
  }
}
#endif

Matrix &ArrayOfMatrices::get(uint i) {
  HARD_ASSERT(i < this->_Matrices.size());
  return this->_Matrices[i];
}

void quaternion_to_matrix(Matrix &matrix, double w, double x, double y, double z, double tx, double ty, double tz) {
  double xsq = x * x;
  double ysq = y * y;
  double zsq = z * z;
  double s = 1.0 / (w * w + x * x + y * y + z * z);
  matrix.at(0, 0) = 1.0 - 2.0 * s * (ysq + zsq);
  matrix.at(0, 1) = 2.0 * s * (x * y - w * z);
  matrix.at(0, 2) = 2.0 * s * (x * z + w * y);
  matrix.at(0, 3) = tx;
  matrix.at(1, 0) = 2.0 * s * (x * y + w * z);
  matrix.at(1, 1) = 1.0 - 2.0 * s * (xsq + zsq);
  matrix.at(1, 2) = 2.0 * s * (y * z - w * x);
  matrix.at(1, 3) = ty;
  matrix.at(2, 0) = 2.0 * s * (x * z - w * y);
  matrix.at(2, 1) = 2.0 * s * (y * z + w * x);
  matrix.at(2, 2) = 1.0 - 2.0 * s * (xsq + ysq);
  matrix.at(2, 3) = tz;
  matrix.at(3, 0) = 0.0;
  matrix.at(3, 1) = 0.0;
  matrix.at(3, 2) = 0.0;
  matrix.at(3, 3) = 1.0;
}

// Only use this with normalized quaternions
void normalized_quaternion_to_matrix(Matrix &matrix, double w, double x, double y, double z, double tx, double ty, double tz) {
  double xsq = x * x;
  double ysq = y * y;
  double zsq = z * z;
  matrix.at(0, 0) = 1.0 - 2.0 * (ysq + zsq);
  matrix.at(0, 1) = 2.0 * (x * y - w * z);
  matrix.at(0, 2) = 2.0 * (x * z + w * y);
  matrix.at(0, 3) = tx;
  matrix.at(1, 0) = 2.0 * (x * y + w * z);
  matrix.at(1, 1) = 1.0 - 2.0 * (xsq + zsq);
  matrix.at(1, 2) = 2.0 * (y * z - w * x);
  matrix.at(1, 3) = ty;
  matrix.at(2, 0) = 2.0 * (x * z - w * y);
  matrix.at(2, 1) = 2.0 * (y * z + w * x);
  matrix.at(2, 2) = 1.0 - 2.0 * (xsq + ysq);
  matrix.at(2, 3) = tz;
  matrix.at(3, 0) = 0.0;
  matrix.at(3, 1) = 0.0;
  matrix.at(3, 2) = 0.0;
  matrix.at(3, 3) = 1.0;
}

void rotation_matrix_to_quaternion(double &w, double &x, double &y, double &z, const Matrix &m) {
  w = sqrt(1.0 + m.at(0, 0) + m.at(1, 1) + m.at(2, 2)) / 2.0;
  ASSERT(w != 0.0);
  double wo4 = 1.0 / (4.0 * w);
  x = (m.at(2, 1) - m.at(1, 2)) * wo4;
  y = (m.at(0, 2) - m.at(2, 0)) * wo4;
  z = (m.at(1, 0) - m.at(0, 1)) * wo4;
};

void rotor3_translation_to_matrix(Matrix &matrix, double scalar, double xy, double yz, double zx, double tx, double ty, double tz) {
  double w = scalar;
  double x = -yz;
  double y = -zx;
  double z = -xy;
  double xsq = x * x;
  double ysq = y * y;
  double zsq = z * z;
  double s = 1.0 / (w * w + x * x + y * y + z * z);
  matrix.at(0, 0) = 1.0 - 2.0 * s * (ysq + zsq);
  matrix.at(0, 1) = 2.0 * s * (x * y - w * z);
  matrix.at(0, 2) = 2.0 * s * (x * z + w * y);
  matrix.at(0, 3) = tx;
  matrix.at(1, 0) = 2.0 * s * (x * y + w * z);
  matrix.at(1, 1) = 1.0 - 2.0 * s * (xsq + zsq);
  matrix.at(1, 2) = 2.0 * s * (y * z - w * x);
  matrix.at(1, 3) = ty;
  matrix.at(2, 0) = 2.0 * s * (x * z - w * y);
  matrix.at(2, 1) = 2.0 * s * (y * z + w * x);
  matrix.at(2, 2) = 1.0 - 2.0 * s * (xsq + ysq);
  matrix.at(2, 3) = tz;
  matrix.at(3, 0) = 0.0;
  matrix.at(3, 1) = 0.0;
  matrix.at(3, 2) = 0.0;
  matrix.at(3, 3) = 1.0;
}

// Only use this with normalized quaternions
void normalized_rotor3_to_matrix(Matrix &matrix, double scalar, double xy, double yz, double zx, double tx, double ty, double tz) {
  double w = -scalar;
  double x = -yz;
  double y = -zx;
  double z = -xy;
  double xsq = x * x;
  double ysq = y * y;
  double zsq = z * z;
  matrix.at(0, 0) = 1.0 - 2.0 * (ysq + zsq);
  matrix.at(0, 1) = 2.0 * (x * y - w * z);
  matrix.at(0, 2) = 2.0 * (x * z + w * y);
  matrix.at(0, 3) = tx;
  matrix.at(1, 0) = 2.0 * (x * y + w * z);
  matrix.at(1, 1) = 1.0 - 2.0 * (xsq + zsq);
  matrix.at(1, 2) = 2.0 * (y * z - w * x);
  matrix.at(1, 3) = ty;
  matrix.at(2, 0) = 2.0 * (x * z - w * y);
  matrix.at(2, 1) = 2.0 * (y * z + w * x);
  matrix.at(2, 2) = 1.0 - 2.0 * (xsq + ysq);
  matrix.at(2, 3) = tz;
  matrix.at(3, 0) = 0.0;
  matrix.at(3, 1) = 0.0;
  matrix.at(3, 2) = 0.0;
  matrix.at(3, 3) = 1.0;
}

void rotation_matrix_to_rotor3(double &scalar, double &xy, double &yz, double &zx, const Matrix &m) {
  double w = sqrt(1.0 + m.at(0, 0) + m.at(1, 1) + m.at(2, 2)) / 2.0;
  ASSERT(w != 0.0);
  double wo4 = 1.0 / (4.0 * w);
  double x = (m.at(2, 1) - m.at(1, 2)) * wo4;
  double y = (m.at(0, 2) - m.at(2, 0)) * wo4;
  double z = (m.at(1, 0) - m.at(0, 1)) * wo4;
  scalar = w;
  xy = -z;
  yz = -x;
  zx = -y;
};


#if 0
#define VEC_LOG(...)                                                                                                               \
  core::clasp_write_string(fmt::format("{}:{}:{} - ", __FILE__, __LINE__, __FUNCTION__));                                            \
  core::clasp_write_string(msg)
#else
#define VEC_LOG(...)
#endif


//
// All code for internal/external coordinate conversions
//
//
/* The definition of torsion angle with sign convention
torsion angle
Also contains definitions of: anticlinal, antiperiplanar, clinal, periplanar, synclinal, synperiplanar
    https://doi.org/10.1351/goldbook.T06406
    In a chain of atoms A-B-C-D, the dihedral angle between the plane containing
    the atoms A,B,C and that containing B,C,D. In a Newman projection the torsion
    angle is the angle (having an absolute value between 0° and 180°) between bonds
    to two specified (fiducial) groups, one from the atom nearer (proximal) to the
    observer and the other from the further (distal) atom. The torsion angle between
    groups A and D is then considered to be positive if the bond A-B is rotated in
    a clockwise direction through less than 180° in order that it may eclipse the
    bond C-D: a negative torsion angle requires rotation in the opposite sense.

The code for calculating dihedral angles and generating externals from internals
must be consistent - so I gathered it all together below.

x^ = x unit vector.

stub = [x^,y^,z^,translate] columns
x^ = (b-a)/|b-a|
z^ = Cross[x^,(a-c)]/|Cross[x^,(a-c)]|  (problem with linear angles denominator will go to zero)
y^ = Cross[z^,x^]

dihedral is defined as a->b->c->d with right hand convention.
a->b->c curl of fingers then d points parallel to thumb is positive dihedral
d = point being calculated
c = parent of d
b = parent of c (grand parent of d)
a = parent of b (great grand parent of d)

Calculating d from dist,theta,phi(dihedral)
This calculates the x,y,z using polar coordinates.
dist = distance to origin
theta = angle between x^ and d-origin
phi = dihedral angle between a-b and d-c,
      positive is clockwise from a-b to d-c looking down b-c
Positive dihedral angle
        c-----b
       /       \
      d(out)    a
  d is coming out of the plane, a,b,c are in the plane.

dx = dist*cos(theta)
dy = dist*sin(theta)*cos(phi)
dz = -dist*sin(theta)*sin(phi) ; the negative sign is VERY important here
d = stub*d   stub transforms coordinate system at origin along x^,y^,z^ to
             one where the origin is on c, x-axis is b-c and a-b is in the x,y plane.


*/
namespace geom {

DOCGROUP(cando);
CL_DOCSTRING(R"doc(Return the angle in radians between (1.0,0.0) and (dx,dy).
The result will be positive or negative depending on the quadrant of dx,dy.
dx>0,dy>0 negative;
dx<0,dy>0 negative;
dx>0,dy<0 positive;
dx<0,dy<0 positive)doc");

CL_DEFUN double geom__planeVectorAngle(double dx, double dy)
{
  double dlen = std::sqrt(dx*dx+dy*dy);
  if (dlen<SMALL_NUMBER) return 0.0;
  dx = dx/dlen;
  dy = dy/dlen;
  if (std::fabs(dx)<SMALL_NUMBER) {
    return (MY_PI/2.0) * (dy<0.0? 1.0 : -1.0);
  } else if (dy<0.0) {
    return safe_acos(dx);
  } else {
    return -safe_acos(dx);
  }
}


void internalCoordinatesFromPointAndStub(const Vector3 &D,
                                         const Matrix &coordinateSystem,
                                         double &distance,
                                         double &theta,
                                         double &phi) {
  VEC_LOG("stub = \n{}\n", coordinateSystem.asString());
  Vector3 x = coordinateSystem.colX();
  Vector3 y = coordinateSystem.colY();
  Vector3 z = coordinateSystem.colZ();
  VEC_LOG("x = {}\n", x.asString());
  VEC_LOG("y = {}\n", y.asString());
  VEC_LOG("z = {}\n", z.asString());
  Vector3 C = coordinateSystem.getTranslation();
  Vector3 CD = D - C;
  double lengthCD = CD.length();
  distance = lengthCD;
  if (lengthCD < SMALL_NUMBER)
    SIMPLE_ERROR("About to divide by zero");
  Vector3 d = CD * (1.0 / lengthCD);
  VEC_LOG("d = {}\n", d.asString());
  double dx = d.dotProduct(x);
  double dy = d.dotProduct(y);
  double dz = d.dotProduct(z);
  VEC_LOG("dx = {}  dy = {}  dz = {}\n", dx, dy, dz);
  phi = geom::geom__planeVectorAngle(dy, dz);
  if (std::isnan(phi)) core::lisp_nan_error();
  VEC_LOG("  dy = {}   dz = {}\n", dy, dz);
  VEC_LOG("_Phi = {} deg\n", (phi / 0.0174533));
  Vector3 dox(1.0, 0.0, 0.0);
  Vector3 dop(dx, dy, dz);
  double f = dop.dotProduct(dox);
  VEC_LOG("dop = {}\n", dop.asString());
  VEC_LOG("dop.dotProduct(dox) = {}\n", dop.dotProduct(dox));
  if (f > (1.0 - SMALL_NUMBER)) {
    theta = 0.0;
    return;
  }
  if (f < (-1.0 + SMALL_NUMBER)) {
    theta = M_PI;
    return;
  }
  Vector3 doz = dox.crossProduct(dop);
  doz = doz.normalized();
  VEC_LOG("doz = {}\n", doz.asString());
  Vector3 doy = doz.crossProduct(dox);
  VEC_LOG("doy = {}\n", doy.asString());
  double eox = dop.dotProduct(dox);
  double eoy = dop.dotProduct(doy);
  VEC_LOG("eox = {}  eoy = {}\n", eox, eoy);
  //  double eoz = dop.dotProduct(doz); // Must be 0.0
  theta = - geom::geom__planeVectorAngle(eox, eoy);
  if (std::isnan(theta)) core::lisp_nan_error();
  if (theta<0.0) {
    SIMPLE_WARN("Theta {} should never be less than zero", theta);
  }
  VEC_LOG("    theta = {} deg\n", (theta / 0.0174533));
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
void stubFromThreePoints(
    Matrix& transform,
    Vector3 const& parent,
    Vector3 const& grand_parent,
    Vector3 const& great_grand_parent )
{
  const Vector3& c = parent;
  const Vector3& b = grand_parent;
  const Vector3& a = great_grand_parent;
  LOG(("a = %s") , a.asString());
  LOG(("b = %s") , b.asString());
  LOG(("c = %s") , c.asString());
  Vector3 e1(b-c);
  e1 = e1.normalized();
  LOG(("e1 = (b-c).normalized : %s") , e1.asString() );
  Vector3 e3(e1.crossProduct(a-c));
  e3 = e3.normalized();
  LOG(("e3 = (e1.crossProduct(a-c)).normalized : %s") , e3.asString() );
  Vector3 e2(e3.crossProduct(e1));
  LOG(("e2 = (e3.crossProduct(e1)): %s") , e2.asString() );
  LOG(("Stub before being set: %s") , this->_Transform.asStringFormatted());
  transform.colX(e1).colY(e2).colZ(e3).setTranslate(c);
  LOG(("Stub = \n%s") , this->_Transform.asStringFormatted());
} 


DOCGROUP(cando);
CL_DEFUN Vector3 geom__pointFromStubAndInternalCoordinates(const Matrix &coordinateSystem, double distance, double angle,
                                                             double dihedral) {
  Vector3 d2;
  return pointFromStubAndInternalCoordinates(coordinateSystem, distance, angle, dihedral, d2);
};

DOCGROUP(cando);
CL_DEFUN core::DoubleFloat_mv geom__internalCoordinatesFromPointAndStub(const Vector3 &D,
                                                                        const Matrix &coordinateSystem) {
  double distance;
  double theta;
  double phi;
  internalCoordinatesFromPointAndStub(D, coordinateSystem, distance, theta, phi);
  return Values(core::clasp_make_double_float(distance), core::clasp_make_double_float(theta), core::clasp_make_double_float(phi));
}

/*! Build Vector3 external coordinate using a coordinate system in a matrix
and a distance, angle and dihedral angle.
The matrix has the columns (nx ny nz origin)
where nx, ny, nz are the orthogonal unit vectors of the coordinate system
and origin is the origin.
This uses the technique described in J Comput Chem 26: 1063–1068, 2005
but it fixes an error in the paper.  The angle here is wrt the nx axis and not
the preceeding atoms D C B.
 */
Vector3 pointFromStubAndInternalCoordinates(const Matrix &coordinateSystem, double distance, double angle, double dihedral,
                                              Vector3 &d2) {
  double cosTheta = std::cos(angle);
  double sinTheta = std::sin(angle);
  double cosPhi = std::cos(dihedral);
  double sinPhi = std::sin(dihedral);
  d2.set(distance * cosTheta,
         distance * sinTheta * cosPhi ,
         distance * sinTheta * (-sinPhi) ); 
  Vector3 d2p = coordinateSystem * d2;
  return d2p;
}


CL_DOCSTRING(R"dx(Return a vector (0 0 0))dx");
DOCGROUP(cando);
CL_DEFUN Vector3 geom__build_origin()
{
  return Vector3(0.0,0.0,0.0);
}


CL_DOCSTRING(R"dx(Return a vector along the x axis distance away from vb)dx");
DOCGROUP(cando);
CL_DEFUN Vector3 geom__build_using_bond( double distance, const Vector3& vb )
{
  Vector3	vTarget;
  vTarget = Vector3(distance,0.0,0.0);
  return vTarget+vb;
}

//! Build a vector at distance from vb and angle from va
DOCGROUP(cando);
CL_DEFUN Vector3 geom__build_using_bond_angle( double distance, const Vector3& vb,
                                               double angle, const Vector3& va)
{
  Vector3 vd, vdn, vr;
  double	ca, sa;
  vd = va-vb;
  vdn = vd.normalized();
  ca = cos(angle);
  sa = sin(angle);
  vr = Vector3(vdn.getX()*cos(angle)+vdn.getY()*sin(angle),
               -vdn.getX()*sin(angle)+vdn.getY()*cos(angle), 0.0);
  vr = vr.multiplyByScalar(distance);
  vr = vr+vb;
  return vr;
}


DOCGROUP(cando);
/* There is a regression test for this function to make sure dihedral angle sign is correct */
CL_DEFUN Vector3 geom__build_using_bond_angle_dihedral( double distance, const Vector3& vc,
                                                        double angle, const Vector3& vb,
                                                        double dihedral, const Vector3& va)
{
  Vector3 bcDir = vb-vc;
  if ( bcDir.length() == 0.0 ) return Vector3(0.0,0.0,0.0);
  Vector3 bcDirNorm = bcDir.normalized();
  Vector3 dPosDist = bcDirNorm.multiplyByScalar(distance);
    //
    // Now find the axis around which to rotate the bond angle
    //
  Vector3	abDir = va-vb;
  if ( abDir.length() == 0.0 ) return Vector3(0.0,0.0,0.0);
  Vector3 abDirNorm = abDir.normalizedOrZero();
  Vector3 angleAxis = bcDirNorm.crossProduct(abDirNorm);
  if ( angleAxis.length() == 0.0 ) return Vector3(0.0,0.0,0.0);
  Vector3 angleAxisNorm = angleAxis.normalizedOrZero();
  Matrix angleRotation;
  angleRotation.rotationAxis(angle,&angleAxisNorm);
  Vector3 dPosAngle = angleRotation.multiplyByVector3(dPosDist);

	    //
	    // Now rotate around the dihedral bond
	    //
  Matrix dihedralRotation;
  dihedralRotation.rotationAxis(-dihedral,&bcDirNorm);
  Vector3 dPosDihedral = dihedralRotation.multiplyByVector3(dPosAngle);
	    //
	    // Now translate it to atom C
	    //
  Vector3 dPos = dPosDihedral.add(vc);
  return dPos;
}

/*! Return the dihedral in radians
 */
DOCGROUP(cando);
CL_DEFUN double calculateDihedral( const Vector3& va,
                                   const Vector3& vb,
                                   const Vector3& vc,
                                   const Vector3& vd)
{
  Vector3 vab = (va - vb);
  Vector3 vcb = (vc - vb);
  Vector3 vdc = (vd - vc);
  Vector3 vacCross = (vab.crossProduct(vcb)).normalized();
  LOG("vacCross = {},{},{}" , (vacCross.getX()) , (vacCross.getY()) , (vacCross.getZ() ) );
  Vector3 vdcCross = (vdc.crossProduct(vcb)).normalized();
  LOG("vdcCross = {},{},{}" , (vdcCross.getX()) , (vdcCross.getY()) , (vdcCross.getZ() ) );
  Vector3 vCross = vacCross.crossProduct(vdcCross);
  LOG("vCross = {},{},{}" , (vCross.getX()) , (vCross.getY()) , (vCross.getZ() ) );
  double dot = vacCross.dotProduct(vdcCross);
  double dih = safe_acos(dot);
  if (isnan(dih)) {
    SIMPLE_ERROR("dihedral is NAN");
  }
  LOG("dih = {}" , (dih ) );
  double sgn = (vCross.dotProduct(vcb))<0.0?-1.0:+1.0;
//    if ( enantiomer ) return -dih*sgn;
  // The sign is applied correctly here.
  // There is a regression test for this: sign-dihedral
  return dih*sgn;
}

CL_DEFUN double calculateDihedralArray( size_t iva,
                                        size_t ivb,
                                        size_t ivc,
                                        size_t ivd,
                                        core::Array_sp array) {
  size_t len = array->length();
  int64_t max_index = len-2;
  if (max_index<=0) {
    SIMPLE_ERROR("The array does not contain 3D vectors");
  }
  if (gc::IsA<core::SimpleVector_float_sp>(array)) {
    auto sfa = gc::As<core::SimpleVector_float_sp>(array);
    Vector3 va(sfa,iva);
    Vector3 vb(sfa,ivb);
    Vector3 vc(sfa,ivc);
    Vector3 vd(sfa,ivd);
    return calculateDihedral(va,vb,vc,vd);
  }
  SIMPLE_ERROR("Handle array type {}", _rep_(array));
}

/*! Return the angle in radians
 */
DOCGROUP(cando);
CL_DEFUN double calculateAngle( const Vector3& va,
                                const Vector3& vb,
                                const Vector3& vc )
{
  Vector3	vab = (va-vb).normalized();
  Vector3	vcb = (vc-vb).normalized();
  double ang = safe_acos(vab.dotProduct(vcb));
  return ang;
}

};
