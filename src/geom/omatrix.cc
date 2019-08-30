/*
    File: omatrix.cc
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
#include <clasp/core/lispStream.h>
#include <cando/geom/omatrix.h>
#include <cando/geom/vector3.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/array.h>
//#include "core/archiveNode.h"
#include <clasp/core/wrappers.h>


namespace geom {


#ifdef XML_ARCHIVE
void OMatrix_O::archiveBase(core::ArchiveP node)
{
  node->archivePlainObject<Matrix>( "matrix",
                                    "Matrix", this->_Value);
}
#endif

CL_NAME("MAKE-MATRIX");
CL_DEFUN OMatrix_sp OMatrix_O::make(bool ident)
{
  GC_ALLOCATE_VARIADIC(OMatrix_O,om,ident);
  return om;
};
  

string	OMatrix_O::__repr__() const
{
  stringstream ss;
  ss << this->_Value.asString();
  return ss.str();
}

CL_LISPIFY_NAME("setFromString");
CL_DEFMETHOD     void	OMatrix_O::setFromString(const string& s)
{
  this->_Value.setFromString(s);
}

#if 0
void	OMatrix_O::setFromStringFast(const string& s)
{
  this->_Value.setFromStringFast(s);
}
#endif

CL_DEFMETHOD void OMatrix_O::setFromVector(core::Array_sp array)
{
  if (gc::IsA<core::SimpleVector_double_sp>(array)) {
    core::SimpleVector_double_sp svd = gc::As_unsafe<core::SimpleVector_double_sp>(array);
    size_t elements = MIN(16,svd->length());
    for ( size_t i=0; i<elements; ++i ) {
      this->_Value[i] = (*svd)[i];
    }
  } else {
    size_t elements = MIN(16,array->arrayTotalSize());
    for ( size_t i=0; i<elements; ++i ) {
      core::T_sp el = array->rowMajorAref(i);
      core::Number_sp num = gc::As<core::Number_sp>(el);
      double val = core::clasp_to_double(num);
      this->_Value[i] = val;
    }
  }
}
    
void	OMatrix_O::setAll(const Matrix& m)
{
  this->_Value.setAll(m);
}

//string OMatrix_O::asStringFast() { return this->_Value.asStringFast(); }


CL_LISPIFY_NAME("Matrix-asString");
CL_DEFMETHOD     string OMatrix_O::asString() { return this->_Value.asString(); }




void OMatrix_O::setToIdentity() { this->_Value.setToIdentity(); }
    /*! Set the matrix element index as the value for testing */
void OMatrix_O::setToIndexAsValue() { this->_Value.setToIndexAsValue();};
void OMatrix_O::translate(Vector3 pos) { this->_Value.translate(pos); }
void OMatrix_O::rotationX(double radians) { this->_Value.rotationX(radians); }
void OMatrix_O::rotationY(double radians) { this->_Value.rotationY(radians); }
void OMatrix_O::rotationZ(double radians) { this->_Value.rotationZ(radians); }



CL_LISPIFY_NAME("make-matrix-identity");
CL_DEFUN OMatrix_sp make_matrix_identity()
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.setToIdentity();
  return res;
}

CL_LISPIFY_NAME("make-m4-translate");
CL_DEFUN OMatrix_sp make_m4_translate(const Vector3& pos)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.translate(pos);
  return res;
}

CL_LISPIFY_NAME("make-m4-rotate-x");
CL_DEFUN OMatrix_sp make_m4_rotate_x(double radians)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.rotationX(radians);
  return res;
}

CL_LISPIFY_NAME("make-m4-rotate-y");
CL_DEFUN OMatrix_sp make_m4_rotate_y(double radians)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.rotationY(radians);
  return res;
}

CL_LISPIFY_NAME("make-m4-rotate-z");
CL_DEFUN OMatrix_sp make_m4_rotate_z(double radians)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.rotationZ(radians);
  return res;
}

OMatrix_sp OMatrix_O::clone() const
{
  GC_COPY(OMatrix_O,clone,*this);
  return clone;
}

CL_NAME("M*M");
CL_DEFMETHOD Matrix OMatrix_O::m_TIMES_m(const Matrix& other) const
{
  Matrix res = this->_Value*other;
  return res;
};

CL_NAME("M*V")
CL_DEFMETHOD Vector3 OMatrix_O::m_TIMES_v(const Vector3& vec) const
{
  Vector3 res = this->_Value*vec;
  return res;
}


core::List_sp OMatrix_O::encode() {
  core::Vector_sp v = core::core__make_vector(cl::_sym_DoubleFloat_O,16);
  for ( size_t i(0); i<16; ++i ) {
    v->rowMajorAset(i,core::clasp_make_double_float(this->_Value[i]));
  }
  return core::Cons_O::create(core::Cons_O::create(INTERN_(kw,m),v),_Nil<core::T_O>());
}

void OMatrix_O::decode(core::List_sp c) {
  core::Vector_sp v = gc::As<core::Vector_sp>(oCdr(oCar(c)));
  for ( size_t i(0); i<16; ++i ) {
    this->_Value[i] = core::clasp_to_double(v->rowMajorAref(i));
  }
}

CL_DOCSTRING("Return the (values EIGEN-VALUES EIGEN-VECTOR-MATRIX) of the MATRIX.");
CL_DEFUN core::T_mv geom__eigen_system(OMatrix_sp matrix)
{
  Vector4 evals;
  Matrix evecs;
  matrix->_Value.eigenSystem(evals,evecs);
  core::List_sp lvals = core::Cons_O::createList(core::clasp_make_double_float(evals.getW()),
                                                 core::clasp_make_double_float(evals.getX()),
                                                 core::clasp_make_double_float(evals.getY()),
                                                 core::clasp_make_double_float(evals.getZ()));
  geom::OMatrix_sp mm = geom::OMatrix_O::create();
  mm->setAll(evecs);
  return Values(mm,lvals);
};

};
