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
DOCGROUP(cando);
CL_DEFUN OMatrix_sp OMatrix_O::make(bool ident)
{
  auto om = gctools::GC<OMatrix_O>::allocate(ident);
  return om;
};

OMatrix_sp OMatrix_O::m4(core::List_sp args)
{
  if (core::cl__length(args)>16) {
    SIMPLE_ERROR("The list of arguments {} is too long to fit in a 16 element homogeneous matrix", _rep_(args));
  }
  int index = 0;
  auto om = gctools::GC<OMatrix_O>::allocate(true);
  for ( auto cur : args ) {
    core::T_sp obj = CONS_CAR(cur);
    if (gc::IsA<core::Number_sp>(obj)) {
      core::Number_sp num = gc::As_unsafe<core::Number_sp>(obj);
      double val = core::clasp_to_double(num);
      om->_Value[index] = val;
      core::clasp_write_string(fmt::format("{} ", val ) );
    } else {
      TYPE_ERROR(obj,cl::_sym_Number_O);
    }
    index = index + 1;
  }
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

CL_DEFMETHOD void OMatrix_O::setTranslation(const Vector3& pos)
{
  this->_Value[3] = pos.getX();
  this->_Value[7] = pos.getY();
  this->_Value[11] = pos.getZ();
}

CL_DEFMETHOD Vector3 OMatrix_O::getTranslation()
{
  Vector3 pos;
  pos.set( this->_Value[3], this->_Value[7], this->_Value[11]);
  return pos;
}
/*
   0  1  2  3
   4  5  6  7
   8  9 10 11
*/
CL_DEFMETHOD Vector3 OMatrix_O::getX_Column() const
{
  Vector3 pos;
  pos.set( this->_Value[0], this->_Value[4], this->_Value[8]);
  return pos;
}

CL_DEFMETHOD Vector3 OMatrix_O::getY_Column() const
{
  Vector3 pos;
  pos.set( this->_Value[0+1], this->_Value[4+1], this->_Value[8+1]);
  return pos;
}

CL_DEFMETHOD Vector3 OMatrix_O::getZ_Column() const
{
  Vector3 pos;
  pos.set( this->_Value[0+2], this->_Value[4+2], this->_Value[8+2]);
  return pos;
}

CL_DEFMETHOD Matrix OMatrix_O::flipXY() const {
  Matrix flipped;
  flipped = this->_Value.flipXY();
  return flipped;
}

CL_LISPIFY_NAME("geom:at");
CL_DEFUN_SETF double setf_at(double v, OMatrix_sp obj, int row, int col) {
  return (obj->_Value.at(row, col) = v);
};

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
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_matrix_identity()
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.setToIdentity();
  return res;
}

CL_LISPIFY_NAME("make-m4-identity");
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_m4_identity()
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.setToIdentity();
  return res;
}


CL_LISPIFY_NAME("make-matrix-scale");
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_matrix_scale(double scale)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.scale(scale);
  return res;
}

CL_LISPIFY_NAME("make-m4-translate");
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_m4_translate(const Vector3& pos)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.translate(pos);
  return res;
}

CL_LISPIFY_NAME("make-m4-rotation-rows");
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_m4_rotation_rows(const Vector3& row0, const Vector3& row1, const Vector3& row2)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.setRow(0
                     ,row0.getX()
                     ,row0.getY()
                     ,row0.getZ()
                     ,0.0);
  res->_Value.setRow(1
                     ,row1.getX()
                     ,row1.getY()
                     ,row1.getZ()
                     ,0.0);
  res->_Value.setRow(2
                     ,row2.getX()
                     ,row2.getY()
                     ,row2.getZ()
                     ,0.0);
  res->_Value.setRow(3
                     , 0.0
                     , 0.0
                     , 0.0
                     , 1.0 );
  return res;
}

CL_LISPIFY_NAME("make-m4-rows");
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_m4_rows(const Vector3& row0, const Vector3& row1, const Vector3& row2, const Vector3& toOriginTranslation)
{
  Matrix toOrigin;
  toOrigin.translate(toOriginTranslation);
  Matrix rotate;
  rotate.setRow(0
                ,row0.getX()
                ,row0.getY()
                ,row0.getZ()
                ,0.0);
  rotate.setRow(1
                ,row1.getX()
                ,row1.getY()
                ,row1.getZ()
                ,0.0);
  rotate.setRow(2
                ,row2.getX()
                ,row2.getY()
                ,row2.getZ()
                ,0.0);
  rotate.setRow(3
                , 0.0
                , 0.0
                , 0.0
                , 1.0 );
  OMatrix_sp res = OMatrix_O::make(false);
  res->_Value = rotate.multiplyByMatrix(toOrigin);
  return res;
}

CL_DOCSTRING(R"doc(Decompose the MATRIX into a T_translation x T_rotation pair)doc");
CL_DEFUN OMatrix_mv geom__decompose_into_translation_rotation(OMatrix_sp matrix) {
  Matrix mrot = matrix->_Value.extract3x3();
  Matrix mtrans;
  mtrans.setToIdentity();
  mtrans.setTranslate(matrix->_Value);
  OMatrix_sp omrot = OMatrix_O::create();
  omrot->_Value = mrot;
  OMatrix_sp omtrans = OMatrix_O::create();
  omtrans->_Value = mtrans;
  return Values(omtrans,omrot);
}


CL_LISPIFY_NAME("make-m4-rotation-columns");
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_m4_rotation_columns(const Vector3& col0, const Vector3& col1, const Vector3& col2)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.setCol(0
                     ,col0.getX()
                     ,col0.getY()
                     ,col0.getZ()
                     ,0.0);
  res->_Value.setCol(1
                     ,col1.getX()
                     ,col1.getY()
                     ,col1.getZ()
                     ,0.0);
  res->_Value.setCol(2
                     ,col2.getX()
                     ,col2.getY()
                     ,col2.getZ()
                     ,0.0);
  res->_Value.setCol(3
                     , 0.0
                     , 0.0
                     , 0.0
                     , 1.0 );
  return res;
}


CL_LISPIFY_NAME("make-m4-columns");
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_m4_columns(const Vector3& col0, const Vector3& col1, const Vector3& col2, const Vector3& col3)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.setCol(0
                     ,col0.getX()
                     ,col0.getY()
                     ,col0.getZ()
                     ,0.0);
  res->_Value.setCol(1
                     ,col1.getX()
                     ,col1.getY()
                     ,col1.getZ()
                     ,0.0);
  res->_Value.setCol(2
                     ,col2.getX()
                     ,col2.getY()
                     ,col2.getZ()
                     ,0.0);
  res->_Value.setCol(3
                     ,col3.getX()
                     ,col3.getY()
                     ,col3.getZ()
                     ,1.0);
  return res;
}

CL_LISPIFY_NAME("make-m4");
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_m4(core::T_sp obj) {
  OMatrix_sp res = OMatrix_O::create();
  if (gc::IsA<core::SimpleVector_sp>(obj)) {
    auto svobj = gc::As_unsafe<core::SimpleVector_sp>(obj);
    if (svobj->length()!=16) goto ERROR;
    for ( size_t idx = 0; idx<16; idx++ ) {
      res->_Value.elements[idx] = core::clasp_to_double((*svobj)[idx]);
    }
    return res;
  }
 ERROR:
  TYPE_ERROR(obj,core::Cons_O::createList(core::Cons_O::createList(cl::_sym_simple_vector),core::make_fixnum(16)));
}



CL_LISPIFY_NAME("make-m4-rotate-x");
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_m4_rotate_x(double radians)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.rotationX(radians);
  return res;
}

CL_LISPIFY_NAME("make-m4-rotate-y");
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_m4_rotate_y(double radians)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.rotationY(radians);
  return res;
}

CL_LISPIFY_NAME("make-m4-rotate-z");
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_m4_rotate_z(double radians)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.rotationZ(radians);
  return res;
}

CL_LISPIFY_NAME("make-m4-rotate-axis");
DOCGROUP(cando);
CL_DEFUN OMatrix_sp make_m4_rotate_axis(double radians,Vector3 axis)
{
  OMatrix_sp res = OMatrix_O::create();
  res->_Value.rotationAxis(radians,&axis);
  return res;
}

CL_LISPIFY_NAME("copy-matrix");
DOCGROUP(cando);
CL_DEFMETHOD OMatrix_sp OMatrix_O::copyMatrix() const
{
  auto clone = gctools::GC<OMatrix_O>::copy(*this);
  return clone;
}

CL_NAME("M+");
CL_DEFMETHOD Matrix OMatrix_O::m_PLUS_(const Matrix& other) const
{
  return _Value + other;
};

CL_NAME("M+!");
CL_DEFMETHOD void OMatrix_O::m_PLUS__BANG_(const Matrix& other)
{
  _Value += other;
};

CL_NAME("M*M");
CL_DEFMETHOD Matrix OMatrix_O::m_TIMES_m(const Matrix& other) const
{
  Matrix res = this->_Value*other;
  return res;
};

CL_NAME("M*V");
CL_DEFMETHOD Vector3 OMatrix_O::m_TIMES_v(const Vector3& vec) const
{
  Vector3 res = this->_Value*vec;
  return res;
}


core::List_sp OMatrix_O::encode() {
  core::Vector_sp v = core::core__make_vector(cl::_sym_DoubleFloat_O,16);
  for ( size_t i(0); i<12; ++i ) {
    v->rowMajorAset(i,core::clasp_make_double_float(this->_Value[i]));
  }
  return core::Cons_O::create(core::Cons_O::create(INTERN_(kw,m),v),nil<core::T_O>());
}

void OMatrix_O::decode(core::List_sp c) {
  core::Vector_sp v = gc::As<core::Vector_sp>(oCdr(oCar(c)));
  for ( size_t i(0); i<12; ++i ) {
    this->_Value[i] = core::clasp_to_double(v->rowMajorAref(i));
  }
  this->_Value[12] = 0.0;
  this->_Value[13] = 0.0;
  this->_Value[14] = 0.0;
  this->_Value[15] = 1.0;
}

CL_DOCSTRING(R"dx(Return the (values EIGEN-VALUES EIGEN-VECTOR-MATRIX) of the MATRIX.)dx");
DOCGROUP(cando);
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


CL_DOCSTRING(R"dx(Convert a homogeneous matrix into a quaternion/translation. This takes the matrix and returns seven numbers, the w,x,y,z,tx,ty,tz where w,x,y,z are the quaternion and the tx,ty,tz is the translation.)dx");
DOCGROUP(cando);
CL_DEFUN core::T_mv geom__matrix_to_quaternion_translation(OMatrix_sp om)
{
  Matrix& m = om->_Value;
  double tx = m.elements[0+3];
  double ty = m.elements[4+3];
  double tz = m.elements[8+3];
  double qw,qx,qy,qz;
  rotation_matrix_to_quaternion(qw,qx,qy,qz,m);
  return Values(core::clasp_make_double_float(qw),
                core::clasp_make_double_float(qx),
                core::clasp_make_double_float(qy),
                core::clasp_make_double_float(qz),
                core::clasp_make_double_float(tx),
                core::clasp_make_double_float(ty),
                core::clasp_make_double_float(tz));
};


CL_DOCSTRING(R"dx(Set the matrix from a quaternion/translation - a set of seven numbers qw,qx,qy,qz,tx,ty,tz)dx");
CL_DEFMETHOD void OMatrix_O::set_from_quaternion(double qw, double qx, double qy, double qz, double tx, double ty, double tz) {
  quaternion_to_matrix(this->_Value,qw,qx,qy,qz,tx,ty,tz);
}

CL_DEFMETHOD void OMatrix_O::set_from_normalized_quaternion(double qw, double qx, double qy, double qz, double tx, double ty, double tz) {
  normalized_quaternion_to_matrix(this->_Value,qw,qx,qy,qz,tx,ty,tz);
}

CL_DEFMETHOD core::T_mv OMatrix_O::rotation_to_quaternion() const {
  double qw, qx, qy, qz;
  rotation_matrix_to_quaternion(qw,qx,qy,qz,this->_Value);
  return Values(core::DoubleFloat_O::create(qw),
                core::DoubleFloat_O::create(qx),
                core::DoubleFloat_O::create(qy),
                core::DoubleFloat_O::create(qz));
};

CL_DEFMETHOD void OMatrix_O::set_from_normalized_rotor3(double rscalar, double rxy, double ryz, double rzx, double tx, double ty, double tz) {
  normalized_rotor3_to_matrix(this->_Value,rscalar,rxy,ryz,rzx,tx,ty,tz);
}

CL_DEFMETHOD core::T_mv OMatrix_O::rotation_to_rotor3() const {
  double rscalar, rxy, ryz, rzx;
  rotation_matrix_to_rotor3(rscalar,rxy,ryz,rzx,this->_Value);
  return Values(core::DoubleFloat_O::create(rscalar),
                core::DoubleFloat_O::create(rxy),
                core::DoubleFloat_O::create(ryz),
                core::DoubleFloat_O::create(rzx));
};


CL_DOCSTRING(R"dx(Convert a homogeneous matrix into a rotor3/translation. This takes the matrix and returns seven numbers, the scalar,xy,yz,zx,tx,ty,tz where w,x,y,z are the rotor3 components and the tx,ty,tz is the translation.)dx");
DOCGROUP(cando);
CL_DEFUN core::T_mv geom__matrix_to_rotor3_translation(OMatrix_sp om)
{
  Matrix& m = om->_Value;
  double tx = m.elements[0+3];
  double ty = m.elements[4+3];
  double tz = m.elements[8+3];
  double rscalar,rxy, ryz, rzx;
  rotation_matrix_to_rotor3(rscalar,rxy,ryz,rzx,m);
  return Values(core::clasp_make_double_float(rscalar),
                core::clasp_make_double_float(rxy),
                core::clasp_make_double_float(ryz),
                core::clasp_make_double_float(rzx),
                core::clasp_make_double_float(tx),
                core::clasp_make_double_float(ty),
                core::clasp_make_double_float(tz));
};

};
